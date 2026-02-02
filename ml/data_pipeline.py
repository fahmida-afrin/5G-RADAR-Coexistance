"""
Data Pipeline for 5G-Radar Coexistence ML Models.

Converts MATLAB simulation outputs (.mat / .csv) into PyTorch Geometric
graph objects suitable for GNN training, and flat tensors for DQN training.

Simulation parameters (from NRwithRADAR.m):
    - NumRBs = 25, SCS = 15 kHz, 4 UEs
    - TDD: 2 DL slots + 8 DL sym, 4 UL sym + 2 UL slots (5 ms period)
    - Radar: PRF 1 kHz, 70 μs pulse, 5 MHz BW, 3 MHz offset, 1 kW @ 5 km
"""

import os
import glob
import numpy as np
import pandas as pd
from pathlib import Path

try:
    import scipy.io as sio
except ImportError:
    sio = None

try:
    import torch
    from torch.utils.data import Dataset, DataLoader
except ImportError:
    torch = None

try:
    from torch_geometric.data import Data, InMemoryDataset
except ImportError:
    Data = None
    InMemoryDataset = None

# ---------------------------------------------------------------------------
# Constants matching NRwithRADAR.m
# ---------------------------------------------------------------------------
NUM_RBS = 25
NUM_UES = 4
SCS_KHZ = 15
SLOTS_PER_FRAME = 10  # for 15 kHz SCS
SYMBOLS_PER_SLOT = 14
DL_CARRIER_FREQ = 3.595e9
RADAR_FREQ_OFFSET = 3e6
RADAR_PRF = 1e3
RADAR_PULSE_WIDTH = 70e-6
RADAR_BW = 5e6
RADAR_POWER = 1000  # W
RADAR_DISTANCE = 5000  # m


# ---------------------------------------------------------------------------
# 1. Loading helpers
# ---------------------------------------------------------------------------

def load_matlab_data(filepath: str) -> dict:
    """Load a .mat or .csv file exported from MATLAB simulation.

    Parameters
    ----------
    filepath : str
        Path to .mat or .csv file.

    Returns
    -------
    dict
        Keys depend on source file; typical keys include
        ``scheduling_logs``, ``rlc_logs``, ``cqi``, ``bler``, ``throughput``.
    """
    ext = Path(filepath).suffix.lower()
    if ext == ".mat":
        if sio is None:
            raise ImportError("scipy is required to read .mat files")
        raw = sio.loadmat(filepath, squeeze_me=True, struct_as_record=True)
        # Drop MATLAB meta keys
        return {k: v for k, v in raw.items() if not k.startswith("__")}
    elif ext == ".csv":
        df = pd.read_csv(filepath)
        return {"dataframe": df}
    else:
        raise ValueError(f"Unsupported file type: {ext}")


def load_all_scenario_data(data_dir: str) -> dict:
    """Load all scenario results from a directory.

    Expected layout::

        data_dir/
            no_radar/        *.mat | *.csv
            radar_rr/        ...
            radar_pf/        ...
            radar_bestcqi/   ...

    Returns
    -------
    dict[str, dict]
        Scenario name → loaded data dict.
    """
    scenarios = {}
    base = Path(data_dir)
    for scenario_dir in sorted(base.iterdir()):
        if scenario_dir.is_dir():
            files = list(scenario_dir.glob("*.mat")) + list(scenario_dir.glob("*.csv"))
            if files:
                scenario_data = {}
                for f in files:
                    scenario_data[f.stem] = load_matlab_data(str(f))
                scenarios[scenario_dir.name] = scenario_data
    return scenarios


# ---------------------------------------------------------------------------
# 2. Graph construction (PRB-level)
# ---------------------------------------------------------------------------

def create_prb_graph(
    cqi_matrix: np.ndarray,
    radar_active: np.ndarray,
    bler_per_ue: np.ndarray | None = None,
    throughput_per_ue: np.ndarray | None = None,
    allocation: np.ndarray | None = None,
) -> "Data":
    """Build a PyTorch-Geometric ``Data`` object from one scheduling snapshot.

    Graph topology
    --------------
    - **Nodes**: one per PRB (``NUM_RBS`` nodes).
    - **Edges**: each PRB is connected to its immediate neighbours (chain)
      *plus* edges to all PRBs that share the same UE allocation (if known).

    Node features  (per PRB)
    -------------------------
    0. PRB index (normalised)
    1-4. CQI reported by each UE for this PRB
    5. Binary radar-active flag for this PRB
    6-9. BLER per UE (if available, else 0)

    Parameters
    ----------
    cqi_matrix : (NUM_UES, NUM_RBS) float
    radar_active : (NUM_RBS,) binary – 1 if radar overlaps this PRB
    bler_per_ue : (NUM_UES,) float, optional
    throughput_per_ue : (NUM_UES,) float, optional
    allocation : (NUM_RBS,) int, optional – UE index allocated per PRB (-1 = free)

    Returns
    -------
    torch_geometric.data.Data
    """
    if Data is None:
        raise ImportError("torch_geometric is required")

    num_prbs = cqi_matrix.shape[1]
    # -- Node features --
    prb_idx = np.arange(num_prbs, dtype=np.float32) / num_prbs
    cqi = cqi_matrix.astype(np.float32) / 15.0  # normalise CQI [0-15] → [0-1]
    radar = radar_active.astype(np.float32).reshape(1, -1)

    if bler_per_ue is None:
        bler_per_ue = np.zeros(NUM_UES, dtype=np.float32)
    bler_feat = np.tile(bler_per_ue.reshape(-1, 1), (1, num_prbs)).astype(np.float32)

    # Stack: (10, NUM_RBS)  →  transpose to (NUM_RBS, 10)
    x = np.vstack([prb_idx.reshape(1, -1), cqi, radar, bler_feat]).T

    # -- Edges (chain + allocation) --
    src, dst = [], []
    for i in range(num_prbs - 1):
        src += [i, i + 1]
        dst += [i + 1, i]

    if allocation is not None:
        for ue_id in range(NUM_UES):
            prbs_for_ue = np.where(allocation == ue_id)[0]
            for a in prbs_for_ue:
                for b in prbs_for_ue:
                    if a != b:
                        src.append(a)
                        dst.append(b)

    edge_index = np.array([src, dst], dtype=np.int64)

    # -- Labels: optimal allocation vector (NUM_RBS,) --
    if allocation is not None:
        y = torch.tensor(allocation, dtype=torch.long)
    else:
        y = torch.zeros(num_prbs, dtype=torch.long)

    data = Data(
        x=torch.tensor(x, dtype=torch.float),
        edge_index=torch.tensor(edge_index, dtype=torch.long),
        y=y,
    )
    # Attach extra metadata
    if throughput_per_ue is not None:
        data.throughput = torch.tensor(throughput_per_ue, dtype=torch.float)

    return data


# ---------------------------------------------------------------------------
# 3. Label generation
# ---------------------------------------------------------------------------

def generate_labels(
    cqi_matrix: np.ndarray,
    radar_active: np.ndarray,
    qos_weights: np.ndarray | None = None,
) -> np.ndarray:
    """Generate *oracle* allocation labels by greedy best-CQI with radar avoidance.

    For each PRB, assign to the UE with highest CQI **unless** radar is
    active on that PRB – in which case mark as unallocated (-1).

    Parameters
    ----------
    cqi_matrix : (NUM_UES, NUM_RBS)
    radar_active : (NUM_RBS,)
    qos_weights : (NUM_UES,), optional – per-UE priority weights

    Returns
    -------
    allocation : (NUM_RBS,) int  – UE index or -1
    """
    num_prbs = cqi_matrix.shape[1]
    if qos_weights is None:
        qos_weights = np.ones(NUM_UES)

    weighted_cqi = cqi_matrix * qos_weights.reshape(-1, 1)
    allocation = np.argmax(weighted_cqi, axis=0).astype(np.int64)
    allocation[radar_active.astype(bool)] = -1
    return allocation


# ---------------------------------------------------------------------------
# 4. Radar overlap helper
# ---------------------------------------------------------------------------

def compute_radar_overlap(
    slot_index: int,
    symbol_index: int,
    num_prbs: int = NUM_RBS,
    radar_prf: float = RADAR_PRF,
    radar_pw: float = RADAR_PULSE_WIDTH,
    scs_khz: float = SCS_KHZ,
    radar_bw: float = RADAR_BW,
    radar_freq_offset: float = RADAR_FREQ_OFFSET,
    nr_bw: float = 5e6,
) -> np.ndarray:
    """Return binary vector (num_prbs,) indicating which PRBs overlap the radar pulse.

    Calculates time & frequency overlap of a radar pulse at the given
    slot/symbol position.
    """
    symbol_duration = 1.0 / (scs_khz * 1e3)  # ~66.7 μs for 15 kHz
    slot_duration = SYMBOLS_PER_SLOT * symbol_duration
    abs_time = slot_index * slot_duration + symbol_index * symbol_duration

    # Radar pulse timing
    pri = 1.0 / radar_prf
    time_in_pri = abs_time % pri
    pulse_active = time_in_pri < radar_pw

    if not pulse_active:
        return np.zeros(num_prbs, dtype=np.float32)

    # Frequency overlap: each PRB is 12 subcarriers × SCS
    prb_bw = 12 * scs_khz * 1e3
    prb_centers = np.arange(num_prbs) * prb_bw + prb_bw / 2

    # Radar band relative to NR carrier centre
    radar_lo = radar_freq_offset - radar_bw / 2
    radar_hi = radar_freq_offset + radar_bw / 2

    # NR PRBs are centred at 0; shift to absolute
    nr_lo = -nr_bw / 2
    prb_abs = nr_lo + prb_centers

    overlap = ((prb_abs + prb_bw / 2) > radar_lo) & ((prb_abs - prb_bw / 2) < radar_hi)
    return overlap.astype(np.float32)


# ---------------------------------------------------------------------------
# 5. Dataset builder
# ---------------------------------------------------------------------------

def build_dataset_from_logs(
    scheduling_logs: list[dict],
    cqi_logs: list[np.ndarray],
    bler_logs: list[np.ndarray] | None = None,
    throughput_logs: list[np.ndarray] | None = None,
) -> list:
    """Convert a sequence of per-slot log entries into a list of PyG Data objects."""
    dataset = []
    for t, (sched, cqi) in enumerate(zip(scheduling_logs, cqi_logs)):
        slot_idx = sched.get("slot", t)
        sym_idx = sched.get("symbol", 0)
        radar = compute_radar_overlap(slot_idx, sym_idx)
        alloc = sched.get("allocation", None)
        if alloc is None:
            alloc = generate_labels(cqi, radar)
        bler = bler_logs[t] if bler_logs else None
        tp = throughput_logs[t] if throughput_logs else None
        data = create_prb_graph(cqi, radar, bler, tp, alloc)
        data.slot = slot_idx
        dataset.append(data)
    return dataset


def split_dataset(data_list: list, train_ratio: float = 0.8, seed: int = 42):
    """Split list of Data objects into train / test."""
    rng = np.random.RandomState(seed)
    idx = rng.permutation(len(data_list))
    split = int(len(data_list) * train_ratio)
    train = [data_list[i] for i in idx[:split]]
    test = [data_list[i] for i in idx[split:]]
    return train, test


# ---------------------------------------------------------------------------
# 6. Synthetic data generator (for testing without MATLAB)
# ---------------------------------------------------------------------------

def generate_synthetic_dataset(
    num_samples: int = 1000,
    num_ues: int = NUM_UES,
    num_prbs: int = NUM_RBS,
    radar_duty_cycle: float = 0.07,
    seed: int = 42,
) -> list:
    """Generate synthetic training data mimicking MATLAB simulation outputs.

    Useful for verifying ML pipeline before real data is available.
    """
    rng = np.random.RandomState(seed)
    dataset = []

    for i in range(num_samples):
        # Random CQI per UE per PRB  [1..15]
        cqi = rng.randint(1, 16, size=(num_ues, num_prbs)).astype(np.float32)

        # Radar active on contiguous block with given duty cycle
        radar = np.zeros(num_prbs, dtype=np.float32)
        if rng.random() < radar_duty_cycle * 10:  # scale up for per-slot
            start = rng.randint(0, num_prbs)
            width = rng.randint(1, max(2, int(num_prbs * 0.4)))
            radar[start: min(start + width, num_prbs)] = 1.0

        alloc = generate_labels(cqi, radar)
        bler = rng.uniform(0, 0.3, size=num_ues).astype(np.float32)
        tp = rng.uniform(0, 50e6, size=num_ues).astype(np.float32)

        data = create_prb_graph(cqi, radar, bler, tp, alloc)
        data.slot = i
        dataset.append(data)

    return dataset


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="5G-Radar data pipeline")
    sub = parser.add_subparsers(dest="cmd")

    p_syn = sub.add_parser("synthetic", help="Generate synthetic dataset")
    p_syn.add_argument("--samples", type=int, default=2000)
    p_syn.add_argument("--output", type=str, default="data/processed/synthetic.pt")

    p_load = sub.add_parser("load", help="Load MATLAB data")
    p_load.add_argument("--dir", type=str, default="data/raw")
    p_load.add_argument("--output", type=str, default="data/processed/dataset.pt")

    args = parser.parse_args()

    if args.cmd == "synthetic":
        print(f"Generating {args.samples} synthetic samples...")
        ds = generate_synthetic_dataset(num_samples=args.samples)
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        if torch is not None:
            torch.save(ds, args.output)
            print(f"Saved to {args.output}")
        else:
            print("PyTorch not available; cannot save .pt file")

    elif args.cmd == "load":
        scenarios = load_all_scenario_data(args.dir)
        print(f"Loaded {len(scenarios)} scenarios: {list(scenarios.keys())}")
    else:
        parser.print_help()
