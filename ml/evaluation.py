"""
Evaluation Framework for 5G-Radar Spectrum Management Models.

Metrics:
  1. Spectrum Efficiency (bits/s/Hz per PRB)
  2. Radar Protection (interference power leaked to radar band)
  3. Evacuation Time (steps to clear radar-affected PRBs)
  4. Throughput (total and per-UE)
  5. Goodput (successful delivery rate)
  6. BLER (block error rate under different conditions)
  7. Fairness (Jain's fairness index)
"""

import numpy as np
from typing import Any

try:
    import torch
    from torch_geometric.data import DataLoader as PyGDataLoader
except ImportError:
    torch = None

from ml.data_pipeline import NUM_RBS, NUM_UES, compute_radar_overlap


# ---------------------------------------------------------------------------
# CQI → Spectral Efficiency mapping (3GPP TS 38.214 Table 5.2.2.1-3)
# ---------------------------------------------------------------------------
CQI_TO_SE = {
    0: 0.0,
    1: 0.1523,
    2: 0.2344,
    3: 0.3770,
    4: 0.6016,
    5: 0.8770,
    6: 1.1758,
    7: 1.4766,
    8: 1.9141,
    9: 2.4063,
    10: 2.7305,
    11: 3.3223,
    12: 3.9023,
    13: 4.5234,
    14: 5.1152,
    15: 5.5547,
}


# ---------------------------------------------------------------------------
# Metric computation
# ---------------------------------------------------------------------------

def compute_spectrum_efficiency(
    allocation: np.ndarray,
    cqi_matrix: np.ndarray,
) -> dict:
    """Compute spectral efficiency for an allocation.

    Parameters
    ----------
    allocation : (NUM_RBS,) int – UE index per PRB (-1 or NUM_UES = unallocated)
    cqi_matrix : (NUM_UES, NUM_RBS) float

    Returns
    -------
    dict with total_se, per_ue_se, avg_se_per_prb
    """
    per_ue = np.zeros(NUM_UES)
    prb_count = np.zeros(NUM_UES)

    for prb in range(len(allocation)):
        ue = allocation[prb]
        if 0 <= ue < NUM_UES:
            cqi = int(np.clip(cqi_matrix[ue, prb], 0, 15))
            per_ue[ue] += CQI_TO_SE.get(cqi, 0.0)
            prb_count[ue] += 1

    total = per_ue.sum()
    allocated_prbs = prb_count.sum()
    avg = total / max(allocated_prbs, 1)

    return {
        "total_se": float(total),
        "per_ue_se": per_ue.tolist(),
        "avg_se_per_prb": float(avg),
        "allocated_prbs": int(allocated_prbs),
    }


def compute_radar_protection(
    allocation: np.ndarray,
    radar_active: np.ndarray,
) -> dict:
    """Measure how well the allocation avoids radar-active PRBs.

    Returns
    -------
    dict with radar_prbs_allocated, total_radar_prbs, protection_ratio
    """
    radar_mask = radar_active.astype(bool)
    total_radar = radar_mask.sum()
    allocated_on_radar = 0

    for prb in range(len(allocation)):
        if radar_mask[prb] and 0 <= allocation[prb] < NUM_UES:
            allocated_on_radar += 1

    protection = 1.0 - (allocated_on_radar / max(total_radar, 1))

    return {
        "radar_prbs_allocated": int(allocated_on_radar),
        "total_radar_prbs": int(total_radar),
        "protection_ratio": float(protection),
    }


def compute_evacuation_time(
    allocation_sequence: list[np.ndarray],
    radar_onset_step: int,
    radar_active: np.ndarray,
) -> int:
    """Steps from radar onset until no UEs allocated on radar PRBs.

    Parameters
    ----------
    allocation_sequence : list of (NUM_RBS,) arrays over time
    radar_onset_step : int – step when radar starts
    radar_active : (NUM_RBS,) binary

    Returns
    -------
    int – number of steps to fully evacuate, or -1 if never evacuated
    """
    radar_mask = radar_active.astype(bool)
    for t in range(radar_onset_step, len(allocation_sequence)):
        alloc = allocation_sequence[t]
        conflict = False
        for prb in range(len(alloc)):
            if radar_mask[prb] and 0 <= alloc[prb] < NUM_UES:
                conflict = True
                break
        if not conflict:
            return t - radar_onset_step
    return -1


def compute_throughput(
    allocation: np.ndarray,
    cqi_matrix: np.ndarray,
    scs_khz: float = 15.0,
    symbol_duration_ms: float = None,
) -> dict:
    """Estimate throughput from allocation and CQI.

    Uses CQI → SE mapping and PRB bandwidth (12 × SCS).
    """
    prb_bw_hz = 12 * scs_khz * 1e3  # 180 kHz for 15 kHz SCS
    if symbol_duration_ms is None:
        symbol_duration_ms = 1.0 / (scs_khz * 1e3) * 1e3  # ~0.0667 ms

    per_ue_tp = np.zeros(NUM_UES)
    for prb in range(len(allocation)):
        ue = allocation[prb]
        if 0 <= ue < NUM_UES:
            cqi = int(np.clip(cqi_matrix[ue, prb], 0, 15))
            se = CQI_TO_SE.get(cqi, 0.0)
            per_ue_tp[ue] += se * prb_bw_hz  # bits/s

    return {
        "total_throughput_bps": float(per_ue_tp.sum()),
        "per_ue_throughput_bps": per_ue_tp.tolist(),
    }


def compute_fairness(per_ue_metric: np.ndarray) -> float:
    """Jain's fairness index."""
    n = len(per_ue_metric)
    s = per_ue_metric.sum()
    sq = (per_ue_metric ** 2).sum()
    if sq == 0:
        return 0.0
    return float(s ** 2 / (n * sq))


def compute_bler_estimate(
    allocation: np.ndarray,
    radar_active: np.ndarray,
    base_bler: float = 0.01,
    radar_bler_boost: float = 0.3,
) -> dict:
    """Estimate BLER given allocation and radar activity.

    Simple model: BLER = base_bler for clean PRBs, base_bler + radar_bler_boost
    for radar-affected PRBs.
    """
    per_ue_bler = np.zeros(NUM_UES)
    per_ue_count = np.zeros(NUM_UES)

    for prb in range(len(allocation)):
        ue = allocation[prb]
        if 0 <= ue < NUM_UES:
            bler = base_bler
            if radar_active[prb] > 0.5:
                bler += radar_bler_boost
            per_ue_bler[ue] += bler
            per_ue_count[ue] += 1

    for ue in range(NUM_UES):
        if per_ue_count[ue] > 0:
            per_ue_bler[ue] /= per_ue_count[ue]

    return {
        "avg_bler": float(per_ue_bler.mean()),
        "per_ue_bler": per_ue_bler.tolist(),
    }


# ---------------------------------------------------------------------------
# Full evaluation pipeline
# ---------------------------------------------------------------------------

def evaluate_allocation(
    allocation: np.ndarray,
    cqi_matrix: np.ndarray,
    radar_active: np.ndarray,
) -> dict:
    """Run all metrics on a single allocation snapshot."""
    return {
        "spectrum_efficiency": compute_spectrum_efficiency(allocation, cqi_matrix),
        "radar_protection": compute_radar_protection(allocation, radar_active),
        "throughput": compute_throughput(allocation, cqi_matrix),
        "fairness": compute_fairness(
            np.array(compute_throughput(allocation, cqi_matrix)["per_ue_throughput_bps"])
        ),
        "bler": compute_bler_estimate(allocation, radar_active),
    }


def evaluate_model(
    model,
    test_data: list,
    device: str = "cpu",
    model_type: str = "gnn",
) -> dict:
    """Evaluate a trained model on test dataset.

    Parameters
    ----------
    model : nn.Module (GNN) or agent (DQN/Hybrid)
    test_data : list of PyG Data objects (for GNN) or state dicts (for DQN)
    device : str
    model_type : one of 'gnn', 'dqn', 'hybrid'

    Returns
    -------
    dict with aggregated metrics across all test samples
    """
    all_metrics = []

    if model_type == "gnn":
        model.eval()
        loader = PyGDataLoader(test_data, batch_size=1, shuffle=False)
        with torch.no_grad():
            for batch in loader:
                batch = batch.to(device)
                out = model(batch)
                pred = out.argmax(dim=1).cpu().numpy()

                # Reconstruct CQI and radar from node features
                x = batch.x.cpu().numpy()
                cqi = x[:, 1:5].T * 15  # denorm
                radar = x[:, 5]

                m = evaluate_allocation(pred, cqi, radar)
                all_metrics.append(m)

    elif model_type in ("dqn", "hybrid"):
        for snap in test_data:
            state = snap if isinstance(snap, np.ndarray) else snap["state"]
            action = model.select_action(state)
            cqi = state[: NUM_UES * NUM_RBS].reshape(NUM_UES, NUM_RBS)
            radar = state[NUM_UES * NUM_RBS: NUM_UES * NUM_RBS + NUM_RBS]
            m = evaluate_allocation(action, cqi, radar)
            all_metrics.append(m)

    return _aggregate_metrics(all_metrics)


def _aggregate_metrics(metrics_list: list[dict]) -> dict:
    """Aggregate metrics across samples."""
    if not metrics_list:
        return {}

    keys = {
        "spectrum_efficiency.total_se",
        "spectrum_efficiency.avg_se_per_prb",
        "radar_protection.protection_ratio",
        "throughput.total_throughput_bps",
        "fairness",
        "bler.avg_bler",
    }

    agg = {}
    for key in keys:
        vals = []
        for m in metrics_list:
            parts = key.split(".")
            v = m
            for p in parts:
                if isinstance(v, dict):
                    v = v.get(p, None)
                else:
                    v = None
                    break
            if v is not None:
                vals.append(float(v))

        if vals:
            agg[key] = {
                "mean": float(np.mean(vals)),
                "std": float(np.std(vals)),
                "min": float(np.min(vals)),
                "max": float(np.max(vals)),
            }

    return agg


# ---------------------------------------------------------------------------
# Baseline evaluators (no ML – direct scheduler simulation)
# ---------------------------------------------------------------------------

def evaluate_round_robin(cqi_matrix: np.ndarray, radar_active: np.ndarray) -> dict:
    """Simulate Round-Robin allocation and evaluate."""
    alloc = np.zeros(NUM_RBS, dtype=np.int64)
    for prb in range(NUM_RBS):
        alloc[prb] = prb % NUM_UES
    return evaluate_allocation(alloc, cqi_matrix, radar_active)


def evaluate_best_cqi(cqi_matrix: np.ndarray, radar_active: np.ndarray) -> dict:
    """Simulate Best-CQI allocation and evaluate."""
    alloc = np.argmax(cqi_matrix, axis=0).astype(np.int64)
    return evaluate_allocation(alloc, cqi_matrix, radar_active)


def evaluate_proportional_fair(
    cqi_matrix: np.ndarray,
    radar_active: np.ndarray,
    past_throughput: np.ndarray | None = None,
) -> dict:
    """Simulate Proportional Fair allocation and evaluate."""
    if past_throughput is None:
        past_throughput = np.ones(NUM_UES)
    pf_metric = cqi_matrix / past_throughput.reshape(-1, 1)
    alloc = np.argmax(pf_metric, axis=0).astype(np.int64)
    return evaluate_allocation(alloc, cqi_matrix, radar_active)


def evaluate_all_baselines(cqi_matrix: np.ndarray, radar_active: np.ndarray) -> dict:
    """Run all baseline evaluations."""
    return {
        "round_robin": evaluate_round_robin(cqi_matrix, radar_active),
        "best_cqi": evaluate_best_cqi(cqi_matrix, radar_active),
        "proportional_fair": evaluate_proportional_fair(cqi_matrix, radar_active),
    }
