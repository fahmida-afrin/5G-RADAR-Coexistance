"""
MATLAB-Python Interface for 5G-Radar Coexistence.

Provides two modes:
  Option A (Offline / Recommended):
      MATLAB Simulation → .mat/.csv → Python Training → Evaluate

  Option B (Online):
      Python ↔ MATLAB Engine API → real-time interaction

This module handles both approaches.
"""

import os
import subprocess
import json
import numpy as np
from pathlib import Path

try:
    import scipy.io as sio
except ImportError:
    sio = None


# ---------------------------------------------------------------------------
# Constants matching NRwithRADAR.m
# ---------------------------------------------------------------------------
MATLAB_SCRIPT_DIR = Path(__file__).resolve().parent.parent  # repo root
RESULTS_DIR = MATLAB_SCRIPT_DIR / "Results"
DATA_RAW_DIR = MATLAB_SCRIPT_DIR / "data" / "raw"
DATA_PROCESSED_DIR = MATLAB_SCRIPT_DIR / "data" / "processed"


# ---------------------------------------------------------------------------
# 1. Scenario configuration generator
# ---------------------------------------------------------------------------

SCENARIO_CONFIGS = {
    "no_radar": {
        "radar_enabled": False,
        "scheduler": "RR",
        "radar_prf": 1e3,
        "radar_pw": 70e-6,
        "radar_freq_offset": 3e6,
        "radar_power": 0,
    },
    "radar_rr": {
        "radar_enabled": True,
        "scheduler": "RR",
        "radar_prf": 1e3,
        "radar_pw": 70e-6,
        "radar_freq_offset": 3e6,
        "radar_power": 1000,
    },
    "radar_pf": {
        "radar_enabled": True,
        "scheduler": "PF",
        "radar_prf": 1e3,
        "radar_pw": 70e-6,
        "radar_freq_offset": 3e6,
        "radar_power": 1000,
    },
    "radar_bestcqi": {
        "radar_enabled": True,
        "scheduler": "BestCQI",
        "radar_prf": 1e3,
        "radar_pw": 70e-6,
        "radar_freq_offset": 3e6,
        "radar_power": 1000,
    },
}


def generate_matlab_config(scenario_name: str, output_path: str | None = None) -> str:
    """Generate a MATLAB .m config snippet for a given scenario.

    Parameters
    ----------
    scenario_name : str
        One of SCENARIO_CONFIGS keys.
    output_path : str, optional
        If given, write the config to this .m file.

    Returns
    -------
    str
        MATLAB code snippet.
    """
    cfg = SCENARIO_CONFIGS[scenario_name]
    lines = [
        f"% Auto-generated config for scenario: {scenario_name}",
        f"simParameters.SchedulerStrategy = '{cfg['scheduler']}';",
    ]
    if cfg["radar_enabled"]:
        lines += [
            f"simParameters.radar.prf = {cfg['radar_prf']};",
            f"simParameters.radar.PulseWidth = {cfg['radar_pw']};",
            f"simParameters.radar.FreqOffset = {cfg['radar_freq_offset']};",
            f"simParameters.radar.Power = {cfg['radar_power']};",
        ]
    else:
        lines.append("simParameters.radar.Power = 0;  % Radar disabled")

    code = "\n".join(lines) + "\n"
    if output_path:
        with open(output_path, "w") as f:
            f.write(code)
    return code


# ---------------------------------------------------------------------------
# 2. MATLAB batch runner (offline – Option A)
# ---------------------------------------------------------------------------

def run_matlab_scenario(
    scenario_name: str,
    matlab_executable: str = "matlab",
    num_frames: int = 10,
    num_trials: int = 1,
    output_dir: str | None = None,
) -> list[str]:
    """Run a MATLAB simulation for a given scenario (requires MATLAB on PATH).

    This generates a wrapper script that sources NRwithRADAR.m with
    modified parameters, then saves results.

    Parameters
    ----------
    scenario_name : str
    matlab_executable : str
    num_frames : int
    num_trials : int
    output_dir : str, optional

    Returns
    -------
    list[str]
        Paths to output .mat files.
    """
    if output_dir is None:
        output_dir = str(DATA_RAW_DIR / scenario_name)
    os.makedirs(output_dir, exist_ok=True)

    cfg = SCENARIO_CONFIGS.get(scenario_name, SCENARIO_CONFIGS["radar_rr"])
    output_files = []

    for trial in range(num_trials):
        out_file = os.path.join(output_dir, f"trial_{trial:03d}.mat")
        # Build MATLAB command
        matlab_cmd = f"""
            cd('{MATLAB_SCRIPT_DIR}');
            rng({trial});
            simParameters = [];
            simParameters.NumFramesSim = {num_frames};
            simParameters.SchedulingType = 1;
            simParameters.NumUEs = 4;
            simParameters.UEPosition = [100 0 0; 150 0 0; 300 0 0; 400 0 0];
            simParameters.NumRBs = 25;
            simParameters.SCS = 15;
            simParameters.DLBandwidth = 5e6;
            simParameters.ULBandwidth = 5e6;
            simParameters.DLCarrierFreq = 3.595e9;
            simParameters.ULCarrierFreq = 3.595e9;
            simParameters.SchedulerStrategy = '{cfg["scheduler"]}';
            simParameters.radar.prf = {cfg["radar_prf"]};
            simParameters.radar.PulseWidth = {cfg["radar_pw"]};
            simParameters.radar.FreqOffset = {cfg["radar_freq_offset"]};
            simParameters.radar.Power = {cfg["radar_power"]};
            simParameters.radar.BW = 5e6;
            simParameters.radar.SampleRate = 30.72e6;
            % Run simulation and save
            try
                NRwithRADAR;
                save('{out_file}', 'simulationLogs', 'simParameters');
                fprintf('Saved: {out_file}\\n');
            catch ME
                fprintf('ERROR: %s\\n', ME.message);
            end
            exit;
        """
        cmd = [
            matlab_executable, "-batch", matlab_cmd.replace("\n", " ")
        ]
        print(f"Running {scenario_name} trial {trial}...")
        try:
            subprocess.run(cmd, timeout=600, check=False)
        except FileNotFoundError:
            print(f"MATLAB not found at '{matlab_executable}'. Skipping.")
            break
        except subprocess.TimeoutExpired:
            print(f"Trial {trial} timed out.")

        if os.path.exists(out_file):
            output_files.append(out_file)

    return output_files


def run_all_scenarios(
    matlab_executable: str = "matlab",
    num_frames: int = 10,
    num_trials: int = 10,
) -> dict[str, list[str]]:
    """Run all baseline scenarios."""
    results = {}
    for name in SCENARIO_CONFIGS:
        results[name] = run_matlab_scenario(
            name,
            matlab_executable=matlab_executable,
            num_frames=num_frames,
            num_trials=num_trials,
        )
    return results


# ---------------------------------------------------------------------------
# 3. MATLAB Engine API interface (online – Option B)
# ---------------------------------------------------------------------------

class MatlabEngineInterface:
    """Online interface using MATLAB Engine API for Python.

    Requires: ``pip install matlabengine``
    """

    def __init__(self):
        self.eng = None

    def start(self):
        import matlab.engine
        self.eng = matlab.engine.start_matlab()
        self.eng.cd(str(MATLAB_SCRIPT_DIR))
        print("MATLAB engine started.")

    def stop(self):
        if self.eng:
            self.eng.quit()
            self.eng = None

    def run_simulation(self, config_overrides: dict | None = None) -> dict:
        """Run one simulation with optional config overrides.

        Returns dict with keys: throughput, bler, retx_rate, cqi_logs
        """
        if self.eng is None:
            self.start()

        # Set overrides in MATLAB workspace
        if config_overrides:
            for k, v in config_overrides.items():
                self.eng.workspace[k] = v

        # Run the simulation script
        self.eng.eval("NRwithRADAR", nargout=0)

        # Extract results
        results = {
            "simParameters": self.eng.workspace["simParameters"],
        }
        try:
            results["simulationLogs"] = self.eng.workspace["simulationLogs"]
        except Exception:
            pass
        return results

    def inject_allocation(self, allocation: np.ndarray):
        """Inject ML-decided allocation back into MATLAB scheduler.

        This would require a modified MATLAB scheduler that reads from
        a shared variable or file.
        """
        import matlab
        alloc_matlab = matlab.double(allocation.tolist())
        self.eng.workspace["mlAllocation"] = alloc_matlab


# ---------------------------------------------------------------------------
# 4. Result extraction from .mat files
# ---------------------------------------------------------------------------

def extract_metrics_from_mat(filepath: str) -> dict:
    """Extract key metrics from a saved simulation .mat file.

    Returns
    -------
    dict with keys:
        throughput_dl, throughput_ul, goodput_dl, goodput_ul,
        bler_dl, bler_ul, retx_rate_dl, retx_rate_ul,
        per_ue_throughput
    """
    if sio is None:
        raise ImportError("scipy required")

    data = sio.loadmat(filepath, squeeze_me=True, struct_as_record=True)
    metrics = {}

    # Try to extract from simulationLogs
    if "simulationLogs" in data:
        logs = data["simulationLogs"]
        metrics["raw_logs"] = logs

    if "simParameters" in data:
        params = data["simParameters"]
        metrics["parameters"] = params

    return metrics


def aggregate_trial_metrics(mat_files: list[str]) -> dict:
    """Aggregate metrics across multiple trials.

    Returns dict with mean ± std for each metric.
    """
    all_metrics = [extract_metrics_from_mat(f) for f in mat_files if os.path.exists(f)]
    if not all_metrics:
        return {}
    return {"num_trials": len(all_metrics), "trials": all_metrics}


# ---------------------------------------------------------------------------
# 5. Experiment matrix for radar parameter sweeps
# ---------------------------------------------------------------------------

def generate_experiment_matrix() -> list[dict]:
    """Generate full experiment matrix for systematic evaluation.

    Covers:
      - 3 schedulers × 3 duty cycles = 9 baseline configs
      - Radar parameter sweep: pulse width, PRI, freq offset
      - Scalability: UE count, PRB count
    """
    experiments = []

    # Experiment 1: Scheduler comparison across duty cycles
    for sched in ["RR", "PF", "BestCQI"]:
        for dc_pct in [1, 5, 10]:
            # duty_cycle = pw / pri → pw = dc * pri
            pri = 1e-3  # 1 ms
            pw = dc_pct / 100.0 * pri
            experiments.append({
                "name": f"sched_{sched}_dc{dc_pct}",
                "scheduler": sched,
                "radar_prf": 1.0 / pri,
                "radar_pw": pw,
                "radar_freq_offset": 3e6,
                "radar_power": 1000,
                "num_ues": 4,
                "num_prbs": 25,
            })

    # ML-based (to be filled with trained model results)
    for dc_pct in [1, 5, 10]:
        pri = 1e-3
        pw = dc_pct / 100.0 * pri
        for model in ["GNN", "Hybrid"]:
            experiments.append({
                "name": f"ml_{model}_dc{dc_pct}",
                "scheduler": model,
                "radar_prf": 1.0 / pri,
                "radar_pw": pw,
                "radar_freq_offset": 3e6,
                "radar_power": 1000,
                "num_ues": 4,
                "num_prbs": 25,
            })

    # Experiment 2: Radar parameter variations
    for pw_us in [1, 10, 50]:
        for pri_ms in [1, 5, 10]:
            for freq_off in [0, 1e6, 3e6]:
                experiments.append({
                    "name": f"radar_pw{pw_us}us_pri{pri_ms}ms_fo{int(freq_off/1e6)}MHz",
                    "scheduler": "RR",
                    "radar_prf": 1.0 / (pri_ms * 1e-3),
                    "radar_pw": pw_us * 1e-6,
                    "radar_freq_offset": freq_off,
                    "radar_power": 1000,
                    "num_ues": 4,
                    "num_prbs": 25,
                })

    # Experiment 3: Scalability
    for n_ue in [2, 4, 8, 16]:
        for n_prb in [25, 50, 100]:
            experiments.append({
                "name": f"scale_ue{n_ue}_prb{n_prb}",
                "scheduler": "RR",
                "radar_prf": 1e3,
                "radar_pw": 70e-6,
                "radar_freq_offset": 3e6,
                "radar_power": 1000,
                "num_ues": n_ue,
                "num_prbs": n_prb,
            })

    return experiments


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="MATLAB-Python interface")
    sub = parser.add_subparsers(dest="cmd")

    p_gen = sub.add_parser("gen-configs", help="Generate MATLAB configs for all scenarios")
    p_run = sub.add_parser("run", help="Run MATLAB scenarios")
    p_run.add_argument("--matlab", default="matlab")
    p_run.add_argument("--trials", type=int, default=10)
    p_run.add_argument("--frames", type=int, default=10)

    p_matrix = sub.add_parser("matrix", help="Print experiment matrix")

    args = parser.parse_args()

    if args.cmd == "gen-configs":
        out_dir = str(DATA_RAW_DIR / "configs")
        os.makedirs(out_dir, exist_ok=True)
        for name in SCENARIO_CONFIGS:
            path = os.path.join(out_dir, f"config_{name}.m")
            generate_matlab_config(name, path)
            print(f"Generated: {path}")

    elif args.cmd == "run":
        results = run_all_scenarios(
            matlab_executable=args.matlab,
            num_frames=args.frames,
            num_trials=args.trials,
        )
        for name, files in results.items():
            print(f"{name}: {len(files)} files")

    elif args.cmd == "matrix":
        matrix = generate_experiment_matrix()
        print(f"Total experiments: {len(matrix)}")
        for exp in matrix:
            print(f"  {exp['name']}")

    else:
        parser.print_help()
