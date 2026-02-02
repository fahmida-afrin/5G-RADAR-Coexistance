# 5G-RADAR-Coexistance
# 5G NR – Radar Coexistence with GNN-DQN Spectrum Management
**Symbol-based NR TDD scheduling with waveform-level radar interference injection + ML-based intelligent PRB allocation**

This repository contains a MATLAB-based 5G NR (TDD) network simulation framework built on the **Communications Toolbox Wireless Network Simulation Library (WNSL)** and **5G Toolbox**, extended to support **symbol-level scheduling** and **waveform-level pulsed radar interference injection**. It includes a **Python-based ML pipeline** with GNN, DQN, and Hybrid GNN-DQN models for intelligent spectrum management under radar interference.

## Project Overview
The goal of this framework is to study **5G NR mid-band coexistence with pulsed radar**, focusing on:
- Symbol-level effects of interference (e.g., which OFDM symbols are most vulnerable)
- Cross-layer consequences (PHY decoding → HARQ → goodput degradation)
- Benchmarking standard schedulers (RR / PF / Best-CQI) under interference
- **GNN-based spectrum management** that learns radar interference patterns
- **Hybrid GNN-DQN** reinforcement learning for adaptive PRB allocation
- Generating high-fidelity logs suitable for ML/GNN-based spectrum management research

---

## Repository Structure

```
5G-RADAR-Coexistance/
├── NRwithRADAR.m                    # Main MATLAB simulation script
├── hAddImpairmentWithRadar.m        # Radar waveform injection (channel model)
├── hNRScheduler.m                   # Base scheduler class
├── hNRSchedulerRoundRobin.m         # Round Robin scheduler
├── hNRSchedulerProportionalFair.m   # Proportional Fair scheduler
├── hNRSchedulerBestCQI.m            # Best-CQI scheduler
├── hNRGNBPhy.m / hNRUEPhy.m        # PHY layer (gNB and UE)
├── hNRGNBMAC.m / hNRUEMAC.m        # MAC layer
├── hWirelessNetworkSimulator.m      # Event-driven network simulator
├── ...                              # Other MATLAB support files
│
├── ml/                              # Python ML pipeline
│   ├── __init__.py
│   ├── data_pipeline.py             # MATLAB → Python data conversion & graph construction
│   ├── gnn_spectrum.py              # GCN, GAT, GraphSAGE models for PRB allocation
│   ├── dqn_spectrum.py              # DQN & Dueling DQN with PER
│   ├── hybrid_gnn_dqn.py           # Hybrid GNN-DQN architecture
│   ├── matlab_interface.py          # MATLAB-Python bridge (offline & online)
│   ├── evaluation.py                # Evaluation metrics (throughput, BLER, fairness, etc.)
│   ├── analysis.py                  # Statistical analysis (t-test, ANOVA, tables)
│   ├── train_gnn.py                 # GNN training script
│   ├── train_hybrid.py              # DQN/Hybrid training script
│   └── generate_figures.py          # Publication-quality figure generation
│
├── data/
│   ├── raw/                         # MATLAB simulation outputs
│   ├── processed/                   # Cleaned data for ML training
│   └── results/                     # Experiment results
│
├── paper/
│   ├── main.tex                     # IEEE-format paper draft
│   ├── references.bib               # Bibliography
│   └── figures/                     # Generated figures
│
└── requirements.txt                 # Python dependencies
```

---

## Prerequisites

### MATLAB Products
- **5G Toolbox**
- **Communications Toolbox**
- **Communications Toolbox Wireless Network Simulation Library (WNSL)** support package

### Python Dependencies
```bash
pip install -r requirements.txt
```

Required packages: PyTorch, PyTorch Geometric, NumPy, SciPy, Pandas, Matplotlib

---

## Quick Start

### 1. Run MATLAB Simulation
```matlab
% In MATLAB, run the main simulation:
NRwithRADAR
```

### 2. Generate Synthetic Training Data
```bash
python -m ml.data_pipeline synthetic --samples 5000 --output data/processed/synthetic.pt
```

### 3. Train GNN Model
```bash
python -m ml.train_gnn --model gcn --epochs 200 --lr 1e-3
python -m ml.train_gnn --model gat --epochs 200
python -m ml.train_gnn --model sage --epochs 200
```

### 4. Train Hybrid GNN-DQN
```bash
python -m ml.train_hybrid --agent hybrid --episodes 1000
python -m ml.train_hybrid --agent dqn --episodes 1000
```

### 5. Generate Figures
```bash
python -m ml.generate_figures --demo
```

### 6. Run Full Experiment Matrix
```bash
python -m ml.matlab_interface matrix        # View all experiment configs
python -m ml.matlab_interface gen-configs    # Generate MATLAB config files
```

---

## ML Models

| Model | Description | Architecture |
|-------|-------------|-------------|
| **GCN** | Graph Convolutional Network | 3-layer GCN, 64 hidden, BatchNorm |
| **GAT** | Graph Attention Network | 3-layer GAT, 4 heads, 64 hidden |
| **GraphSAGE** | Inductive graph learning | 3-layer SAGE, 64 hidden |
| **DQN** | Dueling Double DQN + PER | 256 hidden, per-PRB Q-heads |
| **Hybrid** | GNN feature extractor + DQN | GCN→128-d embedding + Dueling DQN |

---

## Simulation Parameters

| Parameter | Value |
|-----------|-------|
| Carrier Frequency | 3.595 GHz |
| Bandwidth | 5 MHz |
| Subcarrier Spacing | 15 kHz |
| PRBs | 25 |
| UEs | 4 (at 100, 150, 300, 400 m) |
| TDD Period | 5 ms |
| Radar PRF | 1 kHz |
| Radar Pulse Width | 70 μs |
| Radar BW | 5 MHz |
| Radar Freq Offset | 3 MHz |
| Radar Power | 1 kW @ 5 km |

---

## Evaluation Metrics

- **Throughput**: Total and per-UE (bits/s)
- **Spectral Efficiency**: bits/s/Hz per PRB
- **BLER**: Block error rate under radar interference
- **Radar Protection**: % of radar-active PRBs left unallocated
- **Evacuation Time**: Slots to clear radar band after onset
- **Fairness**: Jain's fairness index

---

## Notes
- The scripts call `wirelessnetworkSupportPackageCheck` to verify required MATLAB support packages
- Python ML pipeline can operate in offline mode (recommended) using pre-generated MATLAB data
- Online MATLAB Engine API interface available for real-time interaction
