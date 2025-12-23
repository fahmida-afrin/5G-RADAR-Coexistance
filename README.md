# 5G-RADAR-Coexistance
# 5G NR – Radar Coexistence (MATLAB)  
**Symbol-based NR TDD scheduling with waveform-level radar interference injection**

This repository contains a MATLAB-based 5G NR (TDD) network simulation framework built on the **Communications Toolbox Wireless Network Simulation Library (WNSL)** and **5G Toolbox**, extended to support **symbol-level scheduling** and **waveform-level pulsed radar interference injection**. The codebase is structured to enable reproducible experiments and detailed logging of **throughput, goodput, HARQ retransmissions, and scheduling behavior** under radar interference.


## Project Overview
The goal of this simulation framework is to study **5G NR mid-band coexistence with pulsed radar**, focusing on:
- Symbol-level effects of interference (e.g., which OFDM symbols are most vulnerable)
- Cross-layer consequences (PHY decoding → HARQ → goodput degradation)
- Benchmarking standard schedulers (RR / PF / Best-CQI) under interference
- Generating high-fidelity logs suitable for later ML/GNN-based spectrum management work

---

## Prerequisites
### Required MATLAB Products
- **5G Toolbox**
- **Communications Toolbox**
- **Communications Toolbox Wireless Network Simulation Library (WNSL)** support package

### Notes
- The scripts call `wirelessnetworkSupportPackageCheck` to verify required support packages are installed.

---


