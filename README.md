# Sign-Indefinite Lyapunov

## Overview
This repository contains the code developed during my internship project on the stability analysis and control synthesis of **discrete-time systems with input saturation**. The implementation follows the reference framework in [1] and is further detailed in my internship report (see “Notes” below).

The project develops a discrete-time formulation for systems with input saturation, using **sign-indefinite quadratic Lyapunov functions**. By combining this approach with the **S-procedure** and **Linear Matrix Inequalities (LMIs)**, the framework provides conditions to analyze both **global** and **regional stability**, and to design stabilizing controllers in each case. This formulation reduces conservatism compared to standard quadratic Lyapunov methods, while ensuring asymptotic stability.

---

## Repository Structure

The code is organized into four main cases (each case has the same internal structure):

1. **Global Analysis**  
   - `Analysis_Global/` : MATLAB (LMILAB) solver implementation  
   - `YALMIP_Analysis_Global/` : YALMIP-based solver  
   - `global_analysis_discret.m` : example script demonstrating usage  

2. **Regional Analysis**  
   - `Analysis_Regional/`  
   - `YALMIP_Analysis_Regional/`  
   - `regional_analysis_discret.m`  

3. **Global Synthesis**  
   - `Synthesis_Global/`  
   - `YALMIP_Synthesis_Global/`  
   - `global_synthesis_discret.m`  

4. **Regional Synthesis**  
   - `Synthesis_Regional/`  
   - `YALMIP_Synthesis_Regional/`  
   - `regional_synthesis_discret.m`  

---

## Auxiliary functions

The folder `Auxiliary_function/` contains helper scripts used across all analysis and synthesis cases:

- `Discretize.m` – converts continuous-time models to discrete-time representation.  
- `ItssChurCohn.m` – checks the Schur–Cohn criterion for discrete-time stability.  
- `PlotGlobal.m` – plots results of the global stability analysis.  
- `PlotPhasePlan.m` – plots the phase plane for a given system.  
- `PlotPhasePlanGlobal.m` – plots the phase plane specifically for the global case.  
- `PlotRegional.m` – plots results of the regional stability analysis.  
- `Saturation.m` – defines the symmetric saturation and dead-zone functions.  

---

## Requirements
- MATLAB (with Control System / Robust Control toolboxes for LMIs via LMILAB)  
- OR: MATLAB + **YALMIP** (add to path) and an SDP solver compatible with YALMIP (e.g., **SeDuMi**, **SDPT3**, **MOSEK**)  

---

## Quick Start
1. Clone the repo and open MATLAB in the repository root.  
2. Add folders to the path (example):  
   ```matlab
   addpath(genpath('Auxiliary_function'));


[1] Isabelle Queinnec et al. "Design of Saturating State Feedback With Sign-Indefinite Quadratic Forms." *IEEE Transactions on Automatic Control*, vol. 67, no. 7, 2022, pp. 3507–3520. doi: [10.1109/TAC.2021.3106878](https://doi.org/10.1109/TAC.2021.3106878)
