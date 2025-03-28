# Finite-Key Analysis of Quantum Key Distribution with Characterized Devices Using Entropy Accumulation

This is a public version of the code used in *Finite-Key Analysis of Quantum Key Distribution with Characterized Devices Using Entropy Accumulation* \[[arXiv](https://arxiv.org/abs/2203.06554)\]. This repository pre-dates the open QKD security software.

This code repository contains code to generate all the key rate plots in the above paper. Those examples are used to demonstrate our finite-key security analysis method based on entropy accumulation. Our method works for entanglement-based discrete-variable quantum key distribution protocols with characterized devices.

> [!CAUTION]
> This code uses parallel computing to reduce run time, however CVX is known to be unstable when run in parallel. If you encounter instabilities when running the code, we recommend changing the `parfor` loops to `for` loops in the main function code.

## Installation instructions
> [!CAUTION]
> This repository is for archival and transparency purposes.

This code was designed and tested with Matlab 2023b, though any other recent edition should work. In addition to Matlab, this code requires the following additional resources:
 - [CVX](https://cvxr.com/cvx/download/) v2.2, a library for solving convex optimization problems in MATLAB.
 - The following Mathworks's toolboxes that can be installed with Matlab:
   - Parallel Computing
   - Optimization
   - Symbolic Math
 - For plotting, a Python 3 installation with the additional packages:
   - Scipy
   - Numpy
   - Matplotlib

1. Download the latest release on the side bar or clone with git
```
git clone https://github.com/Optical-Quantum-Communication-Theory/Finite-Key-Analysis-of-Quantum-Key-Distribution-with-Characterized-Devices
```
2. Unzip in your preferred directory and add the folders `EATSolver` and `subroutines` to your Matlab path.
3. Each figure is generated using a combination of Matlab scripts for data, and Python functions for plotting. Here is the list of all assications. Run the Matlab scripts first, then the appropriate Python plotting function.

   | Figure  | Python Plotting Script   | Curves and Matlab Scripts                                                                                                                    | Aditional Notes |
   | ------- | ------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------- | --------------- |
   | Fig. 2  | `plotFig2.py`            | <ul><li> Algorithm 1: `fig2curve1.m` </li><li> Algorithm 2: `qubitBB84v3QScript.m` </li></ul>                                                |                 |
   | Fig. 3  | `plotFig3.py`            | <ul><li> Theorem 5: `qubitBB84v3QScript.m` </li><li> Theorem 15: `qubitBB84v2QScript.m` </li></ul>                                           |                 |
   | Fig. 4a | `plotFig4a.py`           | All curves: `qubitBB84_explore_t_Script.m`                                                                                                   |                 |
   | Fig. 4b | `plotFig4b.py`           | All curves: `qubitBB84_explore_t_Script.m`                                                                                                   | Data is post-processed in `plotFig4b.py`. |
   | Fig. 5  | `plotFig5.py`            | <ul><li> Algorithm 1: `sixfouralg1Script.m` </li><li> Algorithm 2: `sixfourv3QScript.m` </li></ul>                                           |                 |
   | Fig. 6  | `plotFig6.py`            | <ul><li> Theorem 5: `sixfourv3QScript.m` </li><li> Theorem 15: `sixfourv2QScript.m` </li></ul>                                               |                 |
   | Fig. 7a | `plotFig7a.py`           | <ul><li> d=2, EAT: `hd2mubScriptd2v3.m` </li><li> d=3, EAT: `hd2mubScriptd3v3.m` </li></ul>                                                  | Data for d=2 and d=3 postselection curves are included in `ps.mat`. Asymptotic key rates are evaluated directly in `plotfig7a.py` using known analytical formula. |
   | Fig. 7b | `plotFig7b.py`           | <ul><li> d=5, EAT: `hd2mubScriptd5v3.m` </li><li> d=7, EAT: `hd2mubScriptd7v3.m` </li></ul>                                                  | Data for d=5 and d=7 postselection curves are included in `ps.mat`. Asymptotic key rates are evaluated directly in `plotfig7b.py` using known analytical formula. |
   | Fig. 8  | `plotFig8.py`            | <ul><li> L = 0 km: `opticalBB84L0Script.m` </li><li> L = 10 km: `opticalBB84L10Script.m` </li><li> L = 20 km: `opticalBB84L20Script.m` </li> |                 |
