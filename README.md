# Calculating mutual information rates in a multi-state receptor system
This repository contains MATLAB code relating to Fig.2b from the paper "State- versus Reaction-Based Information Processing in Biochemical Networks" by Moor et al. (https://doi.org/10.1103/nfk6-8x5s). The code was used with MATLAB version 25.1.0.2943329 (R2025a) (The MathWorks Inc., Natick, Massachusetts).

## Content
The repository contains the following files:
1. `RunMSReceptorSimulations.m`: Main script to run and plot analyses.  
2. `GenerateReceptorModel.m`: Defines the stoichiometry and rates of the multi-state receptor model.
3. `GenerateLNAODEFilter.m`: Automatically generates differential equations for the required conditional covariances based on the Linear Noise Approximation.
4. `simulation_ms_receptor.mat`: Data file storing the calculated results.

## Reproducing the figures of the paper 
To reproduce Fig.2b from the paper, please run RunMSReceptorSimulations.m and refer to documentation inside the script. For plotting the results only, set `runSimulations=0`. For recomputing the resuls, set `runSimulations=1`, which will overwrite the existing file `simulation_ms_receptor.mat`. Note that the computation can take several minutes to complete.
