SCR-antibiotics
================

This repository includes scripts for generating data according to the described data generating mechanism, calculating estimands, and running a numerical example, as described in the paper "Causal effects on non-terminal event time with application to antibiotic usage and future resistance" by Zehavi, Obolski, Chowers and Nevo, [Causal effects on non-terminal event time with application to antibiotic usage and future resistance](https://arxiv.org/abs/2506.09624).


This repository includes three main folders and an additional folder with simulation results files. All folders listed below are subfolders of the main folder "Estimands Numerical Example."

(1) Data_generation:
-----------------

Contains R scripts named `SimDataWeibFrail`, used to generate datasets based on our simulation settings.

(2) Estimand_calculations:
----------------------------------

Implements the core logic for calculating the presented estimands (SACE, FICE, and AICE). Includes a script that provides helper functions and procedures.

(3) Main_run:
----------------------------------

Two main scripts (for the main text and the supplementary materials) that run the full pipeline:
- Generating data using the `SimDataWeibFrail` script, according to parameter values given in `GetScenarioParams_paper` script.
- Applying estimand computations from `CalcTrueCausalParams` script.
- Producing summary results and the figures shown in the paper.

(4) Example_results:
----------------------------------

A folder containing Rdata files with the results of the numerical example. This folder will be provided on demand by the first author.



## Note
Ensure the path is modified to reflect the configuration of your local system.

  
