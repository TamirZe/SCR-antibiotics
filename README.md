SCR-antibiotics
================

This repository includes scripts for generating data according to the described data generating mechanism, calculating estimands, and running a numerical example, as described in the paper "Causal effects on non-terminal event time with application to antibiotic usage and future resistance" by Zehavi, Obolski, Chowers and Nevo, [Causal effects on non-terminal event time with application to antibiotic usage and future resistance](https://arxiv.org/abs/2506.09624).

This repository includes two main folders, `Estimands Numerical Example` and `Simulaion studies`.

The folders listed below are subfolders of the first main folder `Estimands Numerical Example`.

(1) Data_generation:
-----------------

Contains R scripts named `SimDataWeibFrail`, used to generate datasets based on our simulation settings.

(2) Estimand_calculations:
----------------------------------

Implements the core logic for calculating the presented estimands (SACE, FICE, and AICE). Includes a script that provides helper functions and procedures.

Note that the files in this folder are used for the original version of the numerical example (V0).

(3) Main_run:
----------------------------------

For V0, in the folder `V0`, there are two main scripts (for the main text and the supplementary materials) 
that run the full pipeline:
- Generating data using the `SimDataWeibFrail` script, according to parameter values given in `GetScenarioParams_paper` script.
- Applying estimand computations from `CalcTrueCausalParams` script.
- Producing summary results and the figures shown in the paper.

For the current version V1, the script `numerical_example_estimands_and_bounds_V1` runs the full pipeline,
for both the main text and the supplementary materials (using `main_bool=TRUE` and `main_bool=FALSE`, respectively):
- Generating data using the `SimDataWeibFrail` script, according to parameter values given in `GetScenarioParams_paper` script.
- Applying estimand and bounds computations.
- Producing summary results and the figures shown in the paper.

(4) Example_results:
----------------------------------

A folder containing Rdata files with the results of the numerical example. This folder will be provided on demand by the first author.


The folders listed below are subfolders of the first main folder `Simulaion studies`.


(1) Main_run_simulations:
----------------------------------


(1) Scripts_estimation_for_simulations:
----------------------------------
This folder containד the main **R scripts used for estimation procedures in the simulation studies**.
This folder holds the **core estimation scripts** used after the datasets have already been generated (using the file `generate_datasets_for_sim` described below).

Inside this folder, the following files are shown:

- `EMcausalSC_GH_current_FICE_est_time`  
  An R script that appears to implement the **EM-based estimation procedure** for the causal model, likely including the current version of the **FICE estimator** and timing-related components.

- `CalcRMST_FICE_est`  
  An R script that performs te monte carlo simulations, using the estimates obtained form the EM algorithm, from `EMcausalSC_GH_current_FICE_est_time`.


(3) Summary_after_boot:
----------------------------------


## Note
Ensure the path is modified to reflect the configuration of your local system.

  
