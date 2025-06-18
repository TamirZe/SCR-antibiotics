SCR-antibiotics
================

# Estimands Numerical Example

This repository includes scripts for generating data according to the described data generating mechanism, calculating estimands, and running a numerical example, as described in the paper "Causal effects on non-terminal event time with application to antibiotic usage and future resistance"
by Zehavi, Obolski, Chowers and Nevo.

This repository includes three main folders, and an additional folder with simulation results files.



## Folder Structure

### `Data_generation/`
Contains R scripts used to generate datasets based on our simulation settings.

### `Estimand_calculations/`
Implements the core logic for calculating the presented estimands (SACE, FICE, and AICE). Includes a script that provides helper functions and procedures.

### `Main_run/`
Main scripts, for the main text and for the supplementary materials, separately, that run the full pipeline:
- Loads data from `Data_generation`
- Applies estimand computations from `Estimand_calculations`
- Produce summary results and figures shown in the paper.

## How to Run
1. Begin by modifying the path to reflect your local environment

2. Clone the repository:
   ```bash
   git clone https://github.com/your-username/Estimands_numerical_example.git
   cd Estimands_numerical_example
  
