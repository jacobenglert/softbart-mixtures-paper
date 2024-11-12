---
---
---

# softbart-mixtures-paper

<!-- badges: start -->

<!-- badges: end -->

This repository contains code used to reproduce the simulation results in the accompanying manuscript entitled "Modeling Joint Health Effects of Environmental Exposure Mixtures with Bayesian Additive Regression Trees". The data used in the application is not able to be shared, so only the simulation study results may be reproduced. However, the design of the simulation study is based on the observed data, and the code used to fit and summarize the models in the real data analysis is identical to the code used in the simulation study.

### Repository Contents:

-   `software/`

    -   `install.R`: R script for installing necessary packages from CRAN and GitHub.

    -   `run-mcmc.R`: R script containing the main MCMC algorithm for fitting models.

    -   `helper-functions.R`: R script containing assorted helper functions used throughout the repository.

-   `sim-study/`

    -   `HPC/`

        -   `run-sim.sh`: Shell script for running the simulation for a given set of parameters. Calls `R/run-sim.R`.

        -   `run-sim-study.sh`: Shell script for looping through parameter settings in the simulation study, calling `run-sim.sh` for each setting. Upon completion this calls `summarize-sims.sh` to compile results.

        -   `summarize-sims.sh`: Shell script for summarizing all simulation results by calling `R/summarize-sims.R`.

    -   `Params/`

        -   `params.rds`: Table of parameters for each setting in the simulation study.

        -   `sim-study-data.RData`: Information from the real data analysis used to inform the simulation.

    -   `R/`

        -   `set-params.R`: R script for specifying the BART settings to be explored in the simulation study. The output is stored in `Params/params.rds`.

        -   `run-sim.R`: R script for running a simulation with a given key and seed (see below for more details).

        -   `summarize-sims.R`: R script for summarizing simulation results (e.g., calculating Monte Carlo standard errors). The output is stored in `Results/sim-study-summary.rds`.

        -   `analyze-sim-study.R`: R script for analyzing the results of the simulation study and producing the figures and tables for the manuscript.

    -   `Results/`

        -   `sim-study-summary.rds`: File containing the summarized results for each setting in the simulation. This has been included in the repository so that the entire simulation need not be run again.

        -   `Figures/`: Directory containing simulation study figures used in the manuscript.

### Getting Started

The first thing to do is to ensure all package listed in the software/install.R file are installed. This includes several which are installable from CRAN, and 3 (`jrpg`, `pdpd`, and `SoftBart`) which may be installed from GitHub. `jrpg` is used to sample the latent Pólya-gamma random variables in the MCMC algorithm, and `pdpd` is used to estimate accumulated local effects throughout the posterior distribution. `SoftBart` is also available on CRAN, however if you would like to store model fits to make predictions from later, the GitHub fork will allow for that.

### Reproducing Simulation Study Results

Individuals models in the simulation study can take several hours to run, and so the analysis was originally performed using a SLURM scheduler on the High Performance Computing (HPC) cluster at the Rollins School of Public Health at Emory University. Recreating the entire simulation is possible by reconfiguring the shell scripts located in `sim-study/HPC/` to work within an a cluster environment. To go this route, you can simply run `HPC/run-sim-study.sh` while in the `sim-study` subfolder. This script loops through all parameter settings considered in the simulation and upon completion summarizes the simulations and saves the result to a single file located at `sim-study/Results/sim-study-summary.rds`.

Running everything will take a considerable amount of time depending on the resources available to you. To simply see how one iteration of the simulation works one might follow these steps:

1.  Open `sim-study/R/run-sim.R` in RStudio or your IDE of choice.

2.  Comment out and modify the following lines:

    ``` r
    # args <- commandArgs(trailingOnly = TRUE)
    key <- 8 #as.numeric(args[1])
    seed <- 1 #as.numeric(args[2])
    ```

3.  Proceed to run the program. This will run one simulate of the setting with 25 sparse and soft trees. To examine other settings, look up the key in the simulation parameters table located at `sim-study/Params/params.rds`.

Alternatively, the complete summarized simulation results are stored at `sim-study/Results/sim-study-summary.rds`. Running `sim-study/R/analyze-sim-study.R` reads in this file and creates the tables and figures included in the manuscript. The figures are saved to `sim-study/Results/Figures/`.
