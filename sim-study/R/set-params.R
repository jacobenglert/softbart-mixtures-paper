# Program Name: set-params.R
# Description:  Set parameters for simulation study

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)


# Specify Parameters ------------------------------------------------------

# Model Parameters
nT        <- 300
num_tree  <- c(10, 25, 50, 100)
k         <- 2
base      <- 0.95
power     <- 2
sparse    <- c(TRUE, FALSE)
soft      <- c(TRUE, FALSE)

# Compile
params <- crossing(nT, num_tree, k, base, power, sparse, soft) |>
  mutate(key = row_number()) |>
  select(key, everything())


# Export
saveRDS(params, here('sim-study','Params','params.rds'))
