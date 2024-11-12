# Program Name: summarize-sims.R
# Description:  Summarize and compile simulation study results.

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(here)


# Identify files and keys -------------------------------------------------
all_files <- list.files(here('sim-study','Results','temp'), full.names = T)
all_keys <- unique(as.numeric(substr(basename(all_files), 1, 3)))


# Summarize results -------------------------------------------------------

# Helpful functions
prop_se <- function (x) sqrt(mean(x) * (1 - mean(x)) / length(x))
cont_se <- function (x) sd(x) / sqrt(length(x))

# Summary function
get_summmary <- function (key) {
  
  # Read key files
  key_idx <- which(as.numeric(substr(basename(all_files), 1, 3)) == key)
  key_results <- lapply(all_files[key_idx], readRDS)
  
  # Summarize bias, coverage, and RMSE with Monte Carlo standard error
  stats <- lapply(key_results, '[[', 'stats') |>
    bind_rows() |>
    group_by(param) |>
    reframe(
      est_mean = ifelse(param %in% c('G','nu'), NA, mean(est)),
      est_se = ifelse(param %in% c('G','nu'), NA, cont_se(est)),
      bias_mean = mean(bias),
      bias_se = cont_se(bias),
      coverage_mean = mean(coverage),
      coverage_se = prop_se(coverage),
      rmse_mean = ifelse(param %in% c('G','nu'), mean(rmse), NA),
      rmse_se = ifelse(param %in% c('G','nu'), cont_se(rmse), NA)
    ) |>
    distinct() |>
    ungroup() |>
    mutate(est_lower = est_mean - 1.96 * est_se, 
           est_upper = est_mean + 1.96 * est_se,
           bias_lower = bias_mean - 1.96 * bias_se, 
           bias_upper = bias_mean + 1.96 * bias_se,
           coverage_lower = coverage_mean - 1.96 * coverage_se, 
           coverage_upper = coverage_mean + 1.96 * coverage_se,
           rmse_lower = rmse_mean - 1.96 * rmse_se, 
           rmse_upper = rmse_mean + 1.96 * rmse_se)
  
  # Summarize first order ALE
  ale1 <- lapply(key_results, '[[', 'ale1') |>
    bind_rows() |>
    mutate(est_se = cont_se(est), .by = c(x, var)) |>
    summarise(across(c(est, lcl, ucl, truth, est_se), mean),
              .by = c(x, var))
  
  # Summarize second order ALE
  ale2 <- lapply(key_results, '[[', 'ale2') |>
    bind_rows() |>
    mutate(est_se = cont_se(est), .by = c(x1, x2, var1, var2)) |>
    summarise(across(c(est, lcl, ucl, truth, est_se), mean),
              .by = c(x1, x2, w1, w2, h1, h2, var1, var2)) |>
    arrange(var1, var2, x1, x2)
  
  # Summarize first + second order ALE
  ale3 <- lapply(key_results, '[[', 'ale3') |>
    bind_rows() |>
    mutate(est_se = cont_se(est), .by = c(x1, x2, var1, var2)) |>
    summarise(across(c(est, lcl, ucl, truth, est_se), mean),
              .by = c(x1, x2, w1, w2, h1, h2, var1, var2)) |>
    arrange(var1, var2, x1, x2)
  
  stats$key <- key
  ale1$key <- key
  ale2$key <- key
  ale3$key <- key
  
  return (list(stats = stats, ale1 = ale1, ale2 = ale2, ale3 = ale3))

}

# Summarize simulation results for each key
results <- lapply(all_keys, get_summmary)


# Export results ----------------------------------------------------------
saveRDS(results, here('sim-study','Results','sim-study-summary.rds'))
