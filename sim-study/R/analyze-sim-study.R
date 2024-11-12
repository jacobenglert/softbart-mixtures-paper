# Program Name: analyze-sim-study.R
# Description:  Create figures and tables for the simulation study.

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(patchwork)
library(here)
library(knitr)
library(kableExtra)


# Load Simulation Statistics and Parameter Settings -----------------------

results <- readRDS(here('sim-study','Results','sim-study-summary.rds'))
params <- readRDS(here('sim-study','Params','params.rds'))


# Plotting Function -------------------------------------------------------

stats_plot <- function (data) {
  
  data |>
    ggplot(aes(x = factor(num_tree), y = mean, ymin = lower, ymax = upper, 
               color = soft, lty = sparse, shape = sparse)) +
    geom_hline(aes(yintercept = ref), col = 'gray', lty = 2) +
    geom_pointrange(size = 0.3, position = position_dodge2(0.5)) +
    facet_wrap(~type, scales = 'free') +
    theme_bw() +
    theme(legend.position = 'bottom',
          text = element_text(size = 14)) +
    labs(x = 'Number of Trees',
         y = 'Simulation Mean (95% Monte Carlo Interval)',
         color = 'Soft Trees?',
         lty = 'Sparse Trees?',
         shape = 'Sparse Trees?')
  
}


# Basic Simulation Statistics ---------------------------------------------

# All simulation statistics
sim_stats <- lapply(results, '[[', 'stats') |>
  bind_rows() |>
  left_join(params, by = join_by('key' == 'key'))

# Prepare for plotting
sim_stats_transform <- sim_stats |>
  select(key, param, any_of(colnames(params)), starts_with(c('bias','coverage','rmse'))) |>
  pivot_longer(cols = starts_with(c('bias','coverage','rmse')),
               names_to = c(".value", "value_type"),
               names_pattern = "([a-z]+)_([a-z]+)") |>
  pivot_longer(cols = c('bias','coverage','rmse'), names_to = 'type') |>
  pivot_wider(id_cols = c(key, param, any_of(colnames(params)), type), 
              names_from = value_type, values_from = value) |>
  mutate(ref = case_when(type == 'bias' ~ 0, 
                         type == 'coverage' ~ 0.95))

# Plot bias and coverage for global parameters
sim_stats_transform |>
  filter(param %in% c('alpha','rho','tau2','xi') & type %in% c('bias','coverage')) |>
  stats_plot() +
  facet_grid(type ~ param, scales = 'free_y')

# ggsave(here('sim-study','Results','Figures','sim-global-stats.png'), height = 5, width = 7)

# Figure B1: Plot bias and coverage for confounders
sim_stats_transform |>
  filter(param %in% c('beta1','beta2','beta3','beta4') & type %in% c('bias','coverage')) |>
  mutate(param = gsub('beta','X', param)) |>
  stats_plot() +
  facet_grid(type ~ param, scales = 'free_y')

ggsave(here('sim-study','Results','Figures','sim-beta-stats.png'), height = 5, width = 7)

# Random effects accuracy
sim_stats_transform |>
  filter(param == 'nu') |>
  mutate(type = factor(type, levels = c('bias','coverage','rmse'),
                       labels = c('Average Bias','Average 95% CrI Coverage','RMSE'))) |>
  stats_plot()

# BART prediction accuracy
sim_stats_transform |>
  filter(param == 'G') |>
  mutate(type = factor(type, levels = c('bias','coverage','rmse'),
                       labels = c('Average Bias','Average 95% CrI Coverage','RMSE'))) |>
  stats_plot()
# ggsave(here('sim-study','Results','Figures','sim-bart-stats.png'), height = 5, width = 8)

# Table 1: Simulation statistics for BART predictions
bart_stats_tbl <- sim_stats_transform |>
  mutate(mean = ifelse(type %in% c('bias','rmse'), mean * 10, mean),
         se = ifelse(type %in% c('bias','rmse'), se * 10, se)) |>
  mutate(pt_mcse = paste0(format(round(mean, 2), nsmall = 2), ' (', 
                          format(round(se, 4), nsmall = 4), ')')) |>
  pivot_wider(id_cols = c(param, sparse, soft, num_tree), 
              names_from = type, values_from = pt_mcse) |>
  filter(param == 'G') |>
  arrange(num_tree, soft, sparse) |>
  select(num_tree, soft, sparse, bias, coverage, rmse)

bart_stats_kbl <- bart_stats_tbl |>
  mutate(soft = ifelse(soft, '\\checkmark', ''),
         sparse = ifelse(sparse, '\\checkmark', '')) |>
  kbl(format = 'latex', booktabs = TRUE, linesep = c('','','','\\addlinespace'),
      align = c('c'), label = 'sim-bart-stats', escape = FALSE,
      col.names = c('$T$\\footnotemark[1]',
                    'Soft\\footnotemark[2]', 
                    'Sparse\\footnotemark[3]',
                    'Bias $\\times 10$ (MCSE)',
                    'Coverage (MCSE)',
                    'RMSE $\\times 10$ (MCSE)'),
                    # 'Bias\\footnotemark[4] $\\times 10$ (MCSE)',
                    # 'Coverage\\footnotemark[5] (MCSE)',
                    # 'RMSE\\footnotemark[6] $\\times 10$ (MCSE)'),
      caption = 'Average bias, average 95\\% CrI coverage,
      and RMSE for BART predictions across simulations.') |>
  footnote(general_title = '', escape = FALSE,
           general = c('CrI: Bayesian posterior credible interval.',
                       'MCSE: Monte Carlo Standard Error.'),
           symbol = c('Number of trees.',
                      'Soft BART used (Linero and Yang, 2018).',
                      'Sparse branching process used (Linero, 2018).'))
                      # '$\\\\text{Bias} = \\\\frac{1}{N} \\\\sum_{i=1}^N \\\\left[ \\\\hat{f}(\\\\bz_i) - f(\\\\bz_i) \\\\right]$ for each simulation.',
                      # '$\\\\text{Coverage} = \\\\frac{1}{N} \\\\sum_{i=1}^N \\\\mathbb{I}\\\\left[ \\\\hat{f}_{0.025}(\\\\bz_i) \\\\le f(\\\\bz_i) \\\\le \\\\hat{f}_{0.975}(\\\\bz_i) \\\\right]$ for each simulation.',
                      # '$\\\\text{RMSE} = \\\\frac{1}{N} \\\\sum_{i=1}^N \\\\left[ \\\\hat{f}(\\\\bz_i) - f(\\\\bz_i) \\\\right]^2$ for each simulation.'))
  
bart_stats_kbl

# Table B1: simultion statistics for global parameters
all_stats_tbl <- sim_stats_transform |>
  # mutate(mean = ifelse(type %in% c('bias','rmse'), mean * 10, mean),
  #        se = ifelse(type %in% c('bias','rmse'), se * 10, se)) |>
  mutate(pt_mcse = paste0(format(round(mean, 3), nsmall = 3), ' (', 
                          format(round(se, 3), nsmall = 3), ')')) |>
  pivot_wider(id_cols = c(param, sparse, soft, num_tree), 
              names_from = type, values_from = pt_mcse) |>
  filter(!(param %in% c('G','nu','alpha'))) |>
  arrange(param, num_tree, soft, sparse) |>
  select(param, num_tree, soft, sparse, bias, coverage)
  
global_stats_kbl <- all_stats_tbl |>
  filter(param %in% c('rho','tau2','xi')) |>
  pivot_wider(id_cols = num_tree:sparse, names_from = param, values_from = bias:coverage) |>
  select(num_tree:sparse, ends_with('rho'), ends_with('tau2'), ends_with('xi')) |>
  mutate(soft = ifelse(soft, '\\checkmark', ''),
         sparse = ifelse(sparse, '\\checkmark', '')) |>
  kbl(format = 'latex', booktabs = TRUE, linesep = c('','','','\\addlinespace'),
      align = c('c'), label = 'sim-global-stats', escape = FALSE,
      col.names = c('$T$\\footnotemark[1]',
                    'Soft\\footnotemark[2]', 
                    'Sparse\\footnotemark[3]',
                    rep(c('Bias (MCSE)','Coverage (MCSE)'), times = 3)),
      caption = 'Bias and 95\\% CrI coverage for global parameters in simulation study.') |>
  add_header_above(c(' ' = 3, '$rho$' = 2, '$tau2$' = 2, '$xi$' = 2)) |>
  footnote(general_title = '', escape = FALSE,
           general = c('CrI: Bayesian posterior credible interval.',
                       'MCSE: Monte Carlo Standard Error.'),
           symbol = c('Number of trees.',
                      'Soft BART used (Linero and Yang, 2018).',
                      'Sparse branching process used (Linero and Yang, 2018).'))

global_stats_kbl


# Summarize BART with ALE -------------------------------------------------

# For ALE, we will select a single setting

key <- 8

# Figure 1: First order ALE
results[[key]]$ale1 |>
  filter(!(x %in% c(min(x), max(x))), .by = var) |>
  ggplot(aes(x = x, ymin = lcl, ymax = ucl)) +
  # ggplot(aes(x = x, ymin = est - 1.96 * est_se, ymax = est + 1.96 * est_se)) +
  geom_hline(yintercept = 0, lty = 2, color = 'gray') +
  geom_ribbon(alpha = 0.3) +
  geom_line(aes(y = est, col = 'Estimate')) +
  geom_line(aes(y = truth, col = 'Truth'), lty = 2) +
  facet_wrap(~factor(var, levels = unique(var)), ncol = 5) +
  theme_bw() +
  theme(legend.position = 'bottom',
        text = element_text(size = 11)) +
  scale_color_manual(values = c("Truth" = "red", "Estimate" = "black")) +
  scale_x_continuous(limits = c(0.15, 0.85), breaks = c(0.2, 0.5, 0.8)) +
  labs(#title = 'ALE Plots of First-Order Main Effects',
    x = 'Exposure Value',
    y = 'ALE Main Effect (95% Credible Interval)',
    color = '') +
  guides(color = guide_legend(override.aes = list(lty = c(1, 2))))

fname <- paste0('sim-ale-main-effects-key-', key, '.png')
ggsave(here('sim-study','Results','Figures', fname), height = 3.5, width = 5.5)

# Second order ALE plotting template
gg_second_order_ale <- function (data) {
  ggplot(data = data, mapping = aes(xmin = x1 - w1, xmax = x1 + w2, ymin = x2 - h1, ymax = x2 + h2, fill = est)) +
    geom_rect() +
    theme_bw() +
    coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
    scale_x_continuous(expand = c(0,0), breaks = c(0.25, 0.50, 0.75)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0.25, 0.50, 0.75)) +
    facet_grid(var2 ~ var1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 270),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 12)) +
    labs(x = 'Exposure 1 Value',
         y = 'Exposure 2 Value',
         fill = 'ALE')
}

# Figure B2: Second-order ALE only
ale2 <- results[[key]]$ale2 |>
  filter(!(x1 %in% c(min(x1), max(x1)) | x2 %in% c(min(x2), max(x2))), .by = c(var1, var2))
limits <- c(min(ale2$truth, ale2$est), max(ale2$truth, ale2$est))

est2 <- ale2 |>
  gg_second_order_ale() +
  scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
  labs(title = '(a) Posterior Mean') +
  coord_fixed()

true2 <- ale2 |>
  gg_second_order_ale() +
  geom_rect(aes(fill = truth)) +
  scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
  labs(title = '(b) Truth') +
  coord_fixed()

ale2plt <- est2 + true2 + plot_layout(guides = 'collect')

fname <- paste0('sim-ale-interaction-only-key-', key, '.png')
ggsave(filename = here('sim-study','Results','Figures', fname),
       plot = ale2plt, height = 4, width = 7)


# Figure B3: First + Second order ALE
ale3 <- results[[key]]$ale3 |>
  filter(!(x1 %in% c(min(x1), max(x1)) | x2 %in% c(min(x2), max(x2))), .by = c(var1, var2))
limits <- c(min(ale3$truth, ale3$est), max(ale3$truth, ale3$est))

est3 <- ale3 |>
  gg_second_order_ale() +
  scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
  labs(title = '(a) Posterior Mean') +
  coord_fixed()

true3 <- ale3 |>
  gg_second_order_ale() +
  geom_rect(aes(fill = truth)) +
  scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
  labs(title = '(b) Truth') +
  coord_fixed()

ale3plt <- est3 + true3 + plot_layout(guides = 'collect')

fname <- paste0('sim-ale-main-effects-and-interaction-key-', key, '.png')
ggsave(filename = here('sim-study','Results','Figures', fname), 
       plot = ale3plt, height = 4, width = 7)
