# Program Name: run-sim.R
# Purpose:      Simulate overdispersed counts from a complicated mean function
#               and attempt to recover using proposed model.

# Load Packages -----------------------------------------------------------
library(sf)
library(pdpd)
library(mvtnorm)
library(spdep)
library(here)
source(here('software','run-mcmc.R'))


# Get Parameters ----------------------------------------------------------

# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)
key <- as.numeric(args[1])
seed <- as.numeric(args[2])
id <- paste0(sprintf("%03d", key), '-', sprintf("%03d", seed))

# Read parameters from file
params <- readRDS(here('sim-study','Params','params.rds'))
for (param in colnames(params)) {
  assign(param, unlist(params[key, param, drop = TRUE]))
}

# Number of quantile intervals to use for ALE
K <- 40


# Import Reference Data ---------------------------------------------------

load(here('sim-study','Params','sim-study-data.RData'))


# Simulate Data -----------------------------------------------------------

set.seed(seed)
nT  <- nT                 # Length of time-series
nS  <- length(population) # Number of locations
n   <- nT*nS              # Total number of observations

# Confounders
nx <- 4
x <- t(replicate(n, runif(nx, 0, 1)))
colnames(x) <- paste0('X', 1:nx)

# Exposures
nz <- 10
nza <- 5
nzb <- 5
zcov <- rbind(cbind(R, matrix(0, nzb, nzb)),
              cbind(matrix(0, nza, nza), diag(nzb)))

set.seed(1) # Use the same seed to ensure exposure quantiles do not change
z <- rmvnorm(n, sigma = zcov)
z <- apply(z, 2, \(x) (x - min(x)) / (max(x) - min(x))) #+ 1e-4
colnames(z) <- paste0('Z', 1:nz)

# Fixed parameters
set.seed(seed)  # reset seed
alpha <- -10    # fixed intercept
beta  <- matrix(c(-2, -1, 1, 2), ncol = 1) # fixed effects

# Spatial random effects
# Spatial weight matrices
W <- poly2nb(geometry) |> nb2mat(style = 'B')
D <- diag(rowSums(W))

tau2  <- 0.3    # marginal variance
rho   <- 0.9    # correlation
nu    <- rmvnorm(1, sigma = tau2 * solve(D - rho * W))[1,] # effects
x_nu  <- kronecker(matrix(1, nrow = nT, ncol = 1), diag(nS)) # design matrix

# BART component (Friedman 2001)
f <- function (x) {
  (10*sin(pi*x[,1]*x[,2]) + 20*(x[,3] - .5)^2 + 10*x[,4] + 5*x[,5]) / 5
}

population <- rep(population, times = nT) # Repeat population vector
geometry <- rep(geometry, times = nT)   # Repeat geometry vector
offset <- log(population)                      # Create offset

# Linear predictor
eta <- offset + alpha + x %*% beta + f(z) + x_nu %*% nu

# Dispersion parameter
xi <- 1

# Simulate counts
y <- rnbinom(n, size = xi, prob = 1 / (1 + exp(eta)))


# Fit Model ---------------------------------------------------------------

# Space indicators
s <- rep(1:nS, times = nT)

set.seed(seed)
fit <- run_mcmc(x = x, y = y, z = z, offset = offset,
                s = s, geometry = geometry,
                num_tree = num_tree, k = k, base = base, power = power,
                num_burn = 5000, num_thin = 5, num_save = 1000,
                sparse = sparse, soft = soft)


# Compute ALE -------------------------------------------------------------

# Define prediction function
pred_fun <- function (newdata) {
  sweep(fit$bart$predict_all(as.matrix(newdata)), 2, log(fit$xi), '+')
}

# Compute first-order ALE for each predictor
z_vars <- colnames(z)
ale1 <- lapply(z_vars, \(j) bayes_ale(z, pred_fun, vars = j, k = K, f = f))

ale1 <- do.call(rbind, ale1)
ale1$id <- id

# Compute second-order ALEs for first 5 predictors
pairs <- combn(z_vars[1:5], 2, simplify = FALSE)
ale2 <- lapply(pairs, \(j) bayes_ale(z, pred_fun, vars = j, k = K, f = f, 
                                     include_main_effects = FALSE))

ale2 <- do.call(rbind, ale2)
ale2$id <- id

# Compute second-order + first-order ALEs for first 5 predictors
ale3 <- lapply(pairs, \(j) bayes_ale(z, pred_fun, vars = j, k = K, f = f, 
                                     include_main_effects = TRUE))

ale3 <- do.call(rbind, ale3)
ale3$id <- id


# Compute Simulation Statistics -------------------------------------------

# Implicit intercept
alpha_true      <- alpha + mean(f(z)) + mean(nu) + log(xi)
alpha_est       <- median(fit$logmean)
alpha_bias      <- alpha_est - alpha_true
alpha_lower     <- quantile(fit$alpha, 0.025)
alpha_upper     <- quantile(fit$alpha, 0.975)
alpha_coverage  <- as.numeric(alpha_true >= alpha_lower & alpha_true <= alpha_upper)

# Fixed effects
beta_est    <- apply(fit$beta, 2, \(x) quantile(x, 0.500))
beta_bias   <- beta_est - beta
beta_lower  <- apply(fit$beta, 2, \(x) quantile(x, 0.025))
beta_upper  <- apply(fit$beta, 2, \(x) quantile(x, 0.975))
beta_coverage <- as.numeric(beta >= beta_lower & beta <= beta_upper)

# Dispersion parameter
xi_est      <- median(fit$xi)
xi_bias     <- xi_est - xi
xi_lower    <- quantile(fit$xi, 0.025)
xi_upper    <- quantile(fit$xi, 0.975)
xi_coverage <- as.numeric(xi >= xi_lower & xi <= xi_upper)

# Spatial correlation
rho_est       <- median(fit$rho)
rho_bias      <- rho_est - rho
rho_lower     <- quantile(fit$rho, 0.025)
rho_upper     <- quantile(fit$rho, 0.975)
rho_coverage  <- as.numeric(rho >= rho_lower & rho <= rho_upper)

# Spatial marginal variance
tau2_est      <- median(fit$tau2)
tau2_bias     <- tau2_est - tau2
tau2_lower    <- quantile(fit$tau2, 0.025)
tau2_upper    <- quantile(fit$tau2, 0.975)
tau2_coverage <- as.numeric(tau2 >= tau2_lower & tau2 <= tau2_upper)

# Spatial random effects
nu_true     <- nu - mean(nu)
nu_est      <- apply(fit$nu, 2, median)
nu_bias     <- mean(nu_est - nu_true)
nu_lower    <- apply(fit$nu, 2, quantile, probs = 0.025)
nu_upper    <- apply(fit$nu, 2, quantile, probs = 0.975)
nu_coverage <- mean(as.numeric(nu_true >= nu_lower & nu_true <= nu_upper))
nu_rmse     <- sqrt(mean((nu_est - nu_true)^2))

# BART component
bart_true     <- f(z) + alpha + mean(nu) + log(xi)
bart_mean     <- apply(pred_fun(z), 1, mean)
bart_median   <- apply(pred_fun(z), 1, median)
bart_lower    <- apply(pred_fun(z), 1, quantile, probs = 0.025)
bart_upper    <- apply(pred_fun(z), 1, quantile, probs = 0.975)
bart_bias     <- mean(bart_median - bart_true)
bart_coverage <- mean(as.numeric(bart_true >= bart_lower & bart_true <= bart_upper))
bart_rmse     <- sqrt(mean((bart_mean - bart_true)^2))


# Compile Results ---------------------------------------------------------

stats <- data.frame(
  param = c('alpha', paste0('beta', 1:nx), 'xi', 'rho', 'tau2', 'nu', 'G'),
  est = c(alpha_est, beta_est, xi_est, rho_est, tau2_est, NA, NA),
  lower = c(alpha_lower, beta_lower, xi_lower, rho_lower, tau2_lower, NA, NA),
  upper = c(alpha_upper, beta_upper, xi_upper, rho_upper, tau2_upper, NA, NA),
  bias = c(alpha_bias, beta_bias, xi_bias, rho_bias, tau2_bias, nu_bias, bart_bias),
  coverage = c(alpha_coverage, beta_coverage, xi_coverage, rho_coverage,
               tau2_coverage, nu_coverage, bart_coverage),
  rmse = c(NA, rep(NA, nx), NA, NA, NA, nu_rmse, bart_rmse)
)
stats$id <- id

results <- list(stats = stats, ale1 = ale1, ale2 = ale2, ale3 = ale3)


# Output Results ----------------------------------------------------------

fname <- paste0(id, '.rds')
saveRDS(results, here('sim-study','Results','temp', fname))
