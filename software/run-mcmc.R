
source(here::here('software','helper-functions.R'))

run_mcmc <- function (x, y, z, offset = NULL, 
                      s, geometry,
                      num_tree = 20, k = 2, base = 0.95, power = 2,
                      num_burn = 5000, num_thin = 10, num_save = 1000,
                      soft = TRUE, sparse = TRUE
) {
  
  x <- as.matrix(x)
  z <- as.matrix(z)
  
  if (is.null(offset)) offset <- rep(0, nrow(x))
  
  # Sample sizes and dimensions
  n   <- length(y)
  ns  <- length(unique(s))
  p_x <- ncol(x)
  p_z <- ncol(z)
  
  # Set up design matrix for spatial random effects
  s_unique <- match(unique(s), s)
  x_nu <- spam::spam(0, nrow = n, ncol = ns)
  x_nu[cbind(1:n, match(s, unique(s)))] <- 1
  
  # Determine neighborhood adjaceny matrix and number of neighbors
  W <- spdep::nb2mat(spdep::poly2nb(geometry[s_unique]), style = 'B')
  D <- diag(rowSums(W))

  # Initialize parameters for MCMC sampler
  G     <- numeric(n)           # BART predictor
  beta  <- numeric(p_x)         # Confounder regression coefficients
  xi    <- 1                    # Dispersion parameter
  nu    <- stats::rnorm(ns)     # Spatial random effects
  tau2  <- 1 / stats::rgamma(1, 0.1 + ns / 2, 0.1 + (1 / 2))  # Spatial variance
  rho   <- 0.9                                            # Spatial correlation
  
  # Linear predictor
  fixeff <- as.numeric(x %*% beta)
  raneff <- as.numeric(x_nu %*% nu)
  eta <- offset + fixeff + G + raneff
  
  # Pre-calculate discrete prior distribution for rho
  lambda    <- eigen(solve(D) %*% W, only.values = TRUE)$values
  rho_vals  <- stats::qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
  rho_ll0   <- sapply(rho_vals, \(x) 0.5 * sum(log(1 - x * lambda)), simplify = TRUE)
  
  # Specify fixed effect prior distribution
  b <- rep(0, p_x)
  B <- diag(p_x) * 1e4
  B_inv <- diag(1 / diag(B))
 
  # Tuning parameter for xi proposal distribution (MH option)
  # s_xi <- 0.1
  
  # Create BART objects
  bart_hypers <- SoftBart::Hypers(X = z, Y = G, sigma_hat = 1,
                                  num_tree = num_tree, 
                                  beta = power, gamma = base, k = k)
  bart_opts <- SoftBart::Opts(update_sigma = FALSE,
                              num_burn = num_burn,
                              num_save = num_save * num_thin,
                              num_print = 50000,
                              update_s = TRUE, update_alpha = TRUE)
  
  if (!sparse) {
    bart_opts$update_s <- FALSE
    bart_opts$update_alpha <- FALSE
  }
  
  if (!soft) {
    bart_hypers$width <- 1e-4
    bart_opts$update_tau <- FALSE
  }
  
  sampler <- SoftBart::MakeForest(hypers = bart_hypers, opts = bart_opts, warn = FALSE)
  
  # Allocate posterior storage
  post <- list(
    alpha = numeric(num_save),
    beta = matrix(nrow = num_save, ncol = p_x, dimnames = list(NULL, colnames(x))),
    nu = matrix(nrow = num_save, ncol = ns),
    xi = numeric(num_save),
    rho = numeric(num_save),
    tau2 = numeric(num_save),
    var_counts = matrix(0, nrow = num_save, ncol = p_z, dimnames = list(NULL, colnames(z))),
    logmean = numeric(num_save),
    loglik = numeric(num_save)
  )
  
  # Allocate running storage for WAIC calculation
  ll <- ll2 <- ell <- numeric(n)
  
  # MCMC
  xi_acc <- 0
  pb <- progress::progress_bar$new(
    format = "[:bar] Burning in :current/:total. Total time elapsed: :elapsedfull",
    total = num_burn, clear = FALSE, width = 100)
  for (k in seq_len(num_burn)) {
    
    # Step 1) Sample latent Polya-Gamma random variables
    omega <- jrpg::jrpg(y + xi, eta)[,1]
    
    # Convert to Gaussian form
    y_star <- (y - xi) / (2 * omega)  # y_star ~ N(eta, diag(1 / omega))
    
    # Update spatial weight matrices
    D_rho_W <- spam::as.spam(D - rho * W)
      
    # Step 2) Update spatial random effects
    r_nu <- y_star - offset - fixeff - G
    nu_Sigma <- spam::solve(spam::crossprod.spam(x_nu * sqrt(omega)) + (1 / tau2) * (D_rho_W))
    nu_mu <- nu_Sigma %*% spam::crossprod.spam(x_nu, omega * r_nu)
    nu <- mvtnorm::rmvnorm(n = 1, mean = nu_mu, sigma = nu_Sigma)[1,]
    nu <- nu - mean(nu)
    raneff <- as.numeric(x_nu %*% nu)
    
    # Step 3) Update spatial random effects variance
    tau2 <- 1 / stats::rgamma(1, 0.1 + ns / 2, 0.1 + (nu %*% (D_rho_W) %*% nu) / 2)
    
    # Step 4) Update spatial random effects correlation
    rho_ll <- rho_ll0 + rho_vals / (2 * tau2) * as.numeric(nu %*% W %*% nu)
    rho <- sample(rho_vals, size = 1, prob = exp(rho_ll - max(rho_ll)))
    
    # Step 5) Update fixed effects
    r_beta <- y_star - offset - G - raneff
    beta_Sigma <- solve(B_inv + crossprod(x * sqrt(omega)))
    beta_mu <- beta_Sigma %*% (B_inv %*% b + crossprod(x, omega * r_beta))
    beta <- mvtnorm::rmvnorm(n = 1, mean = beta_mu, sigma = beta_Sigma)[1,]
    fixeff <- as.numeric(x %*% beta)
    
    # Step 6) Update BART
    r_G <- y_star - offset - fixeff - raneff
    G <- sampler$do_gibbs_weighted(z, r_G, omega, z, 1)[1,]
    
    # Update linear predictor
    eta <- offset + fixeff + G + raneff
    
    # Step 7) Update dispersion parameter
    # # (Metropolis-Hastings proposed from centered truncated normal)
    # q <- pmin(0.9999, 1 / (1 + exp(eta))) # 1 - Pr(success)
    # xi_prop <- rltnorm(mean = xi, sd = s_xi, lower = 0)
    # r_xi <- sum(dnbinom(y, size = xi_prop, prob = q, log = TRUE)) -
    #   sum(dnbinom(y, size = xi, prob = q, log = TRUE)) +
    #   dltnorm(xi, mean = xi_prop, sd = s_xi, lower = 0) -
    #   dltnorm(xi_prop, mean = xi, sd = s_xi, lower = 0)
    # 
    # if (log(stats::runif(1)) < r_xi) {
    #   xi <- xi_prop
    #   xi_acc <- xi_acc + 1
    # }
    # 
    # # Update tuning parameter
    # if (k %% 100 == 0) {
    #   xi_acc_rate <- xi_acc / k
    #   if(xi_acc_rate > 0.6) s_xi <- 1.1 * s_xi
    #   if(xi_acc_rate < 0.2) s_xi <- 0.8 * s_xi
    #   cat(s_xi)
    # }
    
    # (Gibbs)
    l <- sapply(1:n, function (i) sum(rbinom(y[i], 1, round(xi / (xi + 1:y[i] - 1), 6)))) # Could try to avoid loop; in rounding avoids numerical stability
    
    # Update r from conjugate gamma distribution given l and psi
    q <- pmin(0.9999, 1 / (1 + exp(eta))) # 1 - Pr(success)
    xi <- rgamma(1, 0.01 + sum(l), 0.01 - sum(log(q)))
    
    
    pb$tick()
  }
  
  # Saved iterations
  pb <- progress::progress_bar$new(
    format = "[:bar] Sampling from posterior :current/:total. Total time elapsed: :elapsedfull",
    total = num_save, clear = FALSE, width = 100)
  for (k in seq_len(num_save)) {
    
    for (k2 in seq_len(num_thin)) {
      
      # Step 1) Sample latent Polya-Gamma random variables
      omega <- jrpg::jrpg(y + xi, eta)[,1]
      
      # Convert to Gaussian form
      y_star <- (y - xi) / (2 * omega)  # y_star ~ N(eta, diag(1 / omega))
      
      # Update spatial weight matrices
      D_rho_W <- spam::as.spam(D - rho * W)
      
      # Step 2) Update spatial random effects
      r_nu <- y_star - offset - fixeff - G
      nu_Sigma <- spam::solve(spam::crossprod.spam(x_nu * sqrt(omega)) + (1 / tau2) * (D_rho_W))
      nu_mu <- nu_Sigma %*% spam::crossprod.spam(x_nu, omega * r_nu)
      nu <- mvtnorm::rmvnorm(n = 1, mean = nu_mu, sigma = nu_Sigma)[1,]
      nu <- nu - mean(nu)
      raneff <- as.numeric(x_nu %*% nu)
      
      # Step 3) Update spatial random effects variance
      tau2 <- 1 / stats::rgamma(1, 0.1 + ns / 2, 0.1 + (nu %*% (D_rho_W) %*% nu) / 2)
      
      # Step 4) Update spatial random effects correlation
      rho_ll <- rho_ll0 + rho_vals / (2 * tau2) * as.numeric(nu %*% W %*% nu)
      rho <- sample(rho_vals, size = 1, prob = exp(rho_ll - max(rho_ll)))
      
      # Step 5) Update fixed effects
      r_beta <- y_star - offset - G - raneff
      beta_Sigma <- solve(B_inv + crossprod(x * sqrt(omega)))
      beta_mu <- beta_Sigma %*% (B_inv %*% b + crossprod(x, omega * r_beta))
      beta <- mvtnorm::rmvnorm(n = 1, mean = beta_mu, sigma = beta_Sigma)[1,]
      fixeff <- as.numeric(x %*% beta)
      
      # Step 6) Update BART
      r_G <- y_star - offset - fixeff - raneff
      G <- sampler$do_gibbs_weighted(z, r_G, omega, z, 1)[1,]
      
      # Update linear predictor
      eta <- offset + fixeff + G + raneff
      
      # Step 7) Update dispersion parameter
      # # (Metropolis-Hastings proposed from centered truncated normal)
      # q <- pmin(0.9999, 1 / (1 + exp(eta))) # 1 - Pr(success)
      # xi_prop <- rltnorm(mean = xi, sd = s_xi, lower = 0)
      # r_xi <- sum(dnbinom(y, size = xi_prop, prob = q, log = TRUE)) -
      #   sum(dnbinom(y, size = xi, prob = q, log = TRUE)) +
      #   dltnorm(xi, mean = xi_prop, sd = s_xi, lower = 0) -
      #   dltnorm(xi_prop, mean = xi, sd = s_xi, lower = 0)
      # 
      # if (log(stats::runif(1)) < r_xi) {
      #   xi <- xi_prop
      # }
      
      # (Gibbs)
      l <- sapply(1:n, function (i) sum(stats::rbinom(y[i], 1, round(xi / (xi + 1:y[i] - 1), 6)))) # Could try to avoid loop; in rounding avoids numerical stability
      
      # Update r from conjugate gamma distribution given l and psi
      q <- pmin(0.9999, 1 / (1 + exp(eta))) # 1 - Pr(success)
      xi <- stats::rgamma(1, 0.01 + sum(l), 0.01 - sum(log(q)))
    }
    
    # Store current values
    post$alpha[k]        <- mean(G)
    post$beta[k,]        <- beta
    post$nu[k,]          <- nu
    post$xi[k]           <- xi
    post$rho[k]          <- rho
    post$tau2[k]         <- tau2
    post$var_counts[k,]  <- sampler$get_counts()[,1]
    post$logmean[k]      <- log(xi) + post$alpha[k]
    
    # WAIC
    ll_curr <- stats::dnbinom(y, size = xi, prob = q, log = TRUE)
    ll  <- ll + ll_curr
    ll2 <- ll2 + ll_curr^2
    ell <- ell + exp(ll_curr)
    
    post$loglik[k] <- sum(ll_curr)
    
    pb$tick()
  }
  
  # Manually thin BART sampler
  keep_iters <- seq(from = num_burn + num_thin, 
                    to = num_burn + num_save * num_thin, 
                    by = num_thin)
  sampler <- SoftBart::MakeForest(
    hypers = sampler$get_hypers(),
    opts = sampler$get_opts(),
    saved_forests = sampler$get_saved_forests()[keep_iters],
    warn = FALSE
  )
  
  post$bart <- sampler
  
  # Update WAIC metrics
  pWAIC       <- sum(ll2 / num_save - (ll / num_save)^2)
  lppd        <- sum(log(ell / num_save))
  post$pWAIC  <- pWAIC
  post$lppd   <- lppd
  post$WAIC   <- -2*lppd +2*pWAIC
  
  return (post)
}


