
# Save a softbart object
save_fit <- function (model, file) {
  saved_fit <- list(fit = model,
                    hypers = model$bart$get_hypers(),
                    opts = model$bart$get_opts(),
                    forests = model$bart$get_saved_forests())
  saveRDS(saved_fit, file)
}

# Load a softbart object
read_fit <- function (file) {
  saved_fit <- readRDS(file)
  model <- saved_fit$fit
  model$bart <- SoftBart::MakeForest(hypers = saved_fit$hypers, 
                                     opts = saved_fit$opts, 
                                     saved_forests = saved_fit$forests, 
                                     warn = FALSE)
  return (model)
}

# Random truncated normal
rltnorm <- function (mean = 0, sd = 1, lower = 0) {
  x <- -1
  while (x < lower) x <- stats::rnorm(n = 1, mean = mean, sd = sd)
  return (x)
}

# Log density of truncated normal
dltnorm <- function (x, mean = 0, sd = 1, lower = 0) {
  if (x <= lower) return (0)
  log_pdf_normal <- dnorm(x, mean = mean, sd = sd, log = TRUE)
  cdf_trunc <- pnorm(lower, mean = mean, sd = sd)
  log_cdf_trunc <- log1p(-cdf_trunc)
  
  return (log_pdf_normal - log_cdf_trunc)
}


# Rescale x values back to original scale
rescale_ale1 <- function (z, ale) {
  var <- unique(ale$var)
  ale$x <- quantile(z[, var], seq(0, 1, length.out = max(ale$k)))
  return (ale)
}

# Rescale x values and plotting coordinates back to original scale
rescale_ale2 <- function (z, ale) {
  
  # Determine original scale of first predictor
  var1 <- unique(ale$var1)
  x1old <- unique(ale$x1)
  x1new <- quantile(z[, var1], seq(0, 1, length.out = max(ale$k1)))
  
  # Determine original scale of second predictor
  var2 <- unique(ale$var2)
  x2old <- unique(ale$x2)
  x2new <- quantile(z[, var2], seq(0, 1, length.out = max(ale$k2)))
  
  # Determine original scale of plotting boundaries
  w1new <- c(diff(x1new)[1], diff(x1new)) / 2
  w2new <- c(rev(abs(diff(rev(x1new)))), rev(abs(diff(x1new)))[1]) / 2
  h1new <- c(diff(x2new)[1], diff(x2new)) / 2
  h2new <- c(rev(abs(diff(rev(x2new)))), rev(abs(diff(x2new)))[1]) / 2
  
  # Update ALE object
  index <- 1
  for (i in 1:length(x1new)) {
    for (j in 1:length(x2new)) {
      ale$x1[index] <- x1new[i]
      ale$x2[index] <- x2new[j]
      ale$w1[index] <- w1new[i]
      ale$w2[index] <- w2new[i]
      ale$h1[index] <- h1new[j]
      ale$h2[index] <- h2new[j]
      index <- index + 1
    }
  }
  return (ale)
}
