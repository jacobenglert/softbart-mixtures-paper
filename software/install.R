# Program Name: install.R
# Purpose:      install R packages needed to run the code in this repo.


# Install from CRAN -------------------------------------------------------
install.packages(c('remotes',
                   'here',
                   'tidyverse',
                   'sf',
                   'mvtnorm',
                   'spdep',
                   'spam',
                   'progress'))



# Install from GitHub -----------------------------------------------------

library(remotes)

# SoftBart is available on CRAN, but this fork contains a modification
# that allows for preserving model fits across sessions.
install_github('jacobenglert/SoftBART')

# Pólya-Gamma random variable generator
install_github('jacobenglert/jrpg', ref = 'v0.1.0', build_vignettes = TRUE)

# Partial dependence for posterior distributions
install_github('jacobenglert/pdpd', ref = 'v0.1.1', build_vignettes = TRUE)