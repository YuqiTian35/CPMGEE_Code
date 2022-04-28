#### the primary setting (alpha=0.7) ####
# section 4.1

library(diagonals)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(repolr)
library(rms)
library(Matrix)
library(tidyverse)
library(rlist)
library(SimCorMultRes)
library(geepack)
library(cpmgee)

sample_size <- 1000
cluster_size <- 6

# beta
beta_coefficients <- c(1, 1)
# covariate
x <- rep(rnorm(sample_size), each = cluster_size) # time-invariant
t <- rep(seq(0, 1, 0.2), sample_size)

# latent correlation matrix for the NORTA method
latent_correlation_matrix <- matrix(0.7, nrow=cluster_size, ncol=cluster_size)
diag(latent_correlation_matrix) <- 1

data <- simulated_cont_data$simdata %>% mutate(y = qchisq(pnorm(y/2), 5)) %>% mcar
data <- data %>% binning(k=300)

# simulation of ordinal responses
simulated_cont_data <- rcont.clm(clsize = cluster_size, 
                                 betas = beta_coefficients, 
                                 xformula = ~ x + t, 
                                 cor.matrix = latent_correlation_matrix)

# fully-iterated repolr
mod_repolr_full <- repolr_full(formula = yc ~ x + t,
                               data = data,
                               categories = length(unique(data$y)),
                               subjects = 'id',
                               initial = 'orm',
                               times = 1:6,
                               corr.mod = 'uniform',
                               alpha = 0.5)

# new data
new.data <- data.frame(x=c(0, 1), t=0.2)

# conditional metrics
cond_mean_full <- mean.cpmgee(mod_repolr_full, data$y, new.data)
cond_median_full <- quantile.cpmgee(mod_repolr_full, data$y, new.data, probs = 0.5)
cond_cdf_full <- cdf.cpmgee(mod_repolr_full, data$y, new.data, at.y=5)

