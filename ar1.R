####  AR1 structure ####
# Web appendix C - 3.4

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
latent_correlation_matrix <- diag(cluster_size)  
for(i in 1:(cluster_size-1)){
  for(j in (i+1):cluster_size){
    latent_correlation_matrix[i,j] <- latent_correlation_matrix[j,i] <-
      0.7^(j - i)
  }
}

data <- simulated_cont_data$simdata %>% mutate(y = qchisq(pnorm(y/2), 5)) %>% mcar
data <- data %>% binning(k=300)

# simulation of ordinal responses
simulated_cont_data <- rcont.clm(clsize = cluster_size, 
                                 betas = beta_coefficients, 
                                 xformula = ~ x + t, 
                                 cor.matrix = latent_correlation_matrix)
# CPM exchangeable
mod_cpm_ar1 <- cpmgee(formula = yc ~ x + t,
                     data = data,
                     categories = length(unique(data$y)),
                     subjects = 'id',
                     initial = 'orm',
                     times = 1:6,
                     corr.mod = 'ar1',
                     alpha = 0.5)

# CPM exchangeable
mod_cpm_ex <- cpmgee(formula = yc ~ x + t,
                     data = data,
                     categories = length(unique(data$y)),
                     subjects = 'id',
                     initial = 'orm',
                     times = 1:6,
                     corr.mod = 'uniform',
                     alpha = 0.5)

# CPM independence
mod_cpm_ind <- orm(y_cont ~ x + t, data = data, x=T, y=T)
mod_cpm_ind <- robcov_fast(fit = mod_cpm_ind, cluster= data$id)

# GEE - exchangeable
mod_gee_ex <- geeglm(2*qnorm(pchisq(y_cont, 5)) ~ x + t, 
                     data = data, id = data$id, corstr = 'exchangeable')


# new data
new.data <- data.frame(x=c(0, 1), t=0.2)

# conditional metrics for CPM exchangeable
cond_mean_ar1 <- mean.cpmgee(mod_cpm_ar1, data$y, new.data)
cond_median_ar1 <- quantile.cpmgee(mod_cpm_ar1, data$y, new.data, probs = 0.5)
cond_cdf_ar1 <- cdf.cpmgee(mod_cpm_ar1, data$y, new.data, at.y=5)

# conditional metrics for CPM exchangeable
cond_mean_ex <- mean.cpmgee(mod_cpm_ex, data$y, new.data)
cond_median_ex <- quantile.cpmgee(mod_cpm_ex, data$y, new.data, probs = 0.5)
cond_cdf_ex <- cdf.cpmgee(mod_cpm_ex, data$y, new.data, at.y=5)

# conditional metrics for CPM independence
cond_mean_ind <- mean.cpmgee(mod_cpm_ind, data$y, new.data)
cond_median_ind <- quantile.cpmgee(mod_cpm_ind, data$y, new.data, probs = 0.5)
cond_cdf_ind <- cdf.cpmgee(mod_cpm_ind, data$y,new.data, at.y = 5)

