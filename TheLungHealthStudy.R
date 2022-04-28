#### Application - The Lung Health Study #### 
source('conditional_function.R')

# load data
load('lhs4fev2.rdata')
dat <- lhs4fev2
dat$visit <- dat$visit - 1
# remove missing data
dat <- dat %>% filter(id != 377418)
dat <- dat %>% mutate(site = as.factor(site))

# CPM independence
mod_cpm_ind <- orm(fev ~ visit * rs12194741 + bmi + age.10 + bmi0.5 + cigs0.10 + pack0.20 + site,
                 data = dat, x=T, y=T)
mod_cpm_ind <- robcov(fit = mod_cpm_ind, cluster= dat$id)

# new data for conditional metrics
p <- length(mod_cpm_ind$coefficients) - (length(mod_cpm_ind$yunique) - 1)

new.data.0 <- new.data.1 <- matrix(NA, ncol = p, nrow=length(seq(0, 4)))
colnames(new.data.0) <- colnames(new.data.1) <-
  names(mod_cpm_ind$coefficients[length(mod_cpm_ind$yunique):length(mod_cpm_ind$coefficients)])

new <- seq(0, 4, 1)
new.data.0[,1] <- new.data.1[,1] <- new
new.data.0[,2] <- 0; new.data.1[,2] <- 1 # rs12194741
new.data.0[,3] <- new.data.1[,3] <- median(dat$bmi) # bmi
new.data.0[,4] <- new.data.1[,4] <- median(dat$age.10)
new.data.0[,5] <- new.data.1[,5] <- median(dat$bmi0.5)
new.data.0[,6] <- new.data.1[,6] <- median(dat$cigs0.10)
new.data.0[,7] <- new.data.1[,7] <- median(dat$pack0.20)
new.data.0[,8] <- new.data.1[,8] <- 0 # site 2
new.data.0[,9] <- new.data.1[,9] <- 0 # site 3
new.data.0[,10] <- new.data.1[,10] <- 0 # site 4
new.data.0[,11] <- new.data.1[,11] <- 0 # site 5
new.data.0[,12] <- new.data.1[,12] <- 0 # site 6
new.data.0[,13] <- new.data.1[,13] <- 0 # site 7
new.data.0[,14] <- new.data.1[,14] <- 1 # site 8
new.data.0[,15] <- new.data.1[,15] <- 0 # site 9
new.data.0[,16] <- new.data.1[,16] <- 0 # site 10
new.data.0[,17] <- 0; new.data.1[,17] <- seq(0, 4, 1) # visit*rs12194741

mean_cpm_0 <- mean.orm(mod_cpm_ind, new.data = new.data.0, se=TRUE)
mean_cpm_1 <- mean.orm(mod_cpm_ind, new.data = new.data.1, se=TRUE)
median_cpm_0 <- quantile.orm(mod_cpm_ind, new.data = new.data.0, probs = 0.5, se=TRUE)
median_cpm_1 <- quantile.orm(mod_cpm_ind, new.data = new.data.1, probs = 0.5, se=TRUE)
cdf_cpm_0 <- cdf.orm(mod_cpm_ind, new.data.0, at.y=3, se = TRUE)
cdf_cpm_1 <- cdf.orm(mod_cpm_ind, new.data.1, at.y=3, se = TRUE)


#### CPM AR1 ####
mod_cpm_ar1 <- cpmgee(formula = fev_ord ~ visit * rs12194741 + bmi + age.10 + bmi0.5 + cigs0.10 + pack0.20 + site, 
                      data = dat0,
                      categories = length(unique(dat0$fev_ord)),
                      subjects = 'id',
                      initial = 'orm',
                      times = 1:5,
                      corr.mod = 'ar1',
                      alpha = 0.5)

mean_ar1_0 <- mean.repolr(mod_cpm_ar1, dat0$fev, new.data.0)
mean_ar1_1 <- mean.repolr(mod_cpm_ar1, dat0$fev, new.data.1)
median_ar1_0 <- quantile.repolr(mod_cpm_ar1, dat0$fev, new.data.0, probs = 0.5)
median_ar1_1 <- quantile.repolr(mod_cpm_ar1, dat0$fev, new.data.1, probs = 0.5)
cdf_ar1_0 <- cdf.repolr(mod_cpm_ar1, dat0$fev, new.data.0, at.y=3)
cdf_ar1_1 <- cdf.repolr(mod_cpm_ar1, dat0$fev, new.data.1, at.y=3)
