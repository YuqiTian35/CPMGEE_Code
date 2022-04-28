#### Application - CD4:CD8 Ratio ####

source('conditional_function.R')

# load data
hgb.d<-read.csv("cd4_cd8_base_hgb_28jul2015.csv") 
ncds.d<-read.csv("cd4cd8rat_ncds_04jun2015.csv")
demo.d<-read.csv("cd4cd8rat_pts_04jun2015.csv")
labs.d<-read.csv("cd4cd8ratio_labs_04jun2015.csv")
extra.d<-read.csv("extra_variables_26oct2015.csv")
unk.vl.d<-read.csv("pts_unk_ud_vl_duration_22sep2015.csv")

# recode data
d<-demo.d
d$cd4<-d$BASELINE_CD4_COUNT
d$cd8<-d$BASELINE_CD8_COUNT
d$cd4.cd8.ratio<-with(d, cd4/cd8)
d$y<-d$cd4.cd8.ratio

d<-d[!is.na(d$y),]
d<-d[!is.na(d$BASELINE_VL),]
d<-d[d$BASELINE_VL<400,]

d$race<-with(d, ifelse(raceeth=="African American", "African American",
                       ifelse(raceeth=="Caucasian", "Caucasian",
                              ifelse(raceeth=="Hispanic", "Hispanic", "Other"))))
d$female<-ifelse(d$sex=="F",1,0)
d$age<-d$AGE_AT_T0
d$route<-with(d, ifelse(RISK_FACTOR_1=="Other"|RISK_FACTOR_1=="Unknown","Other/Unknown",
                        as.character(RISK_FACTOR_1)))
d$hcv<-ifelse(is.na(d$HCV_DX_PRIOR_T0),0,d$HCV_DX_PRIOR_T0)
d$hbv<-ifelse(is.na(d$HBV_DX_PRIOR_T0),0,d$HBV_DX_PRIOR_T0)
d$year<-d$YEAR_OF_T0
d$y1<-with(d,ifelse(y<=0.1,.1,
                    ifelse(y>=2,2,y)))
d$y2<-with(d,ifelse(y<=0.2,.2,
                    ifelse(y>=1.5,1.5,y)))

# recode
cd4.d1<-labs.d[labs.d$testName=="CD4 COUNT",c("AGE_AT_RESULT_DATE","CFAR_PID","RESULT_NUMERIC")]
cd4.d1$age<-cd4.d1$AGE_AT_RESULT_DATE
cd4.d1$cd4<-cd4.d1$RESULT_NUMERIC
cd4.d<-cd4.d1[,c("CFAR_PID","age","cd4")]
cd8.d1<-labs.d[labs.d$testName=="CD8 COUNT",]
cd8.d1$age<-cd8.d1$AGE_AT_RESULT_DATE
cd8.d1$cd8<-cd8.d1$RESULT_NUMERIC
cd8.d<-cd8.d1[,c("CFAR_PID","age","cd8")]

d1<-merge(cd4.d, cd8.d, by=c("CFAR_PID","age"), all=FALSE)
d1$ratio<-with(d1, cd4/cd8)

d$age.b<-d$age
d$ratio.b<-d$cd4.cd8.ratio
d.basic<-d[,c("CFAR_PID","age.b","ratio.b","race","female","route","hcv","hbv","year")]

d.long<-merge(d1,d.basic,by="CFAR_PID",all.x=TRUE)
ord<-order(d.long$CFAR_PID,d.long$age)
d.long<-d.long[ord,]
d.long<-d.long[!is.na(d.long$age.b),]   #### removing those not meeting inclusion criteria

d.long$years<-d.long$age-d.long$age.b

d.long5<-d.long[d.long$years<5,]
d.long1<-d.long[d.long$years<1,]

d.long1$CFAR_PID <- as.character(d.long1$CFAR_PID)

# Select people with 1-7 observations (99% subjects) 
dat <- d.long1 %>% 
  filter(CFAR_PID %in% names(table(d.long1$CFAR_PID))[table(d.long1$CFAR_PID) <= 7])

data.frame(number_of_records = 1:7,
           n = sapply(1:7 , function(x) sum((as.numeric(table(dat$CFAR_PID))  == x))),
           pct = sapply(1:7 , function(x) round(100 * mean((as.numeric(table(dat$CFAR_PID))  == x)), 2))) %>% 
  kable %>% kable_styling(full_width = F)

## add visit variable for each subjects
dat$visit <- lapply(table(dat$CFAR_PID), function(x) seq(1, x)) %>% unlist
# recode baseline year
dat$year_adjusted <- dat$year - 1999

#### CPM independnece ####
dd1 <- datadist(dat)
options(datadist='dd1')
mod_cpm <- orm(ratio ~ years + year_adjusted + race + age.b+ female + route + hcv + hbv,
               data = dat,
               x=T, y=T)
mod_cpm <- robcov(fit = mod_cpm, cluster= dat$CFAR_PID)

# conditional metrics
## new data
p <- length(mod_cpm$coefficients) - (length(mod_cpm$yunique) - 1)
new.data <- matrix(NA, ncol = p, nrow=length(seq(0.05, 0.95, 0.05)))
colnames(new.data) <- names(mod_cpm$coefficients[length(mod_cpm$yunique):length(mod_cpm$coefficients)])

new <- seq(0.05, 0.95, 0.05)
new.data[,1] <- new
new.data[,2] <- 10 # year_adjusted
new.data[,3] <- 1 # Caucasian
new.data[,4] <- 0 # Caucasian 
new.data[,5] <- 0 # Caucasian 
new.data[,6] <- median(dat$age.b) # age.b
new.data[,7] <- 0 # male
new.data[,8] <- 0 # MSM
new.data[,9] <- 1 # MSM
new.data[,10] <- 0 # MSM
new.data[,11] <- 0 # HCV
new.data[,12] <- 0 # HBV

mean_cpm <- mean.orm(mod_cpm, new.data = new.data, se=TRUE)
median_cpm <- quantile.orm(mod_cpm, new.data = new.data, probs = 0.5, se=TRUE)
cdf_est <- cdf.orm(mod_cpm, new.data, at.y=1, se = TRUE)


# convert id to numeric
dat$id <- as.numeric(factor(dat$CFAR_PID, 
                            levels=unique(dat$CFAR_PID)))

#### CPM exchangeable - rounding 2 decimal place #### 
dat2 <- dat
dat2$ratio <- round(dat$ratio, 2)
dat2$ratio_ord <- as.numeric(factor(dat2$ratio,
                                    levels=sort(unique(dat2$ratio))))

mod_repolr_round2 <- cpmgee(formula = ratio_ord ~ years + year_adjusted + race + age.b+ female + route + hcv + hbv, 
                                       data = dat2,
                                       categories = length(unique(dat2$ratio)),
                                       subjects = 'id',
                                       initial = 'orm',
                                       times = 1:7,
                                       corr.mod = 'uniform',
                                       alpha = 0.5)

mean_round <- mean.repolr(mod_repolr_round2, dat2$ratio, new.data)
median_round <- quantile.repolr(mod_repolr_round2, dat2$ratio, new.data, probs = 0.5)
cdf_round <- cdf.repolr(mod_repolr_round2, dat2$ratio, new.data, at.y=1)

#### CPM exchangeable - binning Mb=1000 ####

binning <- function(data, k){
  N <- length(data$ratio)
  q <- N %/% k
  r <- N %% k
  random_list <- sample(c(rep(q, k - r), rep(q + 1, r)))
  # median for each bin
  bin_median <- data %>% 
    arrange(ratio) %>% 
    mutate(ratio_b_ord = unlist(lapply(1:k, function(i) return(rep(i, random_list[i]))))) %>% 
    group_by(ratio_b_ord) %>% 
    summarise(ratio_b = median(ratio)) 
  # add both binning and binning ordinal
  data <- data %>% 
    arrange(ratio) %>% 
    mutate(ratio_b_ord = unlist(lapply(1:k, function(i) return(rep(i, random_list[i]))))) %>% 
    left_join(bin_median, by='ratio_b_ord') 
  # reorder the ordinal because some medians are the same
  data$ratio_b_ord <- dense_rank(data$ratio_b)
  data$ratio_cont <- data$ratio
  data$ratio_1000 <- data$ratio_b_ord
  data$ratio <- data$ratio_b
  data <- data %>% 
    dplyr::select(-c(ratio_b_ord, ratio_b)) %>% 
    arrange(CFAR_PID, visit)
  return(data)
}


dat <- binning(dat, k=1000)

mod_repolr_1000 <- cpmgee(formula = ratio_1000 ~ years + year_adjusted + race + age.b+ female + route + hcv + hbv,
                          data = dat,
                          categories = length(unique(dat$ratio_1000)),
                          subjects = 'id',
                          initial = 'orm',
                          times = 1:5,
                          corr.mod = 'uniform',
                          alpha = 0.5)

mean_bin <- mean.repolr(mod_repolr_1000, dat$ratio, new.data)
median_bin <- quantile.repolr(mod_repolr_1000, dat$ratio, new.data, probs = 0.5)
cdf_bin <- cdf.repolr(mod_repolr_1000, dat$ratio, new.data, at.y=1)


#### GEE exchangeable ####
mod_gee_ex <- geeglm(log(ratio) ~ years + year_adjusted + race + age.b + female + route + hcv + hbv,
                     data = dat,
                     id = dat$CFAR_PID,
                     corstr = "exchangeable")
new.data.gee <- cbind(1, new.data)
mean_gee <- new.data.gee %*% coef(mod_gee_ex)
se_gee <- sqrt(diag(new.data.gee %*%tcrossprod(mod_gee_ex$geese$vbeta, new.data.gee)))
ci_lower_gee <- mean_gee - qnorm(0.975) * se_gee
ci_upper_gee <- mean_gee + qnorm(0.975) * se_gee



