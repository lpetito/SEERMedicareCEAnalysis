###############################
# Do a per-protocol analysis
# for erlotinib trial emulation
###############################

rm(list=ls())

# Load required libraries
require('survival')
require('rms')
#require('RcmdrPlugin.KMggplot2')
require('doParallel')
require('foreach')
require('DescTools')
require('reshape2')

#When on don, be in "proj/hiv_causal/programs"
sysinf <- Sys.info()
if(sysinf['sysname']=="Windows"){
  datpath <- 'W:/seermed2/data/lpetito/erlotinib/baseline.csv'
  savedir <- "W:/seermed2/programs/lpetito/gemzar/plots/eccp/"
} else{
  datpath <- '/proj/seermed2/data/lpetito/erlotinib/baseline.csv'
  savedir <- "/proj/seermed2/programs/lpetito/gemzar/plots/eccp/"
}

# Formulas to use in analysis for "outcome regression"
hrform <- Tevent ~  t + I(t*t) + trt_clone + stage_four + agedx + I(agedx*agedx) + ervisit_b + cci_b + perf_b + cholang_b + pneumon_b
ixform <- Tevent ~ t + I(t*t) + trt_clone + trt_clonet + trt_clonet2 + stage_four + agedx + I(agedx*agedx) + ervisit_b + cci_b + perf_b + cholang_b + pneumon_b

# Formulas to use in analysis for weights
wtformN <- A ~ t + I(t*t) + trt_b + Alag + stage_four + agedx + I(agedx*agedx) + ervisit_b + cci_b + perf_b + cholang_b + pneumon_b
wtformD <- A ~ t + I(t*t) + trt_b + Alag + stage_four + agedx + I(agedx*agedx) + ervisit_b + cci_b + perf_b + cholang_b + pneumon_b +
  ervisit + cci +cholang + pneumon

# Load and clean data, call it "df"
# This data should already be in long form, cloned, and expanded from SAS
# It needs to have the following variables to use in this analysis:

### Structural variables ###
# matchid - the ID that identifies individuals not clones
# newid - the ID that identifies each clone 
# t - time
# Tevent - binary, only 1 when event occurs
# A - indicator for whether this individual has ever used erlotinib in the past (or currently)
# Alag - A from time t-1 (lagged)
# trt_b - indicator for whether this individual initiated erlotinib before the end of the grace period (12 weeks)
# trt_clone - indicator for the artificial assignment of an individual to one arm or the other
# trt_clonet - product term between t and trt_clone
# trt_clonet2 - product term between t, t, and trt_clone
# touseweight - indicator for whether this observation should be used in weight creation
# touse - indicator for whether this observation should be used in the outcome regression
# Note: the completely expanded data is only used in creation of the adjusted (for baseline covariates) survival curves

### Baseline covariates ###
# stage_four - 1 if DX at stage 4, 0 otherwise
# agedx - age at DX
# ervisit_b - # of ED visits in year before DX
# cci_b - Charlson comorbidity score in year before DX
# perf_b - Davidoff score in year before DX (binary)
# cholang_b - Cholangitis or biliary tract obstruction in year before DX
# pneumon_b - pneumonia in year before DX

### Time-varying covariates ###
# ervisit - ER visit from t-1 to t
# cci - Charlson using data from t-11 to t
# cholang - cholangitis or biliarty tract obstruction from t-1 to t
# pneumon - pneumonia from t-1 to t

df <- read.csv(file=datpath, header=T)

# Weight creation ---------------------------------------------------------

create.weights <- function(dflong, wtformN, wtformD, q=0.995){
  ### Inputs
  # dflong - cleaned data from SAS
  # wtformN - numerator formula
  # wtformD - denominator formula
  # q - upper truncation limit - default is 99.5th percentile
  
  # Weights numerator
  wtmodelN <- glm(wtformN, data=dflong[dflong$touseweight==1 & dflong$t != 0,], family = quasibinomial())
  dflong$pnum <- predict(wtmodelN, newdata = dflong, type='response') # Pr A = 1 | baseline
  # Weights denominator
  wtmodelD <- glm(wtformD, data=dflong[dflong$touseweight==1 & dflong$t != 0,], family = quasibinomial())
  dflong$pdenom <- predict(wtmodelD, newdata = dflong, type='response') # Pr A = 1 | baseline, TV
  
  dflong$numCont <- with(dflong, A*pnum + (1-A)*(1-pnum))
  dflong$denCont <- with(dflong, A*pdenom + (1-A)*(1-pdenom))
  
  # Set contribution at baseline to 1
  dflong$numCont[dflong$t==0] <- 1
  dflong$denCont[dflong$t==0] <- 1
  
  # Numerator
  dflong$k1_0 <- with(dflong, ave(numCont, newid, FUN=cumprod))
  # Denominator
  dflong$k1_w <- with(dflong, ave(denCont, newid, FUN=cumprod))
  
  # Create both stabilized and unstabilized weights
  dflong$stabw <- with(dflong, k1_0 / k1_w)
  dflong$unstabw <- with(dflong, 1 / k1_w)
  
  tau <- quantile(dflong$stabw, q)
  dflong$stabw[dflong$stabw > tau] <- tau
  
  return(dflong)
}

df <- create.weights(df, wtformN, wtformD, q=0.99)


# Model fitting and summaries ---------------------------------------------

# Weight models - recreate to save
wtmodelN <- glm(wtformN, data=df[df$touseweight==1,], family = quasibinomial())
wtmodelD <- glm(wtformD, data=df[df==1,], family = quasibinomial())

# Unadjusted unweighted
plrfit <- glm(Tevent ~ t + I(t*t) + trt_clone, data = df[df$touse==1,], family=quasibinomial())
# summary(plrfit); exp(coef(plrfit)) # HR: 1.075

# Adjusted unweighted
adjplrfit <- glm(hrform, 
                 data = df[df$touse==1,], family=binomial())
# summary(adjplrfit); exp(coef(adjplrfit)) # HR: 1.032

# Adjusted weighted
wadjplrfit <- glm(hrform, 
                  data = df[df$touse==1,], weights = stabw, family=quasibinomial())
# summary(wadjplrfit); exp(coef(wadjplrfit)) # HR: 1.05

# Adjusted weighted ix
wadjplrixfit <- glm(ixform, 
                    data = df[df$touse==1,], weights=stabw, family=quasibinomial())
#summary(wadjplrixfit); exp(coef(wadjplrixfit)) # HR: 1.05

# Create model summaries to save 
splrfit <- summary(plrfit)$coefficients
sadjplrfit <- summary(adjplrfit)$coefficients
swadjplrfit <- summary(wadjplrfit)$coefficients
swadjplrixfit <- summary(wadjplrixfit)$coefficients
swtmodelN <- summary(wtmodelN)$coefficients
swtmodelD <- summary(wtmodelD)$coefficients
save(splrfit, sadjplrfit, swadjplrfit, swadjplrixfit, swtmodelN, swtmodelD,
     file=paste0(savedir, "model_summary.RData"))


# Creating data needed to produce adjusted survival curves ---------------------------------------

# In tarceva clones
tarceva <- df[df$trt_clone==1,]
tarceva$p <- 1 - predict(wadjplrixfit, newdata=tarceva, type='response')
tarceva$s <- with(tarceva, ave(p, matchid, FUN=cumprod))

# Now create gemzar only df - simulated
gemzar <- df[df$trt_clone==0,]
gemzar$p <- 1 - predict(wadjplrixfit, newdata=gemzar, type='response')
gemzar$s <- with(gemzar, ave(p, matchid, FUN=cumprod))

both <- rbind(tarceva, gemzar)
both <- both[,c('s', 'trt_clone', 't')]

# Calculate the mean survival at each visit within each treatment arm
results <- aggregate(s ~ t + trt_clone, FUN=mean, data=both)

# Add a variable that treats randomization as a factor
results$treatedf <- factor(results$trt_clone, labels = c("Gemcitabine alone", "Gemcitabine + Erlotinib"))

# Save results for later use
save(results,
     file = paste0(savedir, "surv_curves.RData"))


# Code to facilitate bootstrapping for survival curves --------------------

# Set the seed first
set.seed(123)
B=500 # Number of bootstrap reps
ids <- unique(df$matchid) # Unique IDs
nind <- length(ids) # Number of individuals
bootids1 <- matrix(ids[sample(nind, nind*B, replace=T)], nrow=nind, ncol=B) # This is a matrix that has B columns, with the n resampled IDs already selected

library('plyr')

bootids2 <- apply(bootids1, 2, FUN=function(b){
  temp <- data.frame(matchid = b)
  return(plyr::count(temp, vars = 'matchid'))
}) # This returns a matrix with 2 variables: the ID for each unique individual ('matchid'), and the frequency they appeared in a particular bootstrap sample ('freq')


# A slightly modified function to create bootstrapped weights - add frequency weights
create.weights.boot <- function(dflong, wtformN, wtformD, q=0.99){
  fitN <- glm(wtformN, data=dflong[dflong$touseweight==1,], family = quasibinomial(), weights = freq)
  
  dflong$pnum <- predict(glm(wtformN, data=dflong[dflong$touseweight==1,], family = quasibinomial(), weights = freq),
                         newdata = dflong, type='response') # Pr A = 1 | baseline
  # Weights denominator
  dflong$pdenom <- predict(glm(wtformD, data=dflong[dflong$touseweight==1,], family = quasibinomial(), weights = freq),
                           newdata = dflong, type='response') # Pr A = 1 | baseline, TV
  
  dflong$numCont <- with(dflong, A*pnum + (1-A)*(1-pnum))
  dflong$denCont <- with(dflong, A*pdenom + (1-A)*(1-pdenom))
  
  # Set contribution at baseline to 1
  dflong$numCont[dflong$t %in% c(0)] <- 1
  dflong$denCont[dflong$t %in% c(0)] <- 1
  
  # Numerator
  dflong$k1_0 <- with(dflong, ave(numCont, newid, FUN=cumprod))
  # Denominator
  dflong$k1_w <- with(dflong, ave(denCont, newid, FUN=cumprod))
  
  # Create both stabilized and unstabilized weights
  dflong$stabw <- with(dflong, k1_0 / k1_w)
  #dflong$unstabw <- with(dflong, 1 / k1_w)
  
  tau <- quantile(dflong$stabw, q, na.rm=T)
  dflong$stabw[dflong$stabw > tau] <- tau
  
  return(dflong)
}

system.time(bootres <- foreach(i=1:B, 
                               # .packages('data.table'), 
                               .combine=rbind) %do% {
                                 # Load libraries
                                 require(reshape2)
                                 require(plyr)
                                 require(dplyr)
                                 require(rms)
                                 
                                 # Add "freq" variable to bootdf
                                 bootdf <- left_join(df, bootids2[[i]], by = "matchid")
                                 bootdf$freq[is.na(bootdf$freq)] <- 0
                                 
                                 # Order data by newid to sort by clones instead of IDs
                                 bootdf <- bootdf[order(bootdf$matchid, bootdf$cloneid, bootdf$t),]
                                 
                                 # Create weights
                                 bootdf <- create.weights.boot(bootdf, wtformN, wtformD, q=0.99)
                                 
                                 if(sum(is.na(bootdf$stabw)) != 0){ # Something went horribly wrong
                                   return(data.frame(t = NA,
                                                                                     treatedf = NA,
                                                                                     s = NA,
                                                                                     boot = i))
                                 } else {
                                   bootdf$p <- 1 - predict(glm(ixform, 
                                                               data = bootdf[bootdf$touse==1,], weights = stabw*freq, family=quasibinomial()),
                                                           newdata=bootdf, type='response')
                                   #bootdf.expand$s <- with(bootdf.expand[bootdf.expand$freq != 0], ave(p, bid, FUN=cumprod))
                                   bootdf$s <- NA # Initialize cumulative survival probability
                                   bootdf$s[bootdf$freq != 0] <- with(bootdf[bootdf$freq != 0,], ave(p, newid, FUN=cumprod)) # Calculate 
                                   
                                   # Calculate the mean survival at each visit within each treatment arm
                                   results <- plyr::ddply(bootdf[,c("t", "trt_clone", "s", "freq")], .(t, trt_clone), summarise, s = weighted.mean(s, freq))
                                   
                                   # Add a variable that treats randomization as a factor
                                   results$treatedf <- factor(results$trt_clone, labels = c("Gemcitabine", "Erlotinib"))
                                   
                                   results <- results[,c("t", "treatedf", 's')]
                                   results$boot <- i
                                   return(results)
                                 } 
                                
                               })

save(bootres, file=paste0(savedir,"ppboot500.RData"))


# Code to calculate hazard ratio ------------------
basecoef <- coef(adjplrfit) # Estimates from adjusted unweighted hr model
tvcoef <- coef(wadjplrfit) # Estimates from adjusted weighted model
# Baseline dataset - keeps original people (n=940) only - 
baseline <- df[df$t==0 & df$touseweight==1,] 
baseline$agedx2 <- baseline$agedx*baseline$agedx
# IMPORTANT to keep variables in order from model fits
baseline <- baseline[,c("stage_four", "agedx", "agedx2", "ervisit_b", "cci_b", "perf_b", "cholang_b", "pneumon_b")]
# Only keep required things to keep memory use down
rm(list=setdiff(ls(), c("basecoef", 'tvcoef', 'baseline', 'savedir', 
                        'wtformN', 'wtformD', 'hrform', 'df',
                        'create.weights.boot')))

# Write function to do MC simulation to calculate marginal hazard ratios
MCsimHR <- function(i=1, Nsim = 100000, baseline, maxT=72, coeffit, dfSetA0, dfSetA1){
  
  # First select 100,000 individuals from uncloned data - only keep baseline variables (without A)
  set.seed(seeds[i])
  ids <- sample_n(baseline, size=Nsim, replace=T) # This is the data to sample from
  
  # Get predictions when A=1 
  predsA1 <- t(apply(ids, 1, FUN=function(x){
    ps <- 1 - 1 / (1 + exp( -1 * as.matrix(suppressWarnings(cbind(1, dfSetA1, t(x)))) %*% coeffit ))
    ps[1] <- 1 # Set probability of dying at time 0 to 1
    ss <- cumprod(ps)
    rbinom(maxT, 1, ss)
  }))
  # Get predictions when A=0 
  predsA0 <- t(apply(ids, 1, FUN=function(x){
    ps <- 1 - 1 / (1 + exp( -1 * as.matrix(suppressWarnings(cbind(1, dfSetA0, t(x)))) %*% coeffit ))
    ps[1] <- 1 # Set probability of dying at time 0 to 1
    ss <- cumprod(ps)
    rbinom(maxT, 1, ss)
  }))
  
  # Identify end of followup for each individual
  firstZeroA1 <- apply(predsA1, 1, FUN=function(x){
    if (sum(x) == maxT) { return(maxT+1) } else {
      return(which.min(x))
    }
  })
  
  firstZeroA0 <- apply(predsA0, 1, FUN=function(x){
    if (sum(x) == maxT) { return(maxT+1) } else {
      return(which.min(x))
    }
  })
  
  # Create data to fit model in
  survdf <- data.frame(time = c(firstZeroA1, firstZeroA0),
                       A = rep(c(1, 0), each=Nsim))
  survdf$event <- as.numeric(survdf$time!= (maxT+1)) # Those with maxT+1 survived to the end of follow-up
  survdf$time[survdf$time > maxT] <- maxT # Reset those who survived so their maxT is end of follow-up
  baseA <- exp(coef(coxph(Surv(time, event) ~ A, data=survdf)))
  return(baseA)
}

# Code to do MC simulation to calculate marginal hazard ratios + bootstrap CIs
set.seed(1)
B <- 501
maxT <- 72
seeds <- sample.int(100000000, size=B)

# Set up treatment, time, and interaction variables
t <- 0:maxT
dfSetA0 <- data.frame(t=t, t2 = t^2, A=0)
dfSetA1 <- data.frame(t=t, t2 = t^2, A=1)

# Run MC simulation
margHR <- MCsimHR(i=1, Nsim = 100000, baseline, maxT=72, coeffit=basecoef, dfSetA0, dfSetA1)
tvmargHR <- MCsimHR(i=1, Nsim = 100000, baseline, maxT=72, coeffit=tvcoef, dfSetA0, dfSetA1)

save(margHR, tvmargHR, file=paste0(savedir, "hrests.RData"))

# Code to facilitate bootstrapping for hazard ratios ----------------------

library(doParallel)
library(foreach)

registerDoParallel(cores=4)
getDoParWorkers()

# Set up bootstrap samples
ids <- unique(df$matchid)
nind <- length(ids)
bootids1 <- matrix(ids[sample.int(nind, nind*B, replace=T)], nrow=nind, ncol=B)

bootids2 <- apply(bootids1, 2, FUN=function(b){
  temp <- data.frame(matchid = b)
  return(plyr::count(temp, vars = 'matchid'))
}) # This returns a matrix with 2 variables: the ID for each unique individual ('matchid'), and the frequency they appeared in a particular bootstrap sample ('freq')


bootres1 <- foreach(i=2:B, .combine=rbind) %dopar% {
  # Load libraries
  require(reshape2)
  require(plyr)
  require(dplyr)
  require(rms)
  
  # Add "freq" variable to bootdf
  bootdf <- left_join(df, bootids2[[i]], by = "matchid")
  bootdf$freq[is.na(bootdf$freq)] <- 0
  
  # Order data by newid to sort by clones instead of IDs
  bootdf <- bootdf[order(bootdf$matchid, bootdf$cloneid, bootdf$t),]
  
  # Create weights
  bootdf <- create.weights.boot(bootdf, wtformN, wtformD, q=0.99)
  
  # Get model coefficients - WEIGHTED
  wcoef <- coef(glm(hrform, 
                    data = bootdf[bootdf$touse==1,], 
                    weights = stabw, family=quasibinomial()) )
  
  # Get model coefficients - UNWEIGHTED
  bcoef <- coef(glm(hrform, 
                    data = bootdf[bootdf$touse==1,], 
                    family=quasibinomial()))
  
  # Create baseline dataset from bootstrap sample
  baseline <- bootdf[bootdf$t==0 & bootdf$touseweight==1, c('stage_four', 'agedx', 'ervisit_b', 'cci_b', 'perf_b', 'cholang_b', 'pneumon_b')]
  baseline$agedx2 <- baseline$agedx*baseline$agedx
  baseline <- baseline[,c('stage_four', 'agedx', 'agedx2',  'ervisit_b', 'cci_b', 'perf_b', 'cholang_b', 'pneumon_b')]
  
  # Run MC simulation in bootstrapped data
  margHRb <- MCsimHR(i=i, Nsim = 100000, baseline, maxT=72, coeffit=bcoef, dfSetA0, dfSetA1)
  tvmargHRb <- MCsimHR(i=i, Nsim = 100000, baseline, maxT=72, coeffit=wcoef, dfSetA0, dfSetA1)
  
  # Return data frame
  return(data.frame(baseA = margHRb,
                    tvA = tvmargHRb))
}

save(bootres1, file=paste0(savedir, 'hrbootstraps.RData'))
