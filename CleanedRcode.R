# setting work directory
setwd("C:/Users/Ziming Chen/Documents/Ziming")

# reading dataset
master1 <- read.csv("master3.csv")

# packages
library(rjags)
library(R2jags)
library(extraDistr)
library(dplyr)
library(brms)
library(INLA)
library(brinla)
library(remotes)

# setting priors for INLA

#----half t prior-----------------------------------
#https://becarioprecario.bitbucket.io/inla-gitbook/ch-priors.html#sec:priors
HT.prior = "expression:
sigma = exp(-theta/2);
nu = 3;
log_dens = 0 - 0.5 * log(nu * pi) - (-0.1207822);
log_dens = log_dens - 0.5 * (nu + 1) * log(1 + sigma * sigma);
log_dens = log_dens - log(2) - theta / 2;
return(log_dens);
"

#--------prior list------------------
prior.list = list(
  normal1 = list(prec = list(prior = "normal", param = c(0, 0.1))),
  normal2 = list(prec = list(prior = "normal", param = c(0, 0.1))),
  h.t1 = list(prec = list(prior = HT.prior)),
  h.t2 = list(prec = list(prior = HT.prior)),
  h.t3 = list(prec = list(prior = HT.prior))
) 
h.t = list(prec = list(prior = HT.prior))

############################################
######      Primary Outcome      ###########
######   organ support free days  ##########
############################################

#recoding the organ support free days into 10 categories 
#since INLA can only fit proportional odds model for up to 
#10 categories
master1$OSFD3 <- case_when(master1$OSFD == 1 ~ 1,
                           master1$OSFD == 2 ~ 2,
                           master1$OSFD >= 3 & master1$OSFD <= 12 ~ 3,
                           master1$OSFD >= 13 & master1$OSFD <= 15 ~ 4,
                           master1$OSFD >= 16 & master1$OSFD <= 17 ~ 5,
                           master1$OSFD >= 18 & master1$OSFD <= 19 ~ 6,
                           master1$OSFD == 20 ~ 7,
                           master1$OSFD == 21 ~ 8,
                           master1$OSFD >= 22 & master1$OSFD <= 23 ~ 9,
                           master1$OSFD == 24 ~ 10
)

#-------------------------------------------
#---------------JAGS------------------------
#-------------------------------------------


# making the parameters in JAGS code
Nsites <- length(levels(as.factor(master1$site)))
Nage <- length(levels(as.factor(master1$age)))
Ntime <- length(levels(as.factor(master1$time)))
NOSFD <- length(levels(as.factor(master1$OSFD)))

# making the data structure for JAGS model
jagsdata_1 <- with(master1, list(TxARM = as.numeric(as.factor(TxARM)) - 1,
                                   site = site,
                                   age = age,
                                   time = time,
                                   gender = gender,
                                   N = nrow(master1), 
                                   Nsites = Nsites,
                                   Nage = Nage,
                                   Ntime = Ntime,
                                   y = OSFD3))

lm_jags_1 <- function(){
  # Likelihood:
  for (i in 1:N){
    
    ## Calculate categorical likelihood of the outcome:
    y[i] ~ dcat(pr[i,1:10]) #10 is the number of response catergories (OSFD -1 ~ 22 recoded)
    
    ## The odds is fixed between observation categories
    ## This is the same as any linear regressions except that the intercept isn't included here:
    odds[i] <- beta[1]*gender[i]+beta[2]*TxARM[i]+randomeffect1[site[i]]+randomeffect2[age[i]]+randomeffect3[time[i]]
    
    ## forming the linear predictor using a logit link
    for(n in 1:9){
      logit(g[i,n]) <- alpha[n] - odds[i]
    }
    g[i,10] <- 1
    
    ## Then calculate the non-cumulative probabilities
    pr[i,1] <- g[i,1]
    for(k in 2:9){
      pr[i,k] <- g[i,k] - g[i,k-1]
    }
    pr[i,10] <- 1 - g[i,9]
  }
  
  # Priors:
  ## Calculate the Categories-1 independent intercepts
  for(q in 1:9){
    alpha[q] ~ dnorm(0, 0.1)#0.1 precision
  }
  for (h in 1:2){
    beta[h] ~ dnorm(0, 0.1)#0.1 precision
  }
  
  ## These lines give the prior distributions for the 'random' (hierarchical/nested) parameters to be estimated:
  for(w in 1:3){
    sigma[w] ~ dt(0, 2.5, 3)%_%T(0,)# half t distribution
  }
  tau1 <- 1 / (sigma[1] * sigma[1])# convert to precision
  tau2 <- 1 / (sigma[2] * sigma[2])# convert to precision
  tau3 <- 1 / (sigma[3] * sigma[3])# convert to precision
  for(m in 1:Nsites){
    randomeffect1[m] ~ dnorm(0, tau1)
  }#tau is precision
  for(k in 1:Nage){
    randomeffect2[k] ~ dnorm(0, tau2)
  }#tau is precision  
  for(p in 1:Ntime){
    randomeffect3[p] ~ dnorm(0, tau3)
  }#tau is precision  
}

#initial values for mcmc runs
init_values_1 <- function(){
  list(alpha = sort(rnorm(9)), sigma = rht(3,3,2.5), beta = rnorm(2))
}

#specify the parameters of interest
params_1 <- c("alpha", "beta", "sigma")

#fitting the JAGS runs
fit_lm_1 <- jags(data = jagsdata_1, inits = init_values_1_2, parameters.to.save = params_1_2, model.file = lm_jags_1_2,
                   n.chains = 4, n.iter = 12000, n.burnin = 2000, n.thin = 20, DIC = F)
fit_lm_1


#-------------------------------------------
#---------------STAN------------------------
#-------------------------------------------
prior1_1 = get_prior(OSFD3 ~ gender + TxARM + (1|age) + (1|site) + (1|time), data = master1, family=cumulative("logit"))
prior1_1
brmmd1_1 <- brm(OSFD3 ~ gender + TxARM + (1|age) + (1|site) + (1|time), data = master1, family=cumulative("logit"), 
                prior=prior1_1, cores = getOption("mc.cores", 4),
                chains = 4, iter = 12000, warmup = 2000, thin = 20,control = list(adapt_delta = 0.99, max_treedepth = 12))
summary(brmmd1_1)

#-------------------------------------------
#---------------INLA------------------------
#-------------------------------------------
data.inla <- inla(OSFD3 ~ gender + TxARM + f(age,model = "iid", hyper = HT.prior) + f(site,model = "iid", hyper = HT.prior) + f(time,model = "iid", hyper = HT.prior), family='pom',
                   data = master1, 
                   control.fixed = list(mean = list(gender = 0, TxARM = 0), prec = 0.1),
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config = TRUE))
summary(data.inla)



###################################################
######      Secondary Outcome 1     ###############
######   survival without organ support  ##########
###################################################

#-------------------------------------------
#---------------JAGS------------------------
#-------------------------------------------
jagsdata2 <- with(master1, list(TxARM = as.numeric(as.factor(TxARM)) - 1,
                                site = site,
                                age = age,
                                time = time,
                                gender = gender,
                                N = nrow(master1), 
                                Nsites = Nsites,
                                Nage = Nage,
                                Ntime = Ntime,
                                y = survival_with_no_organ_support))
lm_jags2 <- function(){
  # Likelihood:
  for (i in 1:N){
    
    ## Calculate categorical likelihood of the outcome:
    y[i] ~ dbern(p[i])
    #Note that in order to prevent arithmetic overflows (particularly with the clog-log model,
    #I am going to constrain the estimated linear predictor to between -15 and 15.
    logit(p[i]) <- max(-15, min(15,-(alpha + beta[1]*gender[i]+beta[2]*TxARM[i]+randomeffect1[site[i]]+randomeffect2[age[i]]+randomeffect3[time[i]])))
  }
  # Priors:
  for (h in 1:2){
    beta[h] ~ dnorm(0, 0.1)#0.1 precision
  }
  alpha ~ dnorm(0, 0.1)
  ## These lines give the prior distributions for the 'random' (hierarchical/nested) parameters to be estimated:
  for(w in 1:3){
    sigma[w] ~ dt(0, 2.5, 3)%_%T(0,)# standard deviation of random effect (variance between sites) -- adjust the scale -- maybe unif(0,10)
  }
  tau1 <- 1 / (sigma[1] * sigma[1])# convert to precision
  tau2 <- 1 / (sigma[2] * sigma[2])# convert to precision
  tau3 <- 1 / (sigma[3] * sigma[3])# convert to precision
  for(m in 1:Nsites){
    randomeffect1[m] ~ dnorm(0, tau1)
  }#tau is precision
  for(k in 1:Nage){
    randomeffect2[k] ~ dnorm(0, tau2)
  }#tau is precision  
  for(p in 1:Ntime){
    randomeffect3[p] ~ dnorm(0, tau3)
  }#tau is precision  
}

#initial values for mcmc runs
init_values2 <- function(){
  list(alpha = rnorm(1), sigma = rht(3,3,2.5), beta = rnorm(2))
}

# parameters of interest
params2 <- c("alpha","beta", "sigma")

# fit the JAGS model
fit_lm2 <- jags(data = jagsdata2, inits = init_values2, parameters.to.save = params2, model.file = lm_jags2,
                n.chains = 4, n.iter = 12000, n.burnin = 2000, n.thin = 20, DIC = F)
fit_lm2


#-------------------------------------------
#---------------STAN------------------------
#-------------------------------------------
#getting the priors
prior2 = get_prior(survival_with_no_organ_support ~ gender + TxARM + (1|age) + (1|site) + (1|time), data = master1, family=bernoulli("logit"))
prior2

# fitting the stan model
brmmd2 <- brm(survival_with_no_organ_support ~ gender + TxARM + (1|age) + (1|site) + (1|time), data = master1, family=bernoulli("logit"),
              prior = prior2, cores = getOption("mc.cores", 4),
              chains = 4, iter = 12000, warmup = 2000, thin = 20,control = list(adapt_delta = 0.99, max_treedepth = 12))
summary(brmmd2)



#-------------------------------------------
#---------------INLA------------------------
#-------------------------------------------
dat.inla1 <- inla(survival_with_no_organ_support ~ gender + TxARM + f(age,model = "iid", hyper = h.t) + f(site,model = "iid", hyper = h.t) + f(time,model = "iid", hyper = h.t), family='binomial',
                  data = master1,
                  control.fixed = list(mean = list(gender = 0, TxARM = 0), prec = 0.1),
                  control.family=list(link='logit'),
                  control.predictor=list(link=1, compute=TRUE),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config = TRUE))
summary(dat.inla1)



###################################################
######      Secondary Outcome 2     ###############
#######   length of hospital stay  ################
###################################################

#-------------------------------------------
#---------------JAGS------------------------
#-------------------------------------------
censored <- master1$censor == -1
is.censored <- censored*1
t.to.death <- master1$length_of_hospital_stay
t.to.death[censored] <- NA
t.to.death
length(unique(t.to.death))#45
t <- c(na.omit(unique(t.to.death)), 90)
#censored length of hospital stay
t.cen <- rep(0, times=length(censored))
t.cen[censored] <- master1$length_of_hospital_stay[censored] 
t.cen
eps = 0.000001; # used to guard against numerical imprecision in step function
master1$censor2 <- ifelse(master1$censor == 1, 1, 0)
#recoding censoring status -- 0 for censored and 1 for uncensored
jagsdata3_1 <- with(master1, list(TxARM = as.numeric(as.factor(TxARM)) - 1,
                                  site = site,
                                  age = age,
                                  time = time,
                                  gender = gender,
                                  N = nrow(master1), # number of patients
                                  Nsites = Nsites,
                                  Nage = Nage,
                                  Ntime = Ntime,
                                  obs.t = length_of_hospital_stay, # observed failure or censoring time for each patient
                                  T_ = length(unique(t.to.death))-1, #number of unique failure times
                                  fail = censor2,
                                  t = sort(t),
                                  eps = 0.000001))# used to guard against numerical imprecision in step function#
#step function in BUGS equivalent in R
step <- function(x){
  if(x >= 0){
    return(1)
  }
  else{
    return(0)
  }
}
Y <- matrix(0, nrow = jagsdata3_1$N, ncol = jagsdata3_1$T_)
dN <- matrix(0, nrow = jagsdata3_1$N, ncol = jagsdata3_1$T_)
for(i in 1:jagsdata3_1$N){
  for(j in 1:jagsdata3_1$T_){
    #at risk Y[i,j] = 1 if obs.t[i] >= t[j]
    Y[i,j] <- step(jagsdata3_1$obs.t[i] - jagsdata3_1$t[j] + jagsdata3_1$eps)
    #counting process
    dN[i,j] <- Y[i,j]*step(jagsdata3_1$t[j+1]-jagsdata3_1$obs.t[i]-jagsdata3_1$eps)*jagsdata3_1$fail[i]
  }
}
#----------jags survival official ---------------
jagsdata3 <- with(master1, list(TxARM = as.numeric(as.factor(TxARM)) - 1,
                                site = site,
                                age = age,
                                time = time,
                                gender = gender,
                                N = nrow(master1), # number of patients
                                Nsites = Nsites,
                                Nage = Nage,
                                Ntime = Ntime,
                                T_ = length(unique(t.to.death))-1, #number of unique failure times
                                Y = Y,
                                t = sort(t),
                                dN = dN))
lm_jags3 <- function(){
  for(j in 1:T_){
    for(i in 1:N){
      #poisson likelihood trick
      dN[i,j] ~ dpois(Idt[i,j])
      #intensity
      Idt[i,j] <- Y[i,j]*exp(beta[1]*gender[i]+beta[2]*TxARM[i]+randomeffect1[site[i]]+randomeffect2[age[i]]+randomeffect3[time[i]])*dL0[j]
    }
    dL0[j] ~ dgamma(mu[j], c)
    mu[j] <- dL0.star[j] * c
  }
  c <- 0.001 
  r <- 0.1
  for (j in 1 : T_) {
    dL0.star[j] <- r * (t[j+1]-t[j])
  }
  for (h in 1:2){
    beta[h] ~ dnorm(0, 0.1)#0.1 variance correction
  }
  ## These lines give the prior distributions for the 'random' (hierarchical/nested) parameters to be estimated:
  for(w in 1:3){
    sigma[w] ~ dt(0, 2.5, 3)%_%T(0,)# standard deviation of random effect (variance between sites) -- adjust the scale -- maybe unif(0,10)
  }
  tau1 <- 1 / (sigma[1] * sigma[1])# convert to precision
  tau2 <- 1 / (sigma[2] * sigma[2])# convert to precision
  tau3 <- 1 / (sigma[3] * sigma[3])# convert to precision
  for(m in 1:Nsites){
    randomeffect1[m] ~ dnorm(0, tau1)
  }#tau is precision
  for(k in 1:Nage){
    randomeffect2[k] ~ dnorm(0, tau2)
  }#tau is precision  
  for(p in 1:Ntime){
    randomeffect3[p] ~ dnorm(0, tau3)
  }
}

# setting initial values
init_values3 <- function(){
  list(sigma = rht(3,3,2.5), beta = rnorm(2))
}

# specifying parameter of interest
params3 <- c("beta", "sigma")

# fitting the model
fit_lm3 <- jags(data = jagsdata3, inits = init_values3, parameters.to.save = params3, model.file = lm_jags3,
                n.chains = 4, n.iter = 12000, n.burnin = 2000, n.thin = 20, DIC = F)
fit_lm3

#-------------------------------------------
#---------------STAN------------------------
#-------------------------------------------
#With the exception of categorical, ordinal, and mixture families, left, right, and interval censoring can be modeled through y | cens(censored) ~ predictors. 
#The censoring variable (named censored in this example) should contain the values 'left', 'none', 'right', and 'interval' 
#(or equivalently -1, 0, 1, and 2) to indicate that the corresponding observation is left censored, not censored, right censored, or interval censored.
master2 <- data.frame(X = master1$X,
                      PATID = master1$PATID,
                      TxARM = master1$TxARM,
                      site = as.numeric(as.factor(master1$site)),
                      age = as.numeric(master1$age)+1, # recoding to be greater than 1 to fit stan
                      time = as.numeric(master1$time),
                      OSFD = master1$OSFD + 2,#this is to make sure when fitting a logistic regression(ordinal) the values can be greater than 1
                      gender = master1$gender,
                      platform = master1$platform,
                      length_of_hospital_stay = master1$length_of_hospital_stay + 1,
                      censor = master1$censor)
master2$censor2 <- ifelse(master2$censor == 1, 0, 1)
prior7 = get_prior(length_of_hospital_stay | cens(censor2) ~ gender + TxARM + (1|age) + (1|site) + (1|time), data = master2, family=cox(link = "log", bhaz = NULL))
prior7
brmmd7 <- brm(length_of_hospital_stay | cens(censor2) ~ gender + TxARM + (1|age) + (1|site) + (1|time), data = master2, family=cox(link = "log", bhaz = NULL),
              prior = prior7, cores = getOption("mc.cores", 4),
              chains = 4, iter = 12000, warmup = 2000, thin = 20,control = list(adapt_delta = 0.99, max_treedepth = 12))
summary(brmmd7)

#-------------------------------------------
#---------------INLA------------------------
#-------------------------------------------
library(survival)
master2 <- data.frame(X = master1$X,
                      PATID = master1$PATID,
                      TxARM = master1$TxARM,
                      site = as.numeric(as.factor(master1$site)),
                      age = as.numeric(master1$age)+1, # recoding to be greater than 1 to fit stan
                      time = as.numeric(master1$time),
                      OSFD = master1$OSFD + 2,#this is to make sure when fitting a logistic regression(ordinal) the values can be greater than 1
                      gender = master1$gender,
                      platform = master1$platform,
                      length_of_hospital_stay = master1$length_of_hospital_stay + 1,
                      censor = master1$censor)

master2$censor3 <- ifelse(master2$censor == 1, 1, 0)
#recoding censoring status -- 0 for censored and 1 for uncensored

surv.master <- inla.surv(master2$length_of_hospital_stay, master2$censor3)
#Note how censored observations are marked with a + to indicate that the time to event 
#has been censored and it is larger than the value shown.
cox.master <- inla(surv.master ~ gender + TxARM + f(age,model = "iid", hyper = h.t) + f(site,model = "iid", hyper = h.t) + f(time,model = "iid", hyper = h.t), data = master2, 
                   family = "coxph",
                   control.fixed = list(mean = list(gender = 0, TxARM = 0), prec = 0.1),
                   control.hazard = list(hyper = list(prec = list(param = c(0.001, 0.001)))),
                   control.predictor=list(link=1, compute=TRUE),
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config = TRUE))
summary(cox.master)












