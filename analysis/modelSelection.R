#### MODEL SELECTION -litter, temp, and precip as variables
#input.data <- STER
n.iter = 30000

# input is STER, all time points
  y <- input.data$ratio
  input.data$litterDepth_mean[is.nan(input.data$litterDepth_mean)] <- NA
  if (climate == "daymet"){
    z <-
      cbind(rep(1, length(y)), # use daymet data
            input.data$min_temp.C_avg,
            input.data$precip.mm_avg,
            input.data$litterDepth_mean)
    colnames(z) <- c("betaIntercept", "betaTmin", "betaPrecip","betaLitter")
  } else {
    z <-
      cbind(rep(1, length(y)), # use means from ERA5 data
            input.data$temp.mean,
            input.data$precip.mean,
            input.data$litterDepth_mean)
    colnames(z) <- c("betaIntercept", "betaTmin", "betaPrecip","betaLitter")
  }
  data <-
    list(
      OBS = log(y),
      n = length(y),
      x_ic = 0,
      tau_ic = 0.00001,
      a_obs = 0.1,
      r_obs = 0.1,
      a_add = 0.1,
      r_add = 0.1
    )
  data[["Z"]] <- z
  
  ourmodel <- "model {
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
  
  # #### Random Effects
  # tau_alpha~dgamma(0.1,0.1)
  # for(t in 1:n){
  # alpha[t]~dnorm(0,tau_alpha)
  # }
  
  #### Fixed Effects
  beta_IC~dnorm(0,0.001)
  betaIntercept~dnorm(0,0.001)
  betaTmin~dnorm(0,0.001)
  betaPrecip~dnorm(0,0.001)
  betaLitter~dnorm(0,0.001)
  muTmin~dnorm(0,0.001)
  muPrecip~dnorm(0,0.001)
  muLitter~dnorm(0,0.001)
  tauTmin~dgamma(0.01,0.01)
  tauPrecip~dgamma(0.01,0.01)
  tauLitter~dgamma(0.01,0.01)
  
  #### Data Model
  for(t in 1:n){
  OBS[t] ~ dnorm(x[t],tau_obs)
  Z[t,2] ~ dnorm(muTmin,tauTmin)
  Z[t,3] ~ dnorm(muPrecip,tauPrecip)
  Z[t,4] ~ dnorm(muLitter,tauLitter)
  }
  
  #### Process Model
  for(t in 2:n){
  mu[t] <- beta_IC*x[t-1]  + betaIntercept*Z[t,1] + betaTmin*Z[t,2] + betaPrecip*Z[t,3] + betaLitter*Z[t,4]
  
  x[t]~dnorm(mu[t],tau_add)
  }
}"

  j.model   <- jags.model (file = textConnection(ourmodel),
                           data = data,n.adapt = 2000,
                           n.chains = 3)
  
  jags.out   <- coda.samples (
    model = j.model,
    variable.names = c(
      "x",
      "tau_add",
      "tau_obs",
      "beta_IC",
      "betaIntercept",
      "betaTmin",
      "betaPrecip",
      "betaLitter"
    ),
    n.iter = n.iter
  )
  
  # check psrf scores; we want these to be under 1.1
  if (any(gelman.diag(jags.out)[[1]] > 1.12)) {
    cat("PSRF values above 1.1. Re-assess model inputs and convergence.")
  } else
    cat("No PSRF values are above 1.1. Model is sufficiently converged.")
  
DIC1 <- dic.samples(j.model, 
                    variable.names = c(
                      "x",
                      "tau_add",
                      "tau_obs",
                      "beta_IC",
                      "betaIntercept",
                      "betaTmin",
                      "betaPrecip",
                      "betaLitter"
                    ),
                    n.iter = n.iter)

DIC1







y <- input.data$ratio

if (climate == "daymet"){
  z <-
    cbind(rep(1, length(y)), # use daymet data
          input.data$min_temp.C_avg,
          input.data$precip.mm_avg)
  colnames(z) <- c("betaIntercept", "betaTmin", "betaPrecip")
} else {
  z <-
    cbind(rep(1, length(y)), # use means from ERA5 data
          input.data$temp.mean,
          input.data$precip.mean)
  colnames(z) <- c("betaIntercept", "betaTmin", "betaPrecip")
}
data <-
  list(
    OBS = log(y),
    n = length(y),
    x_ic = 0,
    tau_ic = 0.00001,
    a_obs = 0.1,
    r_obs = 0.1,
    a_add = 0.1,
    r_add = 0.1
  )
data[["Z"]] <- z

ourmodel <- "model {

#### Priors
x[1] ~ dnorm(x_ic,tau_ic)
tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)

# #### Random Effects
# tau_alpha~dgamma(0.1,0.1)
# for(t in 1:n){
# alpha[t]~dnorm(0,tau_alpha)
# }

#### Fixed Effects
beta_IC~dnorm(0,0.001)
betaIntercept~dnorm(0,0.001)
betaTmin~dnorm(0,0.001)
betaPrecip~dnorm(0,0.001)
muTmin~dnorm(0,0.001)
muPrecip~dnorm(0,0.001)
tauTmin~dgamma(0.01,0.01)
tauPrecip~dgamma(0.01,0.01)

#### Data Model
for(t in 1:n){
OBS[t] ~ dnorm(x[t],tau_obs)
Z[t,2] ~ dnorm(muTmin,tauTmin)
Z[t,3] ~ dnorm(muPrecip,tauPrecip)
}

#### Process Model
for(t in 2:n){
mu[t] <- beta_IC*x[t-1]  + betaIntercept*Z[t,1] + betaTmin*Z[t,2] + betaPrecip*Z[t,3]

x[t]~dnorm(mu[t],tau_add)
}
}"

j.model   <- jags.model (file = textConnection(ourmodel),
                         data = data,n.adapt = 2000,
                         n.chains = 3)

jags.out   <- coda.samples (
  model = j.model,
  variable.names = c(
    "x",
    "tau_add",
    "tau_obs",
    "beta_IC",
    "betaIntercept",
    "betaTmin",
    "betaPrecip"
  ),
  n.iter = n.iter
)

# check psrf scores; we want these to be under 1.1
if (any(gelman.diag(jags.out)[[1]] > 1.12)) {
  cat("PSRF values above 1.1. Re-assess model inputs and convergence.")
} else
  cat("No PSRF values are above 1.1. Model is sufficiently converged.")


DIC2 <- dic.samples(j.model, 
                    variable.names = c(
                      "x",
                      "tau_add",
                      "tau_obs",
                      "beta_IC",
                      "betaIntercept",
                      "betaTmin",
                      "betaPrecip"),
                    n.iter = n.iter)

DIC2 # tmin, precip

DIC1 # tmin, precip, litter















y <- input.data$ratio

if (climate == "daymet"){
  z <-
    cbind(rep(1, length(y)), # use daymet data
          input.data$min_temp.C_avg)
  colnames(z) <- c("betaIntercept", "betaTmin")
} else {
  z <-
    cbind(rep(1, length(y)), # use means from ERA5 data
          input.data$temp.mean)
  colnames(z) <- c("betaIntercept", "betaTmin")
}
data <-
  list(
    OBS = log(y),
    n = length(y),
    x_ic = 0,
    tau_ic = 0.00001,
    a_obs = 0.1,
    r_obs = 0.1,
    a_add = 0.1,
    r_add = 0.1
  )
data[["Z"]] <- z

ourmodel <- "model {

#### Priors
x[1] ~ dnorm(x_ic,tau_ic)
tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)

# #### Random Effects
# tau_alpha~dgamma(0.1,0.1)
# for(t in 1:n){
# alpha[t]~dnorm(0,tau_alpha)
# }

#### Fixed Effects
beta_IC~dnorm(0,0.001)
betaIntercept~dnorm(0,0.001)
betaTmin~dnorm(0,0.001)
muTmin~dnorm(0,0.001)
tauTmin~dgamma(0.01,0.01)

#### Data Model
for(t in 1:n){
OBS[t] ~ dnorm(x[t],tau_obs)
Z[t,2] ~ dnorm(muTmin,tauTmin)
}

#### Process Model
for(t in 2:n){
mu[t] <- beta_IC*x[t-1]  + betaIntercept*Z[t,1] + betaTmin*Z[t,2] 

x[t]~dnorm(mu[t],tau_add)
}
}"

j.model   <- jags.model (file = textConnection(ourmodel),
                         data = data,n.adapt = 2000,
                         n.chains = 3)

jags.out   <- coda.samples (
  model = j.model,
  variable.names = c(
    "x",
    "tau_add",
    "tau_obs",
    "beta_IC",
    "betaIntercept",
    "betaTmin"),
  n.iter = n.iter
)

# check psrf scores; we want these to be under 1.1
if (any(gelman.diag(jags.out)[[1]] > 1.12)) {
  cat("PSRF values above 1.1. Re-assess model inputs and convergence.")
} else
  cat("No PSRF values are above 1.1. Model is sufficiently converged.")


DIC3 <- dic.samples(j.model, 
                    variable.names = c(
                      "x",
                      "tau_add",
                      "tau_obs",
                      "beta_IC",
                      "betaIntercept",
                      "betaTmin"),
                    n.iter = n.iter)

DIC3 ## tmin

DIC2 # tmin, precip

DIC1 # tmin, precip, litter













y <- input.data$ratio

if (climate == "daymet"){
  z <-
    cbind(rep(1, length(y)), # use daymet data
          input.data$precip.mm_avg)
  colnames(z) <- c("betaIntercept", "betaPrecip")
} else {
  z <-
    cbind(rep(1, length(y)), # use means from ERA5 data
          input.data$precip.mean)
  colnames(z) <- c("betaIntercept", "betaPrecip")
}
data <-
  list(
    OBS = log(y),
    n = length(y),
    x_ic = 0,
    tau_ic = 0.00001,
    a_obs = 0.1,
    r_obs = 0.1,
    a_add = 0.1,
    r_add = 0.1
  )
data[["Z"]] <- z

ourmodel <- "model {

#### Priors
x[1] ~ dnorm(x_ic,tau_ic)
tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)

# #### Random Effects
# tau_alpha~dgamma(0.1,0.1)
# for(t in 1:n){
# alpha[t]~dnorm(0,tau_alpha)
# }

#### Fixed Effects
beta_IC~dnorm(0,0.001)
betaIntercept~dnorm(0,0.001)
betaPrecip~dnorm(0,0.001)
muPrecip~dnorm(0,0.001)
tauPrecip~dgamma(0.01,0.01)

#### Data Model
for(t in 1:n){
OBS[t] ~ dnorm(x[t],tau_obs)
Z[t,2] ~ dnorm(muPrecip,tauPrecip)
}

#### Process Model
for(t in 2:n){
mu[t] <- beta_IC*x[t-1]  + betaIntercept*Z[t,1] + betaPrecip*Z[t,2]

x[t]~dnorm(mu[t],tau_add)
}
}"

j.model   <- jags.model (file = textConnection(ourmodel),
                         data = data,n.adapt = 2000,
                         n.chains = 3)

jags.out   <- coda.samples (
  model = j.model,
  variable.names = c(
    "x",
    "tau_add",
    "tau_obs",
    "beta_IC",
    "betaIntercept",
    "betaPrecip"
  ),
  n.iter = n.iter
)

# check psrf scores; we want these to be under 1.1
if (any(gelman.diag(jags.out)[[1]] > 1.12)) {
  cat("PSRF values above 1.1. Re-assess model inputs and convergence.")
} else
  cat("No PSRF values are above 1.1. Model is sufficiently converged.")


DIC4 <- dic.samples(j.model, 
                    variable.names = c(
                      "x",
                      "tau_add",
                      "tau_obs",
                      "beta_IC",
                      "betaIntercept",
                      "betaPrecip"),
                    n.iter = n.iter)

# MODEL SELECTION using rjags::dic.samples()

DIC4 ## precipitation
DIC3 ## temperature
DIC2 ## temperature, precipitation
DIC1 ## temperature, precipitation, litterDepth




