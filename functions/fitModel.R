# function for fitting all time-series data

input.data <- STER
n.iter = 30000
fitModel <- function(input.data, n.iter = 30000, climate = "era5") {
  y <- input.data$ratio
  
  if (climate == "daymet"){
  z <-
    cbind(rep(1, length(y)), # use daymet data
          input.data$precip.mm_avg)
  colnames(z) <- c("betaIntercept", 
                   "betaPrecip")
  } else {
  z <-
    cbind(rep(1, length(y)), # use means from ERA5 data
          input.data$precip.mean)
  colnames(z) <- c("betaIntercept", 
                   "betaPrecip")
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
  
  # jags.out   <- coda.samples (model = j.model,
  #                             variable.names = c("tau_add","tau_obs", "beta_IC", "betaIntercept", "betaPrecip"),
  #                             n.iter = n.iter)
  # view trace plots
  #plot(jags.out)

  
  # take samples of all of our parameters
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
  
  # check psrf scores; these should be under 1.1
  if (any(gelman.diag(jags.out)[[1]] > 1.1)) {
    cat("PSRF values above 1.1. Re-assess model inputs and convergence.")
  } else
    cat("No PSRF values are above 1.1. Model is sufficiently converged.")
  
  mc3.out <- jags.out
  out = list(
    params = NULL,
    predict = NULL,
    model = ourmodel,
    data = data
  )
  mfit = as.matrix(mc3.out, chains = TRUE)
  pred.cols = union(grep("x[", colnames(mfit), fixed = TRUE),
                    grep("mu[", colnames(mfit), fixed = TRUE))
  chain.col = which(colnames(mfit) == "CHAIN")
  out$predict = mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
  out$params   = mat2mcmc.list(mfit[, -pred.cols])
  return(out)
}
