# function for fitting all time-series data
fitModel <- function(input.data, n.iter=30000) {
y <- input.data$ratio
z <- cbind(rep(1,length(y)), input.data$min_temp.C_avg, input.data$precip.mm_avg)
colnames(z) <- c("betaIntercept", "betaTmin", "betaPrecip")
data <- list(OBS=log(y),n=length(y), x_ic = 0,tau_ic = 0.00001,a_obs=0.1,
               r_obs=0.1,a_add=0.1,r_add=0.1)
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
                         data = data,
                         # inits = init,
                         n.chains = 3)

# jags.out   <- coda.samples (model = j.model,
#                             variable.names = c("tau_add","tau_obs", "beta_IC", "betaIntercept", "betaTmin", "betaPrecip"),
#                             n.iter = n.iter)
# view trace plots
#plot(jags.out)

# check psrf scores; we want these to be under 1.1
if (any(gelman.diag(jags.out)[[1]] > 1.1)) { 
  cat("PSRF values above 1.1. Re-assess model inputs and convergence.")
  } else cat("No PSRF values are above 1.1. Model is sufficiently converged.")

# plot correlation matrix of all of our parameters
out <- as.matrix(jags.out)
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add",
                                               "tau_obs", "beta_IC", "betaIntercept", 
                                               "betaTmin", "betaPrecip"),
                            n.iter = n.iter)


mc3.out <- jags.out
out = list(params=NULL,predict=NULL,model=ourmodel,data=data)
mfit = as.matrix(mc3.out,chains=TRUE)
pred.cols = union(grep("x[",colnames(mfit),fixed=TRUE),grep("mu[",colnames(mfit),fixed=TRUE))
chain.col = which(colnames(mfit)=="CHAIN")
out$predict = mat2mcmc.list(mfit[,c(chain.col,pred.cols)])
out$params   = mat2mcmc.list(mfit[,-pred.cols])
return(out)
}