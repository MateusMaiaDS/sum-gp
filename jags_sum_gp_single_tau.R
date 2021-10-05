# Header ------------------------------------------------------------------

# A sum of K Gaussian process models in JAGS
# Andrew Parnell

# This file fits a sum of Gaussian Process (GP) regression models to data in JAGS, and produces predictions/forecasts
# It is not meant to converge as the K GPs are non-identifiable

# Some boiler plate code to clear the workspace and load in required packages
#rm(list = ls())
library(R2jags)
library(MASS) # Useful for mvrnorm function
library(dplyr)
library(tidyverse)

# Simulate data -----------------------------------------------------------

# Some R code to simulate data
N <- 100 # can take to N = 100 but fitting gets really slow ...
true_mu <- 0
true_tau <- 10 # Tau^-1 results the nugget term
true_phi <- 0.1
true_nu <- 0.1
# source("gpbart/motivation_simulation_examples.R")


# # Creating a 1-D motivating example
sim_dataset <- motivating_stump_single_tree(n = N,
                                            mu = true_mu,
                                            tau = true_tau, # The variance it should be 12.1 (almost)
                                            nu = true_nu,
                                            phi = true_phi,
                                            seed = 42,seq_boolean = TRUE)



# Getting the x and y

# x <- sim_dataset[,1]
# y <- sim_dataset[,2]

# Setting the seed
set.seed(123)

# Simulating from a sin function
x <- seq(from = -pi,to = pi, length.out =  N)
y <- sqrt((true_tau^-1)*(true_nu^-1))*sin(x)+rnorm(n = N,sd = sqrt(true_tau^-1))

set.seed(NULL)
# Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data
model_code <- "
model
{
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dnorm(sum(h[i,1:K]), tau)
  }

  
  # Gaussian Processes
  for (k in 1:K) {
    Omega[k,1:N,1:N] <- inverse(Sigma[k,1:N,1:N])
    h[1:N,k] ~ dmnorm(m[1:N,k], Omega[k,1:N,1:N])


    for (i in 1:N) {
      m[i,k] <- mu[k]
      Sigma[k,i,i] <- ((nu[k])^(-1)) + 1e-12
      for(j in (i+1):N) {
        Sigma[k,i,j] <- ((nu[k])^(-1)) * exp( - 1/(2*pow(phi[k],2)) * pow(x[i] - x[j], 2) )
        Sigma[k,j,i] <- Sigma[k,i,j]
      }
    }
  }

  # Priors
  for (k in 1:K) {
    mu[k] ~ dnorm(0, tau_mu )
    phi[k] ~ dunif(0,10)
  }
  
  
  tau ~ dgamma(a_tau, d_tau)

}
"

# Set up the data
K <- 5
prior_factor <- 0.001
# sd_hat <- naive_sigma(x = x, y = y)^(-2)
K_bart <- 2
model_data <- list(N = N, y = y,
                   x = x, K = K,
                   a_tau = true_tau*prior_factor,
                   d_tau = prior_factor,
                   tau_mu = (K_bart^2*K)/((max(y)-min(y))^2 ),
                   nu = rep(1,K))


# tau =  K*true_tau)
# , phi = rep(1, 10), nu = rep(1, 10)

# Choose the parameters to watch
model_parameters <- c( "mu", "h","phi", "tau")

# Run the model - can be slow
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code),
  n.chains = 1,n.iter = 2000
) # Number of different starting positions. Only need 1 here as it won't converge
# saveRDS(model_run, file = paste0('model_run_',K,'.rds'))
#model_run <- readRDS(file = 'model_run_10.rds')

# Simulated results -------------------------------------------------------

# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
print(model_run)

# I'm just going to use the mean values of the parameters
tau <- model_run$BUGSoutput$mean$tau
mu <- model_run$BUGSoutput$mean$mu
h <- model_run$BUGSoutput$mean$h
nu <- model_data$nu
phi <- model_run$BUGSoutput$mean$phi #model_run$BUGSoutput$mean$phi

# Now create predictions
N_new <- 100
x_new <- seq(-pi, pi, length = N_new)
JAGS_Sigma_new <- array(NA, dim = c(model_data$K, model_data$N, N_new))
JAGS_Sigma_star <- JAGS_cov_h_new <- array(NA, dim = c(model_data$K, N_new, N_new))
JAGS_Sigma <- array(NA, dim = c(model_data$K, model_data$N, model_data$N))
JAGS_h_new <- matrix(NA, ncol = model_data$K, nrow = N_new)
JAGS_cov_h <- array(NA, dim = c(model_data$K, N, N))

for (k in 1:model_data$K) {
  JAGS_Sigma_new[k, , ] <- c((nu[k]^-1)) * exp(-1 / (2 * phi[k]^2) * outer(x, x_new, "-")^2)
  JAGS_Sigma_star[k, , ] <- c((nu[k]^-1)) * exp(-1 / (2 * phi[k]^2) * outer(x_new, x_new, "-")^2) #+ diag((1/tau), N_new)
  JAGS_Sigma[k, , ] <- c((nu[k]^-1)) *  exp(-1 / (2 * phi[k]^2) * outer(x, x, "-")^2) + diag((1/tau), N)
  
  JAGS_h_new[, k] <- t(JAGS_Sigma_new[k, , ]) %*% solve(JAGS_Sigma[k, , ], h[, k] - mu[k]) + mu[k]
  JAGS_cov_h_new[k, , ] <- JAGS_Sigma_star[k, , ] - t(JAGS_Sigma_new[k, , ]) %*% solve(JAGS_Sigma[k, , ], JAGS_Sigma_new[k, , ])
}


JAGS_pred_mean <- rowSums(JAGS_h_new)
JAGS_pred_sd <- sqrt(1/tau)
JAGS_pred_sd_train <- sqrt(rowSums(apply(JAGS_Sigma, 1, diag)))

# Plotting each tree dataset.
plot_each_tree <- JAGS_h_new %>% as.data.frame() %>% add_column(x = x_new)
colnames(plot_each_tree) <- c(paste0("tree_",1:K),"x")
plot_each_tree <- plot_each_tree %>% pivot_longer(cols = starts_with("tree"),names_to = "tree") %>% rename(pred = value)

# Plot output
# Plotting the data
ggplot()+
  geom_line(data = plot_each_tree, mapping = aes(x = x,
                                                 y = pred, col = tree),show.legend = FALSE,
            alpha = 0.5)+
  geom_ribbon(data = data.frame(x = x_new,
                                ci_up = JAGS_pred_mean + 1.96*JAGS_pred_sd,
                                ci_low = JAGS_pred_mean - 1.96*JAGS_pred_sd),
              mapping = aes(x = x,
                            ymin = ci_up,
                            ymax = ci_low), alpha = 0.2, fill = "red")+
  geom_line(data = data.frame(x = x_new,
                              y = JAGS_pred_mean),
            mapping = aes( x = x,
                           y = y), col = "red")+
  geom_point(data = data.frame(x = x,
                               y = y), 
             mapping = aes(x = x, y = y), )+
  # geom_line(data = data.frame(x = x_new,
  #                             y = sqrt(true_nu^(-1)*(true_tau^-1))*sin(x_new)),
  #           mapping = aes(x = x, y = y), col = "blue")+
  # ylim(c(-2,2))+
  ggtitle(paste("JAGS"))+
  theme_classic()

# Pretty good!


# Printing the posterior mean
print(nu)
print(mu)
print(phi)
print(tau)

