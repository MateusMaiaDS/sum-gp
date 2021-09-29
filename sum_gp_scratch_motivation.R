# Loading all the files
rm(list = ls())
library(tidyverse)

# Set the gp-sum as working directory


# Verify the working directories
source("rbase_files/sum_gp_scratch.R")
source("fast_gp_multiple_tau.R")

# Generating the synthetic data
# Some R code to simulate data
N <- 25 # can take to N = 100 but fitting gets really slow ...
mu <- 0 # Mu value
tau <- true_tau <- 10 # Tau^-1 results the nugget term
phi <- 0.1 # Phi main value
nu <- true_nu <- 0.1 # Nu main value
n_iter <- 1000 # Number of MCMC samples
K <- 10 # Number of trees
# Setting other parameters
prior_factor <- 0.01 # This prior factor here is used to control the variance over \tau
a_tau <- (K^3)*true_tau * prior_factor
d_tau <- prior_factor
nu_vector <- rep(nu, K)
phi_vector <- rep(phi, K)

# Setting c
c <- K

# Setting the predictions
predictions <- matrix(0,
  nrow = N,
  ncol = K
)

# Getting the data
# Simulating from a sin function
set.seed(123)
x <- seq(from = -pi, to = pi, length.out = N) %>% as.matrix()
y <- sqrt((true_tau^-1) * (true_nu^-1)) * sin(x) + rnorm(n = N, sd = sqrt(true_tau^-1)) %>% as.matrix()

# plot(x,y,pch=20)

# Now create predictions
N_new <- 100
x_grid <- seq(-pi, pi, length = N_new) %>% as.matrix()
colnames(x) <- colnames(x_grid) <- "x"


# Generating the model
gp_sum_mod <- gp_sum(
  x = x,
  y = y,
  K = K,
  tau = tau,
  mu = mu,
  nu_vector = rep(0.1, K),
  phi_vector = rep(0.1, K),
  predictions = predictions,
  n_iter = n_iter,
  a_tau = a_tau,
  d_tau = d_tau,
  c = c
)


# Gathering all final results now
# I'm just going to use the mean values of the parameters
tau <- gp_sum_mod$tau_store %>% apply(2, mean)
mu <- gp_sum_mod$mu_store %>% apply(2, mean)
h <- gp_sum_mod$prediction_store %>% apply(c(2, 3), mean)
nu <- gp_sum_mod$nu_store %>% apply(2, mean)
phi <- gp_sum_mod$phi_store %>% apply(2, mean)


# Now create predictions
N_new <- 100
x_new <- seq(-pi, pi, length = N_new)


# Now create predictions
N_new <- 100
x_new <- seq(-pi, pi, length = N_new)
JAGS_Sigma_new <- array(NA, dim = c(K, N, N_new))
JAGS_Sigma_star <- JAGS_cov_h_new <- array(NA, dim = c(K, N_new, N_new))
JAGS_Sigma <- array(NA, dim = c(K, N, N))
JAGS_h_new <- matrix(NA, ncol = K, nrow = N_new)
JAGS_cov_h <- array(NA, dim = c(K, N, N))

for (k in 1:K) {
  JAGS_Sigma_new[k, , ] <- c((nu[k]^-1) * tau[k]^-1) * exp(-1 / (2 * phi[k]^2) * outer(c(x), x_new, "-")^2)
  JAGS_Sigma_star[k, , ] <- c((nu[k]^-1) * tau[k]^-1) * exp(-1 / (2 * phi[k]^2) * outer(x_new, x_new, "-")^2) + diag(1 / tau[k], N_new)
  JAGS_Sigma[k, , ] <- c((nu[k]^-1) * (tau[k]^-1)) * exp(-1 / (2 * phi[k]^2) * outer(c(x), c(x), "-")^2) + diag(1 / tau[k], N)

  JAGS_h_new[, k] <- t(JAGS_Sigma_new[k, , ]) %*% solve(JAGS_Sigma[k, , ], h[, k] - mu[k]) + mu[k]
  JAGS_cov_h_new[k, , ] <- JAGS_Sigma_star[k, , ] - t(JAGS_Sigma_new[k, , ]) %*% solve(JAGS_Sigma[k, , ], JAGS_Sigma_new[k, , ])
}


JAGS_pred_mean <- rowSums(JAGS_h_new)
JAGS_pred_sd <- sqrt(rowSums(apply(JAGS_cov_h_new, 1, diag)))
JAGS_pred_sd_train <- sqrt(rowSums(apply(JAGS_Sigma, 1, diag)))

# Plotting each tree dataset.
plot_each_tree <- JAGS_h_new %>%
  as.data.frame() %>%
  add_column(x = x_new)
colnames(plot_each_tree) <- c(paste0("tree_", 1:K), "x")

plot_each_tree <- plot_each_tree %>%
  pivot_longer(cols = starts_with("tree"), names_to = "tree") %>%
  rename(pred = value)

# Plot output
# Plotting the data
ggplot() +
  geom_line(
    data = plot_each_tree, mapping = aes(
      x = x,
      y = pred, col = tree
    ), show.legend = FALSE,
    alpha = 0.5
  ) +
  geom_ribbon(
    data = data.frame(
      x = x_new,
      ci_up = JAGS_pred_mean + 1.96 * JAGS_pred_sd,
      ci_low = JAGS_pred_mean - 1.96 * JAGS_pred_sd
    ),
    mapping = aes(
      x = x,
      ymin = ci_up,
      ymax = ci_low
    ), alpha = 0.2, fill = "red"
  ) +
  geom_line(
    data = data.frame(
      x = x_new,
      y = JAGS_pred_mean
    ),
    mapping = aes(
      x = x,
      y = y
    ), col = "red"
  ) +
  geom_point(
    data = data.frame(
      x = x,
      y = y
    ),
    mapping = aes(x = x, y = y),
  ) +
  geom_line(
    data = data.frame(
      x = x_new,
      y = sqrt((true_nu^(-1)) * (true_tau^-1)) * sin(x_new)
    ),
    mapping = aes(x = x, y = y), col = "blue"
  ) +
  ggtitle(paste("GPBART From Scratch")) +
  theme_classic()




# Veryfing all the parameters
print(nu)
print(mu)
print(phi)
print(tau)
