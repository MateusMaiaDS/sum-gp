# Some R code to simulate data
N <- 40 # can take to N = 100 but fitting gets really slow ...
set.seed(123)
mu <- 0
tau <- 10 # Tau^-1 results the nugget term
phi <- 0.1
nu <- 0.1
n_iter <- 1000
K <- 2
prior_factor <- 0.01
a_tau = (K^3)*tau*prior_factor
d_tau = prior_factor
nu_vector <- rep(nu, K)
phi_vector <- rep(phi, K)

c <- K


predictions <- matrix(0,
                      nrow = N,
                      ncol = K)

# Setting other parameters
prior_factor <- 100


# Getting the data
# Simulating from a sin function
x <- seq(from = -pi,to = pi, length.out =  N) %>% as.matrix
# y <- (tau^-1)*((nu^-1)*sin(x)+rnorm(n = N,sd = 1)) %>% as.matrix
y <- sqrt((tau^-1)*(nu^-1))*sin(x)+rnorm(n = N,sd = sqrt(tau^-1)) %>% as.matrix

# Now create predictions
N_new <- 100
x_grid <- seq(-pi, pi, length = N_new) %>% as.matrix()
colnames(x) <- colnames(x_grid) <- "x"



