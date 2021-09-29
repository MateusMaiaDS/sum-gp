# Function to calculate the sum of guassian processes from the scratch
gp_sum <- function(
         x,
         y,
         K,
         tau,
         mu,
         nu_vector,
         phi_vector,
         predictions,
         n_iter,
         a_tau,
         d_tau,
         c
         ){
  
  # Calculating the residuals
  
  # Creating the stores for each parameter
  nu_store <- matrix(NA, nrow = n_iter, ncol = K)
  phi_store <- matrix(NA, nrow = n_iter, ncol = K)
  mu_store <- matrix(NA, nrow = n_iter, ncol = K)
  tau_store <- matrix(NA, nrow = n_iter, ncol = K)
  
  
  # Auxiliar vectors 
  mu_vector <- rep(mu,K)
  tau_vector <- rep(tau,K)
  
  # Saving the prediction store
  prediction_store <- array(0, dim = c(n_iter,
                                        length(y),
                                        K))
  # Setting the progress bar
  progress_bar <- utils::txtProgressBar(
    min = 1, max = n_iter,
    style = 3, width = 50,
    label = "Running GP-Sum-Sampler..."
  )
  
  
  # Initializing over the number of iterations
  for(i in 1:n_iter){
    
    # Counting ....
    utils::setTxtProgressBar(progress_bar, i)
    
    
    # Iterating over the number of trees
    for(j in 1:K){
      
      # Calculating the residuals
      current_partial_residuals <- y - rowSums(predictions[,-j,drop = FALSE])
      
      # Calculating all the Omegas, and Omegas + I to not be always repeating it
      
      # Calculating the Omega 
      omega <- kernel_function(distance_matrix = symm_distance_matrix(m1 = x),nu = nu,phi = phi)
      
      # Setting the n
      n <- dim(omega)[1]
      
      # Calculating S
      S <- (sum(solve(omega+diag(n))))+1/c
      
      # Caculating R^{-1}%*%(I+\Omega)^{-1}%*%R
      R_omega_I_R <- crossprod(current_partial_residuals,
                               solve((diag(n)+omega),current_partial_residuals))
      
      
      # Caculating R^{-1}%*%(I+\Omega)^{-1}%*%1
      R_omega_I_one <- rowSums(crossprod(current_partial_residuals,solve(diag(n)+omega)))
      
      
      # Update tau values
      tau_vector[j] <- update_tau(current_partial_residuals = current_partial_residuals,
                                  x = x,
                                  y = y,
                                  a_tau = a_tau,
                                  d_tau = d_tau,
                                  phi = phi_vector[j],
                                  nu = nu_vector[j],
                                  S = S,
                                  R_omega_I_R = R_omega_I_R,
                                  R_omega_I_one = R_omega_I_one)
      
      # Update mu values 
      mu_vector[j] <- update_mu(current_partial_residuals = current_partial_residuals,
                                x = x,
                                y = y,
                                tau = tau_vector[j],
                                nu = nu_vector[j],
                                phi = phi_vector[j],
                                S = S,
                                R_omega_I_one = R_omega_I_one)
      
      # Updating the values for the h
      predictions[,j] <- update_h(current_partial_residuals = current_partial_residuals,
                                  x = x,
                                  y = y,
                                  tau = tau_vector[j],
                                  mu = mu_vector[j],
                                  nu = nu_vector[j],
                                  phi = phi_vector[j],
                                  omega = omega)
      
      # Updating \nu and \phi
      phi_vector[j] <- update_phi(current_partial_residuals = current_partial_residuals,
                               x = x,
                               y = y,
                               tau  = tau_vector[j],
                               mu  = mu_vector[j],
                               nu = nu_vector[j],
                               phi = phi_vector[j],
                               omega)
        
    }# Finishing the iteration over the trees
     
    # Saving the results 
    nu_store[i,] <- nu_vector
    tau_store[i, ] <- tau_vector
    mu_store[i,] <- mu_vector
    phi_store[i,] <- phi_vector
    prediction_store[i, , ] <- predictions
    
  }
  
  # Saving and restoring all the results
  return(list(nu_store = nu_store,
              tau_store = tau_store,
              mu_store = mu_store,
              phi_store = phi_store,
              prediction_store = prediction_store))
    
     
}


# Function to update_tau 

update_tau <- function(
           current_partial_residuals,
           x,
           y,
           a_tau,
           d_tau,
           phi,
           nu,
           S,
           R_omega_I_one,
           R_omega_I_R
           ){

  
  # Calculating the n
  n <- length(y)

  # Calculating the shape
  tau_shape <- 0.5 * n + a_tau
  tau_rate <- 0.5*(R_omega_I_R)-
              0.5*((S^-1))*((R_omega_I_one)^2)+d_tau
  
  new_gamma <- rgamma(n = 1,
                      shape = tau_shape,
                      rate = tau_rate)
  
  # Finish the function
  return(new_gamma)
  
}

# Function to update mu
update_mu <- function(current_partial_residuals,
          x,
          y,
          tau,
          nu,
          phi,
          S,
          R_omega_I_one){
  
  # Calculating the mean and sd for mu
  mu_mean <- (S^-1)*R_omega_I_one
  
  # Calculating the sd
  mu_sd <- sqrt(((tau*S)^-1))
  
  # Getting the new sample
  new_mu <- rnorm(n = 1,mean = mu_mean,sd = mu_sd)
  
  return(new_mu)
}


# Update h from the GP
update_h <- function(current_partial_residuals,
                     x,
                     y,
                     nu,
                     phi,
                     tau,
                     mu,
                     omega){
  
  # Calculating omega
  n <- nrow(omega)
  
  # Calculating g_mean 
  g_mean <- mu + crossprod(omega, solve( (diag(n) + omega ), 
                                         (current_partial_residuals-mu)
                                         )
                                       )
  
  # Return the new predicted values
  return(g_mean)
}


# Updating h from the GP
update_phi <- function(current_partial_residuals,
           x ,
           y ,
           tau,
           mu ,
           nu ,
           phi,
           omega){
  
  # Since phi only have a uniform distribution 
  #I only need to calculate the distribution of the posterior
  phi_proposal <- runif(1,min = 0,max = 10)
  
  # Calculating all old quantities 
  omega_old <- omega
  n <- dim(omega)[1]
  
  # Calculating g_mean_old
  g_mean_old <- mu + crossprod(omega_old, solve( (diag(n) + omega_old ), 
                                         (current_partial_residuals-mu)))
  
  g_var_old <- (tau^-1)*(omega_old-crossprod(omega_old, solve( (diag(n) + omega_old ), 
                                                       (omega_old))))
  # Calculating the g_mean_new 
  log_likelihood_old <- mvtnorm::dmvnorm(x = c(current_partial_residuals),mean = g_mean_old,
                                         sigma = g_var_old,log = TRUE)
  
  
  # ========
  # ======== Calculating new quantities
  # ========
  
  # Calculating all new quantities 
  omega_new <- kernel_function(distance_matrix = as.matrix(dist(x)),nu = nu,phi = phi_proposal)

  # Calculating g_mean_new
  g_mean_new <- mu + crossprod(omega_new, solve( (diag(n) + omega_new ), 
                                             (current_partial_residuals-mu)))
  
  g_var_new <- (tau^-1)*(omega_new-crossprod(omega_new, solve( (diag(n) + omega_new ), 
                                                       (omega_new))))
  # Calculating the g_mean_new 
  log_likelihood_new <- mvtnorm::dmvnorm(x = c(current_partial_residuals),mean = g_mean_new,
                                         sigma = g_var_new,log = TRUE)
  
  # Calculating the likelihood from the old phi we would have
  acceptance_phi <- exp(log_likelihood_new-log_likelihood_old)
  
  if(log_likelihood_new == -Inf){
    return(phi)
  }
  
  # Accepting or not
  if (runif(1) < acceptance_phi) { #
    return(phi_proposal)
  } else {
    return(phi)
  }
  
}


# Read one of the trees as example 
rMVN2 <- function(b, Q) {
  p   <- ncol(Q)
  U   <- chol(Q)
  z   <- rnorm(p)
  backsolve(U,z,transpose = FALSE, k=p)+b
}



