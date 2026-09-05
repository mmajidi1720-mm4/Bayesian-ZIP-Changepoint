# Load required libraries
library(MASS)
library(coda)
library(ggplot2)


# Generate covariates
generate_covariates <- function(n, q = 3, r = 3) {
  Z <- matrix(0, n, q)
  W <- matrix(0, n, r)
  
  Z[, 1] <- 1  # Intercept
  W[, 1] <- 1  # Intercept
  
  for (j in 2:q) {
    Z[, j] <- rnorm(n)
  }
  
  for (j in 2:r) {
    W[, j] <- rnorm(n)
  }
  
  return(list(Z = Z, W = W))
}

# Generate ZIP data
generate_zip_data <- function(n, m, beta1, beta2, eta1, eta2, Z, W) {
  Y <- numeric(n)
  
  for (t in 1:n) {
    if (t <= m) {
      beta_t <- beta1
      eta_t <- eta1
    } else {
      beta_t <- beta2
      eta_t <- eta2
    }
    
    z_t <- Z[t, ]
    w_t <- W[t, ]
    
    lambda_t <- exp(sum(z_t * beta_t))
    p_t <- 1 / (1 + exp(-sum(w_t * eta_t)))
    
    B_t <- rbinom(1, 1, p_t)
    
    if (B_t == 1) {
      Y[t] <- 0
    } else {
      Y[t] <- rpois(1, lambda_t)
    }
  }
  
  return(Y)
}

# Data augmentation
data_augmentation <- function(Y, beta, eta, Z, W, indices) {
  n_seg <- length(indices)
  V <- numeric(n_seg)
  B <- numeric(n_seg)
  
  for (i in 1:n_seg) {
    t <- indices[i]
    
    z_t <- if(is.matrix(Z)) Z[t, ] else matrix(Z[t, ], nrow=1)
    w_t <- if(is.matrix(W)) W[t, ] else matrix(W[t, ], nrow=1)
    
    lambda_t <- exp(sum(z_t * beta))
    p_t <- 1 / (1 + exp(-sum(w_t * eta)))
    
    if (Y[t] == 0) {
      prob_B1 <- p_t / (p_t + (1 - p_t) * exp(-lambda_t))
      B[i] <- rbinom(1, 1, prob_B1)
      
      if (B[i] == 1) {
        V[i] <- rpois(1, lambda_t)
      } else {
        V[i] <- 0
      }
    } else {
      B[i] <- 0
      V[i] <- Y[t]
    }
  }
  
  return(list(V = V, B = B))
}

# Update beta
update_beta <- function(V, B, eta, Z, W, indices, beta_0, sigma_beta_sq) {
  
  neg_log_post <- function(beta) {
    ll <- 0
    for (i in 1:length(indices)) {
      t <- indices[i]
      z_t <- if(is.matrix(Z)) Z[t, ] else matrix(Z[t, ], nrow=1)
      lambda_t <- exp(sum(z_t * beta))
      ll <- ll + V[i] * sum(z_t * beta) - lambda_t
    }
    
    prior <- -0.5 * sum((beta - beta_0)^2) / sigma_beta_sq
    return(-(ll + prior))
  }
  
  init_beta <- rep(0, length(beta_0))
  opt_result <- tryCatch({
    optim(init_beta, neg_log_post, method = "BFGS", hessian = TRUE)
  }, error = function(e) {
    optim(init_beta, neg_log_post, method = "Nelder-Mead", hessian = FALSE)
  })
  
  beta_max <- opt_result$par
  
  if(!is.null(opt_result$hessian)) {
    Sigma_max <- tryCatch({
      solve(opt_result$hessian)
    }, error = function(e) {
      diag(0.1, length(beta_0))
    })
  } else {
    Sigma_max <- diag(0.1, length(beta_0))
  }
  
  beta_new <- mvrnorm(1, beta_max, Sigma_max)
  
  return(beta_new)
}

# Update eta
update_eta <- function(V, B, beta, Z, W, indices, eta_0, sigma_eta_sq) {
  
  neg_log_post <- function(eta) {
    ll <- 0
    for (i in 1:length(indices)) {
      t <- indices[i]
      w_t <- if(is.matrix(W)) W[t, ] else matrix(W[t, ], nrow=1)
      ll <- ll + B[i] * sum(w_t * eta) - log(1 + exp(sum(w_t * eta)))
    }
    
    prior <- -0.5 * sum((eta - eta_0)^2) / sigma_eta_sq
    return(-(ll + prior))
  }
  
  init_eta <- rep(0, length(eta_0))
  opt_result <- tryCatch({
    optim(init_eta, neg_log_post, method = "BFGS", hessian = TRUE)
  }, error = function(e) {
    optim(init_eta, neg_log_post, method = "Nelder-Mead", hessian = FALSE)
  })
  
  eta_max <- opt_result$par
  
  if(!is.null(opt_result$hessian)) {
    Sigma_max <- tryCatch({
      solve(opt_result$hessian)
    }, error = function(e) {
      diag(0.1, length(eta_0))
    })
  } else {
    Sigma_max <- diag(0.1, length(eta_0))
  }
  
  eta_new <- mvrnorm(1, eta_max, Sigma_max)
  
  return(eta_new)
}

# Fixed update_changepoint with numerical stability
update_changepoint <- function(m, Y, beta1, beta2, eta1, eta2, Z, W, n) {
  
  # Hybrid scheme
  if (runif(1) < 0.5) {
    delta <- sample(c(-1, 1), 1)
    m_prop <- m + delta
  } else {
    m_prop <- sample(2:(n-1), 1)
  }
  
  # Check bounds
  if (m_prop < 2 || m_prop > n - 1) {
    return(m)
  }
  
  # Compute log posterior
  log_post_curr <- compute_log_posterior_m(m, Y, beta1, beta2, eta1, eta2, Z, W, n)
  log_post_prop <- compute_log_posterior_m(m_prop, Y, beta1, beta2, eta1, eta2, Z, W, n)
  
  # Check for NA or Inf
  if (is.na(log_post_curr) || is.na(log_post_prop) || 
      is.infinite(log_post_curr) || is.infinite(log_post_prop)) {
    return(m)
  }
  
  log_alpha <- log_post_prop - log_post_curr
  
  # Additional check
  if (is.na(log_alpha) || is.infinite(log_alpha)) {
    return(m)
  }
  
  if (log(runif(1)) < log_alpha) {
    return(m_prop)
  } else {
    return(m)
  }
}

# Fixed compute_log_posterior_m with numerical stability
compute_log_posterior_m <- function(m, Y, beta1, beta2, eta1, eta2, Z, W, n) {
  
  log_lik <- 0
  
  # Segment 1
  for (t in 1:m) {
    z_t <- Z[t, ]
    w_t <- W[t, ]
    
    lambda_t <- exp(sum(z_t * beta1))
    logit_p <- sum(w_t * eta1)
    p_t <- 1 / (1 + exp(-logit_p))
    
    # Numerical stability
    if (is.infinite(lambda_t) || is.na(lambda_t)) {
      return(-Inf)
    }
    
    if (Y[t] == 0) {
      # Use log-sum-exp trick
      log_term1 <- log(p_t)
      log_term2 <- log(1 - p_t) - lambda_t
      log_lik <- log_lik + log(exp(log_term1) + exp(log_term2))
    } else {
      log_lik <- log_lik + log(1 - p_t) + Y[t] * log(lambda_t) - 
        lambda_t - lfactorial(Y[t])
    }
    
    if (is.na(log_lik) || is.infinite(log_lik)) {
      return(-Inf)
    }
  }
  
  # Segment 2
  for (t in (m+1):n) {
    z_t <- Z[t, ]
    w_t <- W[t, ]
    
    lambda_t <- exp(sum(z_t * beta2))
    logit_p <- sum(w_t * eta2)
    p_t <- 1 / (1 + exp(-logit_p))
    
    if (is.infinite(lambda_t) || is.na(lambda_t)) {
      return(-Inf)
    }
    
    if (Y[t] == 0) {
      log_term1 <- log(p_t)
      log_term2 <- log(1 - p_t) - lambda_t
      log_lik <- log_lik + log(exp(log_term1) + exp(log_term2))
    } else {
      log_lik <- log_lik + log(1 - p_t) + Y[t] * log(lambda_t) - 
        lambda_t - lfactorial(Y[t])
    }
    
    if (is.na(log_lik) || is.infinite(log_lik)) {
      return(-Inf)
    }
  }
  
  log_prior_m <- -log(n - 2)
  
  return(log_lik + log_prior_m)
}

# Fixed between_model_update with numerical stability
between_model_update <- function(model, Y, beta1, beta2, eta1, eta2, m, Z, W, n,
                                 sigma_b, sigma_d, sigma_beta_sq, sigma_eta_sq) {
  
  q <- length(beta1)
  r <- length(eta1)
  
  if (model == 0) {
    # Birth move: M0 -> M00
    
    m_prop <- sample(2:(n-1), 1)
    
    b_tilde <- rnorm(q, 0, sigma_b)
    d_tilde <- rnorm(r, 0, sigma_d)
    
    beta2_prop <- beta1 + b_tilde
    eta2_prop <- eta1 + d_tilde
    
    log_lik_M00 <- compute_log_likelihood_M00(Y, beta1, beta2_prop, eta1, eta2_prop, 
                                              m_prop, Z, W)
    log_lik_M0 <- compute_log_likelihood_M0(Y, beta1, eta1, Z, W, n)
    
    # Check for NA/Inf
    if (is.na(log_lik_M00) || is.na(log_lik_M0) || 
        is.infinite(log_lik_M00) || is.infinite(log_lik_M0)) {
      return(list(model = 0, beta1 = beta1, beta2 = beta2, 
                  eta1 = eta1, eta2 = eta2, m = m))
    }
    
    log_prior_ratio <- 0  # Equal model priors
    log_proposal_ratio <- -log(n - 2)
    
    log_alpha <- log_lik_M00 - log_lik_M0 + log_prior_ratio + log_proposal_ratio
    
    if (is.na(log_alpha) || is.infinite(log_alpha)) {
      return(list(model = 0, beta1 = beta1, beta2 = beta2, 
                  eta1 = eta1, eta2 = eta2, m = m))
    }
    
    if (log(runif(1)) < log_alpha) {
      return(list(model = 1, beta1 = beta1, beta2 = beta2_prop, 
                  eta1 = eta1, eta2 = eta2_prop, m = m_prop))
    } else {
      return(list(model = 0, beta1 = beta1, beta2 = beta2, 
                  eta1 = eta1, eta2 = eta2, m = m))
    }
    
  } else {
    # Death move: M00 -> M0
    
    beta1_prop <- (m * beta1 + (n - m) * beta2) / n
    eta1_prop <- (m * eta1 + (n - m) * eta2) / n
    
    log_lik_M0 <- compute_log_likelihood_M0(Y, beta1_prop, eta1_prop, Z, W, n)
    log_lik_M00 <- compute_log_likelihood_M00(Y, beta1, beta2, eta1, eta2, m, Z, W)
    
    if (is.na(log_lik_M00) || is.na(log_lik_M0) || 
        is.infinite(log_lik_M00) || is.infinite(log_lik_M0)) {
      return(list(model = 1, beta1 = beta1, beta2 = beta2, 
                  eta1 = eta1, eta2 = eta2, m = m))
    }
    
    log_prior_ratio <- 0
    log_proposal_ratio <- log(n - 2)
    
    log_alpha <- log_lik_M0 - log_lik_M00 + log_prior_ratio + log_proposal_ratio
    
    if (is.na(log_alpha) || is.infinite(log_alpha)) {
      return(list(model = 1, beta1 = beta1, beta2 = beta2, 
                  eta1 = eta1, eta2 = eta2, m = m))
    }
    
    if (log(runif(1)) < log_alpha) {
      return(list(model = 0, beta1 = beta1_prop, beta2 = beta2, 
                  eta1 = eta1_prop, eta2 = eta2, m = m))
    } else {
      return(list(model = 1, beta1 = beta1, beta2 = beta2, 
                  eta1 = eta1, eta2 = eta2, m = m))
    }
  }
}

#  compute_log_likelihood_M0
compute_log_likelihood_M0 <- function(Y, beta, eta, Z, W, n) {
  log_lik <- 0
  
  for (t in 1:n) {
    z_t <- Z[t, ]
    w_t <- W[t, ]
    
    lambda_t <- exp(sum(z_t * beta))
    logit_p <- sum(w_t * eta)
    p_t <- 1 / (1 + exp(-logit_p))
    
    if (is.infinite(lambda_t) || is.na(lambda_t)) {
      return(-Inf)
    }
    
    if (Y[t] == 0) {
      log_term1 <- log(p_t)
      log_term2 <- log(1 - p_t) - lambda_t
      log_lik <- log_lik + log(exp(log_term1) + exp(log_term2))
    } else {
      log_lik <- log_lik + log(1 - p_t) + Y[t] * log(lambda_t) - 
        lambda_t - lfactorial(Y[t])
    }
    
    if (is.na(log_lik) || is.infinite(log_lik)) {
      return(-Inf)
    }
  }
  
  return(log_lik)
}

# compute_log_likelihood_M00
compute_log_likelihood_M00 <- function(Y, beta1, beta2, eta1, eta2, m, Z, W) {
  log_lik <- 0
  n <- length(Y)
  
  # Segment 1
  for (t in 1:m) {
    z_t <- Z[t, ]
    w_t <- W[t, ]
    
    lambda_t <- exp(sum(z_t * beta1))
    logit_p <- sum(w_t * eta1)
    p_t <- 1 / (1 + exp(-logit_p))
    
    if (is.infinite(lambda_t) || is.na(lambda_t)) {
      return(-Inf)
    }
    
    if (Y[t] == 0) {
      log_term1 <- log(p_t)
      log_term2 <- log(1 - p_t) - lambda_t
      log_lik <- log_lik + log(exp(log_term1) + exp(log_term2))
    } else {
      log_lik <- log_lik + log(1 - p_t) + Y[t] * log(lambda_t) - 
        lambda_t - lfactorial(Y[t])
    }
    
    if (is.na(log_lik) || is.infinite(log_lik)) {
      return(-Inf)
    }
  }
  
  # Segment 2
  for (t in (m+1):n) {
    z_t <- Z[t, ]
    w_t <- W[t, ]
    
    lambda_t <- exp(sum(z_t * beta2))
    logit_p <- sum(w_t * eta2)
    p_t <- 1 / (1 + exp(-logit_p))
    
    if (is.infinite(lambda_t) || is.na(lambda_t)) {
      return(-Inf)
    }
    
    if (Y[t] == 0) {
      log_term1 <- log(p_t)
      log_term2 <- log(1 - p_t) - lambda_t
      log_lik <- log_lik + log(exp(log_term1) + exp(log_term2))
    } else {
      log_lik <- log_lik + log(1 - p_t) + Y[t] * log(lambda_t) - 
        lambda_t - lfactorial(Y[t])
    }
    
    if (is.na(log_lik) || is.infinite(log_lik)) {
      return(-Inf)
    }
  }
  
  return(log_lik)
}




# Main RJMCMC algorithm
rjmcmc_zip <- function(Y, Z, W, n_iter = 50000, burn_in = 10000,
                       beta_0 = c(0, 0, 0), eta_0 = c(0, 0, 0),
                       sigma_beta_sq = 1000, sigma_eta_sq = 1000,
                       sigma_b = 0.5, sigma_d = 0.5) {
  
  n <- length(Y)
  q <- length(beta_0)
  r <- length(eta_0)
  
  # Initialize
  model <- 1
  m <- round(n / 2)
  beta1 <- beta_0
  beta2 <- beta_0
  eta1 <- eta_0
  eta2 <- eta_0
  
  # Storage
  model_samples <- numeric(n_iter)
  m_samples <- numeric(n_iter)
  beta1_samples <- matrix(0, n_iter, q)
  beta2_samples <- matrix(0, n_iter, q)
  eta1_samples <- matrix(0, n_iter, r)
  eta2_samples <- matrix(0, n_iter, r)
  
  # MCMC loop
  for (iter in 1:n_iter) {
    
    # Data augmentation
    if (model == 0) {
      aug_data <- data_augmentation(Y, beta1, eta1, Z, W, 1:n)
    } else {
      indices1 <- 1:m
      indices2 <- (m+1):n
      
      aug_data1 <- data_augmentation(Y, beta1, eta1, Z, W, indices1)
      aug_data2 <- data_augmentation(Y, beta2, eta2, Z, W, indices2)
      
      aug_data <- list(
        V1 = aug_data1$V,
        B1 = aug_data1$B,
        V2 = aug_data2$V,
        B2 = aug_data2$B
      )
    }
    
    # Within-model updates
    if (model == 1) {
      # Update m
      m <- update_changepoint(m, Y, beta1, beta2, eta1, eta2, Z, W, n)
      
      indices1 <- 1:m
      indices2 <- (m+1):n
      
      # Re-do data augmentation
      aug_data1 <- data_augmentation(Y, beta1, eta1, Z, W, indices1)
      aug_data2 <- data_augmentation(Y, beta2, eta2, Z, W, indices2)
      
      # Update parameters
      beta1 <- update_beta(aug_data1$V, aug_data1$B, eta1, Z, W, indices1,
                           beta_0, sigma_beta_sq)
      
      beta2 <- update_beta(aug_data2$V, aug_data2$B, eta2, Z, W, indices2,
                           beta_0, sigma_beta_sq)
      
      eta1 <- update_eta(aug_data1$V, aug_data1$B, beta1, Z, W, indices1,
                         eta_0, sigma_eta_sq)
      
      eta2 <- update_eta(aug_data2$V, aug_data2$B, beta2, Z, W, indices2,
                         eta_0, sigma_eta_sq)
    } else {
      # Update parameters for M0
      beta1 <- update_beta(aug_data$V, aug_data$B, eta1, Z, W, 1:n,
                           beta_0, sigma_beta_sq)
      
      eta1 <- update_eta(aug_data$V, aug_data$B, beta1, Z, W, 1:n,
                         eta_0, sigma_eta_sq)
    }
    
    # Between-model update
    result <- between_model_update(model, Y, beta1, beta2, eta1, eta2,
                                   m, Z, W, n, sigma_b, sigma_d,
                                   sigma_beta_sq, sigma_eta_sq)
    
    model <- result$model
    beta1 <- result$beta1
    beta2 <- result$beta2
    eta1 <- result$eta1
    eta2 <- result$eta2
    m <- result$m
    
    # Store samples
    model_samples[iter] <- model
    m_samples[iter] <- m
    beta1_samples[iter, ] <- beta1
    beta2_samples[iter, ] <- beta2
    eta1_samples[iter, ] <- eta1
    eta2_samples[iter, ] <- eta2
    
    if (iter %% 1000 == 0) {
      cat("Iteration:", iter, "\n")
    }
  }
  
  # Remove burn-in
  keep <- (burn_in + 1):n_iter
  
  return(list(
    model = model_samples[keep],
    m = m_samples[keep],
    beta1 = beta1_samples[keep, ],
    beta2 = beta2_samples[keep, ],
    eta1 = eta1_samples[keep, ],
    eta2 = eta2_samples[keep, ]
  ))
}

###################################################

############################################################
############################################################

# Example: Generate data with change-point
n <- 200
m_true <- 100

# Generate covariates
cov_data <- generate_covariates(n)
Z <- cov_data$Z
W <- cov_data$W

# True parameters
beta1_true <- c(0.50, 0.80, -1.0)
beta2_true <- c(-0.50, -0.80, 2.0)
eta1_true <- c(-10.0, 0.10, 0.04)
eta2_true <- c(-2.0, 3.0, 6.0)

# Generate data
Y <- generate_zip_data(n, m_true, beta1_true, beta2_true, 
                       eta1_true, eta2_true, Z, W)

# Run RJMCMC
results <- rjmcmc_zip(Y, Z, W, n_iter = 50000, burn_in = 10000)

# Analyze results
cat("Posterior probability of change-point model:", 
    mean(results$model), "\n")

cat("Estimated change-point location:", 
    mean(results$m[results$model == 1]), "\n")

cat("True change-point:", m_true, "\n")

# Plot posterior distribution of m
hist(results$m[results$model == 1], breaks = 50, 
     main = "Posterior Distribution of Change-point Location",
     xlab = "m", col = "lightblue")
abline(v = m_true, col = "red", lwd = 2, lty = 2)

# Parameter estimates
cat("\nBeta1 estimates:\n")
print(colMeans(results$beta1[results$model == 1, ]))
cat("True values:", beta1_true, "\n")

cat("\nBeta2 estimates:\n")
print(colMeans(results$beta2[results$model == 1, ]))
cat("True values:", beta2_true, "\n")
