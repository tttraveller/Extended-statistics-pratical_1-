#task4
#generate the 121 lambda values
log_lambda_grid <- seq(-1, 5, by = 0.05)
lambda_grid <- exp(log_lambda_grid)
n_lambda <- length(lambda_grid)

# Store the results for each penalty value
mbic_values <- rep(NA_real_, n_lambda)
edf_values <- rep(NA_real_, n_lambda)
loglik_values <- rep(NA_real_, n_lambda)
convergence_codes <- rep(NA_integer_, n_lambda)

# Use warm starts to improve optimisation efficiency.
beta_current <- beta_start

for (i in seq_along(lambda_grid)) {
  lambda_current <- lambda_grid[i]
  # Fit the penalised Poisson model(by BFGS)
  fit <- fit_poisson(
    y = y,
    lambda = lambda_current,
    start = beta_current
  )
  
  # Record the convergence and skip unsuccessful fits
  convergence_codes[i] <- fit$convergence
  
  if (fit$convergence != 0L) {
    warning(paste("Optimisation failed at lambda =", lambda_current))
    next
  }
  
  # Extract the estimated parameters and update the parameters
  beta_hat <- fit$par
  beta_current <- beta_hat
  
  # Compute Hessians
  H_lambda <- poission_lasso_hess(
    beta_hat, y, lambda = lambda_current
  )
  
  H_0 <- poission_lasso_hess(
    beta_hat, y, lambda = 0
  )
  
  # Compute the effective degrees of freedom
  # avoid the numerical Error since the Hessian matrix is numerically singular
  edf <- tryCatch(
    sum(diag(solve(H_lambda, H_0))),
    error = function(e) {
      warning(paste("EDF calculation failed at lambda =", lambda_current))
      NA_real_
    }
  )
  
  if (!is.finite(edf)) next
  
  # Evaluate the unpenalised log-likelihood
  loglik <- -poission_lasso_objective(
    beta_hat, y, lambda = 0
  )
  
  # Compute the effective number of changepoints
  m_tilde <- edf - 1
  
  # Modified BIC
  mbic <- 3 * m_tilde * log(length(y)) - 2 * loglik
  
  # Results
  edf_values[i] <- edf
  loglik_values[i] <- loglik
  mbic_values[i] <- mbic
}

# Combine the results
q4_results <- data.frame(
  log_lambda = log_lambda_grid,
  lambda = lambda_grid,
  EDF = edf_values,
  mBIC = mbic_values,
  convergence = convergence_codes
)

# Check that valid results are available
valid_index <- which(is.finite(mbic_values))

if (length(valid_index) == 0L) {
  stop("No valid mBIC values were obtained.")
}

best_index <- valid_index[
  which.min(mbic_values[valid_index])
]

lambda_best <- lambda_grid[best_index]
log_lambda_best <- log_lambda_grid[best_index]
edf_best <- edf_values[best_index]
m_tilde_best <- edf_best - 1

# Display the selected model
cat("Best lambda:", lambda_best, "\n")
cat("Best log(lambda):", log_lambda_best, "\n")
cat("Minimum mBIC:", mbic_values[best_index], "\n")
cat("Effective degrees of freedom:", edf_best, "\n")
cat("Effective number of changepoints:", m_tilde_best, "\n")
cat("Number of invalid mBIC values:",
    sum(!is.finite(mbic_values)), "\n")
cat("Number of failed optimisations:",
    sum(convergence_codes != 0L, na.rm = TRUE), "\n")

# Plot mBIC against log(lambda)
plot(log_lambda_grid, mbic_values,
     type = "l",
     lwd = 2,
     xlab = "log(lambda)",
     ylab = "mBIC",
     main = "Task 4: Selection of lambda")

abline(v = log_lambda_best, col = "red", lty = 2)
points(log_lambda_best, mbic_values[best_index],
       col = "red", pch = 16)


#task5

# Refit the model
fit_best <- fit_poisson(
  y = y,
  lambda = lambda_best,
  start = beta_start
)

# Check the convergence
cat("Task 5 BFGS convergence code:", fit_best$convergence, "\n")

if (fit_best$convergence != 0L) {
  stop("Task 5 optimisation did not converge.")
}

# Extract the estimated parameters
beta_best <- fit_best$par

# Compute the estimated mean function
mu_best <- exp(cumsum(beta_best))

# Estimate the number of changepoints
m_hat <- floor(m_tilde_best)

# Identify the changepoints using the largest absolute coefficients
beta_changes <- abs(beta_best[-1L])

if (m_hat > 0L) {
  sorted_index <- order(beta_changes, decreasing = TRUE)
  tau_hat <- sort(sorted_index[seq_len(m_hat)] + 1L)
} else {
  tau_hat <- integer(0)
}

# Display the estimated changepoints
cat("Selected lambda:", lambda_best, "\n")
cat("Estimated number of changepoints:", m_hat, "\n")
cat("Estimated changepoint locations:", tau_hat, "\n")

# Plot the observations, true mean and estimated mean
plot(time, y,
     pch = 16,
     col = "grey70",
     xlab = "Time",
     ylab = "Count",
     main = "Task 5: Best fitted changepoint model")

lines(time, mu_true, col = "black", lwd = 2, type = "s")
lines(time, mu_best, col = "blue", lwd = 2, type = "s")

legend("topright",
       legend = c("Observations", "True mean", "Estimated mean"),
       col = c("grey70", "black", "blue"),
       pch = c(16, NA, NA),
       lty = c(NA, 1, 1),
       lwd = c(NA, 2, 2),
       bty = "n")
