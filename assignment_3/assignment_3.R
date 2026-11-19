# Practical 3: Poisson Changepoint Model
# GitHub repository: https://github.com/tttraveller/Extended-statistics-pratical_1-.git
#
# Group number: 36
# Group member contributions:
# Lianjie Wu(s2970562): Tasks 1-3
#   Data simulation, smoothed LASSO objective and derivatives,
#   finite-difference verification, and initial BFGS fitting.
# Sicheng Liang(s2956212): Tasks 4-5
#   mBIC-based penalty selection, optimal model refitting,
#   and changepoint estimation.
# Josie Zhuo(s2959889): Tasks 6-7
#   Parametric bootstrap, bias correction, confidence intervals,
#   and final visualisation.
# All members participated in code review, testing, and integration.

#Task 1: Simulate data from a poisson changepoint model
set.seed(12345)
n <- 100
time <- 1:n

tau_true <- c(1, 30, 70, n+1)  ## cpt positions (with tau_0 and tau_{m+1})
m_true <- length(tau_true)-2  ## number of cpts (m)
theta_true <- c(10, 2, 7)  ## segment parameters (theta_1, theta_2, theta_3)

mu_true <- numeric(n)  ## mean function

for(i in seq_along(theta_true))
  mu_true[tau_true[i]:(tau_true[i+1]-1)] <- theta_true[i]
y <- rpois(n, lambda = mu_true)  ## simulate data y

# plot the observations, true mean and changepoint locations
plot(time, y, pch = 16, col = "grey")
lines(mu_true, col="black", lwd=2, type ="s")
abline(v = tau_true[-c(1,length(tau_true))], lty="dashed")

#Task 2: Smoothed LASSO objective and its derivatives

# Compute the penalised negative Poisson log-likelihood
# The cumulative sum of beta gives the log-mean at each time point
# Only beta[2:n] is penalised, since beta[1] is the baseline
poission_lasso_objective <- function(beta, y, lambda = 1, eps = 1e-6) {
  eta <- cumsum(beta)
  mu <- exp(eta)
  neg_loglik <- sum(mu - y * eta + lgamma(y + 1))
  penalty <- lambda * sum(sqrt(beta[-1L]^2 + eps))
  neg_loglik + penalty
}

# Compute the gradient of the penalised objective
# Reverse cumulative sums compute the negative log-likelihood gradient
# The smoothed LASSO derivative is added for beta[2:n]
poission_lasso_grad <- function(beta, y, lambda = 1, eps = 1e-6) {
  mu <- exp(cumsum(beta))
  gradient <- rev(cumsum(rev(mu - y)))
  gradient[-1L] <- gradient[-1L] +
    lambda * beta[-1L] / sqrt(beta[-1L]^2 + eps)
  gradient
}

# Compute the Hessian matrix of the penalised objective
# Each likelihood entry sums mu from max(i,j) to n
poission_lasso_hess <- function(beta, y, lambda = 1, eps = 1e-6) {
  p <- length(beta)
  mu <- exp(cumsum(beta))
  tail_mu <- rev(cumsum(rev(mu)))
  max_index <- outer(seq_len(p), seq_len(p), FUN = pmax)
  H <- matrix(tail_mu[as.vector(max_index)], nrow = p, ncol = p)
  j <- seq.int(2L, p)
  H[cbind(j, j)] <- H[cbind(j, j)] +
    lambda * eps / (beta[j]^2 + eps)^(3 / 2)
  H
}

# Check the analytical derivatives using central finite differences
# Start from a constant-mean model with no initial changepoints
beta_start <- c(log(mean(y)), rep(0, n - 1L))

# Test the baseline and three selected change parameters
check_idx <- c(1L, 30L, 70L, 100L)
h <- 1e-6

# Approximate the gradient using difference in objective values
numeric_grad <- sapply(check_idx, function(j) {
  bp <- beta_start
  bm <- beta_start
  bp[j] <- bp[j] + h
  bm[j] <- bm[j] - h
  (poission_lasso_objective(bp, y) -
      poission_lasso_objective(bm, y)) / (2 * h)
})
analytic_grad <- poission_lasso_grad(beta_start, y)
grad_error <- max(abs(numeric_grad - analytic_grad[check_idx]))

# Approximate Hessian columns using differences in gradients
numeric_hess_cols <- sapply(check_idx, function(j) {
  bp <- beta_start
  bm <- beta_start
  bp[j] <- bp[j] + h
  bm[j] <- bm[j] - h
  (poission_lasso_grad(bp, y) - poission_lasso_grad(bm, y)) / (2 * h)
})
analytic_hess <- poission_lasso_hess(beta_start, y)
hess_error <- max(abs(numeric_hess_cols - analytic_hess[, check_idx]))

# Report the largest errors and stop if they exceed tolerance
cat("Finite-difference gradient max error:", grad_error, "\n")
cat("Finite-difference Hessian max error:", hess_error, "\n")
stopifnot(grad_error < 1e-3, hess_error < 1e-2)


# Task 3: Initial model fitting with lambda=1

# Fit the penalised Poisson model using BFGS
# If optimisation fails to converge, restart from the last estimate
# with an increased iteration limit
fit_poisson <- function(y, lambda, start) {
  fit <- optim(par = start,
               fn = poission_lasso_objective,
               gr = poission_lasso_grad,
               y = y, lambda = lambda, eps = 1e-6,
               method = "BFGS",
               control = list(maxit = 2000, reltol = 1e-9))
  if (fit$convergence != 0L) {
    fit <- optim(par = fit$par,
                 fn = poission_lasso_objective,
                 gr = poission_lasso_grad,
                 y = y, lambda = lambda, eps = 1e-6,
                 method = "BFGS",
                 control = list(maxit = 4000, reltol = 1e-9))
  }
  fit
}

# Fit the model with the initial penalty lambda=1
fit_initial <- fit_poisson(y, lambda = 1, start = beta_start)

# Check convergence and the remaining gradient magnitude
cat("Task 3 BFGS convergence code:", fit_initial$convergence, "\n")
cat("Task 3 maximum absolute gradient:",
    max(abs(poission_lasso_grad(fit_initial$par, y, lambda = 1))), "\n")
if (fit_initial$convergence != 0L) stop("Task 3 did not converge.")

# Recover the estimated mean from the fitted log-mean parameters
mu_initial <- exp(cumsum(fit_initial$par))

# Compare observations with the true and estimated mean functions
plot(
  time, y,
  pch = 16, col = "grey70",
  xlab = "Time", ylab = "Count",
  main = "Task 3: Initial fit (lambda = 1)"
)

lines(time, mu_true, col = "black", lwd = 2, type = "s")
lines(time, mu_initial, col = "blue", lwd = 2, type = "s")

legend(
  "topright",
  legend = c("Observations", "True mean", "Estimated mean"),
  col = c("grey70", "black", "blue"),
  pch = c(16, NA, NA),
  lty = c(NA, 1, 1),
  cex = 0.8,
  bty = "n"
)


#Task 4: Penalty selection via grid search and mBIC

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
  # Handle possible numerical failures when solving the Hessian system.
  edf <- tryCatch(
    #compute the trace(H_lambda^{-1} * H_0)
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
# Warn if any grid points failed to converge or produced invalid mBIC values
if (any(convergence_codes != 0L, na.rm = TRUE) ||
    any(!is.finite(mbic_values))) {
  warning("Some grid search results are invalid. Inspect Task 4.")
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


#task5 Optimal model refitting and changepoint identification

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
cat("Task 5 maximum absolute gradient:",
    max(abs(poission_lasso_grad(beta_best, y, lambda = lambda_best))), "\n")

# Plot the observations, true mean and estimated mean
plot(time, y,
     pch = 16,
     col = "grey70",
     xlab = "Time",
     ylab = "Count",
     main = "Task 5: Best fitted changepoint model"
     )

lines(time, mu_true, col = "black", lwd = 2, type = "s")
lines(time, mu_best, col = "blue", lwd = 2, type = "s")

legend("topright",
       legend = c("Observations", "True mean", "Estimated mean"),
       col = c("grey70", "black", "blue"),
       pch = c(16, NA, NA),
       lty = c(NA, 1, 1),
       lwd = c(NA, 2, 2),
       cex = 0.8,
       bty = "n")


# Task 6: Parametric bootstrap and bias correction

# Use the selected penalty for all 200 bootstrap samples.
B <- 200L
n_obs <- length(y)

stopifnot(
  length(time) == n_obs,
  length(mu_true) == n_obs,
  length(mu_best) == n_obs,
  length(beta_best) == n_obs,
  all(is.finite(mu_best)), all(mu_best > 0),
  is.finite(lambda_best), lambda_best > 0
)

# Fix the seed and store each fitted bootstrap mean function.
set.seed(20261117)
mu_boot <- matrix(NA_real_, nrow = B, ncol = n_obs)
boot_convergence <- integer(B)

# Simulate Poisson data and refit with the selected penalty.

for (b in seq_len(B)) {
  y_sim <- rpois(n_obs, lambda = mu_best)
  fit_sim <- fit_poisson(y_sim, lambda = lambda_best,
                          start = beta_best)
    
  # Verify convergence before accepting the fitted mean.
  boot_convergence[b] <- fit_sim$convergence
  if (fit_sim$convergence != 0L) {
    stop(sprintf("Task 6: bootstrap fit %d did not converge.", b))
  }
    
  mu_sim <- exp(cumsum(fit_sim$par))
  if (any(!is.finite(mu_sim)) || any(mu_sim <= 0)) {
    stop(sprintf("Task 6: invalid fitted mean for sample %d.", b))
  }
  mu_boot[b, ] <- mu_sim
}


# Estimate the bootstrap mean and correct its pointwise bias.
mu_bar <- colMeans(mu_boot)
bootstrap_bias <- mu_bar - mu_best
mu_bc <- mu_best - bootstrap_bias

# Reflect the 95% percentile limits around the original estimate.
percentile_limits <- apply(mu_boot, 2L, quantile,
                           probs = c(0.025, 0.975), names = FALSE)
q_L <- percentile_limits[1L, ]
q_U <- percentile_limits[2L, ]
q_L_bc <- 2 * mu_best - q_U
q_U_bc <- 2 * mu_best - q_L

stopifnot(all(is.finite(mu_bc)), all(is.finite(q_L_bc)),
          all(is.finite(q_U_bc)), all(q_L_bc <= q_U_bc))

cat("Task 6 completed bootstrap fits:",
    sum(boot_convergence == 0L), "of", B, "\n")
# Compare observed counts and the three mean functions.
plot(time, y, pch = 16, col = "grey70",
     ylim = range(c(0, y, mu_true, mu_best, mu_bar)),
     xlab = "Time", ylab = "Count",
     main = "Task 6: Bootstrap mean comparison")
lines(time, mu_true, col = "black", lwd = 2, type = "s")
lines(time, mu_best, col = "blue", lwd = 2, type = "s")
lines(time, mu_bar, col = "darkorange3", lwd = 2, type = "s")
legend("topright",
       legend = c("Observations", "True mean", "Estimated mean",
                  "Bootstrap mean"),
       col = c("grey70", "black", "blue", "darkorange3"),
       pch = c(16, NA, NA, NA),
       lty = c(NA, 1, 1, 1),
       lwd = c(NA, 2, 2, 2),
       cex = 0.8, bty = "n")


# Task 7: Final plot with bias-corrected confidence limits

# Report negative basic-bootstrap limits without altering them.
negative_lower_idx <- which(q_L_bc < 0)
if (length(negative_lower_idx) > 0L) {
  cat("Task 7: negative corrected lower limits at time points:",
      paste(time[negative_lower_idx], collapse = ", "), "\n")
}

# Include all plotted series when choosing the vertical range.
final_ylim <- range(c(0, y, mu_true, mu_best, mu_bc,
                      q_L_bc, q_U_bc))

plot(time, y, pch = 16, col = "grey70", ylim = final_ylim,
     xlab = "Time", ylab = "Count",
     main = "Task 7: Bias-corrected 95% confidence limits")
# Give the lower and upper confidence limits distinct line types.
lines(time, q_L_bc, col = "darkgreen", lwd = 1.5,
      lty = 2, type = "s")
lines(time, q_U_bc, col = "darkgreen", lwd = 1.5,
      lty = 3, type = "s")
lines(time, mu_true, col = "black", lwd = 2, type = "s")
lines(time, mu_best, col = "blue", lwd = 2, type = "s")
lines(time, mu_bc, col = "red3", lwd = 2, type = "s")
legend("topright",
       legend = c("Observations", "True mean", "Estimated mean",
                  "Bias-corrected mean", "95% lower limit",
                  "95% upper limit"),
       col = c("grey70", "black", "blue", "red3",
               "darkgreen", "darkgreen"),
       pch = c(16, NA, NA, NA, NA, NA),
       lty = c(NA, 1, 1, 1, 2, 3),
       lwd = c(NA, 2, 2, 2, 1.5, 1.5),
       cex = 0.8, bty = "n")
