#Task 2: Smoothed LASSO objective and its derivatives

# Compute the penalised negative Poisson log-likelihood
# The cimulative sum of beta gives the log-mean at each time point
# Only beta[2:n]is penalised, since beta[1] is the baseline
poission_lasso_objective <- function(beta, y, lambda = 1, eps = 1e-6) {
  eta <- cumsum(beta)
  mu <- exp(eta)
  neg_loglik <- sum(mu - y * eta + lgamma(y + 1))
  penalty <- lambda * sum(sqrt(beta[-1L]^2 + eps))
  neg_loglik + penalty
}

# Compute the gradient of the penalised objective
# Reverse cumulative sums give the derivatives of the likelihood
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

# Compare observations with the ture and estimated mean functions
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
  cex = 0.85,
  bty = "n"
)