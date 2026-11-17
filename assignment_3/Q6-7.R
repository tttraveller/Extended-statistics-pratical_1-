# Task 6: Parametric bootstrap and bias correction

# Simulate from the best fitted model using a fixed penalty.
B <- 200L
n_obs <- length(y)

# Check the fitted values and required inputs.
stopifnot(
  length(time) == n_obs,
  length(mu_true) == n_obs,
  length(mu_best) == n_obs,
  length(beta_best) == n_obs,
  all(is.finite(mu_best)),
  all(mu_best > 0),
  is.finite(lambda_best),
  lambda_best > 0
)

set.seed(20261117)

# Store one fitted mean function per bootstrap replicate.
mu_boot <- matrix(NA_real_, nrow = B, ncol = n_obs)
boot_convergence <- integer(B)
boot_max_gradient <- numeric(B)

# Simulate Poisson data and refit with the original estimates.
boot_elapsed <- system.time({
  for (b in seq_len(B)) {
    y_sim <- rpois(n_obs, lambda = mu_best)
    fit_sim <- fit_poisson(
      y = y_sim,
      lambda = lambda_best,
      start = beta_best
    )
    
    boot_convergence[b] <- fit_sim$convergence
    if (fit_sim$convergence != 0L) {
      stop(sprintf("Task 6: bootstrap fit %d did not converge.", b))
    }
    
    # Check the gradient after the initial BFGS fit.
    grad_sim <- max(abs(poission_lasso_grad(
      fit_sim$par, y_sim, lambda = lambda_best
    )))
    
    # Refine fits with comparatively large gradients.
    # The threshold is a diagnostic choice, not a convergence rule.
    if (grad_sim > 0.1) {
      fit_refined <- optim(
        par = fit_sim$par,
        fn = poission_lasso_objective,
        gr = poission_lasso_grad,
        y = y_sim,
        lambda = lambda_best,
        eps = 1e-6,
        method = "BFGS",
        control = list(maxit = 4000, reltol = 1e-12)
      )
      
      # Keep the refined fit only if it converged and improved.
      if (fit_refined$convergence == 0L &&
          is.finite(fit_refined$value) &&
          fit_refined$value <= fit_sim$value) {
        fit_sim <- fit_refined
      }
    }
    
    mu_sim <- exp(cumsum(fit_sim$par))
    if (any(!is.finite(mu_sim)) || any(mu_sim <= 0)) {
      stop(sprintf("Task 6: invalid mean at bootstrap fit %d.", b))
    }
    
    mu_boot[b, ] <- mu_sim
    boot_max_gradient[b] <- max(abs(poission_lasso_grad(
      fit_sim$par, y_sim, lambda = lambda_best
    )))
  }
})[["elapsed"]]

# Average fitted means across bootstrap samples.
mu_bar <- colMeans(mu_boot)

# Estimate bootstrap bias and correct the fitted mean.
bootstrap_bias <- mu_bar - mu_best
mu_bc <- mu_best - bootstrap_bias

# Compute pointwise 2.5% and 97.5% bootstrap quantiles.
percentile_limits <- apply(
  mu_boot, 2L, quantile,
  probs = c(0.025, 0.975), names = FALSE
)
q_L <- percentile_limits[1L, ]
q_U <- percentile_limits[2L, ]

# Form basic bootstrap intervals using reflected quantiles.
# Keep the original limits for statistical reporting.
q_L_bc <- 2 * mu_best - q_U
q_U_bc <- 2 * mu_best - q_L

stopifnot(
  all(is.finite(mu_bc)),
  all(is.finite(q_L_bc)),
  all(is.finite(q_U_bc)),
  all(q_L_bc <= q_U_bc)
)

cat("Task 6 completed bootstrap fits:",
    sum(boot_convergence == 0L), "of", B, "\n")
cat("Task 6 bootstrap elapsed time (seconds):",
    round(boot_elapsed, 2), "\n")
cat("Task 6 maximum absolute gradient:",
    signif(max(boot_max_gradient), 4), "\n")

# Compare observations, true mean and fitted mean functions.
plot(time, y,
     pch = 16, col = "grey70",
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

cat("Task 6 plot displayed.\n")

# Task 7: Final plot with bias-corrected confidence limits

# Means from Poisson models cannot be negative. The basic bootstrap
# interval may nevertheless have a negative lower limit near a change.
# Clip lower limits only for plotting; keep q_L_bc unchanged.
negative_lower_idx <- which(q_L_bc < 0)
q_L_bc_plot <- pmax(0, q_L_bc)

if (length(negative_lower_idx) > 0L) {
  cat("Task 7 negative lower limits:",
      length(negative_lower_idx), "\n")
  cat("Task 7 affected time points:",
      paste(time[negative_lower_idx], collapse = ", "), "\n")
  cat("Task 7 lower limits clipped to 0 for plotting only.\n")
}

# Include all displayed curves when setting the y-axis range.
final_ylim <- range(c(0, y, mu_true, mu_best, mu_bc,
                      q_L_bc_plot, q_U_bc))

plot(time, y,
     pch = 16, col = "grey70", ylim = final_ylim,
     xlab = "Time", ylab = "Count",
     main = "Task 7: Bias-corrected 95% confidence limits")

# Draw interval bounds before the fitted mean curves.
lines(time, q_L_bc_plot, col = "darkgreen", lwd = 1.5,
      lty = 3, type = "s")
lines(time, q_U_bc, col = "darkgreen", lwd = 1.5,
      lty = 3, type = "s")
lines(time, mu_true, col = "black", lwd = 2, type = "s")
lines(time, mu_best, col = "blue", lwd = 2, type = "s")
lines(time, mu_bc, col = "red3", lwd = 2, type = "s")

# Disclose the plotting-only truncation in the legend.
limits_label <- if (length(negative_lower_idx) > 0L) {
  "95% limits (lower clipped at 0)"
} else {
  "95% corrected limits"
}

legend("topright",
       legend = c("Observations", "True mean", "Estimated mean",
                  "Bias-corrected mean", limits_label),
       col = c("grey70", "black", "blue", "red3", "darkgreen"),
       pch = c(16, NA, NA, NA, NA),
       lty = c(NA, 1, 1, 1, 3),
       lwd = c(NA, 2, 2, 2, 1.5),
       cex = 0.76, bty = "n")

cat("Task 7 plot displayed.\n")
