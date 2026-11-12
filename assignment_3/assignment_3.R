set.seed(12345)
n <- 100
time <- 1:n
tau_true <- c(1, 30, 70, n+1)
m_true <- length(tau_true)-2
theta_true <- c(10, 2, 7)
mu_true <- numeric(n)
## cpt positions (with tau_0 and tau_{m+1})
## number of cpts (m)
## segment parameters (theta_1, theta_2, theta_3)
## mean function
for(i in seq_along(theta_true))
  mu_true[tau_true[i]:(tau_true[i+1]-1)] <- theta_true[i]
y <- rpois(n, lambda = mu_true)
## simulate data y
plot(time, y, pch = 16, col = "grey")
lines(mu_true, col="black", lwd=2, type ="s")
abline(v = tau_true[-c(1,length(tau_true))], lty="dashed")