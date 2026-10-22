# Question 6: Visualisation of duel simulation results
set.seed(202606)

# ------------------------------------------------------------
# Part 1: Distribution of duel duration
# ------------------------------------------------------------
z <- sim.duels(
  c(1, 0, 0),
  c(0, 0, 0),
  0.5,
  10000
)

tt <- z[, 2]
md <- median(tt)
F <- cumsum(tabulate(tt)) / length(tt)
png("Q6_duration.png", width = 1000, height = 700, res = 120)

plot(
  0:length(F),
  c(0, F),
  type = "s",
  xlab = "Number of tosses",
  ylab = "Estimated P(T <= t)",
  main = "Duel duration: 100 versus 000, p = 0.5"
)

lines(
  c(md, md),
  c(0, 1),
  col = 2,
  lty = 2
)

lines(
  c(0, max(tt)),
  c(0.5, 0.5),
  col = "grey50",
  lty = 3
)

text(
  md,
  0.15,
  paste("Median =", md),
  pos = 4,
  col = 2
)

dev.off()

# ------------------------------------------------------------
# Part 2: Monte Carlo variability
# ------------------------------------------------------------

est <- rep(NA_real_, 50)

for (i in seq_along(est)) {
  est[i] <- pwin(
    c(0, 1, 1),
    c(1, 1, 1),
    p = 0.5,
    nreps = 2000
  )
}


ref <- 7 / 8
png("Q6_variability.png", width = 1000, height = 700, res = 120)

plot(
  seq_along(est),
  est,
  type = "b",
  pch = 19,
  ylim = range(c(est, ref)) + c(-0.005, 0.005),
  xlab = "Independent run (2000 duels each)",
  ylab = "Estimated P(011 beats 111)",
  main = "Monte Carlo variability: 50 independent estimates"
)

lines(
  c(1, 50),
  c(ref, ref),
  col = 2,
  lty = 2
)

text(
  40,
  ref,
  "7/8 = 0.875",
  pos = 1,
  col = 2
)

dev.off()


cat("\nQuestion 6 results:\n")

cat(
  "Median duel duration:",
  md, "\n"
)

cat(
  "Mean of 50 win-probability estimates:",
  round(mean(est), 4), "\n"
)

cat(
  "Standard deviation of estimates:",
  round(sd(est), 4), "\n"
)

cat(
  "\nPlots saved as Q6_duration.png and Q6_variability.png\n"
)
