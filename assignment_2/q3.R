# estimate the probability that pattern a appears before pattern b
# each replicate uses the shared-stream simulation
pwin <- function(a, b, p = 0.5, nreps = 10000){
  stopifnot(length(a) == length(b), any(a != b))
  
  # Winners are recorded as 1 or 0, so their mean gives the estimated win probability
  mean(sim.duels(a,b,p,nreps)[,1])
}

# estimate the win probability from 10000 independent duels
set.seed(202603)
nr <- 10000
phat <- pwin(c(1,0,0),c(0,0,0),nreps = nr)


#Additionally, calculate the Monte Carlo standard error of the estimated probability
se <- sqrt(phat*(1-phat)/nr)

# calcuate a Wilson 95% confidence interval for the win probability
# qnorm(0.975) gives the 97.5th percentile of the standard normal distribution
z95 <- qnorm(0.975)
den <- 1+z95^2/nr
centre <- (phat+z95^2/(2*nr))/den
half <- z95*sqrt(phat*(1-phat)/nr+z95^2/(4*nr^2))/den

#display the estimate, standard error and the confidence interval
cat("\nQ3: P(100 beats 000), Monte Carlo SE and Wilson 95% interval\n")
print(data.frame(estimate = phat, SE = se, lower = centre-half, upper = centre+half),
      digits = 4, row.names = FALSE)
