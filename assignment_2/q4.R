# calculate the exact expected waiting time for a binary pattern a
#each matching prefix and suffix contributes the inverse probability of the prefix
ewait <- function(a, p = 0.5){
  k <- length(a)
  stopifnot(k > 0, all(a %in% 0:1), p > 0, p < 1)
  ew <- 0
  
  # check all possible overlap lengths
  for(j in seq_len(k)){
    if(all(a[1:j] == a[(k-j+1):k])){
      
      # count the heads in the matching prefix and calculate its probability
      # add the inverse probability to the expected waiting time
      h <- sum(a[1:j])
      ew <- ew + 1/(p^h * (1-p)^(j-h))
    }
  }
  ew
}

# compare exact waiting time with Monte Carlo estimates for several patterns
set.seed(202604)
cat("\nQ4: exact waits and Monte Carlo checks (MC SE shown)\n")
aa <- list(c(1,1),c(0,0),c(1,0),c(0,1),c(1,0,1),c(1,1,1))
pp <- c(0.5,0.5,0.5,0.5,0.25,0.7)
wait.check <- data.frame(pattern = sapply(aa,paste,collapse = ""),
                         p = pp, exact = NA_real_, MC = NA_real_, SE = NA_real_)
for(i in seq_along(aa)){
  
  # using the same pattern twice makes the simulation stop at its first appearance
  # the return toss count is therefore the waiting time for that pattern
  w <- sim.duels(aa[[i]],aa[[i]],pp[i],10000)[,2]
  
  #calculate the exact mean, simulated mean and Monte Carlo standard error
  wait.check[i,3:5] <- c(ewait(aa[[i]],pp[i]),mean(w),sd(w)/sqrt(length(w)))
}

#compare the simulation estimates with the exact values
wait.check$relative.error <- (wait.check$MC-wait.check$exact)/wait.check$exact
wait.check$standardised.error <- (wait.check$MC-wait.check$exact)/wait.check$SE
print(wait.check, digits = 4, row.names = FALSE)