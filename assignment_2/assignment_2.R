# Group members:
# Sicheng Liang (s2956212)
# Lianjie Wu (s2970562)
# Josie Zhuo (s2959889)
#
# Contribution Statement:
# - Josie Zhuo
# - Sicheng Liang 
# - Lianjie Wu
#
#
# GitHub repository:
# https://github.com/tttraveller/Extended-statistics-pratical_1-.git
#
# Practical 2: Sequence duels and the perils of independence
#
# This task compares the performance between different binary patterns in a shared coin toss stream. The process estimates the win probabilities, calculates the first waiting times and finds non-transitive cycles for fair and biased coins.
#
# Functions defined:
#   firsthit()       - finds the first occurrence of a pattern
#   sim.duels()      - simulates several duels using shared streams
#   duel()           - simulates one duel
#   pwin()           - estimates the probability that a wins
#   ewait()          - calculates the exact expected waiting time
#   patterns()       - generates all binary patterns of length k
#   pmat()           - calculates the pairwise win-probability matrix
#   best.reply()     - finds the best reply to each pattern
#   shortest.cycle() - searches for the shortest cycle of clear wins

# Task 1: firsthit()
# Find the first occurence of a binary pattern a in a toss stream x
# Return the position where the pattern ends, or NA if it never appear
firsthit <- function(x, a){
  k <- length(a)
  stopifnot(k > 0)
  
  # a pattern cannot appear if the stream is shorter than pattern
  if(length(x) < k){
    return(NA_integer_)
  }
  
  # check all possible strating positions of a lenth-k pattern
  n <- length(x)-k+1
  match <- rep(TRUE,n)
  
  # compare each pattern elements with the corresponding stream positions
  # A starting position remains TRUE only if all k elements match
  for(j in seq_len(k)){
    match <- match & (x[j:(j+n-1)] == a[j])
  }
  ii <- which(match)
  if(length(ii) == 0){
    return(NA_integer_)
  }
  # Convert the first starting position to the position of the last toss
  ii[1]+k-1
}

# Compare the results with know answers
cat("\nQ1: first completion positions\n")
x <- c(1,0,0,1,0,1,1,1,0)
print(data.frame(
  example = c("101 in example stream", "011 in example stream",
              "absent", "stream too short", "overlap", "one toss"),
  expected = c(6,7,NA,NA,3,2),
  observed = c(firsthit(x,c(1,0,1)), firsthit(x,c(0,1,1)),
               firsthit(c(0,0,0),c(1,1)), firsthit(c(1),c(1,0)),
               firsthit(c(1,1,1,1),c(1,1,1)), firsthit(c(0,1),1))))


# Task 2: sim.duels()
# Simulate nreps duels between two binary patterns using shared toss streams
# Each row represents a duel, with p being the probability of heads "1"
# return a matrix which containing the winner and number of tosses for each duel
sim.duels <- function(a, b, p, nreps){
  k <- length(a)
  stopifnot(k > 0, length(b) == k, all(a %in% 0:1), all(b %in% 0:1),
            p > 0, p < 1, nreps >= 1, nreps == floor(nreps))
  
  # Generate tosses in blocks
  # active tracks unfinished duels and tail keeps across block boundaries
  block <- max(32,k)
  ans <- matrix(NA_real_, nreps, 2)
  colnames(ans) <- c("winner","tosses")
  active <- seq_len(nreps)
  tail <- matrix(numeric(0), nreps, 0)
  nt <- 0
  
  # Continue generate tosses only for duels that have not finished
  while(length(active) > 0){
    n <- length(active)
    x <- cbind(tail, matrix(sample(0:1, n*block, replace = TRUE,
                                   prob = c(1-p,p)), nrow = n))
    # repeat both patterns across rows for comparison with each toss stream
    A <- matrix(a,n,k,byrow = TRUE)
    B <- matrix(b,n,k,byrow = TRUE)
    finished <- rep(FALSE,n)
    
    # check each lenth-k window in order to find the first completed pattern
    for(m in k:ncol(x)){
      window <- x[,(m-k+1):m,drop = FALSE]
      
      # rowSums counts mismatches in each row
      # zero means a full match
      ha <- rowSums(window != A) == 0
      hb <- rowSums(window != B) == 0
      
      # record the first match for each unfinished duel
      # delete the retained tail so that earlier tosses are not  counted twice
      ii <- which(!finished & (ha | hb))
      if(length(ii) > 0){
        ans[active[ii],] <- cbind(as.numeric(ha[ii]), nt+m-ncol(tail))
        finished[ii] <- TRUE
      }
      if(all(finished)){
        break
      }
    }
    
    # keep only unfinished duels and the last k-1 tosses
    # these retained tosses allow patterns to cross block boundaries
    ii <- which(!finished)
    active <- active[ii]
    tail <- x[ii,seq_len(k-1)+ncol(x)-k+1,drop = FALSE]
    nt <- nt+block
  }
  ans
}

# Simulate one duel between two different pattens of the same length
# winner is 1 if a appears first, or 0 if b appears first
duel <- function(a, b, p = 0.5){
  stopifnot(length(a) == length(b), any(a != b))
  as.numeric(sim.duels(a,b,p,1)[1,])
}

# ues a fixed seed to make the example reproducible
set.seed(202602)
cat("\nQ2: one duel, 101 against 011 (winner, tosses)\n")
print(duel(c(1,0,1),c(0,1,1)))



# Task 3: pwin()
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



# Task 4: ewait()
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


#Task 5
#generate all binary patterns of length k in standard counting order. 
patterns <- function(k) {
  stopifnot(k >= 1, k == floor(k))
  a <- matrix(0, 2^k, k)
  for (j in seq_len(k)) {
    a[, j] <- (0:(2^k - 1) %/% 2^(k - j)) %% 2
  }
  a
}

pmat <- function(k, p = 0.5, nreps = 10000) {
  a <- patterns(k)
  #generate labels and reshape the number in row
  lab <- apply(a, 1, paste, collapse = "")
  M <- matrix(NA_real_, nrow(a), nrow(a), dimnames = list(lab, lab))
  #iterate all i<j in matrix
  for (i in 1:(nrow(a) - 1)) {
    for (j in (i + 1):nrow(a)) {
      M[i, j] <- pwin(a[i, ], a[j, ], p, nreps)
      M[j, i] <- 1 - M[i, j]  # Symmetrical complement (reduce the computation consumption)
    }
  }
  M
}

# Task 5(a)
# Identify the best counter-strategy to each pattern.
best.reply <- function(M) {
  ii <- apply(M, 2, which.max)
  data.frame(
    pattern = colnames(M),
    reply   = rownames(M)[ii],
    win     = round(M[cbind(ii, seq_len(ncol(M)))], 4),
    row.names = NULL
  )
}


# Task 5(b): shortest.cycle()
# Use Breadth-First Search (BFS) in Graph thory to find the shortest directed cycle of clear wins.
shortest.cycle <- function(M, cutoff = 0.5) {
  #build the Adjacency Matrix
  A <- !is.na(M) & M > cutoff
  best <- integer(0)
  #iterate through each vertex in turn as the start node to search for cycles.
  for (s in seq_len(nrow(M))) {
    paths <- vector("list", nrow(M))
    paths[[s]] <- s
    queue <- s
    h <- 1
    #Breadth-First Search search from distance = 1,2,3 step by step.
    while (h <= length(queue)) {
      v <- queue[h]
      h <- h + 1
      for (w in which(A[v, ])) {
        #case1: if the lose node is the starting node, it becomes a closed cycle(s-v-...-w-s)
        if (w == s) {
          z <- c(paths[[v]], s)
          #update the best cycle when this is the first or shorter cycle.
          if (length(best) == 0 || length(z) < length(best)) {
            best <- z
          }
        } else if (length(paths[[w]]) == 0) {
          paths[[w]] <- c(paths[[v]], w)
          queue <- c(queue, w)
        }
      }
    }
  }
  #return the best cycle
  best
}
# Execution & Analysis
set.seed(42)
nreps <- 10000
M <- pmat(3, p = 0.5, nreps = nreps)

# --- 5(a) Best Replies ---
br <- best.reply(M)
cat("\nQ5(a): Best replies for each pattern (k = 3, p = 0.5):\n")
print(br, row.names = FALSE)

# Because win probabilities are estimated via Monte Carlo simulation, exact equality 
# (M == 0.5) will not occur due to sampling noise.
# Maximum standard error of a Bernoulli proportion estimate is 1 / (2 * sqrt(nreps)).
# We set tolerance to 3 SEs (~0.015) to identify pairs that are essentially evenly matched.

tol <- 3 / (2 * sqrt(nreps))
ij <- which(upper.tri(M) & abs(M - 0.5) <= tol, arr.ind = TRUE)

cat("\nQ5(b): Approximately evenly matched pairs (|P - 0.5| <=", round(tol, 4), "):\n")
print(data.frame(
  pattern_A = rownames(M)[ij[, 1]],
  pattern_B = colnames(M)[ij[, 2]],
  win_A     = round(M[ij], 4)
), row.names = FALSE)


# Explanation
# There are 10 approximately evenly matched pairs under the fair coin.
# These pairs are not included as clear winning edges.
# No 3-cycle is found among the clear wins.
# The shortest cycle has 4 patterns: 001->011->110->100->001.


cyc <- shortest.cycle(M, cutoff = 0.5 + tol)

if (length(cyc) > 0) {
  cat("\nShortest cycle of clear wins (Length =", length(cyc) - 1, "patterns):\n")
  from <- cyc[-length(cyc)]
  to   <- cyc[-1]
  print(data.frame(
    beats        = rownames(M)[from],
    next_pattern = rownames(M)[to],
    win_prob     = round(M[cbind(from, to)], 4)
  ), row.names = FALSE)
} else {
  cat("\nNo clear cycle found.\n")
}


# --- 5(c) ewait vs Head-to-Head Performance ---
cat("\nQ5(c): Isolated expected waiting time vs Competitive performance:\n")
a3 <- patterns(3)
ew <- apply(a3, 1, ewait, p = 0.5)
ranking <- data.frame(
  pattern   = rownames(M),
  wait_ew   = ew,
  mean_win  = round(rowMeans(M, na.rm = TRUE), 4)
)
print(ranking[order(ew), ], row.names = FALSE)
cat("\nQuantified Counter-Example:\n")
cat("Patterns '001' (TTH) and '100' (HTT) both have identical isolated E(W) = 8 tosses.\n")
cat("Yet in a direct head-to-head duel, estimated P('100' beats '001') =", 
    round(M["100", "001"], 4), "(approx 3/4 to 1/4 advantage).\n")
cat("Thus, isolated waiting times fail to determine competitive head-to-head outcomes.\n")



# Task 6
# Visualisation of duel simulation results

set.seed(202606)

# Task 6(a)
# Part 1: Distribution of duel duration

# simulate 10000 duels between 100 and 000
# record the number of tosses until either pattern appears first
z <- sim.duels(c(1,0,0),c(0,0,0),0.5,10000)

# extract the toss counts and calculate the median
tt <- z[,2]
md <- median(tt)

# calculate the empirical cumulative probability for each toss count
F <- cumsum(tabulate(tt))/length(tt)

# save the first plot as a PNG file
png("Q6_duration.png",width = 1000,height = 700,res = 120)

# plot the empirical cumulative distribution of duel duration
plot(0:length(F),c(0,F),type = "s",
     xlab = "Number of tosses",
     ylab = "Estimated P(T <= t)",
     main = "Duel duration: 100 versus 000, p = 0.5")

# mark the median duel duration with a vertical line
lines(c(md,md),c(0,1),col = 2,lty = 2)

# add a horizontal line at probability 0.5
lines(c(0,max(tt)),c(0.5,0.5),col = "grey50",lty = 3)

# show the median value on the plot
text(md,0.15,paste("Median =",md),pos = 4,col = 2)

dev.off()


# Task 6(b)
# Part 2: Monte Carlo variability

# repeat the win probability estimation 50 times
# each estimate uses 2000 independent duels between 011 and 111
est <- rep(NA_real_,50)

for(i in seq_along(est)){
  est[i] <- pwin(c(0,1,1),c(1,1,1),p = 0.5,nreps = 2000)
}

# theoretical win probability for 011 against 111
ref <- 7/8

# save the second plot as a PNG file
png("Q6_variability.png",width = 1000,height = 700,res = 120)

# plot the 50 estimates to show their simulation variability
plot(seq_along(est),est,type = "b",pch = 19,
     ylim = range(c(est,ref))+c(-0.005,0.005),
     xlab = "Independent run (2000 duels each)",
     ylab = "Estimated P(011 beats 111)",
     main = "Monte Carlo variability: 50 independent estimates")

# compare the estimates with the theoretical probability
lines(c(1,50),c(ref,ref),col = 2,lty = 2)

# label the theoretical probability
text(40,ref,"7/8 = 0.875",pos = 1,col = 2)

dev.off()


# display the median duel duration and the variability of the estimates
cat("\nQuestion 6 results:\n")

cat("Median duel duration:",md,"\n")

# calculate the average of the 50 estimated probabilities
cat("Mean of 50 win-probability estimates:",
    round(mean(est),4),"\n")

# calculate the standard deviation between the 50 estimates
cat("Standard deviation of estimates:",
    round(sd(est),4),"\n")

cat("\nPlots saved as Q6_duration.png and Q6_variability.png\n")




# Task 7
set.seed(42)
for(p in c(0.25,0.7)){
  Mp <- pmat(3,p,nreps)
  bp <- best.reply(Mp)
  bp$fair.reply <- br$reply
  cat("\nQ7: best replies, p =",p,"(fair reply shown for comparison)\n")
  print(bp, row.names = FALSE)
  cp <- shortest.cycle(Mp,0.5+tol)
  if(length(cp) == 0){
    cat("No cycle among clear wins at cutoff",0.5+tol,"\n")
  }else{
    cat("Shortest cycle of clear wins:",length(cp)-1,"patterns\n")
    from <- cp[-length(cp)]
    to <- cp[-1]
    print(data.frame(beats = rownames(Mp)[from], next.pattern = rownames(Mp)[to],
                     win = round(Mp[cbind(from,to)],4)), row.names = FALSE)
  }
}
# At p = 0.25, the best reply to 001(TTH) changes from 100(HTT) to 000(TTT); at p = 0.7,
# the best reply to 110(HHT) changes from 011(THH) to 111(HHH). The fair-coin rule is not universal.

# Non-transitivity survives both biases: clear 3-cycles replace the shortest
# clear 4-cycle at p = 0.5. The observed cycles are 000->001->100->000 (p = 0.25)
# and 011->111->101->011 (p = 0.7).

# Cycle claims use only wins above 0.5 + three maximum Monte Carlo SEs;
# no detected cycle means none in this graph, not proof for near-tied true odds.

# Isolated waiting times describe marginal distributions, but cannot decide the probability of win P(Wa < Wb).

# Independent tosses do not imply independent waiting times. Even separate marginal
# distributions would not determine the joint law, only compared with the means of the waiting time will discards 
#even more critical information.
# Thus, equal means therefore do not imply a fair duel; Q5(c) gives a direct example.

# The hazard is assuming independence when it share the same environment. Focus only 
#on individual means makes us confident, but that will leads to false when predict who 
#will win. In this game, both players watch the same coin tosses, they are not independent.
#Comparing only their average waiting times in isolation leads to completely wrong predictions.