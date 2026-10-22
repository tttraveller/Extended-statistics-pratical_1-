
# Group members:
# Sicheng Liang
# Lianjie Wu
# Josie Zhuo
#
# Contribution Statement:
# Please complete the actual contributions of each member.
# GitHub repository:
#
# Practical 2: Sequence duels and the perils of independence
#
# This script simulates duels between binary patterns in a shared
# coin toss stream. It estimates win probabilities, calculates
# expected waiting times, and investigates non-transitive cycles.
#
# Functions defined:
#   firsthit()       - finds the first occurrence of a pattern
#   sim.duels()      - simulates multiple duels at once
#   duel()           - simulates a single duel
#   pwin()           - estimates the probability that a beats b
#   ewait()          - calculates the exact expected waiting time
#   patterns()       - generates all binary patterns of length k
#   pmat()           - builds the win-probability matrix
#   best.reply()     - finds the best reply to each pattern
#   shortest.cycle() - finds the shortest cycle of clear wins


# Task 1
# firsthit: find the first occurrence of pattern a in toss stream x
# input:
#   x: binary vector of coin tosses
#   a: binary pattern to search for
# output:
#   The ending position of the first match, or NA if not found
# steps:
#   (a) check whether the stream is long enough
#   (b) compare the pattern with all possible windows
#   (c) return the position where the first match ends

firsthit <- function(x, a){
  k <- length(a)
  stopifnot(k > 0)
  
  # (a) A pattern cannot appear if the stream is too short
  if(length(x) < k){
    return(NA_integer_)
  }
  
  n <- length(x)-k+1
  match <- rep(TRUE,n)
  
  # (b) Compare each element of a with the matching stream positions
  for(j in seq_len(k)){
    match <- match & (x[j:(j+n-1)] == a[j])
  }
  
  ii <- which(match)
  
  if(length(ii) == 0){
    return(NA_integer_)
  }
  
  # (c) Convert the starting position to the ending position
  ii[1]+k-1
}


# Task 2
# sim.duels: simulate several duels using shared coin toss streams
# input:
#   a, b: binary patterns of the same length
#   p: probability of heads (1)
#   nreps: number of independent duels
# output:
#   ans: a matrix with winner and number of tosses for each duel
# details:
#   Tosses are generated in blocks to avoid a fixed stream length.
#   Unfinished duels keep the last k-1 tosses for the next block.

sim.duels <- function(a, b, p, nreps){
  k <- length(a)
  
  stopifnot(k > 0, length(b) == k,
            all(a %in% 0:1), all(b %in% 0:1),
            p > 0, p < 1,
            nreps >= 1, nreps == floor(nreps))
  
  # Generate tosses in blocks
  block <- max(32,k)
  ans <- matrix(NA_real_,nreps,2)
  colnames(ans) <- c("winner","tosses")
  
  # active stores unfinished duels; tail keeps tosses across blocks
  active <- seq_len(nreps)
  tail <- matrix(numeric(0),nreps,0)
  nt <- 0
  
  while(length(active) > 0){
    n <- length(active)
    
    # Generate another block only for unfinished duels
    x <- cbind(tail,
               matrix(sample(0:1,n*block,replace = TRUE,
                             prob = c(1-p,p)),nrow = n))
    
    A <- matrix(a,n,k,byrow = TRUE)
    B <- matrix(b,n,k,byrow = TRUE)
    finished <- rep(FALSE,n)
    
    # Check each length-k window in order
    for(m in k:ncol(x)){
      window <- x[,(m-k+1):m,drop = FALSE]
      
      # rowSums counts mismatches; zero means a full match
      ha <- rowSums(window != A) == 0
      hb <- rowSums(window != B) == 0
      
      # Record the first completed pattern in each unfinished duel
      ii <- which(!finished & (ha | hb))
      
      if(length(ii) > 0){
        ans[active[ii],] <- cbind(as.numeric(ha[ii]),
                                  nt+m-ncol(tail))
        finished[ii] <- TRUE
      }
      
      if(all(finished)){
        break
      }
    }
    
    # Keep only unfinished duels
    ii <- which(!finished)
    active <- active[ii]
    
    # Retain the last k-1 tosses for matches across block boundaries
    if(k > 1){
      tail <- x[ii,(ncol(x)-k+2):ncol(x),drop = FALSE]
    }else{
      tail <- matrix(numeric(0),length(ii),0)
    }
    
    nt <- nt+block
  }
  
  return(ans)
}


# duel: simulate one duel between two different patterns
# input:
#   a, b: distinct binary patterns of equal length
#   p: probability of heads
# output:
#   c(winner,tosses), where winner is 1 if a wins and 0 if b wins

duel <- function(a, b, p = 0.5){
  stopifnot(length(a) == length(b),any(a != b))
  
  return(as.numeric(sim.duels(a,b,p,1)[1,]))
}


# Task 3
# pwin: estimate the probability that pattern a appears before b
# input:
#   a, b: distinct binary patterns
#   p: probability of heads
#   nreps: number of independent duels
# output:
#   Estimated win probability of pattern a

pwin <- function(a, b, p = 0.5, nreps = 10000){
  stopifnot(length(a) == length(b),any(a != b))
  
  # Winners are 1 or 0, so their mean gives the win probability
  return(mean(sim.duels(a,b,p,nreps)[,1]))
}


# Task 4
# ewait: calculate the exact expected waiting time of a pattern
# input:
#   a: binary pattern
#   p: probability of heads
# output:
#   ew: expected number of tosses until a first appears
# steps:
#   (a) check each possible overlap length
#   (b) calculate the probability of a matching prefix
#   (c) sum the inverse probabilities

ewait <- function(a, p = 0.5){
  k <- length(a)
  stopifnot(k > 0,all(a %in% 0:1),p > 0,p < 1)
  
  ew <- 0
  
  # (a) Check whether the first j tosses equal the last j tosses
  for(j in seq_len(k)){
    if(all(a[1:j] == a[(k-j+1):k])){
      
      # (b) Count heads in the matching prefix
      h <- sum(a[1:j])
      
      # (c) Add the inverse probability
      ew <- ew+1/(p^h*(1-p)^(j-h))
    }
  }
  
  return(ew)
}


# Task 5
# patterns: generate all binary patterns of length k
# input:
#   k: pattern length
# output:
#   A matrix with 2^k rows, one pattern per row
#   Patterns are arranged in binary counting order.

patterns <- function(k){
  stopifnot(k >= 1,k == floor(k))
  
  a <- matrix(0,2^k,k)
  
  for(j in seq_len(k)){
    a[,j] <- (0:(2^k-1) %/% 2^(k-j)) %% 2
  }
  
  return(a)
}


# pmat: estimate the full matrix of pairwise win probabilities
# input:
#   k: pattern length
#   p: probability of heads
#   nreps: number of simulated duels for each pair
# output:
#   M[i,j]: estimated probability that pattern i beats pattern j
#   The diagonal is NA because a pattern cannot duel itself.

pmat <- function(k, p = 0.5, nreps = 10000){
  a <- patterns(k)
  
  # Create labels for the rows and columns
  lab <- apply(a,1,paste,collapse = "")
  M <- matrix(NA_real_,nrow(a),nrow(a),
              dimnames = list(lab,lab))
  
  # Only simulate one direction for each pair
  for(i in seq_len(nrow(a)-1)){
    for(j in (i+1):nrow(a)){
      M[i,j] <- pwin(a[i,],a[j,],p,nreps)
      
      # The opposite win probability is its complement
      M[j,i] <- 1-M[i,j]
    }
  }
  
  return(M)
}


# Task 5(a)
# best.reply: find the best counter-pattern for each pattern
# input:
#   M: matrix of pairwise win probabilities
# output:
#   A data frame containing each pattern, its best reply,
#   and the estimated probability that the reply wins
# details:
#   Each column represents the pattern being challenged.
#   The largest entry in that column gives its best reply.

best.reply <- function(M){
  ii <- apply(M,2,which.max)
  
  return(data.frame(
    pattern = colnames(M),
    reply = rownames(M)[ii],
    win = round(M[cbind(ii,seq_len(ncol(M)))],4),
    row.names = NULL
  ))
}


# Task 5(b)
# shortest.cycle: find the shortest directed cycle of clear wins
# input:
#   M: matrix of pairwise win probabilities
#   cutoff: minimum probability required for a winning edge
# output:
#   A vector of pattern indices forming the shortest cycle
#   The first index is repeated at the end to close the cycle.
# details:
#   Breadth-first search (BFS) checks the shortest paths first.
#   An edge from i to j means pattern i beats pattern j.

shortest.cycle <- function(M, cutoff = 0.5){
  # Build the adjacency matrix
  A <- !is.na(M) & M > cutoff
  best <- integer(0)
  
  # Search for cycles starting from each pattern
  for(s in seq_len(nrow(M))){
    paths <- vector("list",nrow(M))
    paths[[s]] <- s
    
    queue <- s
    h <- 1
    
    # BFS checks paths in increasing order of length
    while(h <= length(queue)){
      v <- queue[h]
      h <- h+1
      
      for(w in which(A[v,])){
        
        # Case 1: Return to the starting pattern
        if(w == s){
          z <- c(paths[[v]],s)
          
          if(length(best) == 0 || length(z) < length(best)){
            best <- z
          }
          
          # Case 2: Continue to a pattern not yet visited
        }else if(length(paths[[w]]) == 0){
          paths[[w]] <- c(paths[[v]],w)
          queue <- c(queue,w)
        }
      }
    }
  }
  
  return(best)
}


# ========== Main execution ==========


# Task 1: Test firsthit with known examples

cat("\n=== Task 1: First completion positions ===\n")

x <- c(1,0,0,1,0,1,1,1,0)

# Compare the function results with manually checked positions
print(data.frame(
  example = c("101 in example stream","011 in example stream",
              "absent","stream too short","overlap","one toss"),
  expected = c(6,7,NA,NA,3,2),
  observed = c(firsthit(x,c(1,0,1)),
               firsthit(x,c(0,1,1)),
               firsthit(c(0,0,0),c(1,1)),
               firsthit(c(1),c(1,0)),
               firsthit(c(1,1,1,1),c(1,1,1)),
               firsthit(c(0,1),1))
))


# Task 2: Simulate one duel

set.seed(202602)

cat("\n=== Task 2: One duel between 101 and 011 ===\n")

# Return the winner and the number of tosses
print(duel(c(1,0,1),c(0,1,1)))


# Task 3: Estimate a win probability

set.seed(202603)
nr <- 10000

# Estimate the probability that 100 beats 000
phat <- pwin(c(1,0,0),c(0,0,0),nreps = nr)

# Monte Carlo standard error
se <- sqrt(phat*(1-phat)/nr)

# Wilson 95% confidence interval
z95 <- qnorm(0.975)
den <- 1+z95^2/nr
centre <- (phat+z95^2/(2*nr))/den
half <- z95*sqrt(phat*(1-phat)/nr+
                   z95^2/(4*nr^2))/den

cat("\n=== Task 3: Estimated win probability ===\n")

print(data.frame(
  estimate = phat,
  SE = se,
  lower = centre-half,
  upper = centre+half
),digits = 4,row.names = FALSE)


# Task 4: Exact waiting times and Monte Carlo checks

set.seed(202604)

aa <- list(c(1,1),c(0,0),c(1,0),
           c(0,1),c(1,0,1),c(1,1,1))

pp <- c(0.5,0.5,0.5,0.5,0.25,0.7)

wait.check <- data.frame(
  pattern = sapply(aa,paste,collapse = ""),
  p = pp,
  exact = NA_real_,
  MC = NA_real_,
  SE = NA_real_
)

for(i in seq_along(aa)){
  
  # Use the shared simulator to record a single pattern's waiting time.
  # The same pattern is placed on both sides; only tosses are used.
  w <- sim.duels(aa[[i]],aa[[i]],pp[i],10000)[,2]
  
  wait.check[i,3:5] <- c(
    ewait(aa[[i]],pp[i]),
    mean(w),
    sd(w)/sqrt(length(w))
  )
}

# Compare simulated and exact waiting times
wait.check$relative.error <- (wait.check$MC-wait.check$exact)/
  wait.check$exact

wait.check$standardised.error <- (wait.check$MC-wait.check$exact)/
  wait.check$SE

cat("\n=== Task 4: Exact and simulated waiting times ===\n")
print(wait.check,digits = 4,row.names = FALSE)


# Task 5: Win-probability matrix

set.seed(42)
nreps <- 10000

M <- pmat(3,p = 0.5,nreps = nreps)


# Task 5(a): Best replies for all length-3 patterns

br <- best.reply(M)

cat("\n=== Task 5(a): Best replies, p = 0.5 ===\n")
print(br,row.names = FALSE)


# Task 5(b): Evenly matched pairs and non-transitive cycles

# Monte Carlo estimates are not exactly 0.5 due to sampling variation.
# Three maximum standard errors are used as a tolerance.
tol <- 3/(2*sqrt(nreps))

ij <- which(upper.tri(M) & abs(M-0.5) <= tol,
            arr.ind = TRUE)

cat("\n=== Task 5(b): Approximately even pairs ===\n")

print(data.frame(
  pattern_A = rownames(M)[ij[,1]],
  pattern_B = colnames(M)[ij[,2]],
  win_A = round(M[ij],4)
),row.names = FALSE)

# Several pairs are evenly matched under a fair coin.
# The four complementary pairs are not the only even pairs.
# These pairs cannot form clear winning edges.
# The shortest cycle has 4 patterns, and no clear 3-cycle is found.

cyc <- shortest.cycle(M,cutoff = 0.5+tol)

if(length(cyc) > 0){
  cat("\nShortest cycle of clear wins:",
      length(cyc)-1,"patterns\n")
  
  from <- cyc[-length(cyc)]
  to <- cyc[-1]
  
  print(data.frame(
    beats = rownames(M)[from],
    next_pattern = rownames(M)[to],
    win_prob = round(M[cbind(from,to)],4)
  ),row.names = FALSE)
  
}else{
  cat("\nNo clear cycle found.\n")
}


# Task 5(c): Expected waiting times and duel performance

cat("\n=== Task 5(c): Waiting times and win probabilities ===\n")

a3 <- patterns(3)
ew <- apply(a3,1,ewait,p = 0.5)

ranking <- data.frame(
  pattern = rownames(M),
  wait_ew = ew,
  mean_win = round(rowMeans(M,na.rm = TRUE),4)
)

print(ranking[order(ew),],row.names = FALSE)

# 001 and 100 have the same expected waiting time,
# but 100 beats 001 with probability 0.75.
cat("\nPatterns 001 and 100 both have E(W) = 8 tosses.\n")
cat("Estimated P(100 beats 001) =",
    round(M["100","001"],4),"\n")
cat("Theoretical probability = 0.75\n")

# Equal expected waiting times do not imply equal duel probabilities.


# Task 6(a): Distribution of duel duration

set.seed(202606)

# Simulate 10000 duels between 100 and 000
z <- sim.duels(c(1,0,0),c(0,0,0),0.5,10000)

tt <- z[,2]
md <- median(tt)

# Empirical cumulative distribution of duel duration
F <- cumsum(tabulate(tt))/length(tt)

png("Q6_duration.png",width = 1000,height = 700,res = 120)

plot(0:length(F),c(0,F),type = "s",
     xlab = "Number of tosses",
     ylab = "Estimated P(T <= t)",
     main = "Duel duration: 100 versus 000, p = 0.5")

# Mark the median and probability 0.5
lines(c(md,md),c(0,1),col = 2,lty = 2)
lines(c(0,max(tt)),c(0.5,0.5),col = "grey50",lty = 3)
text(md,0.15,paste("Median =",md),pos = 4,col = 2)

dev.off()


# Task 6(b): Monte Carlo variability

# Run 50 independent estimates, with 2000 duels each
est <- rep(NA_real_,50)

for(i in seq_along(est)){
  est[i] <- pwin(c(0,1,1),c(1,1,1),nreps = 2000)
}

# For a fair coin, P(011 beats 111) = 7/8
ref <- 7/8

png("Q6_variability.png",width = 1000,height = 700,res = 120)

plot(seq_along(est),est,type = "b",pch = 19,
     ylim = range(c(est,ref))+c(-0.005,0.005),
     xlab = "Independent run (2000 duels each)",
     ylab = "Estimated P(011 beats 111)",
     main = "Monte Carlo variability: 50 independent estimates")

# Compare the estimates with the theoretical probability
lines(c(1,50),c(ref,ref),col = 2,lty = 2)
text(40,ref,"7/8 = 0.875",pos = 1,col = 2)

dev.off()

cat("\n=== Task 6: Simulation summary ===\n")
cat("Duration median:",md,"\n")
cat("Mean estimated probability:",round(mean(est),4),"\n")
cat("Standard deviation:",round(sd(est),4),"\n")
cat("Plots saved to Q6_duration.png and Q6_variability.png\n")


# Task 7: Effect of coin bias on best replies and cycles

set.seed(42)

for(p in c(0.25,0.7)){
  Mp <- pmat(3,p,nreps)
  bp <- best.reply(Mp)
  
  # Compare the best replies with the fair-coin results
  bp$fair.reply <- br$reply
  
  cat("\n=== Task 7: Best replies, p =",p,"===\n")
  print(bp,row.names = FALSE)
  
  # Find the shortest cycle with clear winning edges
  cp <- shortest.cycle(Mp,0.5+tol)
  
  if(length(cp) == 0){
    cat("No clear cycle found at cutoff",0.5+tol,"\n")
    
  }else{
    cat("Shortest cycle of clear wins:",
        length(cp)-1,"patterns\n")
    
    from <- cp[-length(cp)]
    to <- cp[-1]
    
    print(data.frame(
      beats = rownames(Mp)[from],
      next.pattern = rownames(Mp)[to],
      win = round(Mp[cbind(from,to)],4)
    ),row.names = FALSE)
  }
}


# Task 7: Results and comments

# At p = 0.25, the best reply to 001 changes from 100 to 000.
# At p = 0.7, the best reply to 110 changes from 011 to 111.
# Therefore, the best replies depend on the coin bias.

# Non-transitivity still appears in both biased cases.
# At p = 0.25, one cycle is 000 -> 001 -> 100 -> 000.
# At p = 0.7, one cycle is 011 -> 111 -> 101 -> 011.
# These cycles have three patterns, while the shortest clear
# cycle for the fair coin has four patterns.

# Expected waiting time only describes how long one pattern
# takes to appear on its own. It does not determine which
# pattern is more likely to win a duel.

# Although each coin toss is independent, the waiting times
# of the two patterns are not independent. Both patterns
# are observed in the same toss stream, and their overlaps
# affect which one appears first.

# Comparing only the separate expected waiting times can
# therefore give misleading predictions about duel outcomes.
