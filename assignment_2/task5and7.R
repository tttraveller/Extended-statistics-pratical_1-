#task5
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

#(a)
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
#(b)
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


# ------------------------------------------------------------------------------
# Execution & Analysis
set.seed(42)
nreps <- 10000
M <- pmat(3, p = 0.5, nreps = nreps)

# --- 5(a) Best Replies ---
br <- best.reply(M)
cat("\nQ5(a): Best replies for each pattern (k = 3, p = 0.5):\n")
print(br, row.names = FALSE)


# --- 5(b) Even Pairs & Shortest Non-transitive Cycle ---
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
# Exactly 4 pairs are symmetric complements: (000,111), (001,110), (010,101), (011,100).
# A clear cycle cannot cross an evenly matched pair. Overlapping patterns (000, 111, 010, 101)
# are heavily dominated and act as sources. Thus, no 3-cycle exists, and the shortest
# cycle of clear wins is constrained to have length 4 among the remaining patterns.

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

#task7
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

