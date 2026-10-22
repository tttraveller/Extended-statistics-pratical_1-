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
