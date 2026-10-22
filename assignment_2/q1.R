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
