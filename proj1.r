music_raw <- readLines("palestrina.txt")
length(music_raw)
head(music_raw,50)

preprocess <- function(rawdata) {
  # (a) Remove first 3 header lines
  if(length(rawdata) > 3){
    rawdata <- rawdata[-(1:3)]
  }
  # (b) Remove full comment lines starting with #; keep notes like F#5e
  remove_rows <- grepl("^#", rawdata)
  rawdata <- rawdata[!remove_rows]
  
  # (c) Detect blank lines (end of piece), append " ||" to the previous line
  is_blank <- trimws(rawdata) == ""
  keep_line <- !is_blank
  out_lines <- rawdata[keep_line]
  blank_idx <- which(is_blank)
  target_idx <- blank_idx - 1
  target_idx <- target_idx[target_idx >= 1 & target_idx <= length(out_lines)]
  if(length(target_idx) > 0){
    out_lines[target_idx] <- paste0(out_lines[target_idx], " ||")
  }
  
  # (d) Append " |" to normal lines; skip lines already ending with || or |
  for (i in seq_along(out_lines)) {
    line_tmp <- out_lines[i]
    if (!grepl("\\|\\|$", line_tmp) && !grepl("\\|$", line_tmp)) {
      out_lines[i] <- paste0(line_tmp, " |")
    }
  }
  
  # (e) Split into individual tokens by whitespace, remove empty strings
  split_result <- strsplit(out_lines, "\\s+")
  all_tokens <- unlist(split_result)
  clean_tokens <- all_tokens[all_tokens != ""]
  
  return(clean_tokens)
}
# ========== Main execution ==========
music_raw <- readLines("palestrina.txt")
music_clean <- preprocess(music_raw)
# Task5 Markov chain setup
notes <- unique(music_clean)
tokens <- match(music_clean, notes)
# Safety check for end marker ||
end <- which(notes == "||")
if(length(end)==0){
  warning("End marker || not found in tokens")
}

# ==============================================================================
# Task 6: Construct Lagged Transition Matrix make_matrix
# ==============================================================================
make_matrix <- function(tokens, end, mlag = 4) {
  n <- length(tokens)
  
  # (a) Create an (n - mlag) x (mlag + 1) matrix for adjacent token sequences
  M <- matrix(0L, nrow = n - mlag, ncol = mlag + 1)
  for (j in 1:(mlag + 1)) {
    M[, j] <- tokens[j:(n - mlag + j - 1)]
  }
  
  # (b) Discard rows where the double bar-line token '||' appears in the first mlag columns
  # Sequences spanning across two separate musical pieces must not be used
  rows_with_end <- which(rowSums(M[, 1:mlag, drop = FALSE] == end) > 0)
  if (length(rows_with_end) > 0) {
    M <- M[-rows_with_end, , drop = FALSE]
  }
  
  return(M)
}

cat("Constructing transition matrix M (mlag = 4)...\n")
end_token_id <- which(notes == "||")
M <- make_matrix(tokens = tokens, end = end_token_id, mlag = 4)
cat("Matrix M successfully built! Dimensions:", dim(M)[1], "x", dim(M)[2], "\n\n")

# Task 7: Next-Note Sampling Function next_note
# ==============================================================================
next_note <- function(key, M, w = NULL) {
  # 1. Determine maximum lag mlag from matrix dimensions
  mlag <- ncol(M) - 1
  
  # Default mixture weights to a vector of ones of length mlag
  if (is.null(w)) {
    w <- rep(1, mlag)
  }
  
  # 2. Retain only the trailing mlag tokens if key exceeds mlag
  if (length(key) > mlag) {
    key <- tail(key, mlag)
  }
  
  k_max <- length(key)
  cand_tokens <- integer(0)
  cand_probs <- numeric(0)
  
  # 3. Iterate over matching order k from 1 to k_max (k notes suffix)
  for (k in 1:k_max) {
    # According to formula (1): sum_{i=1}^m w_i P(next | v[i:m])
    # A sequence of length k corresponds to i = mlag - k + 1
    # Thus, the highest order (k = mlag) uses w[1], and lowest (k = 1) uses w[mlag]
    w_idx <- mlag - k + 1
    weight_k <- w[w_idx]
    
    # Skip matching when mixture weight is non-positive
    if (weight_k <= 0) next
    
    # Extract suffix of length k and select corresponding columns in M
    mc <- mlag - k + 1
    sub_key <- key[(k_max - k + 1):k_max]
    
    # Vectorized matching via colSums on transposed sub-matrix
    ii <- colSums(!(t(M[, mc:mlag, drop = FALSE]) == sub_key))
    matched_idx <- which(ii == 0)
    n_u <- length(matched_idx)
    
    # If matches are found, accumulate candidate tokens and their weights
    if (n_u > 0) {
      cand_tokens <- c(cand_tokens, M[matched_idx, mlag + 1])
      cand_probs <- c(cand_probs, rep(weight_k / n_u, n_u))
    }
  }
  
  # 4. Sampling decision rules (at most one call to sample)
  total_prob <- sum(cand_probs)
  
  if (length(cand_tokens) == 0 || total_prob == 0) {
    # Case 1: No match found -> sample from marginal next-token distribution
    return(sample(M[, mlag + 1], size = 1))
  } else if (length(cand_tokens) == 1) {
    # Case 2: Exactly one candidate -> return directly without invoking sample()
    return(cand_tokens[1])
  } else {
    # Case 3: Multiple candidates -> sample according to empirical mixture probabilities
    return(sample(cand_tokens, size = 1, prob = cand_probs))
  }
}

# ==============================================================================
# Task 7 Verification Tests
# ==============================================================================
cat("=== Testing next_note function ===\n")

# Test 1: Predict next note given an initial sequence (e.g., G4w D5w)
test_prefix_tokens <- match(c("G4w", "D5w"), notes)
set.seed(42)
pred_token <- next_note(key = test_prefix_tokens, M = M)
cat("Input sequence: G4w D5w\n")
cat("Predicted next Token ID:", pred_token, "| Corresponding note:", notes[pred_token], "\n\n")

# Test 2: Unseen prefix triggering Case 1 fallback sampling
dummy_key <- c(99999L, 99999L)
fallback_token <- next_note(key = dummy_key, M = M)
cat("Fallback sampled note for non-existent prefix:", notes[fallback_token], "\n\n")

# Test 3: Sequential autoregressive generation of 10 notes
seq_test <- test_prefix_tokens
for (step in 1:10) {
  nxt <- next_note(seq_test, M)
  seq_test <- c(seq_test, nxt)
}
cat("10-note sequence generation demo:\n")
cat(paste(notes[seq_test], collapse = " "), "\n")


#Task 8
generate_music<-function(start, M, notes,w=NULL){
  mlag <- ncol(M) - 1
  if (is.null(w)) {
    w <- rep(1, mlag)}
  end_token <- which(notes == "||")
  bar_token <- which(notes == "|")
  gen_seq <- start
  bar_count <- if (start == bar_token) 1 else 0
  while (length(gen_seq) < 200 && bar_count < 8 && tail(gen_seq, 1) != end_token) {
    nxt <- next_note(key = gen_seq, M = M, w = w)
    gen_seq <- c(gen_seq, nxt)
    if (nxt == bar_token) {
      bar_count <- bar_count + 1
    }
  }
  gen_notes <- notes[gen_seq]
  sim_music <- paste(gen_notes, collapse = " ")
  return(sim_music)
}

#test
set.seed(42)  
sim_music <- generate_music(start = 20, M = M, notes = notes, w = c(1, 1, 1, 1))
cat("=== Task 8===\n")
print(sim_music)
print(sim_music)