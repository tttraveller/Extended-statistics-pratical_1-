################################################################################
# Practical 1: Music Simulator (Tasks 3 - 7)
################################################################################

# ==============================================================================
# Task 3: Read Raw Data File
# ==============================================================================
# Read all lines from the raw text score file
music_raw <- readLines("palestrina.txt")

# ==============================================================================
# Task 4: Data Preprocessing Function preprocess
# ==============================================================================
preprocess <- function(music) {
  # (a) Remove the first 3 lines containing composer and source information
  music <- music[-c(1:3)]
  
  # (b) Remove metadata comment lines starting with '#'
  # Note: sharps within note encodings (e.g., F#4q) are not at line starts
  music <- music[!grepl("^#", music)]
  
  # Strip leading and trailing empty lines to avoid misidentifying boundaries
  while (length(music) > 0 && music[1] == "") {
    music <- music[-1]
  }
  while (length(music) > 0 && music[length(music)] == "") {
    music <- music[-length(music)]
  }
  
  # (c) & (d) Add single bar-lines "|" and double bar-lines "||"
  # Locate empty lines separating different pieces of music
  blank_idx <- which(music == "")
  
  # The last measure of each piece is the line immediately preceding a blank line
  # The final line of the non-empty vector represents the end of the last piece
  piece_end_idx <- c(blank_idx - 1, length(music))
  
  is_piece_end <- logical(length(music))
  is_piece_end[piece_end_idx] <- TRUE
  
  # Append " ||" to the final measure of each piece, and " |" to all other measures
  music[is_piece_end]  <- paste0(music[is_piece_end], " ||")
  music[!is_piece_end] <- paste0(music[!is_piece_end], " |")
  
  # Discard empty lines that only contain whitespace and the appended bar-line
  music <- music[music != " |"]
  
  # (e) Split measures into individual tokens (notes, rests, and bar-line symbols)
  all_notes <- unlist(strsplit(music, "\\s+"))
  all_notes <- all_notes[all_notes != ""]
  
  return(all_notes)
}

# Execute preprocessing
cat("Preprocessing data...\n")
music_clean <- preprocess(music_raw)
cat("Preprocessing complete! Total tokens:", length(music_clean), "\n")
cat("Verification of last 10 tokens (must match course handout benchmark):\n")
print(tail(music_clean, 10))
cat("\n")

# ==============================================================================
# Task 5: Construct Vocabulary (notes) and Token Sequence (tokens)
# ==============================================================================
# (a) Vector of unique note encodings acting as the dictionary
notes <- unique(music_clean)

# (b) Convert raw string tokens into integer indices corresponding to 'notes'
tokens <- match(music_clean, notes)
cat("Vocabulary size (unique notes):", length(notes), "\n\n")

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

# ==============================================================================
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
