#Group members:
# Sicheng Liang (s2956212)
# Lianjie Wu (s2970562)
# Josie Zhuo (s2959889)
# Contribution Statement:
# - Josie Zhuo (33.3%): Implemented preprocess function and tokenization setup(mainly focused on task4-5).
# - Sicheng Liang (33.3%): Implemented make_matrix and generate_music, testing suite
# - Lianjie Wu (33.3%): Implemented next_note sampling engine and code review.
# All members contributed equally to debugging, documentation, code review and Git version control.
#
# Practical 1: Markov Music Simulator
#
# This script implements a simple Markov-chain prediction model which is trained on Giovanni Palestrina's choral scores. It preprocesses the raw text score, builds a lagged transition matrix, and then generates simulated music by sampling from the model.
#
# Functions defined:
#   preprocess()    - cleans raw score text into a flat token vector
#   make_matrix()   - constructs an (n-mlag) x (mlag+1) matrix of token sequences
#   next_note()     - samples the next note given a current context (key)
#   generate_music()- simulates a full piece of music up to stopping conditions


#task4
# preprocess: clean raw score text into a flat token vector
# input:
#   music: raw character vector from "palestrina.txt"
# output:
#   clean_tokens: character vector of notes, rests, and bar-lines
# steps:
#   (a) remove first 3 header lines
#   (b) remove comment lines starting with '#' (keep notes like F#5e)
#   (c) append "||" to the last line of each piece
#   (d) append "|" to the end of each measure
#   (e) split each line into individual tokens by whitespace
preprocess <- function(music) {
  rawdata <- music
  # (a) Remove the first 3 header lines
  if(length(rawdata) > 3){
    rawdata <- rawdata[-(1:3)]
  }
  # (b) Remove full comment lines starting with '#' keep notes like 'F#5e'
  remove_rows <- grepl("^#", rawdata)  #define the line which should be removed
  rawdata <- rawdata[!remove_rows]
  
  # (c) Detect blank lines (end of piece), append " ||" to the previous line
  is_blank <- (trimws(rawdata) == "")
  target_idx <- which(is_blank) - 1
  # Ensure target indices are valid (positive) and ignore the consecutive blank lines.
  target_idx <- target_idx[target_idx >= 1 & !is_blank[target_idx]]
  #append double '||' before removing blank lines
  rawdata[target_idx] <- paste0(rawdata[target_idx], " ||")
  out_lines <- rawdata[!is_blank]
  
  
  # (d) Append " |" to normal lines; skip lines already ending with || or |
  for (i in seq_along(out_lines)) {
    line_tmp <- out_lines[i]
    if (!grepl("\\|\\|$", line_tmp) && !grepl("\\|$", line_tmp)) {
      out_lines[i] <- paste0(line_tmp, " |")
    }
  }
  
  # (e) Split into individual tokens by whitespace, flatten into 1D vector, and remove empty strings 
  split_result <- strsplit(out_lines, "\\s+")
  all_tokens <- unlist(split_result)
  clean_tokens <- all_tokens[all_tokens != ""]
  
  return(clean_tokens)
}
# ========== Main execution ==========
music_raw <- readLines("palestrina.txt")
music_clean <- preprocess(music = music_raw)

# Task5 Markov chain setup
notes <- unique(music_clean)
tokens <- match(music_clean, notes)
# Safety check for end marker ||
end <- which(notes == "||")
if(length(end)==0){
  warning("End marker || not found in tokens")
}

# Task 6: Construct Lagged Transition Matrix make_matrix
#
# make_matrix: build a matrix of adjacent token sequences
# input:
#   tokens: integer vector of note IDs
#   end: integer, token ID for "||"
#   mlag: integer, maximum lag (default 4)
# output:
#   M: integer matrix with (n-mlag) rows and (mlag+1) columns.
#      Each row is one sequence of mlag+1 adjacent tokens.
#      Rows with "||" in the first mlag columns are removed.
make_matrix <- function(tokens, end, mlag = 4) {
  n <- length(tokens)
  
  # (a) Create an (n - mlag) x (mlag + 1) matrix for adjacent token sequences
  M <- matrix(0L, nrow = n - mlag, ncol = mlag + 1) #'0L' represent the 'integer' category of tokens and save memory.
  for (j in 1:(mlag + 1)) {
    M[, j] <- tokens[j:(n - mlag + j - 1)]
  }
  
  # (b) Discard rows where the double bar-line token '||' appears in the first mlag columns
  # Sequences spanning across two separate musical pieces must not be used
  rows_with_end <- which(rowSums(M[, 1:mlag, drop = FALSE] == end) > 0)
  if (length(rows_with_end) > 0) {
    M <- M[-rows_with_end, , drop = FALSE]
  }# use 'if' to avoid the dataset do not have matching rows
  
  return(M)
}

cat("Constructing transition matrix M (mlag = 4)...\n")
end_token_id <- which(notes == "||")
M <- make_matrix(tokens = tokens, end = end_token_id, mlag = 4)
cat("Matrix M successfully built! Dimensions:", dim(M)[1], "x", dim(M)[2], "\n\n")

# Task 7: Next-Note Sampling Function next_note
# next_note: sample the next note token given a current context
# input:
#   key: integer vector of recent token IDs
#   M: transition matrix from make_matrix()
#   w: numeric vector of mixture weights 
# output:
#   A single integer token ID for the next note.
# details:
#   Tries matching at full lag order, then reduces by one note at a time.
#   If no match is found, falls back to marginal sampling.
next_note <- function(key, M, w = NULL) {
  mlag <- ncol(M) - 1
  
  #Default mixture weights to a vector of ones of length mlag
  if (is.null(w)) {
    w <- rep(1, mlag)
  }
  
  #only keep the trailing mlag tokens if key > mlag
  if (length(key) > mlag) {
    key <- tail(key, mlag)
  }
  
  k_max <- length(key)
  cand_tokens <- integer(0)
  cand_probs <- numeric(0)
  
  #Check matched for different order models(k from 1 to k_max)
  for (k in 1:k_max) {
    w_idx <- mlag - k + 1
    weight_k <- w[w_idx]
    
    # Skip matching if mixture weight is non-positive
    if (weight_k <= 0) next
    
    # Extract suffix of length k and select corresponding columns in M
    mc <- mlag - k + 1
    sub_key <- key[(k_max - k + 1):k_max]
    
    # find the matching rows by the formula
    ii <- colSums(!(t(M[, mc:mlag, drop = FALSE]) == sub_key))
    matched_idx <- which(ii == 0)
    n_u <- length(matched_idx)
    
    #record candidate next notes and weights
    if (n_u > 0) {
      cand_tokens <- c(cand_tokens, M[matched_idx, mlag + 1])
      cand_probs <- c(cand_probs, rep(weight_k / n_u, n_u))
    }
  }
  total_prob <- sum(cand_probs)
  if (length(cand_tokens) == 0 || total_prob == 0) {
    # Case 1 No match found:sample randomly from all notes in columns (mlag+1)
    return(sample(M[, mlag + 1], size = 1))
  } else if (length(cand_tokens) == 1) {
    # Case 2 Exactly one candidate: return directly without invoking sample()
    return(cand_tokens[1])
  } else {
    # Case 3 Multiple candidates: sample according to empirical mixture probabilities(mixture weights)
    return(sample(cand_tokens, size = 1, prob = cand_probs))
  }
}

# Task 7 Verification Tests

cat("=== Testing next_note function ===\n")

# Test 1: Predict next note from a short prompt
test_prefix_tokens <- match(c("G4w", "D5w"), notes)
set.seed(42) # set the random seed to make sure the results could be replicated
pred_token <- next_note(key = test_prefix_tokens, M = M)
cat("Input sequence: G4w D5w\n")
cat("Predicted next Token ID:", pred_token, "| Corresponding note:", notes[pred_token], "\n\n")

# Test 2: Non-existent prefix to check the random background sampling
dummy_key <- c(99999L, 99999L)
fallback_token <- next_note(key = dummy_key, M = M)
cat("Fallback sampled note for non-existent prefix:", notes[fallback_token], "\n\n")

# Test 3: Sequential auto-regressive generation of 10 notes
seq_test <- test_prefix_tokens
for (step in 1:10) {
  nxt <- next_note(seq_test, M)
  seq_test <- c(seq_test, nxt)
}
cat("10-note sequence generation demo:\n")
cat(paste(notes[seq_test], collapse = " "), "\n")


#Task 8
# generate_music: simulate a full piece of music
# input: 
#   start: Token ID into initiate simulation
#   M: Transition matrix
#   notes: Vector of note names
#   W: Mixture weights
# output:
#   sim_music: A single character string representing the simulated music
# details:
#   Repeatedly calls next_note() until one of the stopping conditions is met:
#   (1) "||" is generated, (2) 8 bars ("|") have been generated,
#   or (3) 200 notes total.

generate_music<-function(start, M, notes,w=NULL){
  mlag <- ncol(M) - 1
  if (is.null(w)) {
    w <- rep(1, mlag)}
  end_token <- which(notes == "||")
  bar_token <- which(notes == "|")
  gen_seq <- start
  bar_count <- if (start == bar_token) 1 else 0
  
  # Continue generation until reaching '||', 8 bars, or 200 tokens
  while (length(gen_seq) < 200 && bar_count < 8 && tail(gen_seq, 1) != end_token) {
    nxt <- next_note(key = gen_seq, M = M, w = w)
    gen_seq <- c(gen_seq, nxt)
    if (nxt == bar_token) {
      bar_count <- bar_count + 1
    }
  }
  
  # Decode numeric tokens into text string
  gen_notes <- notes[gen_seq]
  sim_music <- paste(gen_notes, collapse = " ")
  return(sim_music)
}



# Model evaluation and weight comparison analysis
# compare different mixture weights:
# 1) w1=c(1,1,1,1)
# completely use transitions across all orders from 1 to 4
# It is the most balanced musical output. 
# The 4-note history maintains longer melodic patterns and musical features. The lower lags guarantee smooth transitions and avoid sudden interruption.
#
# 2) w2=c(1,0,0,0)
# use the 4 most recent notes 
# some parts of it are same with original music piece.
# However, because 4-note sequence frequency can be very low, it frequently fails to match and falls back to random marginal draws. It causes occasional sharp note jumps.
#
# 3) w3=c(0,0,0,1)
# completely rely on the single immediate preceding note. 
# note transitions are well-supported by data and practically zero fallback.
# However, it loses the long term memory of bar structure and phrase development, which makes it wander without clear direction.

set.seed(42)
cat("\n=== Task 8: Generating Music with Different Mixture Weights ===\n")

# Test Weight Configuration 1: Full mixture model
sim_w1 <- generate_music(start = 20, M = M, notes = notes, w = c(1, 1, 1, 1))
cat("\n[w1 = (1, 1, 1, 1)] Balanced Mixture:\n", sim_w1, "\n")

# Test Weight Configuration 2: Strict 4th-order lag model
sim_w2 <- generate_music(start = 20, M = M, notes = notes, w = c(1, 0, 0, 0))
cat("\n[w2 = (1, 0, 0, 0)] Pure 4th-Order:\n", sim_w2, "\n")

# Test Weight Configuration 3: Pure 1st-order simple Markov model
sim_w3 <- generate_music(start = 20, M = M, notes = notes, w = c(0, 0, 0, 1))
cat("\n[w3 = (0, 0, 0, 1)] Pure 1st-Order:\n", sim_w3, "\n")