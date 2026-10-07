# ==============================================================================
# Practical 1: Music Simulator
# Course: Statistical Computing
# 
# Group Members:
# - Member 1 Name (UUN: sXXXXXXX)
# - Member 2 Name (UUN: sXXXXXXX)
# - Member 3 Name (UUN: sXXXXXXX)
#
# Contribution Statement:
# - Member 1 (~33.3%): Implemented preprocess function and tokenization setup.
# - Member 2 (~33.3%): Implemented make_matrix and next_note sampling engine.
# - Member 3 (~33.3%): Implemented generate_music, testing suite, and code review.
# All members contributed equally to debugging, documentation, and Git version control.
#
# Overview:
# This script constructs a variable-order Markov chain model to simulate 
# Renaissance choral music based on Palestrina's compositions. It processes 
# scientific pitch notation text, builds lag transition matrices, and generates
# new musical pieces using empirical mixture probability sampling.
# ==============================================================================

# Ensure working directory is set appropriately before execution
# setwd("put/your/local/repo/location/here") ## Commented out for submission


# ==============================================================================
# Task 4: Preprocessing Raw Music Data
# ==============================================================================
# preprocess:
#   Inputs:
#     music - Character vector containing lines from the raw text file.
#   Outputs:
#     A character vector of individual note, rest, and bar-line tokens.
#   Purpose:
#     Clean header/comment metadata, delineate musical pieces with "||", 
#     mark measure boundaries with "|", and tokenize by whitespace.
# ==============================================================================
preprocess <- function(music) {
  # (a) Remove first 3 header lines containing metadata
  if (length(music) > 3) {
    music <- music[-(1:3)]
  }
  
  # (b) Remove comment lines beginning with '#'
  music <- music[!grepl("^#", music)]
  
  # (c) Identify empty lines separating pieces and append " ||" to the last line of each piece
  is_blank <- (trimws(music) == "")
  target_idx <- which(is_blank) - 1
  target_idx <- target_idx[target_idx >= 1 & !is_blank[target_idx]]
  music[target_idx] <- paste0(music[target_idx], " ||")
  
  # Filter out blank rows
  out_lines <- music[!is_blank]
  
  # Ensure the very last line of the entire file ends with " ||"
  n_lines <- length(out_lines)
  if (n_lines > 0 && !grepl("\\|\\|$", out_lines[n_lines])) {
    out_lines[n_lines] <- paste0(out_lines[n_lines], " ||")
  }
  
  # (d) Vectorized addition of " |" to measures not ending in "||" or "|"
  needs_bar <- !grepl("\\|\\|$", out_lines) & !grepl("\\|$", out_lines)
  out_lines[needs_bar] <- paste0(out_lines[needs_bar], " |")
  
  # (e) Split into individual tokens by whitespace
  split_tokens <- unlist(strsplit(out_lines, "\\s+"))
  clean_tokens <- split_tokens[split_tokens != ""]
  
  return(clean_tokens)
}


# ==============================================================================
# Task 6: Construct Lagged Transition Matrix
# ==============================================================================
# make_matrix:
#   Inputs:
#     tokens - Integer vector of token IDs representing the cleaned corpus.
#     end    - Integer token ID corresponding to the double bar-line "||".
#     mlag   - Maximum order/lag of the model (default: 4).
#   Outputs:
#     An (N x (mlag + 1)) matrix where each row is an observed sequence of 
#     length mlag + 1, excluding sequences crossing piece boundaries.
# ==============================================================================
make_matrix <- function(tokens, end, mlag = 4) {
  n <- length(tokens)
  
  # (a) Construct lagged sequences matrix
  M <- matrix(0L, nrow = n - mlag, ncol = mlag + 1)
  for (j in 1:(mlag + 1)) {
    M[, j] <- tokens[j:(n - mlag + j - 1)]
  }
  
  # (b) Discard rows containing "||" within the prefix (first mlag columns)
  has_end <- rowSums(M[, 1:mlag, drop = FALSE] == end) > 0
  if (any(has_end)) {
    M <- M[!has_end, , drop = FALSE]
  }
  
  return(M)
}


# ==============================================================================
# Task 7: Next-Note Sampling Function
# ==============================================================================
# next_note:
#   Inputs:
#     key - Integer vector of preceding token IDs.
#     M   - Transition matrix constructed by make_matrix.
#     w   - Mixture weight vector for lag models (default: rep(1, mlag)).
#   Outputs:
#     A single integer token ID sampled according to the mixture model.
# ==============================================================================
next_note <- function(key, M, w = NULL) {
  mlag <- ncol(M) - 1
  if (is.null(w)) {
    w <- rep(1, mlag)
  }
  
  # Truncate key to trailing mlag tokens if it exceeds mlag
  if (length(key) > mlag) {
    key <- tail(key, mlag)
  }
  
  k_max <- length(key)
  cand_tokens <- integer(0)
  cand_probs <- numeric(0)
  
  # Iterate from order 1 up to k_max
  for (k in 1:k_max) {
    # Indexing: w[1] matches lag mlag (highest order), w[mlag] matches lag 1
    w_idx <- mlag - k + 1
    weight_k <- w[w_idx]
    if (weight_k <= 0) next
    
    mc <- mlag - k + 1
    sub_key <- key[(k_max - k + 1):k_max]
    
    # Vectorized matching recommended in practical sheet
    ii <- colSums(!(t(M[, mc:mlag, drop = FALSE]) == sub_key))
    matched_idx <- which(ii == 0)
    n_u <- length(matched_idx)
    
    if (n_u > 0) {
      cand_tokens <- c(cand_tokens, M[matched_idx, mlag + 1])
      cand_probs <- c(cand_probs, rep(weight_k / n_u, n_u))
    }
  }
  
  total_prob <- sum(cand_probs)
  
  # Decision rules (ensuring at most one call to sample)
  if (length(cand_tokens) == 0 || total_prob == 0) {
    # Case 1: Fallback marginal sampling
    return(sample(M[, mlag + 1], size = 1))
  } else if (length(cand_tokens) == 1) {
    # Case 2: Exactly one candidate -> return directly
    return(cand_tokens[1])
  } else {
    # Case 3: Multiple candidates -> mixture probability sampling
    return(sample(cand_tokens, size = 1, prob = cand_probs))
  }
}


# ==============================================================================
# Task 8: Music Generation Pipeline
# ==============================================================================
# generate_music:
#   Inputs:
#     start - Token ID of the starting note.
#     M     - Transition matrix from make_matrix.
#     notes - Vector of unique note characters.
#     w     - Mixture weights (default: rep(1, mlag)).
#   Outputs:
#     A single character string representing the generated music.
# ==============================================================================
generate_music <- function(start, M, notes, w = NULL) {
  mlag <- ncol(M) - 1
  if (is.null(w)) {
    w <- rep(1, mlag)
  }
  
  end_token <- which(notes == "||")
  bar_token <- which(notes == "|")
  
  gen_seq <- start
  bar_count <- if (start == bar_token) 1 else 0
  
  # Loop until reaching ||, 8 measures, or 200 notes
  while (length(gen_seq) < 200 && bar_count < 8 && tail(gen_seq, 1) != end_token) {
    nxt <- next_note(key = gen_seq, M = M, w = w)
    gen_seq <- c(gen_seq, nxt)
    if (nxt == bar_token) {
      bar_count <- bar_count + 1
    }
  }
  
  return(paste(notes[gen_seq], collapse = " "))
}


# ==============================================================================
# Execution & Task 8 Analysis
# ==============================================================================
# 1. Load and preprocess corpus
music_raw <- readLines("palestrina.txt")
music_clean <- preprocess(music = music_raw)

# 2. Vocabulary extraction and tokenization
notes <- unique(music_clean)
tokens <- match(music_clean, notes)
end_token_id <- which(notes == "||")

# 3. Construct transition matrix (default mlag = 4)
M <- make_matrix(tokens = tokens, end = end_token_id, mlag = 4)

# 4. Generate music across different weighting configurations
set.seed(42)
start_note_id <- 20 # Arbitrary valid initial note token

# w1 = (1, 1, 1, 1): Full mixture model across all orders (1 to 4)
# Captures both local melodic context and flexible short-term transitions.
sim_w1 <- generate_music(start = start_note_id, M = M, notes = notes, w = c(1, 1, 1, 1))

# w2 = (1, 0, 0, 0): Pure 4th-order lag model
# High local structural coherence, but prone to data sparsity / fallback sampling.
sim_w2 <- generate_music(start = start_note_id, M = M, notes = notes, w = c(1, 0, 0, 0))

# w3 = (0, 0, 0, 1): Pure 1st-order simple Markov model
# Highly fluid note transitions with large support, but lacks multi-bar melody structure.
sim_w3 <- generate_music(start = start_note_id, M = M, notes = notes, w = c(0, 0, 0, 1))