preprocess <- function(rawdata) {
  # (a) Remove first 3 lines of metadata (Composer, Source, etc.)
  if (length(rawdata) > 3) {
    rawdata <- rawdata[-(1:3)]
  }
  
  # (b) Remove full-line comments starting with '#', keep sharps inside notes
  is_comment <- grepl("^\\s*#", rawdata)
  rawdata <- rawdata[!is_comment]
  
  # (c) Mark piece boundaries: append double bar "||" to lines before blank lines
  is_blank <- trimws(rawdata) == ""
  blank_indices <- which(is_blank)
  end_piece_indices <- integer(0)
  
  for (idx in blank_indices) {
    target <- idx - 1
    while (target >= 1 && is_blank[target]) {
      target <- target - 1
    }
    if (target >= 1) {
      end_piece_indices <- c(end_piece_indices, target)
    }
  }
  
  # Mark last non-blank line as piece ending
  last_non_blank <- tail(which(!is_blank), 1)
  if (length(last_non_blank) > 0) {
    end_piece_indices <- unique(c(end_piece_indices, last_non_blank))
  }
  
  if (length(end_piece_indices) > 0) {
    rawdata[end_piece_indices] <- paste(trimws(rawdata[end_piece_indices]), "||")
  }
  
  # Remove blank lines
  non_blank_data <- rawdata[!is_blank]
  
  # (d) Append single bar "|" to normal lines (skip if already ends with | or ||)
  out_lines <- character(length(non_blank_data))
  for (i in seq_along(non_blank_data)) {
    cur_line <- trimws(non_blank_data[i])
    if (endsWith(cur_line, "||") || endsWith(cur_line, "|")) {
      out_lines[i] <- cur_line
    } else {
      out_lines[i] <- paste(cur_line, "|")
    }
  }
  
  # (e) Split into individual token vector by whitespace
  split_tokens <- unlist(strsplit(out_lines, "\\s+"))
  clean_tokens <- split_tokens[split_tokens != ""]
  
  return(clean_tokens)
}

# ------------------------------------------------------------------------------
# Task 6: Construct Lagged Transition Matrix
# ------------------------------------------------------------------------------
# Inputs:
#   tokens - integer vector of token IDs
#   end    - integer ID corresponding to the double bar-line token "||"
#   mlag   - maximum order/lag for Markov dependencies (default = 4)
# Output:
#   M      - integer matrix of dimensions (n - mlag) x (mlag + 1)
make_matrix <- function(tokens, end, mlag = 4) {
  n <- length(tokens)
  if (n <= mlag) {
    stop("Length of token vector must be greater than mlag.")
  }
  
  # (a) Build sliding window sequence matrix M
  M <- matrix(0L, nrow = n - mlag, ncol = mlag + 1)
  for (j in 1:(mlag + 1)) {
    M[, j] <- tokens[j:(n - mlag + j - 1)]
  }
  
  # (b) Remove rows with end token to avoid cross-piece transitions
  if (length(end) > 0 && !is.na(end)) {
    rows_with_end <- which(rowSums(M[, 1:mlag, drop = FALSE] == end) > 0)
    if (length(rows_with_end) > 0) {
      M <- M[-rows_with_end, , drop = FALSE]
    }
  }
  
  return(M)
}

# ------------------------------------------------------------------------------
# Helper: Semitone Height Calculation (for leap suppression)
# ------------------------------------------------------------------------------
get_semitone <- function(note_str) {
  if (is.na(note_str) || note_str %in% c("|", "||") || startsWith(note_str, "R")) {
    return(NA_real_)
  }
  p <- sub("[a-z.]+$", "", note_str)
  step <- substr(p, 1, 1)
  acc <- if (grepl("#", p)) 1 else if (grepl("b", p)) -1 else 0
  oct <- as.numeric(sub("^[A-G][#b]?", "", p))
  
  base_steps <- c(C = 0, D = 2, E = 4, F = 5, G = 7, A = 9, B = 11)
  if (!step %in% names(base_steps) || is.na(oct)) return(NA_real_)
  return((oct + 1) * 12 + base_steps[step] + acc)
}

# ------------------------------------------------------------------------------
# Task 7: Next-Note Mixture Sampling Function (Rhythm & Melody Enhanced)
# ------------------------------------------------------------------------------
# Inputs:
#   key - integer vector representing recent token history
#   M   - transition matrix produced by make_matrix
#   w   - mixture weights vector of length mlag (default = rep(1, mlag))
# Output:
#   A single integer token representing the predicted next note.
next_note <- function(key, M, w = NULL) {
  mlag <- ncol(M) - 1
  
  if (is.null(w)) {
    w <- rep(1, mlag)
  }
  
  if (length(key) > mlag) {
    key <- tail(key, mlag)
  }
  
  k_max <- length(key)
  cand_tokens <- integer(0)
  cand_probs <- numeric(0)
  
  # Iterate over orders for pattern matching
  for (k in 1:k_max) {
    if (w[k] <= 0) next
    
    sub_key <- key[(k_max - k + 1):k_max]
    mc <- mlag - k + 1
    
    # Vectorized pattern matching
    ii <- colSums(!(t(M[, mc:mlag, drop = FALSE]) == sub_key))
    matched_idx <- which(ii == 0)
    n_u <- length(matched_idx)
    
    if (n_u > 0) {
      cand_tokens <- c(cand_tokens, M[matched_idx, mlag + 1])
      cand_probs <- c(cand_probs, rep(w[k] / n_u, n_u))
    }
  }
  
  total_prob <- sum(cand_probs)
  
  # Sampling logic (at most one sample call)
  if (length(cand_tokens) == 0 || total_prob <= 0) {
    return(sample(M[, mlag + 1], size = 1))
  } else if (length(cand_tokens) == 1) {
    return(cand_tokens[1])
  } else {
    # Rhythm & melody optimization for multiple candidates
    if (exists("notes", envir = parent.frame())) {
      dict_notes <- get("notes", envir = parent.frame())
      cand_note_chars <- dict_notes[cand_tokens]
      last_token <- tail(key, 1)
      last_note_char <- dict_notes[last_token]
      
      penalties <- rep(1, length(cand_tokens))
      
      # 1. Rhythm adjustment: boost medium-short durations, suppress long durations
      dur_chars <- sub(".*([a-z])\\.?$", "\\1", cand_note_chars)
      penalties[dur_chars == "b"] <- penalties[dur_chars == "b"] * 0.05
      penalties[dur_chars == "w"] <- penalties[dur_chars == "w"] * 0.30
      penalties[dur_chars == "h"] <- penalties[dur_chars == "h"] * 1.00
      penalties[dur_chars == "q"] <- penalties[dur_chars == "q"] * 2.50
      penalties[dur_chars == "e"] <- penalties[dur_chars == "e"] * 1.80
      
      # 2. Suppress long rests to avoid silence
      is_rest <- startsWith(cand_note_chars, "R")
      penalties[is_rest] <- penalties[is_rest] * 0.05
      
      # 3. Regulate note density per bar
      bars_pos <- which(dict_notes[key] == "|")
      notes_in_cur_bar <- if (length(bars_pos) == 0) length(key) else length(key) - tail(bars_pos, 1)
      is_bar_line <- cand_note_chars == "|"
      if (notes_in_cur_bar < 3) {
        penalties[is_bar_line] <- penalties[is_bar_line] * 0.1
      } else if (notes_in_cur_bar >= 4) {
        penalties[is_bar_line] <- penalties[is_bar_line] * 2.5
      }
      
      # 4. Pitch smoothing: penalize leaps over a perfect fifth (7 semitones)
      last_semi <- get_semitone(last_note_char)
      cand_semis <- vapply(cand_note_chars, get_semitone, numeric(1))
      if (!is.na(last_semi)) {
        dist <- abs(cand_semis - last_semi)
        large_leap <- !is.na(dist) & (dist > 7)
        penalties[large_leap] <- penalties[large_leap] * 0.2
      }
      
      cand_probs <- cand_probs * penalties
      if (sum(cand_probs) <= 0) cand_probs <- rep(1, length(cand_tokens))
    }
    
    return(sample(cand_tokens, size = 1, prob = cand_probs))
  }
}

# ------------------------------------------------------------------------------
# Task 8: Autoregressive Music Simulation Function
# ------------------------------------------------------------------------------
# Inputs:
#   start - single integer token ID (or note character string) to initiate simulation
#   M     - transition matrix produced by make_matrix
#   notes - dictionary character vector of unique notes
#   w     - mixture weights vector
# Output:
#   sim_music - single character string representing the generated music score
generate_music <- function(start, M, notes, w = NULL) {
  mlag <- ncol(M) - 1
  if (is.null(w)) {
    # Default: bias towards higher-order memory for melodic coherence
    w <- c(0.1, 0.3, 0.6, 1.0)
  }
  
  if (is.character(start)) {
    start <- match(start, notes)
    if (is.na(start)) {
      stop("Specified starting note not found in the notes dictionary.")
    }
  }
  
  end_token <- which(notes == "||")
  bar_token <- which(notes == "|")
  
  gen_seq <- start
  bar_count <- if (start == bar_token) 1 else 0
  
  # Generation stops when any condition is met:
  # 1. Sequence reaches 200 tokens
  # 2. 8 bar lines are encountered
  # 3. Double bar line is encountered
  while (length(gen_seq) < 200 && bar_count < 8 && tail(gen_seq, 1) != end_token) {
    nxt <- next_note(key = gen_seq, M = M, w = w)
    gen_seq <- c(gen_seq, nxt)
    if (nxt == bar_token) {
      bar_count <- bar_count + 1
    }
  }
  
  sim_music <- paste(notes[gen_seq], collapse = " ")
  return(sim_music)
}

# ==============================================================================
# Main Execution Pipeline
# ==============================================================================
# Task 3: Read raw score file
music_raw <- readLines("palestrina.txt")
# Task 4: Data preprocessing
music_clean <- preprocess(music_raw)
cat("Preprocessing complete. Total tokens:", length(music_clean), "\n")
cat("Tail 10 tokens:\n")
print(tail(music_clean, n = 10))
# Task 5: Build lexicon and integer token sequence
notes <- unique(music_clean)
tokens <- match(music_clean, notes)
end_token_id <- which(notes == "||")
cat("Lexicon size (unique notes):", length(notes), "\n\n")
# Task 6: Build transition matrix M (default mlag = 4)
cat("Building transition matrix M (mlag = 4)...\n")
M <- make_matrix(tokens = tokens, end = end_token_id, mlag = 4)
cat("Matrix M built successfully. Dimensions:", nrow(M), "x", ncol(M), "\n\n")
# Task 7: Verification test
cat("=== Task 7: next_note Verification Tests ===\n")
set.seed(42)
test_prefix_tokens <- match(c("G4w", "D5w"), notes)
pred_token <- next_note(key = test_prefix_tokens, M = M)
cat("Prefix: G4w D5w -> Predicted next note:", notes[pred_token], "(Token ID:", pred_token, ")\n\n")
# Task 8: Music generation and weight comparison demo
cat("=== Task 8: Simulated Music (Rhythm-Optimized) ===\n")
set.seed(42)
# Start with mid-range main note (D4h), avoid opening rest
start_token <- match("D4h", notes)
if (is.na(start_token)) start_token <- 20
sim_music <- generate_music(start = start_token, M = M, notes = notes, w = c(0.1, 0.3, 0.6, 1.0))
print(sim_music)

# Weight comparison tests:
# w2 = c(1, 0, 0, 0): Pure 1st-order Markov chain, high randomness, no rhythmic memory
# w3 = c(0, 0, 0, 1): Pure 4th-order Markov chain, strictly preserves 4-note phrases
cat("\n=== Testing w2 = c(1, 0, 0, 0) (1st-order Markov) ===\n")
set.seed(42)
print(generate_music(start = start_token, M = M, notes = notes, w = c(1, 0, 0, 0)))
cat("\n=== Testing w3 = c(0, 0, 0, 1) (4th-order Markov) ===\n")
set.seed(42)
print(generate_music(start = start_token, M = M, notes = notes, w = c(0, 0, 0, 1)))

if (file.exists("play_music.r")) {
  source("play_music.r")
  # Play generated score at medium-high rhythmic speed (bpm = 180~200)
  play(str_to_waveform(sim_music, bpm = 180))
}
source("play_music.r")
play(str_to_waveform(sim_music, bpm = 180))

