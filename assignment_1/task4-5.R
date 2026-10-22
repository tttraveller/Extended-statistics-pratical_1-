# Group members:
# Contributions:
# setwd("put/your/local/repo/location/here") ## SUBMIT: keep this commented out
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
