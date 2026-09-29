# Group members:
# Contributions:
# setwd("put/your/local/repo/location/here") ## SUBMIT: keep this commented out
preprocess <- function(rawdata) {
  # (a) 删除前3行（文件头部说明文字）
  if(length(rawdata) >3){
    rawdata <- rawdata[-(1:3)]
  }
  # (b) 删除以#开头的注释行；F#这种音符不受影响
  remove_rows <- grepl("^#", rawdata)
  rawdata <- rawdata[!remove_rows]
  
  # (c)识别空白行，空白行代表乐曲结束，给上一行末尾加上 " ||"
  # 去掉慢速for循环，改用向量运算提速
  is_blank <- trimws(rawdata) == ""
  keep_line <- !is_blank
  out_lines <- rawdata[keep_line]
  blank_idx <- which(is_blank)
  target_idx <- blank_idx - 1
  target_idx <- target_idx[target_idx >=1 & target_idx <= length(out_lines)]
  if(length(target_idx) >0){
    out_lines[target_idx] <- paste0(out_lines[target_idx], " ||")
  }
  
  # (d) 普通行末尾加 " |"；已经有 || 的不加
  for (i in seq_along(out_lines)) {
    line_tmp <- out_lines[i]
    if (!grepl("\\|\\|$", line_tmp) && !grepl("\\|$", line_tmp)) {
      out_lines[i] <- paste0(line_tmp, " |")
    }
  }
  
  # (e) 按空白切分，展开成单个token，删掉空字符串
  split_result <- strsplit(out_lines, "\\s+")
  all_tokens <- unlist(split_result)
  clean_tokens <- all_tokens[all_tokens != ""]
  
  return(clean_tokens)
}

# ========== 主执行代码 ==========
music_raw <- readLines("palestrina.txt")
music_clean <- preprocess(music_raw)
# Task5
notes <- unique(music_clean)
tokens <- match(music_clean, notes)
# 修复：防止找不到||导致end为空
end <- which(notes == "||")
if(length(end)==0){
  warning("没有找到乐曲结束标记 ||")
}
