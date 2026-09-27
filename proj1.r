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
  out_lines <- c()
  for (line in rawdata) {
    trimmed_line <- trimws(line)
    if (trimmed_line == "") {
      if (length(out_lines) > 0) {
        last_pos <- length(out_lines)
        out_lines[last_pos] <- paste0(out_lines[last_pos], " ||")
      }
    } else {
      out_lines <- c(out_lines, line)
    }
  }
  
  # (d) 普通行末尾加 " |"；已经有 || 的不加
  for (i in seq_along(out_lines)) {
    if (!grepl("\\|\\|$", out_lines[i])) {
      out_lines[i] <- paste0(out_lines[i], " |")
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
end <- which(notes == "||")



# 可选：校验打印，提交时可以保留或者删掉
print("前10个token：")
print(head(music_clean,10))
print("最后10个token：")
print(tail(music_clean,10))
print(paste("总token数量：", length(music_clean)))
print(paste("结束标记||对应的整数id：", end))

