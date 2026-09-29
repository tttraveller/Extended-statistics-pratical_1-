music_raw <- readLines("palestrina.txt")
length(music_raw)
head(music_raw,50)

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


print("前10个token：")
print(head(music_clean,10))
print("最后10个token：")
print(tail(music_clean,10))
print(paste("总token数量：", length(music_clean)))
print(paste("结束标记||对应的整数id：", end))


make_matrx<-function(tokens,end,mlag = 4){
  n<-length(tokens)
  # 6.(a) Create an (n- mlag) × (mlag + 1) matrix
  M<-matrix(0,nrow = n-mlag,ncol= mlag+1)
  # slipping windows
  for (i in 0:mlag){
    M[,i+1]<- tokens[(i+1):(i+n-mlag)]
  }
  # 6.(b)define the invalid rows
  invalid_rows<-rowSums(M[,1:mlag,drop = FALSE]==end)>0
  #delect the invalid rows
  M<-M[!invalid_rows, ,drop=FALSE]
  return (M)
}



  
  
