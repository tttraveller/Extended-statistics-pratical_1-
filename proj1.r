music_raw <- readLines("palestrina.txt")
length(music_raw)
head(music_raw,50)

preprocess<-function(music){
  # (a) Remove the first three lines containing details about the composer and the data source.
  music<-music[-(1:3)]
  # (b) All of the lines beginning with # provides supplementary information about the music, such as the
  # voice, time signature, key signature, etc. These lines need to be removed.
  music<-music[which(substr(music,1,1) != '#')]
  return(music)
}
make_matrx<-function(tokens,end,malg = 4){
  n<-lenth(tokens)
  # 6.(a) Create an (n- mlag) × (mlag + 1) matrix
  M<-matrix(0,nrow = n-mlag,ncol= mlag+1)
  # slipping windows
  for (i in 0:mlag){
    M[,i+1]<- tokens[i+1,i+1+nlag]
  }
  # 6.(b)define the invalid rows
  invalid_rows<-rowsum[M[1:mlag,drop = FALSE]==end]>0
  #delect the invalid rows
  M<-M[!invalid_rows, ,drop<-FALSE]
  return (M)
}


  
  
