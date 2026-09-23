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
  # 
}


  
  
