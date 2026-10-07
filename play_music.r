require("audio")
library(tuneR) 

# -------------------------------------------------------------------------


# -------------------------------------------------------------------------


# -------------------------------------------------------------------------


# -------------------------------------------------------------------------


# -------------------------------------------------------------------------


note_to_hz <- function(note, accidental, octave){
  note_data <- c(A = 12, B = 14, C = 3, D = 5, E = 7, F = 8, G = 10, R = NA)
  accidental_data <- c("bb" = -2, "b" = -1, "#" = 1, "##" = 2)
  note_value <- note_data[toupper(note)]
  if(accidental == ""){ acc_val <- 0 } else { acc_val <- accidental_data[accidental] }
  octave_value <- as.numeric(octave)
  note_num = note_value + acc_val + octave_value*12
  freq = 440 * 2^((note_num - 60)/12)
  return(unname(freq))
}

duration_to_num <- function(duration){
  duration_data <- c(b = 8, w = 4, h = 2, q = 1, e = 0.5, s = 0.25)
  if(is.na(duration) || duration == "") return(0)
  dur_char <- substr(duration,1,1)
  dur_value <- duration_data[tolower(dur_char)]
  dot_count <- nchar(gsub("[^.]","",duration))
  dur_multiple <- 1
  if(dot_count >=1) dur_multiple <- dur_multiple + 0.5
  if(dot_count >=2) dur_multiple <- dur_multiple + 0.25
  return(unname(dur_value * dur_multiple))
}

# 一阶低通（截止频率调高，不闷）
lowpass_filter <- function(x, cutoff_hz, sample_rate){
  rc <- 1/(2*pi*cutoff_hz)
  dt <- 1/sample_rate
  alpha <- dt/(rc+dt)
  y <- numeric(length(x)); y[1] <- x[1]
  for(i in 2:length(x)) y[i] <- y[i-1] + alpha*(x[i]-y[i-1])
  return(y)
}

# 大混响：5个延迟，空间感强
simple_reverb <- function(x, sample_rate){
  out <- x
  delays <- c(0.025, 0.045, 0.070, 0.100, 0.140)
  gains  <- c(0.30, 0.25, 0.20, 0.15, 0.10)
  for(i in seq_along(delays)){
    d <- round(delays[i]*sample_rate)
    echo <- c(rep(0,d), x[1:(length(x)-d)])
    out <- out + gains[i]*echo
  }
  # 再加一点点反射尾音
  tail <- c(rep(0, round(0.20*sample_rate)), x[1:(length(x)-round(0.20*sample_rate))])
  out <- out + 0.05*tail
  return(out)
}

make_note_waveform <- function(freq, duration, bpm, sample_rate = 44100) {
  if(duration <= 0 || is.na(duration)) return(numeric(0))
  time_len <- duration / bpm * 60
  time <- seq(0, time_len, by = 1/sample_rate)
  n_samp <- length(time)
  if(is.na(freq)) return(rep(0, n_samp))
  
  # ==== 随机力度（每次演奏强弱不一样，去单调）====
  velocity <- runif(1, 0.85, 1.15)
  
  # ==== 颤音 ====
  vib_depth <- 0.004
  vib_rate  <- 5.5
  freq_mod <- freq * (1 + vib_depth * sin(2*pi*vib_rate*time))
  phase <- cumsum(freq_mod) / sample_rate * 2*pi
  
  # ==== 泛音：高频权重提高，更亮 ====
  harmonics <- list(
    list(mult=1, w=0.55, decay=1.0),
    list(mult=2, w=0.35, decay=0.7),
    list(mult=3, w=0.22, decay=0.45),
    list(mult=4, w=0.12, decay=0.30),
    list(mult=5, w=0.06, decay=0.18)
  )
  wave <- rep(0, n_samp)
  for(h in harmonics){
    env_h <- exp(-time * h$decay * 2.5)
    wave <- wave + h$w * sin(h$mult * phase) * env_h
  }
  
  # ==== ADSR ====
  attack_time   <- 0.012
  decay_time    <- 0.05
  sustain_level <- 0.65
  release_time  <- 0.12
  attack_samples  <- min(round(attack_time*sample_rate), n_samp)
  decay_samples   <- min(round(decay_time*sample_rate),  n_samp - attack_samples)
  release_samples <- min(round(release_time*sample_rate), n_samp)
  sustain_samples <- n_samp - attack_samples - decay_samples - release_samples
  env <- numeric(n_samp)
  env[1:attack_samples] <- seq(0,1,length.out=attack_samples)
  if(decay_samples>0)
    env[(attack_samples+1):(attack_samples+decay_samples)] <- seq(1,sustain_level,length.out=decay_samples)
  if(sustain_samples>0){
    start_sus <- attack_samples+decay_samples+1
    env[start_sus:(n_samp-release_samples)] <- sustain_level
  }
  if(release_samples>0)
    env[(n_samp-release_samples+1):n_samp] <- seq(sustain_level,0,length.out=release_samples)
  wave <- wave * env * velocity
  
  # ==== 低通放开到 9000，不闷 ====
  wave <- lowpass_filter(wave, cutoff_hz = 9000, sample_rate)
  
  if(max(abs(wave))>0) wave <- wave/max(abs(wave))*0.7
  return(wave)
}

str_to_waveform <- function(input_str, bpm = 120){
  tokens <- strsplit(input_str," ")[[1]]
  tokens <- tokens[tokens!="|" & tokens!="||"]
  tokens <- tokens[nchar(tokens)>0]
  
  parse_one_token <- function(token){
    if(grepl("^R",token)){
      dur_part <- substr(token,2,nchar(token))
      return(list(note="R", acc="", oct="4", dur=dur_part))
    }
    m <- regexec("^([A-G])(#|b|bb)?([0-9]+)([whqes]\\.*)$", token)
    reg_res <- regmatches(token,m)[[1]]
    if(length(reg_res)<5) return(NULL)
    note_name <- reg_res[2]
    acc <- reg_res[3]; if(is.na(acc)) acc <- ""
    return(list(note=note_name, acc=acc, oct=reg_res[4], dur=reg_res[5]))
  }
  
  parsed_list <- lapply(tokens, parse_one_token)
  parsed_list <- Filter(Negate(is.null), parsed_list)
  if(length(parsed_list)==0) return(numeric(0))
  
  hz_vec <- c(); dur_vec <- c()
  for(p in parsed_list){
    if(p$note=="R") hz_vec <- c(hz_vec, NA)
    else hz_vec <- c(hz_vec, note_to_hz(p$note, p$acc, p$oct))
    dur_vec <- c(dur_vec, duration_to_num(p$dur))
  }
  
  wave_all <- c()
  overlap <- round(0.01*44100)
  for(i in seq_along(hz_vec)){
    seg <- make_note_waveform(freq=hz_vec[i], duration=dur_vec[i], bpm=bpm)
    if(length(wave_all)>overlap && length(seg)>overlap){
      wave_all[(length(wave_all)-overlap+1):length(wave_all)] <-
        wave_all[(length(wave_all)-overlap+1):length(wave_all)] + seg[1:overlap]
      wave_all <- c(wave_all, seg[(overlap+1):length(seg)])
    } else {
      wave_all <- c(wave_all, seg)
    }
  }
  
  wave_all <- simple_reverb(wave_all, 44100)
  if(max(abs(wave_all))>0) wave_all <- wave_all/max(abs(wave_all))*0.85
  return(wave_all)
}

# ========== 生成 ==========
set.seed(42)  # 固定随机力度，每次听的效果一致
waveform <- str_to_waveform(sim_music, bpm=120)
waveform_int <- round(waveform * 32767)
wave_obj <- Wave(waveform_int, samp.rate = 44100, bit = 16)
writeWave(wave_obj, "markov_output.wav")
cat("✅ markov_output.wav 已生成（亮音色 + 大混响 + 随机力度）\n")
