# ============================================================
# detect_body_movements.R
# and detect real sleep time 
#  detect_body_movements.py  capsule
# ============================================================

detect_body_movements <- function(signal, win_size = 1400, thresh1 = 15) {
  
  n_windows <- floor(length(signal) / win_size)
  window_sd <- numeric(n_windows)
  
  #  calc SD for each window 
  for (i in seq_len(n_windows)) {
    start <- (i - 1) * win_size + 1
    end   <- i * win_size
    seg   <- signal[start:end]
    window_sd[i] <- sqrt(sum((seg - mean(seg))^2) / (length(seg) - 1))
  }
  
  # calc MAD 
  mean_sd  <- mean(window_sd)
  mad_val  <- mean(abs(window_sd - mean_sd))
  thresh2  <- 2 * mad_val
  
  cat("thresh1 (off-bed)  :", thresh1, "\n")
  cat("thresh2 (movement) :", round(thresh2, 2), "\n")
  
  # class each window
  # 1 = sleep   2 = movement        3 = not on bed
  state <- integer(n_windows)
  for (i in seq_len(n_windows)) {
    sd_i <- round(window_sd[i])
    if (sd_i < thresh1) {
      state[i] <- 3
    } else if (sd_i > thresh2) {
      state[i] <- 2
    } else {
      state[i] <- 1
    }
  }
  
  # make mask on samples level
  n_total     <- n_windows * win_size
  sample_flag <- integer(n_total)
  for (i in seq_len(n_windows)) {
    start <- (i - 1) * win_size + 1
    end   <- i * win_size
    sample_flag[start:end] <- state[i]
  }
  
  valid_samples <- sum(sample_flag == 1)
  cat("Windows valid    :", sum(state == 1), "/", n_windows, "\n")
  cat("Samples valid    :", valid_samples, "/", n_total, "\n")
  
  list(
    state        = state,
    sample_flag  = sample_flag,
    window_sd    = window_sd,
    thresh1      = thresh1,
    thresh2      = thresh2,
    win_size     = win_size
  )
}