# ============================================================
# compute_hr.R
#  beat_to_beat.py + compute_vitals.py
# ============================================================

source("/home/vboxuser/project_R/detect_peaks.R")

compute_hr_window <- function(sig_window, time_ms_window, mpd = 1) {
  
  peaks <- detect_peaks_manual(sig_window, mpd = mpd)
  if (length(peaks) < 2) return(NA_real_)
  
  # diff between peaks ms
  peak_times    <- time_ms_window[peaks]
  intervals_ms  <- diff(peak_times)
  if (length(intervals_ms) == 0) return(NA_real_)
  
  # BPM = 60000 / mean diff bet. peaks ms
  bpm <- 60000 / mean(intervals_ms)
  round(bpm, 2)
}

compute_hr_all_windows <- function(cardiac_signal, time_ms,
                                   win_size = 1400, valid_windows) {
  hr_bcg <- c()
  hr_ref_indices <- c()
  
  for (w_idx in valid_windows) {
    t1 <- (w_idx - 1) * win_size + 1
    t2 <- w_idx * win_size
    if (t2 > length(cardiac_signal)) break
    
    hr <- compute_hr_window(
      cardiac_signal[t1:t2],
      time_ms[t1:t2],
      mpd = 1
    )
    
    # if valid
    if (!is.na(hr) && hr >= 30 && hr <= 200) {
      hr_bcg         <- c(hr_bcg, hr)
      hr_ref_indices <- c(hr_ref_indices, w_idx)
    }
  }
  
  cat("Windows value calc HR:", length(hr_bcg), "\n")
  list(hr_bcg = hr_bcg, window_indices = hr_ref_indices)
}