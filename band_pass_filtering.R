# ============================================================
# band_pass_filtering.R
#isolate CG
# band_pass_filtering.py from  capsule
# 2.5 Hz:  5.0 Hz
# ============================================================

bandpass_filter_bcg <- function(signal, fs) {
  
  # lib for signal filter coefficients 
  if (!requireNamespace("signal", quietly = TRUE)) {
    stop("install.packages('signal')")
  }
  
  # Highpass  2.5 Hz — removeـ DC and slowartifacts
  hp <- signal::cheby1(2, 0.5, 2.5 / (fs / 2), type = "high")
  
  # Lowpass عند 5.0 Hz — remove high noise
  lp <- signal::cheby1(4, 0.5, 5.0 / (fs / 2), type = "low")
  
  # apply the two filters 
  sig_hp <- signal::filtfilt(hp, signal)
  sig_bp <- signal::filtfilt(lp, sig_hp)
  
  sig_bp
}