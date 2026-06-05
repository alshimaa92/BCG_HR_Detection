# ============================================================
# detect_peaks.R 
#detect_peaks.py 
# ============================================================

detect_peaks_manual <- function(x, mph = NULL, mpd = 1, threshold = 0) {
  
  n  <- length(x)
  if (n < 3) return(integer(0))
  
  dx <- diff(x)
  
  # fin rising and falling points
  rising  <- c(dx, 0) <= 0
  falling <- c(0, dx) > 0
  peaks   <- which(rising & falling)
  
  # cancel the first and the last
  peaks <- peaks[peaks != 1 & peaks != n]
  if (length(peaks) == 0) return(integer(0))
  
  #filter the min rising
  if (!is.null(mph)) {
    peaks <- peaks[x[peaks] >= mph]
  }
  
  # filter with neighbors
  if (threshold > 0 && length(peaks) > 0) {
    left_diff  <- x[peaks] - x[pmax(peaks - 1, 1)]
    right_diff <- x[peaks] - x[pmin(peaks + 1, n)]
    peaks <- peaks[pmin(left_diff, right_diff) >= threshold]
  }
  
  # filter the min dist between two peaks
  if (mpd > 1 && length(peaks) > 1) {
    peaks <- peaks[order(x[peaks], decreasing = TRUE)]
    keep  <- rep(TRUE, length(peaks))
    for (i in seq_along(peaks)) {
      if (keep[i]) {
        too_close <- which(abs(peaks - peaks[i]) <= mpd & seq_along(peaks) != i)
        keep[too_close] <- FALSE
      }
    }
    peaks <- sort(peaks[keep])
  }
  
  peaks
}