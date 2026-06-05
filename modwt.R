# ============================================================
# modwt.R
# wavelet 
#  modwt_matlab_fft.py + modwt_mra_matlab_fft.py
#  wavelet: bior3.9 —  4
# ============================================================

# filters  bior3.9 (const PyWavelets)
get_bior39_filters <- function() {
  rec_lo <- c(-0.00138, 0.00276, 0.01398, -0.01671, -0.06688,
              0.07937,  0.41700,  0.60637,  0.41700,  0.07937,
              -0.06688, -0.01671,  0.01398,  0.00276, -0.00138)
  rec_hi <- c( 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,
               0.00000, -0.17678,  0.53033, -0.53033,  0.17678,
               0.00000,  0.00000,  0.00000,  0.00000,  0.00000)
  list(lo = rec_lo / sqrt(2), hi = rec_hi / sqrt(2))
}

# MODWT decomposition — div signal into levels
modwt_manual <- function(x, J = 4) {
  filters <- get_bior39_filters()
  N       <- length(x)
  Vhat    <- fft(x)
  lo_fft  <- fft(c(filters$lo, rep(0, N - length(filters$lo))))
  hi_fft  <- fft(c(filters$hi, rep(0, N - length(filters$hi))))
  
  w_list <- vector("list", J + 1)
  
  for (jj in 0:(J - 1)) {
    upfactor <- 2^jj
    idx      <- (upfactor * (0:(N - 1))) %% N + 1
    Gup      <- lo_fft[idx]
    Hup      <- hi_fft[idx]
    w_list[[jj + 1]] <- Re(fft(Hup * Vhat, inverse = TRUE)) / N
    Vhat             <- Gup * Vhat
  }
  w_list[[J + 1]] <- Re(fft(Vhat, inverse = TRUE)) / N
  
  do.call(rbind, w_list)
}

# MODWT MRA — return each level as a separate signal 
modwtmra_manual <- function(w, J = 4) {
  filters <- get_bior39_filters()
  N       <- ncol(w)
  J0      <- nrow(w) - 1
  lo_fft  <- fft(c(filters$lo, rep(0, N - length(filters$lo))))
  hi_fft  <- fft(c(filters$hi, rep(0, N - length(filters$hi))))
  null_v  <- rep(0, N)
  
  imodwtrec <- function(v_in, w_in, lev) {
    upfactor <- 2^lev
    idx  <- (upfactor * (0:(N - 1))) %% N + 1
    Gup  <- Conj(lo_fft[idx])
    Hup  <- Conj(hi_fft[idx])
    Re(fft(Gup * fft(v_in) + Hup * fft(w_in), inverse = TRUE)) / N
  }
  
  mra <- vector("list", J0 + 1)
  
  for (J in J0:1) {
    v   <- null_v
    wco <- w[J, ]
    for (jj in J:1) {
      v   <- imodwtrec(v, wco, jj - 1)
      wco <- null_v
    }
    mra[[J]] <- v
  }
  
  v <- w[J0 + 1, ]
  for (J in (J0 + 1):1) {
    v <- imodwtrec(v, null_v, J - 1)
  }
  mra[[J0 + 1]] <- v
  
  do.call(rbind, mra)
}

# main function return 4th level of cardiac_cycle
extract_cardiac_cycle <- function(signal) {
  w       <- modwt_manual(signal, J = 4)
  mra     <- modwtmra_manual(w, J = 4)
  cardiac <- mra[5, ]
  cardiac
}