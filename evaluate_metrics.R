# ============================================================
# evaluate_metrics.R
#plots
# ============================================================

library(ggplot2)

# ---- MAE ----
calc_mae <- function(predicted, reference) {
  mean(abs(predicted - reference))
}

# ---- RMSE ----
calc_rmse <- function(predicted, reference) {
  sqrt(mean((predicted - reference)^2))
}

# ---- MAPE ----
calc_mape <- function(predicted, reference) {
  mean(abs((predicted - reference) / reference)) * 100
}

# ---- Pearson r manual ----
calc_pearson <- function(x, y) {
  xbar <- mean(x)
  ybar <- mean(y)
  num  <- sum((x - xbar) * (y - ybar))
  den  <- sqrt(sum((x - xbar)^2) * sum((y - ybar)^2))
  if (den == 0) return(NA_real_)
  num / den
}

# ---- P-value  ----
calc_pvalue <- function(r, n) {
  t_stat <- r * sqrt((n - 2) / (1 - r^2))
  2 * pt(-abs(t_stat), df = n - 2)
}

# ---- metrics---
evaluate_all <- function(hr_bcg, hr_ref, label = "results") {
  n    <- min(length(hr_bcg), length(hr_ref))
  pred <- hr_bcg[1:n]
  ref  <- hr_ref[1:n]
  
  mae  <- calc_mae(pred, ref)
  rmse <- calc_rmse(pred, ref)
  mape <- calc_mape(pred, ref)
  r    <- calc_pearson(pred, ref)
  pval <- calc_pvalue(r, n)
  
  cat("\n========", label, "========\n")
  cat("N         :", n, "\n")
  cat("MAE       :", round(mae,  3), "BPM\n")
  cat("RMSE      :", round(rmse, 3), "BPM\n")
  cat("MAPE      :", round(mape, 3), "%\n")
  cat("Pearson r :", round(r,    4), "\n")
  cat("P-value   :", format(pval, scientific = TRUE, digits = 3), "\n")
  
  list(n=n, mae=mae, rmse=rmse, mape=mape, r=r, pval=pval, pred=pred, ref=ref)
}

# ---- Regression Plot ----
plot_regression <- function(m, title = "BCG vs Reference HR") {
  df <- data.frame(BCG = m$pred, Ref = m$ref)
  ggplot(df, aes(x = BCG, y = Ref)) +
    geom_point(alpha = 0.5, colour = "#2563EB", size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, colour = "#1E3A5F") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    annotate("text", x = Inf, y = -Inf,
             label = paste0("r = ", round(m$r, 3),
                            "\np = ", format(m$pval, digits=2, scientific=TRUE)),
             hjust = 1.1, vjust = -0.5, size = 3.5) +
    labs(title = title, x = "BCG HR (BPM)", y = "Reference HR (BPM)") +
    theme_minimal()
}

# ---- Bland-Altman Plot ----
plot_bland_altman <- function(m, title = "Bland-Altman") {
  diff_v <- m$pred - m$ref
  mean_v <- (m$pred + m$ref) / 2
  bias   <- mean(diff_v)
  sd_d   <- sqrt(sum((diff_v - bias)^2) / (length(diff_v) - 1))
  loa_hi <- bias + 1.96 * sd_d
  loa_lo <- bias - 1.96 * sd_d
  
  cat("Bias:", round(bias,2), "| LoA:", round(loa_lo,2), "to", round(loa_hi,2), "\n")
  
  df <- data.frame(mean = mean_v, diff = diff_v)
  ggplot(df, aes(x = mean, y = diff)) +
    geom_point(alpha = 0.5, colour = "#DC2626", size = 1.5) +
    geom_hline(yintercept = bias,   colour = "#1E3A5F", linewidth = 0.8) +
    geom_hline(yintercept = loa_hi, colour = "#1E3A5F", linetype = "dashed") +
    geom_hline(yintercept = loa_lo, colour = "#1E3A5F", linetype = "dashed") +
    labs(title = title,
         x = "Mean of BCG & Reference (BPM)",
         y = "Difference BCG - Reference (BPM)") +
    theme_minimal()
}