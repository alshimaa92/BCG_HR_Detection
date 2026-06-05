
library(writexl)
output_dir <- "./results/metrics5"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

diff_v     <- hr_bcg_all - hr_ref_all
total_n    <- length(hr_bcg_all)
final_bias <- mean(diff_v, na.rm = TRUE)
final_sd_d <- sd(diff_v, na.rm = TRUE)
final_mae  <- mean(abs(diff_v), na.rm = TRUE)
final_rmse <- sqrt(mean(diff_v^2, na.rm = TRUE))
final_mape <- mean(abs(diff_v / hr_ref_all), na.rm = TRUE) * 100
final_r    <- cor(hr_bcg_all, hr_ref_all, method = "pearson", use = "complete.obs")
t_stat     <- final_r * sqrt((total_n - 2) / (1 - final_r^2))
final_pval <- 2 * pt(-abs(t_stat), df = total_n - 2)

summary_df <- data.frame(
  Metric = c("MAE", "RMSE", "MAPE", "Pearson r", "P-value", "Bias (mean diff)", "SD of diff"),
  Value  = c(final_mae, final_rmse, final_mape, final_r, final_pval, final_bias, final_sd_d),
  Unit   = c("BPM", "BPM", "%", "", "", "BPM", "BPM")
)

today_date <- Sys.Date()
csv_file   <- file.path(output_dir, paste0("metrics_summary_", today_date, ".csv"))
excel_file <- file.path(output_dir, paste0("metrics_summary_", today_date, ".xlsx"))
write.csv(summary_df, csv_file, row.names = FALSE)
write_xlsx(summary_df, excel_file)

cat("\n  ------------------------------------------------------------------------\n")
cat(sprintf("  All_Subjects_Pooled  (n = %d paired windows)\n", total_n))
cat("  ------------------------------------------------------------------------\n")
cat(sprintf("  MAE                = %.4f BPM\n", final_mae))
cat(sprintf("  RMSE               = %.4f BPM\n", final_rmse))
cat(sprintf("  MAPE               = %.4f %%\n", final_mape))
cat(sprintf("  Pearson r          = %.4f\n", final_r))
cat(sprintf("  P-value            = %f\n", final_pval))
cat(sprintf("  Bias (mean diff)   = %.4f BPM\n", final_bias))
cat(sprintf("  SD of diff         = %.4f BPM\n\n", final_sd_d))
cat(sprintf("  \u2714 Summary CSV -> %s\n", csv_file))
cat(sprintf("  \u2714 Summary Excel -> %s\n", excel_file))
cat(sprintf("  \u2714 Plots saved in: %s\n\n  Done.\n\n", output_dir))
