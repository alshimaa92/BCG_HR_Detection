# ============================================================
# main.R — Full Pipeline — All Subjects
# ============================================================

source("/home/vboxuser/project_R/load_data.R")
source("/home/vboxuser/project_R/detect_body_movements.R")
source("/home/vboxuser/project_R/band_pass_filtering.R")
source("/home/vboxuser/project_R/modwt.R")
source("/home/vboxuser/project_R/detect_peaks.R")
source("/home/vboxuser/project_R/compute_hr.R")
source("/home/vboxuser/project_R/evaluate_metrics.R")

DATA_ROOT <- "/home/vboxuser/Downloads/dataset(7)/dataset/data"
WIN_SIZE  <- 1400
options(scipen = 999)

cat("\n[ Step 1 ] data load...\n")
file_table <- load_all_subjects(DATA_ROOT)
cat("total nights:", nrow(file_table), "\n")

# total results table
all_results <- data.frame(
  subject  = character(),
  night    = character(),
  n        = integer(),
  mae      = numeric(),
  rmse     = numeric(),
  mape     = numeric(),
  r        = numeric(),
  pval     = numeric(),
  stringsAsFactors = FALSE
)

# ── Loop عon each subject and each night ─────────────────────────
for (row_idx in seq_len(nrow(file_table))) {

  subj_id  <- file_table$subject[row_idx]
  night_id <- file_table$date[row_idx]
  cat("\n========================================\n")
  cat("Subject:", subj_id, "| Night:", night_id, "\n")

  # load data
  bcg <- tryCatch(load_bcg_file(file_table$bcg_file[row_idx]),
                  error = function(e) NULL)
  rr  <- tryCatch(load_rr_file(file_table$rr_file[row_idx]),
                  error = function(e) NULL)
  if (is.null(bcg) || is.null(rr) || nrow(rr) < 10) {
    cat("skip no enough data\n")
    next
  }

  fs <- bcg$fs[1]

  # Step 2: body movements
  movement   <- detect_body_movements(bcg$bcg, win_size = WIN_SIZE)
  valid_wins <- which(movement$state == 1)
  n_trim     <- length(movement$state) * WIN_SIZE
  if (length(valid_wins) < 10) {
    cat("skip no valid  windows \n")
    next
  }

  # Step 3: Bandpass
  sig_filtered <- tryCatch(
    bandpass_filter_bcg(bcg$bcg[1:n_trim], fs),
    error = function(e) NULL
  )
  if (is.null(sig_filtered)) next

  # Step 4: MODWT + HR + Sync
  time_ms     <- bcg$timestamp[1:n_trim]
  hr_bcg_all  <- c()
  hr_ref_all  <- c()

  rr_time_num <- as.numeric(as.POSIXct(
                   rr$timestamp,
                   format = "%Y/%m/%d %H:%M:%S",
                   tz = "UTC")) * 1000 - (8 * 3600 * 1000)

  win_groups <- split(valid_wins, ceiling(seq_along(valid_wins) / 20))

  for (g_idx in seq_along(win_groups)) {
    group <- win_groups[[g_idx]]
    t1    <- (min(group) - 1) * WIN_SIZE + 1
    t2    <- min(max(group) * WIN_SIZE, length(sig_filtered))
    seg   <- sig_filtered[t1:t2]
    if (length(seg) < 100) next

    cardiac <- tryCatch(extract_cardiac_cycle(seg), error = function(e) NULL)
    if (is.null(cardiac)) next

    for (w in seq_along(group)) {
      s1 <- (w - 1) * WIN_SIZE + 1
      s2 <- w * WIN_SIZE
      if (s2 > length(cardiac)) next

      hr <- compute_hr_window(
        cardiac[s1:s2],
        time_ms[(t1 + s1 - 1):(t1 + s2 - 1)],
        mpd = round(fs * 0.5)
      )
      if (is.na(hr) || hr < 30 || hr > 200) next

      mid_idx    <- min(t1 + s1 - 1 + (WIN_SIZE %/% 2), length(time_ms))
      win_mid_ms <- time_ms[mid_idx]
      diffs      <- abs(rr_time_num - win_mid_ms)
      closest    <- which.min(diffs)
      if (diffs[closest] > 60000) next

      hr_bcg_all <- c(hr_bcg_all, hr)
      hr_ref_all <- c(hr_ref_all, rr$hr_ref_bpm[closest])
    }
  }

  if (length(hr_bcg_all) < 10) {
    cat("skip no enough HR estimates \n")
    next
  }

  # Step 5: metrics
  m <- evaluate_all(hr_bcg_all, hr_ref_all,
                    paste("Subject", subj_id, "Night", night_id))

  # save results
  all_results <- rbind(all_results, data.frame(
    subject  = subj_id,
    night    = night_id,
    n        = m$n,
    mae      = round(m$mae,  3),
    rmse     = round(m$rmse, 3),
    mape     = round(m$mape, 3),
    r        = round(m$r,    4),
    pval     = m$pval,
    stringsAsFactors = FALSE
  ))
}

# ── اtotal results────────────────────────────────────────
cat("\n\n========================================\n")
cat("total results for  subjects\n")
cat("========================================\n")
print(all_results)

cat("\n── mean of metrics──\n")
cat("MAE  :", round(mean(all_results$mae,  na.rm=TRUE), 3), "BPM\n")
cat("RMSE :", round(mean(all_results$rmse, na.rm=TRUE), 3), "BPM\n")
cat("MAPE :", round(mean(all_results$mape, na.rm=TRUE), 3), "%\n")
cat("r    :", round(mean(all_results$r,    na.rm=TRUE), 4), "\n")

# ── save results in CSV ───────────────────────────────────
write.csv(all_results,
          "/home/vboxuser/project_R/results_all_subjects.csv",
          row.names = FALSE)
cat("\n save  results_all_subjects.csv\n")

# ── Plots of all data────────────────────────────
cat("\n[ Plots ] plot all results..\n")

# add final BCG & REF from last  subject as ex
p1 <- plot_regression(m)
p2 <- plot_bland_altman(m)
print(p1)
print(p2)

cat("\n[ Done ] !\n")
