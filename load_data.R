# ============================================================
# load_data.R
# BCG & RR from dataset
# ============================================================

library(readr)

load_bcg_file <- function(filepath) {
  
  # اread text
  lines <- readLines(filepath)
  
  # delete header)
  lines <- lines[-1]
  
  # first row bCG,Timestamp,fs
  first  <- strsplit(lines[1], ",")[[1]]
  fs     <- as.numeric(trimws(first[3]))      # 140
  ts_start <- as.numeric(trimws(first[2]))    # timestamp
  
  # the rest row contain  BCG 
  bcg_vals <- as.numeric(trimws(lines))
  bcg_vals <- bcg_vals[!is.na(bcg_vals)]
  
  n        <- length(bcg_vals)
  time_sec <- (seq_len(n) - 1) / fs
  
  # do timestamps first from fs
  # each sample  = 1000/140   ms = 7.14 ms
  timestamp <- ts_start + (seq_len(n) - 1) * (1000 / fs)
  
  cat("fs =", fs, "Hz | samples =", n,
      "| duration =", round(n / fs / 60, 1), "min\n")
  
  data.frame(
    time_sec  = time_sec,
    bcg       = bcg_vals,
    timestamp = timestamp,
    fs        = fs
  )
}
load_rr_file <- function(filepath) {
  
  # file contain header: Timestamp, Heart Rate, RR Interval in seconds
  raw <- read_csv(filepath, col_names = TRUE, show_col_types = FALSE)
  
  colnames(raw) <- c("Timestamp", "HeartRate_bpm", "RR_sec")
  
  # convert numbers
  raw$HeartRate_bpm <- as.numeric(raw$HeartRate_bpm)
  raw$RR_sec        <- as.numeric(raw$RR_sec)
  
  # remove null col and rows
  raw <- raw[!is.na(raw$HeartRate_bpm) & raw$HeartRate_bpm > 0, ]
  
  # remove unvalid values
  raw <- raw[raw$HeartRate_bpm >= 30 & raw$HeartRate_bpm <= 200, ]
  
  cat("RR beats:", nrow(raw),
      "| HR range:", round(min(raw$HeartRate_bpm), 1),
      "-", round(max(raw$HeartRate_bpm), 1), "BPM\n")
  
  data.frame(
    timestamp  = as.character(raw$Timestamp),
    hr_ref_bpm = raw$HeartRate_bpm,
    rr_sec     = raw$RR_sec
  )
}

get_subject_file_pairs <- function(subject_dir) {
  bcg_dir <- file.path(subject_dir, "BCG")
  rr_dir  <- file.path(subject_dir, "Reference", "RR")
  if (!dir.exists(bcg_dir) || !dir.exists(rr_dir)) return(NULL)
  bcg_files <- list.files(bcg_dir, pattern = "_BCG\\.csv$", full.names = TRUE)
  rr_files  <- list.files(rr_dir,  pattern = "_RR\\.csv$",  full.names = TRUE)
  if (length(bcg_files) == 0 || length(rr_files) == 0) return(NULL)
  extract_date <- function(f) strsplit(basename(f), "_")[[1]][2]
  bcg_dates    <- sapply(bcg_files, extract_date)
  rr_dates     <- sapply(rr_files,  extract_date)
  common_dates <- intersect(bcg_dates, rr_dates)
  if (length(common_dates) == 0) return(NULL)
  lapply(common_dates, function(d) {
    list(
      date     = d,
      bcg_file = bcg_files[bcg_dates == d],
      rr_file  = rr_files[rr_dates   == d]
    )
  })
}

load_all_subjects <- function(data_root) {
  subject_dirs <- sort(list.dirs(data_root, recursive = FALSE))
  subject_dirs <- subject_dirs[grepl("/\\d{2}$", subject_dirs)]
  cat("no of subjects:", length(subject_dirs), "\n\n")
  all_pairs <- list()
  for (subj_dir in subject_dirs) {
    subj_id <- basename(subj_dir)
    pairs   <- get_subject_file_pairs(subj_dir)
    if (is.null(pairs)) {
      cat("Subject", subj_id, "-- null\n")
      next
    }
    cat("Subject", subj_id, "--", length(pairs), "night\n")
    for (p in pairs) {
      all_pairs[[length(all_pairs) + 1]] <- data.frame(
        subject  = subj_id,
        date     = p$date,
        bcg_file = p$bcg_file,
        rr_file  = p$rr_file,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, all_pairs)
}
