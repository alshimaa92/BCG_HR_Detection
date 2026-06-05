# BCG Heart Rate Detection — R Implementation

## Project Overview

This project implements a complete pipeline for detecting heart rate (HR)
from Ballistocardiography (BCG) signals using R. It is part of a comparative
analysis between R (Method 1) and Julia (Method 2) implementations.

BCG signals are recorded from a sensor placed under the mattress. The sensor
captures mechanical forces produced by the heart during each beat. This
pipeline extracts heart rate from these signals and compares the results
against a Polar H9 reference device.

| Item        | Detail                                      |
|-------------|---------------------------------------------|
| Course      | Advanced Biostatistics — PhD Coursework     |
| Supervisor  | Dr. Ibrahim Sadek                           |
| Method      | R (Method 1)                                |
| Dataset     | Figshare DOI: 10.6084/m9.figshare.26013157  |
| Reference   | DOI: 10.1038/s41597-024-03950-5             |
| Environment | Ubuntu 24.04 — R terminal                  |

---

## Project Roadmap
Raw BCG Signal (140 Hz, under-mattress sensor)
│
▼
┌─────────────────────────────┐
│  load_data.R                │  Read BCG + RR files, match by date
│  22 subjects · 42 nights    │  Handle non-standard CSV format
└─────────────┬───────────────┘
│
▼
┌─────────────────────────────┐
│  detect_body_movements.R    │  Sliding SD window (500 samples)
│  MANUAL — no sd()           │  States: sleep / movement / off-bed
└─────────────┬───────────────┘
│
▼
┌─────────────────────────────┐
│  band_pass_filtering.R      │  Chebyshev Type I bandpass
│  2.5 Hz – 5.0 Hz            │  Isolates cardiac frequency band
└─────────────┬───────────────┘
│
▼
┌─────────────────────────────┐
│  modwt.R                    │  MODWT wavelet decomposition
│  MANUAL — FFT-based          │  bior3.9 · level 4 · mra[5]
└─────────────┬───────────────┘
│
▼
┌─────────────────────────────┐
│  detect_peaks.R             │  J-peak detection
│  MANUAL — no findpeaks()    │  mpd = round(fs × 0.5) = 70 samples
└─────────────┬───────────────┘
│
▼
┌─────────────────────────────┐
│  compute_hr.R               │  Peak intervals → BPM
│  Timestamp sync (UTC–8h)    │  Match BCG ↔ RR reference
└─────────────┬───────────────┘
│
▼
┌─────────────────────────────┐
│  evaluate_metrics.R         │  MAE · RMSE · MAPE · Pearson r
│  MANUAL — all formulas      │  P-value · Regression · Bland-Altman
└─────────────┬───────────────┘
│
▼
┌─────────────────────────────┐
│  main.R                     │  Orchestrates full pipeline
│  Loops all 42 nights        │  Saves results_all_subjects.csv
└─────────────────────────────┘

---

## File Structure
project_R/
├── main.R                      # Full pipeline orchestrator
├── load_data.R                 # Data loading & BCG/RR file matching
├── detect_body_movements.R     # [MANUAL] Movement artifact removal
├── band_pass_filtering.R       # Chebyshev bandpass 2.5–5.0 Hz
├── modwt.R                     # [MANUAL] MODWT wavelet decomposition
├── detect_peaks.R              # [MANUAL] J-peak detection algorithm
├── compute_hr.R                # [MANUAL] HR estimation + sync
├── evaluate_metrics.R          # [MANUAL] All metrics + ggplot2 plots
├── results_all_subjects.csv    # Output: metrics per night per subject
└── README.md

> **[MANUAL]** means the algorithm is implemented from scratch
> without any black-box functions (no sd(), no findpeaks(), no cor())

---

## Dataset Structure

Download from: https://doi.org/10.6084/m9.figshare.26013157
data/
01/
BCG/          01_YYYYMMDD_BCG.csv    (BCG signal · 140 Hz)
Reference/
RR/         01_YYYYMMDD_RR.csv     (Polar H9 reference HR)
02/ ... 32/

### BCG File Format
- Row 1: `BCG_value, Unix_timestamp_ms_UTC+8, fs`
- Row 2+: `BCG_value` only (non-standard CSV)
- Sampling rate: 140 Hz

### RR File Format
- Columns: `Timestamp | HeartRate_bpm | RR_interval_sec`
- Timestamp format: `yyyy/MM/dd H:mm:ss`
- Event-driven (one row per heartbeat)

### Data Availability
- 32 subjects total
- 22 subjects have matched BCG + RR reference files
- Subjects 21–30: BCG only (no RR reference) — excluded from analysis
- 42 valid nights across 22 subjects

---

## Dependencies

```r
# Install in R terminal
install.packages("readr")    # CSV loading only
install.packages("ggplot2")  # Visualization only
install.packages("signal")   # Filter coefficients only
```

### System Requirements (Ubuntu)

```bash
sudo apt-get update
sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev
```

---

## Usage

```bash
# 1. Clone the repository
git clone https://github.com/alshimaa92/BCG_HR_Detection.git
cd BCG_HR_Detection

# 2. Download and extract the Figshare dataset
# https://doi.org/10.6084/m9.figshare.26013157
# Place extracted folder at: data/

# 3. Install R packages (run once)
sudo R -e "install.packages(c('readr','ggplot2','signal'),
           repos='https://cran.rstudio.com/')"

# 4. Run the full pipeline
R
source("main.R")
```

---

## Manual Implementations

All core algorithmic steps are implemented manually.
Comments in each file indicate where black-box functions were replaced.

### Standard Deviation (detect_body_movements.R)
```r
# Manual SD — replaces sd()
window_sd[i] <- sqrt(sum((seg - mean(seg))^2) / (length(seg) - 1))
```

### MAD — Mean Absolute Deviation (detect_body_movements.R)
```r
# Manual MAD
mad_val <- mean(abs(window_sd - mean(window_sd)))
thresh2 <- 2 * mad_val
```

### MODWT Wavelet (modwt.R)
```r
# Manual FFT-based MODWT — replaces wavelets package
# bior3.9 filters hardcoded from PyWavelets coefficients
# scaled by 1/sqrt(2) for MODWT normalization
```

### Peak Detection (detect_peaks.R)
```r
# Manual peak detection — replaces findpeaks()
dx    <- diff(x)
peaks <- which((c(dx, 0) <= 0) & (c(0, dx) > 0))
# min peak distance: mpd = round(fs * 0.5) = 70 samples
```

### Pearson Correlation (evaluate_metrics.R)
```r
# Manual Pearson r — replaces cor()
num <- sum((x - mean(x)) * (y - mean(y)))
den <- sqrt(sum((x - mean(x))^2) * sum((y - mean(y))^2))
r   <- num / den
```

### P-value (evaluate_metrics.R)
```r
# Manual p-value — replaces cor.test()
t_stat <- r * sqrt((n - 2) / (1 - r^2))
pval   <- 2 * pt(-abs(t_stat), df = n - 2)
```

---

## Results Summary

| Metric    | Mean (22 subjects) | Best Night     | Worst Night     |
|-----------|--------------------|----------------|-----------------|
| MAE       | 23.4 BPM           | 7.4 (Sub 01)   | 36.8 (Sub 06)   |
| RMSE      | 26.6 BPM           | 9.9 (Sub 01)   | 38.9 (Sub 06)   |
| MAPE      | 41.4 %             | 9.0 (Sub 01)   | 75.8 (Sub 02)   |
| Pearson r | -0.032             | 0.200 (Sub 01) | -0.241 (Sub 10) |

**Coverage:** 22 subjects · 42 nights · 331,000+ valid HR windows

### Key Findings
- Subjects with resting HR > 60 BPM → MAE 7–12 BPM (clinically acceptable)
- Subjects with deep-sleep HR 30–50 BPM → MAE 30–36 BPM
- Timestamp synchronization (UTC+8 offset) was the most critical fix
- Peak detection mpd parameter has the largest single impact on accuracy

---

## Critical Bugs Fixed During Development

| Bug | Cause | Fix |
|-----|-------|-----|
| Pearson r = NA | 8h UTC offset between BCG and RR | Subtract 8×3600×1000 ms |
| MAE = 25 BPM bias | mpd=1 caused 2× peak count | mpd = round(fs × 0.5) = 70 |
| BCG reads 1 sample | Non-standard CSV format | Custom readLines() parser |

---

## License

This project is for academic use only.
Dataset © original authors — see reference publication.
