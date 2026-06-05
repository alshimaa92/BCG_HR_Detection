# BCG Heart Rate Detection — R Implementation

## Project Overview
Comparative analysis of heart rate detection from Ballistocardiography (BCG)
signals using R (Method 1) and Julia (Method 2).

**Course:** Advanced Biostatistics  
**Supervisor:** Dr. Ibrahim Sadek  
**Dataset:** Figshare DOI: 10.6084/m9.figshare.26013157  
**Reference:** DOI: 10.1038/s41597-024-03950-5  

---

## Project Roadmap

---

## Dependencies

```r
install.packages("readr")
install.packages("ggplot2")
install.packages("signal")
```

**System requirements (Ubuntu):**
```bash
sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev
```

---

## Usage

```bash
# 1. Clone the repository
git clone https://github.com/YOUR_USERNAME/BCG_HR_Detection.git
cd BCG_HR_Detection

# 2. Download dataset from Figshare
# https://doi.org/10.6084/m9.figshare.26013157
# Extract to: data/

# 3. Run in R terminal
R
source("main.R")
```

---

## File Structure
> **[MANUAL]** = manually implemented — no black-box functions used

---

## Results Summary

| Metric | Mean (All Subjects) | Best Subject |
|--------|-------------------|--------------|
| MAE    | 23.4 BPM          | 7.4 BPM (Sub 01) |
| RMSE   | 26.6 BPM          | 9.9 BPM (Sub 01) |
| MAPE   | 41.4 %            | 9.0 % (Sub 01)   |
| Pearson r | -0.032         | 0.200 (Sub 01)   |

**Dataset:** 22 subjects · 42 nights · 331,000+ valid windows

---

## Manual Implementations

All core algorithms are implemented manually:

- **SD calculation:** `sqrt(sum((x - mean(x))^2) / (n-1))`
- **MAD:** `mean(abs(window_sd - mean(window_sd)))`
- **MODWT:** FFT-based convolution with bior3.9 filters
- **Peak detection:** Sign-change method with min-distance constraint
- **Pearson r:** `sum((x-xbar)(y-ybar)) / sqrt(sum((x-xbar)^2) * sum((y-ybar)^2))`
- **P-value:** `2 * pt(-|t|, df=n-2)`
