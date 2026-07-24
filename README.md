# Hyperspectral Pre-Processing Pipeline (`hyperspecprepro`)

This repository contains the R processing pipeline for handling, filtering, normalizing, and extracting fused **CASI-SASI hyperspectral remote sensing data**. 

---

## 📂 Repository Structure & Scripts
The core processing pipeline is divided into three distinct scripts:

*   **`mosaic.R`** — The complete end-to-end processing pipeline to transform raw flightlines into a master site mosaic. It encompasses all steps from data scaling and spectral band reduction through to filtering, global brightness normalization, and exporting the final mosaics as new raster layers.
*   **`extractplots.R`** — The vector-to-raster extraction pipeline. This script uses plot center coordinates to dynamically crop and mask the master mosaics, exporting the isolated plot data as both high-resolution spatial rasters and tabular data frames.
*   **`SavitzkyGolay_testing.R`** — The development sandbox containing all trial-and-error methodologies. It documents the parameter tuning, combinations matrix evaluations, and the empirical grid search used to determine the optimal window lengths and polynomial orders.

---

## 🛠️ Software & Dependency Environment
The processing architecture is built in **R (v4.5.1)** using optimized geospatial, signal processing, and data manipulation libraries:
* **Spatial & Geospatial Operations:** `terra` (v1.9-11)
* **Signal Processing & Filtering:** `prospectr` (v0.2.8), `signal` (v1.8-1)
* **Data Wrangling & High-Throughput Operations:** `data.table` (v1.18.2.1), `dplyr` (v1.2.1), `tidyr` (v1.3.2)
* **Visualization & Quality Control:** `ggplot2` (v4.0.2)

---

## 🔬 Core Processing Workflow

### 1. Band Selection & Segmentation
* **Band Reduction:** Raw digital numbers are converted to scaled reflectance. Interpolated bands (absorption features) and spectral tails are systematically removed, leaving **78 retained bands** spanning the Visible to SWIR regions.
* **Non-Conterminous Handling:** To account for data gaps from band removal, the 78 bands are partitioned into 5 internally contiguous spectral segments:
  *  **Visible-NIR:** 31 bands
  *  **SWIR 1:** 7 bands
  *  **SWIR 2:** 8 bands
  *  **SWIR 3:** 16 bands
  *  **SWIR 4:** 16 bands

### 2. Hybrid Savitzky-Golay Filtering Architecture
To suppress pixel-level sensor chatter and noise without memory exhaustion from large rasters, a high-throughput hybrid optimization architecture was engineered:
* **Grid Search Tuning:** Moving frame window lengths ($w = 3, 5, 7, 9, 11$) and local polynomial orders ($p = 2, 3, 4$) were tested via an empirical grid search. Mathematical violations ($p \ge w$) were automatically filtered. Visual diagnostics established a **second-order local polynomial ($p = 2$)** with a **5-band moving window ($w = 5$)** as optimal.
* **Core Filtering:** For each segment, `prospectr::savitzkyGolay` processes the core internal bands across the entire pixel matrix simultaneously.
* **Edge Reconstruction:** Because standard moving windows trim terminal boundaries, `signal::sgolayfilt` applies specialized, asymmetric polynomial weights to recover the missing 2 outer bands on each edge. This avoids computationally expensive full-spectrum row loops.

### 3. Global Brightness Normalization & Mosaicking
* **Anomalies Suppressed:** Suppresses cross-track illumination errors, cloud shading effects, and topographic variations.
* **Global Scale Selection:** Evaluated against a segment-wise alternative; global normalization was selected to protect relative spectral shape profiles and prevent artificial discontinuities.
* **Mathematical Calibration:** Calculated as the square root of the sum of squared values for each pixel across all 78 dimensions. Every waveband value is divided by this factor to scale it to a uniform length.
* **Artifact Removal:** Background noise is eliminated by reclassifying any pixel with a maximum reflectance value below 0.001 across all dimensions as `NA`.
* **Compilation:** Cleaned flightlines are compiled into spatial raster collections and combined into a master site mosaic using **mean-average cell aggregation**.

### 4. Vector-to-Raster Plot Extraction Pipeline
* **Target Sampling:** Circular plot boundaries (400-square-meter sampling area) are generated using center coordinates and intersected with the master mosaics.
* **Memory Optimization:** Imagery is dynamically clipped to the immediate bounding area around each individual plot prior to extraction.
* **Bilinear Disaggregation:** To eliminate jagged edges caused by the native 1.25-meter spatial resolution, pixels are disaggregated into a fine **0.25-meter grid** using bilinear interpolation.
* **Data Output:** Mosaics are cropped and masked to the circular boundaries. Outputs are exported both as spatial rasters and as clean tabular data frames with explicit per-pixel tracking coordinates and parent plot identifiers appended.

---
*Note: Large raster datasets (`.bsq`, `.hdr`, `.tif`) are processed using local hardware paths and are excluded from version control tracking via `.gitignore`.*
