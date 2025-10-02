# UniFORM_R_v3
**R implementation of UniFORM (Unified Normalization of Features Observed in Regions of Multiplexed images)**

A comprehensive normalization pipeline for correcting batch effects in multiplexed imaging data, with a focus on prostate cancer (PRAD) tissue analysis.

---
## Features
- ✅ **Histogram-based normalization** - Aligns feature distributions across multiple samples
- ✅ **Automatic registration** - Correlation-based histogram alignment with no manual intervention
- ✅ **Interactive landmark refinement** - Optional web-based interface for manual landmark adjustment
- ✅ **GMM analysis** - Two-component Gaussian Mixture Models to identify cell populations
- ✅ **Comprehensive visualization** - Before/after plots for quality control
- ✅ **Modular design** - Four independent modules for flexible workflow customization
- ✅ **Pure R implementation** - No Python dependencies required (uses base R + graphics)

---
## Background
Multiplexed imaging technologies (e.g., CyCIF, CODEX, IMC) enable simultaneous measurement of dozens of protein markers in tissue samples. However, technical variability between samples creates batch effects that can confound biological interpretation.

**UniFORM** addresses this by:
1. Computing histograms of marker intensities across all cells
2. Identifying optimal alignment between sample distributions
3. Applying intensity corrections to normalize data while preserving biological variation

This R implementation reproduces the UniFORM method, making it accessible to the R community without requiring Python installations.

---
## Installation & Requirements
### Requirements
- **R** ≥ 3.6.0
- **Base R packages only** (graphics, stats, utils)
- No external dependencies required!

### Installation
```bash
# Clone or download this repository
git clone https://github.com/Waasem/UniFORM_R_v3.git
cd UniFORM_R_v3
# No installation needed - just source the files
```

---
## Quick Start
Run the complete pipeline on synthetic PRAD data:
```bash
Rscript PRAD_analysis.R
```

This will:
1. Generate 20 synthetic PRAD samples (20,000 cells each, 20 markers)
2. Compute histograms and fit GMMs
3. Perform automatic registration
4. Generate interactive landmark refinement tools
5. Apply normalization
6. Create visualization plots

**Expected runtime:** ~2-5 minutes  
**Output:** 60+ PNG files + normalized data in `normalized_data_R/`

---
## Detailed Usage
### Complete Workflow (8 Steps)
The `PRAD_analysis.R` script orchestrates the complete pipeline:

#### Step 1: Load Raw Data
```r
# Load pipeline modules
source("uniform_preprocessing.R")
source("uniform_registration.R")
source("uniform_landmark.R")
source("uniform_normalization.R")

# Define your sample names and markers
sample_names <- c('Sample-01', 'Sample-02', 'Sample-03')
markers <- c('DAPI', 'CD3', 'CD45', 'EPCAM')

# Your data should be a list of data frames
# Each data frame contains: cell_id, sample_id, x, y, and marker columns
```

#### Step 2: Histogram Calculation
```r
results <- process_sample_distributions(
  feature_input = feature_data,
  sample_ids = sample_names,
  all_markers = markers,
  markers_to_plot = markers,
  use_normalized = FALSE,
  num_bins = 1024,
  output_figure_path = "histogram_distributions.png"
)
intensity_ranges <- results$intensity_ranges
histograms <- results$histograms
gmm_models <- results$gmm_models
```

#### Step 3: Automatic Registration
```r
registration_results <- automatic_registration(
  histogram_data = histograms,
  all_markers = markers,
  selected_markers = markers,
  sample_ids = sample_names,
  num_bins = 1024
)
shifts_map <- registration_results$shifts_map
chosen_references <- registration_results$chosen_references
implied_landmarks_map <- registration_results$implied_landmarks_map
```

#### Step 4: Landmark Refinement (Optional)
```r
landmark_refinement(
  histogram_data = histograms,
  gmm_models = gmm_models,
  markers = markers,
  sample_ids = sample_names,
  data_source = feature_data,
  implied_landmarks_map = implied_landmarks_map,
  output_directory = "landmark_picking"
)
```

#### Step 5: Manual Landmark Adjustment
1. Open HTML files in `landmark_picking/` folder
2. Click on distribution peaks to set landmark positions
3. Click "Save Landmarks to CSV"
4. Update `landmark_picking/landmark_annotations_filled.csv`

#### Step 6: Re-run Registration with Manual Landmarks
```r
landmark_map <- read_landmark_csv(
  "landmark_picking/landmark_annotations_filled.csv",
  sample_names
)
registration_results <- automatic_registration(
  histogram_data = histograms,
  all_markers = markers,
  selected_markers = markers,
  sample_ids = sample_names,
  landmark_map = landmark_map,
  num_bins = 1024
)
```

#### Step 7: Generate Normalized Data
```r
normalized_data <- generate_normalized_feature(
  feature_input = feature_data,
  sample_ids = sample_names,
  markers = markers,
  intensity_ranges = intensity_ranges,
  shifts_map = shifts_map,
  chosen_references = chosen_references,
  plot_dist = TRUE,
  save_normalized_features = TRUE
)
```

#### Step 8: Verification
```r
normalized_results <- process_sample_distributions(
  feature_input = normalized_data,
  sample_ids = sample_names,
  all_markers = markers,
  use_normalized = TRUE,
  output_figure_path = "distributions_normalized.png"
)
verification <- verify_normalization(
  original_data = feature_data,
  normalized_data = normalized_data,
  markers = markers,
  sample_ids = sample_names
)
```

---
## Pipeline Overview
```
┌─────────────────────────────────────────────────────────────┐
│                    UniFORM_R_v3 Pipeline                    │
└─────────────────────────────────────────────────────────────┘
Input Data (List of Data Frames)
         │
         ├─► uniform_preprocessing.R
         │   ├─ Log-transform intensities
         │   ├─ Compute histograms (1024 bins)
         │   └─ Fit 2-component GMMs
         │
         ├─► uniform_registration.R
         │   ├─ Correlation-based alignment
         │   ├─ Choose reference sample
         │   └─ Compute shift maps
         │
         ├─► uniform_landmark.R (OPTIONAL)
         │   ├─ Generate interactive HTML
         │   ├─ Export landmark CSV
         │   └─ Re-compute shifts with manual landmarks
         │
         └─► uniform_normalization.R
             ├─ Apply shifts to raw intensities
             ├─ Generate normalized columns
             └─ Create before/after visualizations
Output: Normalized Data + QC Plots
```

---
## Input Data Format
Your data must be a **list of data frames**, where each data frame represents one sample:
```r
feature_data <- list(
  # Sample 1
  data.frame(
    cell_id = 1:1000,
    sample_id = "Sample-01",
    x = runif(1000, 0, 1000),
    y = runif(1000, 0, 1000),
    DAPI = rlnorm(1000, 3, 1),
    CD3 = rlnorm(1000, 2.5, 0.8),
    CD45 = rlnorm(1000, 2.8, 0.9)
    # ... more markers
  ),
  # Sample 2
  data.frame(
    cell_id = 1:1000,
    sample_id = "Sample-02",
    # ... same structure
  )
)
```

**Required columns:**
- `cell_id` - Unique cell identifier within sample
- `sample_id` - Sample name (character)
- `x`, `y` - Spatial coordinates (optional, for visualization)
- One column per marker with raw intensity values (numeric)

---
## Output Files
### Visualization Files (PNG)
| File | Description |
|------|-------------|
| `PRAD_feature_histogram_distributions.png` | Grid of all marker distributions (before normalization) |
| `PRAD_feature_distributions_normalized.png` | Grid of all marker distributions (after normalization) |
| `registration_{marker}.png` | Side-by-side before/after for each marker |
| `distribution_{marker}_before.png` | Individual marker distribution (original) |
| `distribution_{marker}_after.png` | Individual marker distribution (normalized) |

### Data Files
| Directory/File | Description |
|----------------|-------------|
| `normalized_data_R/` | Normalized data in RDS format |
| `normalized_data_R/{sample}_normalized.rds` | Per-sample data with `{marker}_normalized` columns |
| `landmark_picking/` | Interactive landmark refinement tools |
| `landmark_picking/group_{n}_landmarks.html` | Interactive plots for sample groups |
| `landmark_picking/landmark_annotations.csv` | CSV template for manual landmarks |
| `landmark_picking/landmark_annotations_filled.csv` | Pre-filled landmarks from automatic registration |

---
## Interactive Landmark Refinement
The optional landmark refinement step creates interactive HTML files powered by Plotly.js:

### Using the Interface
1. **Open HTML file** in any modern web browser:
   ```bash
   open landmark_picking/group_1_landmarks.html
   ```
2. **View distributions** - Each marker shows overlaid histograms for all samples in the group
3. **Set landmarks** - Click on a histogram peak to set the landmark position for that sample
   - Alternatively, enter values manually in the input fields below each plot
4. **Save landmarks** - Click "Save Landmarks to CSV" button
   - Downloads `landmark_annotations.csv` to your downloads folder
5. **Update the filled CSV**:
   ```bash
   cp ~/Downloads/landmark_annotations.csv landmark_picking/landmark_annotations_filled.csv
   ```
6. **Re-run analysis** - Execute `PRAD_analysis.R` again to use your manual landmarks

### Why Manual Refinement?
- **Bimodal distributions** - Some markers have distinct positive/negative populations
- **Reference alignment** - Ensures all samples align at biologically meaningful features
- **Quality control** - Visual inspection catches registration failures

---
## Key Parameters
### Histogram Resolution
```r
num_bins = 1024  # Number of histogram bins
                 # Higher = finer resolution, slower computation
                 # Recommended: 512-2048
```

### Registration Settings
```r
max_shift = 50   # Maximum shift to search during correlation
                 # Increase if samples are very different
```

### Visualization
```r
dpi = 100        # Plot resolution (dots per inch)
plots_per_row = 4  # Grid layout for multi-marker plots
```

### Landmark Refinement
```r
group_size = 4   # Number of samples per HTML file
                 # Smaller = easier to view, more files
```

### Normalization
```r
plot_dist = TRUE              # Generate before/after distribution plots
plot_single_cell_corr = FALSE # Plot single-cell correlations (slow)
gmm_analysis = FALSE          # Include GMM fit plots
save_normalized_features = TRUE  # Save RDS files
```

---
## Example Results
Running on the PRAD dataset (20 samples × 20 markers × 20,000 cells):

### Before Normalization
- Mean intensity varies significantly across samples
- CV (coefficient of variation) of sample means: ~15-25%
- Visible batch effects in PCA/UMAP

### After Normalization
- Sample distributions align at key biological features
- CV of sample means: ~5-10% (typical improvement: 50-80%)
- Batch effects reduced while preserving biological variation

### Console Output Example
```
Processing marker: CD3
  Before normalization:
    Mean range: 2.45 3.12
    SD range: 0.82 0.98
    CV of means: 0.183
  After normalization:
    Mean range: 2.67 2.81
    SD range: 0.85 0.92
    CV of means: 0.042
  Improvement in consistency: 77.1 %
```

---
## Troubleshooting
### Issue: "GMM fitting failed"
**Solution:** This is usually harmless. GMM fits are for visualization only and don't affect normalization.

### Issue: Registration produces poor alignment
**Solutions:**
1. Use interactive landmark refinement to manually set alignment points
2. Increase `num_bins` for finer resolution
3. Check that your data has sufficient dynamic range (not all zeros)
4. Verify log-transformation is appropriate for your data

### Issue: Memory errors with large datasets
**Solutions:**
1. Reduce `num_bins` to 512
2. Process fewer markers at a time
3. Reduce number of cells per sample for testing
4. Disable distribution plotting (`plot_dist = FALSE`)

### Issue: All samples have identical distributions
**Check:**
1. Your input data has real variation between samples
2. You're not accidentally loading the same sample multiple times
3. Marker names match exactly between data and parameter lists

### Issue: HTML files won't open or plots don't display
**Solutions:**
1. Use a modern browser (Chrome, Firefox, Edge)
2. Check browser console for JavaScript errors
3. Files require internet connection to load Plotly.js from CDN
