# UniFORM_R: Universal ImmunoFluorescence nORMalization

![R](https://img.shields.io/badge/R-276DC3?style=flat&logo=r&logoColor=white)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

UniFORM_R is an R implementation of the UniFORM (Universal ImmunoFluorescence nORMalization) pipeline for normalizing multiplex immunofluorescence data. This package provides robust normalization methods to harmonize intensity measurements across different samples, batches, and experimental conditions.

## Features

- Robust normalization methods for multiplex immunofluorescence data
- Harmonization of intensity measurements across samples
- Support for batch effect correction
- Compatible with various experimental conditions

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("Waasem/UniFORM_R_v3")
```

## Usage

### Basic Example

```r
library(UniFORM)

# Create sample immunofluorescence intensity data
set.seed(123)
data <- matrix(rnorm(100), nrow=10, ncol=10)
colnames(data) <- paste0("Marker_", 1:10)

# Normalize using default quantile method
normalized_data <- normalize_uniform(data)

# View results
head(normalized_data)
```

### Batch-Aware Normalization

```r
# Create sample data with batch information
batch <- rep(1:2, each=5)

# Apply batch-aware normalization
normalized_batch <- normalize_uniform(data, method="quantile", batch=batch)
```

### Different Normalization Methods

```r
# Quantile normalization (default)
norm_quantile <- normalize_uniform(data, method="quantile")

# Median normalization
norm_median <- normalize_uniform(data, method="median")

# Z-score normalization
norm_zscore <- normalize_uniform(data, method="zscore")
```

For more examples, see the scripts in `inst/examples/`.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use UniFORM_R in your research, please cite:

```
UniFORM_R: Universal ImmunoFluorescence nORMalization
https://github.com/Waasem/UniFORM_R_v3
```

## Contact

For questions and feedback, please open an issue on GitHub.