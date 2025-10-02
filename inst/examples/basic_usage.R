# UniFORM_R Example Usage
# This script demonstrates how to use the UniFORM package for normalizing
# multiplex immunofluorescence data

# Load the package
library(UniFORM)

# Example 1: Basic normalization
cat("Example 1: Basic normalization\n")
cat("================================\n\n")

# Create sample data (simulating immunofluorescence intensity measurements)
set.seed(123)
n_cells <- 100
n_markers <- 5

# Simulate raw intensity data with different scales
raw_data <- matrix(
  c(
    rnorm(n_cells, mean = 100, sd = 20),  # Marker 1
    rnorm(n_cells, mean = 500, sd = 100), # Marker 2
    rnorm(n_cells, mean = 50, sd = 10),   # Marker 3
    rnorm(n_cells, mean = 200, sd = 40),  # Marker 4
    rnorm(n_cells, mean = 150, sd = 30)   # Marker 5
  ),
  nrow = n_cells,
  ncol = n_markers
)

colnames(raw_data) <- paste0("Marker_", 1:5)
rownames(raw_data) <- paste0("Cell_", 1:n_cells)

cat("Raw data summary (first 5 cells):\n")
print(head(raw_data, 5))
cat("\n")

# Apply quantile normalization
normalized_quantile <- normalize_uniform(raw_data, method = "quantile")

cat("Quantile normalized data summary (first 5 cells):\n")
print(head(normalized_quantile, 5))
cat("\n")

# Example 2: Batch-aware normalization
cat("Example 2: Batch-aware normalization\n")
cat("====================================\n\n")

# Create batch labels (e.g., two different experimental batches)
batch <- rep(1:2, each = n_cells / 2)

cat("Batch information:\n")
cat("  Batch 1:", sum(batch == 1), "cells\n")
cat("  Batch 2:", sum(batch == 2), "cells\n\n")

# Apply batch-aware normalization
normalized_batch <- normalize_uniform(raw_data, method = "quantile", batch = batch)

cat("Batch-normalized data summary (first 5 cells):\n")
print(head(normalized_batch, 5))
cat("\n")

# Example 3: Different normalization methods
cat("Example 3: Comparing normalization methods\n")
cat("==========================================\n\n")

# Create smaller dataset for comparison
small_data <- raw_data[1:10, 1:3]

cat("Original data:\n")
print(small_data)
cat("\n")

# Compare different methods
cat("Quantile normalization:\n")
print(normalize_uniform(small_data, method = "quantile"))
cat("\n")

cat("Median normalization:\n")
print(normalize_uniform(small_data, method = "median"))
cat("\n")

cat("Z-score normalization:\n")
print(normalize_uniform(small_data, method = "zscore"))
cat("\n")

cat("Examples completed successfully!\n")
