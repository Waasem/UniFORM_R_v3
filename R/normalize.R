#' UniFORM Package
#'
#' Universal ImmunoFluorescence nORMalization
#'
#' @description
#' An R implementation of the UniFORM pipeline for normalizing multiplex
#' immunofluorescence data. This package provides robust normalization methods
#' to harmonize intensity measurements across different samples, batches, and
#' experimental conditions.
#'
#' @docType package
#' @name UniFORM-package
#' @aliases UniFORM
#'
#' @keywords internal
"_PACKAGE"

#' Normalize Immunofluorescence Data
#'
#' @description
#' Normalize multiplex immunofluorescence intensity measurements using the
#' UniFORM method. This function harmonizes intensity values across different
#' samples, batches, and experimental conditions.
#'
#' @param data A numeric matrix or data frame containing intensity measurements.
#'   Rows represent features (e.g., cells or regions) and columns represent
#'   markers or channels.
#' @param method Character string specifying the normalization method.
#'   Currently supported: "quantile" (default), "median", "zscore".
#' @param batch Optional vector indicating batch membership for each sample.
#'   If provided, batch-specific normalization will be applied.
#'
#' @return A normalized matrix or data frame with the same dimensions as the
#'   input data.
#'
#' @examples
#' # Create example data
#' set.seed(123)
#' data <- matrix(rnorm(100), nrow=10, ncol=10)
#'
#' # Normalize using default method
#' normalized_data <- normalize_uniform(data)
#'
#' # Normalize with batch information
#' batch <- rep(1:2, each=5)
#' normalized_data_batch <- normalize_uniform(data, batch=batch)
#'
#' @export
normalize_uniform <- function(data, method = "quantile", batch = NULL) {
  
  # Input validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input data must be a matrix or data frame")
  }
  
  # Convert to matrix if data frame
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  
  # Check for numeric data
  if (!is.numeric(data)) {
    stop("Input data must be numeric")
  }
  
  # Perform normalization based on method
  normalized <- switch(
    method,
    "quantile" = normalize_quantile(data, batch),
    "median" = normalize_median(data, batch),
    "zscore" = normalize_zscore(data, batch),
    stop("Unknown normalization method: ", method)
  )
  
  return(normalized)
}

# Internal normalization functions

#' @keywords internal
normalize_quantile <- function(data, batch = NULL) {
  # Placeholder for quantile normalization
  # In a full implementation, this would perform quantile normalization
  
  if (is.null(batch)) {
    # Simple quantile-based normalization
    # Normalize each column to have the same distribution
    normalized <- apply(data, 2, function(x) {
      (x - median(x, na.rm = TRUE)) / mad(x, na.rm = TRUE)
    })
  } else {
    # Batch-aware normalization
    normalized <- data
    for (b in unique(batch)) {
      batch_idx <- which(batch == b)
      normalized[batch_idx, ] <- apply(data[batch_idx, , drop = FALSE], 2, function(x) {
        (x - median(x, na.rm = TRUE)) / mad(x, na.rm = TRUE)
      })
    }
  }
  
  return(normalized)
}

#' @keywords internal
normalize_median <- function(data, batch = NULL) {
  # Median-based normalization
  
  if (is.null(batch)) {
    normalized <- apply(data, 2, function(x) {
      x - median(x, na.rm = TRUE)
    })
  } else {
    normalized <- data
    for (b in unique(batch)) {
      batch_idx <- which(batch == b)
      normalized[batch_idx, ] <- apply(data[batch_idx, , drop = FALSE], 2, function(x) {
        x - median(x, na.rm = TRUE)
      })
    }
  }
  
  return(normalized)
}

#' @keywords internal
normalize_zscore <- function(data, batch = NULL) {
  # Z-score normalization
  
  if (is.null(batch)) {
    normalized <- apply(data, 2, function(x) {
      (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    })
  } else {
    normalized <- data
    for (b in unique(batch)) {
      batch_idx <- which(batch == b)
      normalized[batch_idx, ] <- apply(data[batch_idx, , drop = FALSE], 2, function(x) {
        (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
      })
    }
  }
  
  return(normalized)
}

#' MAD (Median Absolute Deviation)
#'
#' @keywords internal
mad <- function(x, na.rm = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  median(abs(x - median(x)))
}
