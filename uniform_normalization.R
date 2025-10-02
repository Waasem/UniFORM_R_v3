# uniform_normalization.R
# Normalization functions for UniFORM pipeline
# R implementation of Python UniFORM.normalization module

#' Calculate shift in log pixel space
#' @param shift_bins Shift in histogram bins
#' @param bin_width Width of each bin
#' @return Shift in log pixel values
calculate_shift_in_log_pixels <- function(shift_bins, bin_width) {
  return(shift_bins * bin_width)
}

#' Apply normalization shift to intensities
#' @param intensities Original intensity values
#' @param shift Shift value to apply
#' @param log_space Whether values are already in log space
#' @return Normalized intensities
apply_normalization_shift <- function(intensities, shift, log_space = FALSE) {
  if (log_space) {
    # Values already in log space, just add shift
    return(intensities + shift)
  } else {
    # Transform to log, apply shift, transform back
    log_int <- log(pmax(intensities, 1))
    normalized_log <- log_int + shift
    return(exp(normalized_log))
  }
}

#' Plot line histogram
#' @param histograms List of histogram data
#' @param marker Marker name
#' @param sample_ids Sample IDs to plot
#' @param output_path Path to save plot
#' @param plot_density If TRUE plot density, if FALSE plot frequency
plot_line_histogram <- function(histograms, marker, sample_ids = NULL, output_path = NULL, plot_density = TRUE) {
  if (is.null(sample_ids)) {
    sample_ids <- names(histograms[[marker]])
  }

  if (!is.null(output_path)) {
    png(output_path, width = 800, height = 600)
  }

  # Use tab20-like colors
  tab20_colors <- c(
    "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78",
    "#2ca02c", "#98df8a", "#d62728", "#ff9896",
    "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
    "#e377c2", "#f7b6d9", "#7f7f7f", "#c7c7c7",
    "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"
  )
  colors <- rep(tab20_colors, length.out = length(sample_ids))

  # Find range
  all_x <- c()
  all_values <- c()

  for (sample in sample_ids) {
    if (sample %in% names(histograms[[marker]])) {
      hist_data <- histograms[[marker]][[sample]]
      # Use bin_edges if available, otherwise bins
      x_coords <- if (!is.null(hist_data$bin_edges)) hist_data$bin_edges else hist_data$bins
      all_x <- c(all_x, x_coords)

      if (plot_density) {
        bin_width <- if (length(x_coords) > 1) diff(x_coords[1:2]) else 1
        density <- hist_data$counts / sum(hist_data$counts) / bin_width
        all_values <- c(all_values, density)
      } else {
        all_values <- c(all_values, hist_data$counts)
      }
    }
  }

  plot(0, 0, type = "n",
       xlim = range(all_x),
       ylim = c(0, max(all_values) * 1.1),
       xlab = paste0("log(", marker, ")"),
       ylab = if (plot_density) "Density" else "Frequency",
       main = paste(marker, "Distributions"),
       cex.lab = 1.0,   # 10pt
       cex.main = 1.2,  # 12pt
       font.lab = 1,    # normal for label
       font.main = 2)   # bold for title

  for (i in seq_along(sample_ids)) {
    sample <- sample_ids[i]
    if (sample %in% names(histograms[[marker]])) {
      hist_data <- histograms[[marker]][[sample]]
      # Use bin_edges if available, otherwise bins
      x_coords <- if (!is.null(hist_data$bin_edges)) hist_data$bin_edges else hist_data$bins

      if (plot_density) {
        bin_width <- if (length(x_coords) > 1) diff(x_coords[1:2]) else 1
        values <- hist_data$counts / sum(hist_data$counts) / bin_width
      } else {
        values <- hist_data$counts
      }
      lines(x_coords, values, col = colors[i], lwd = 2, type = "l")
    }
  }

  legend("topright", legend = sample_ids, col = colors, lwd = 2, cex = 0.8, bty = "n")

  if (!is.null(output_path)) {
    dev.off()
  }
}

#' Plot correlations and fit line
#' @param x_values X values
#' @param y_values Y values
#' @param title Plot title
#' @param output_path Path to save plot
plot_correlations_and_fit_line <- function(x_values, y_values, title = "", output_path = NULL) {
  if (!is.null(output_path)) {
    png(output_path, width = 600, height = 600)
  }

  # Fit linear model
  fit <- lm(y_values ~ x_values)

  plot(x_values, y_values,
       xlab = "Reference Sample",
       ylab = "Target Sample",
       main = title,
       pch = 16, col = rgb(0, 0, 1, 0.3))

  abline(fit, col = "red", lwd = 2)
  abline(a = 0, b = 1, col = "gray", lty = 2)

  # Add R-squared
  r2 <- summary(fit)$r.squared
  legend("topleft", legend = paste("R² =", round(r2, 3)), bty = "n")

  if (!is.null(output_path)) {
    dev.off()
  }
}

#' Plot GMM components
#' @param gmm_model GMM model from fit_two_component_gmm
#' @param histogram_data Histogram data
#' @param title Plot title
#' @param output_path Path to save plot
plot_gmm <- function(gmm_model, histogram_data, title = "", output_path = NULL) {
  if (!is.null(output_path)) {
    png(output_path, width = 600, height = 500)
  }

  # Normalize histogram to density
  density <- histogram_data$counts / sum(histogram_data$counts) / diff(histogram_data$bins[1:2])

  plot(histogram_data$bins, density,
       type = "l", lwd = 2,
       xlab = "Log(Intensity)",
       ylab = "Density",
       main = title)

  if (!is.null(gmm_model) && !is.na(gmm_model$threshold)) {
    # Plot GMM components
    lines(gmm_model$eval_points, gmm_model$pdf_comp0, col = "blue", lwd = 2, lty = 2)
    lines(gmm_model$eval_points, gmm_model$pdf_comp1, col = "red", lwd = 2, lty = 2)
    lines(gmm_model$eval_points, gmm_model$pdf_comp0 + gmm_model$pdf_comp1,
          col = "green", lwd = 2, lty = 3)

    # Add threshold line
    abline(v = gmm_model$threshold, col = "purple", lwd = 2, lty = 4)

    legend("topright",
           legend = c("Data", "Component 1", "Component 2", "Combined", "Threshold"),
           col = c("black", "blue", "red", "green", "purple"),
           lty = c(1, 2, 2, 3, 4),
           lwd = 2, cex = 0.8)
  }

  if (!is.null(output_path)) {
    dev.off()
  }
}

#' Generate normalized feature data
#' @param feature_input Input data (list of data frames)
#' @param sample_ids Sample IDs
#' @param markers Markers to normalize
#' @param intensity_ranges Intensity ranges from preprocessing
#' @param shifts_map Shifts from registration
#' @param chosen_references Reference samples
#' @param num_bins Number of bins
#' @param dpi DPI for plots
#' @param plot_dist Plot distributions
#' @param plot_single_cell_corr Plot single-cell correlations
#' @param gmm_analysis Perform GMM analysis
#' @param save_normalized_features Save normalized data
#' @return Updated feature_input with normalized values
generate_normalized_feature <- function(
  feature_input,
  sample_ids = NULL,
  markers = NULL,
  intensity_ranges,
  shifts_map,
  chosen_references,
  num_bins = 1024,
  dpi = 100,
  plot_dist = FALSE,
  plot_single_cell_corr = FALSE,
  gmm_analysis = FALSE,
  save_normalized_features = TRUE
) {

  cat("Generating normalized features...\n")

  # Determine input type
  if (is.list(feature_input) && !is.data.frame(feature_input)) {
    data_type <- "list"
    if (is.null(sample_ids)) {
      sample_ids <- paste0("Sample_", seq_along(feature_input))
    }
    if (is.null(markers)) {
      markers <- setdiff(names(feature_input[[1]]), c("cell_id", "sample_id", "x", "y"))
    }
  } else {
    stop("Unsupported input type")
  }

  # Create output directory for normalized data
  if (save_normalized_features) {
    norm_dir <- "normalized_data_R"
    if (!dir.exists(norm_dir)) {
      dir.create(norm_dir)
    }
  }

  # Apply normalization for each marker
  for (marker in markers) {
    if (!(marker %in% names(shifts_map))) {
      cat("  Skipping marker", marker, "(no shift data)\n")
      next
    }

    cat("  Normalizing marker:", marker, "\n")

    shifts <- shifts_map[[marker]]
    ref_sample <- chosen_references[[marker]]

    if (!is.null(intensity_ranges) && marker %in% names(intensity_ranges)) {
      bin_width <- diff(intensity_ranges[[marker]]) / num_bins
    } else {
      bin_width <- 0.01  # Default
    }

    # Apply normalization to each sample
    for (i in seq_along(sample_ids)) {
      sample_id <- sample_ids[i]

      if (sample_id %in% names(shifts)) {
        shift_bins <- shifts[[sample_id]]
        shift_log <- calculate_shift_in_log_pixels(shift_bins, bin_width)

        if (data_type == "list") {
          if (marker %in% names(feature_input[[i]])) {
            # Get original intensities
            original <- feature_input[[i]][[marker]]

            # Apply normalization
            normalized <- apply_normalization_shift(original, shift_log, log_space = FALSE)

            # Add normalized column
            feature_input[[i]][[paste0(marker, "_normalized")]] <- normalized

            cat("    ", sample_id, ": applied shift of", round(shift_log, 3), "\n")
          }
        }
      }
    }

    # Plot distributions if requested
    if (plot_dist) {
      # Use tab20-like colors
      tab20_colors <- c(
        "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78",
        "#2ca02c", "#98df8a", "#d62728", "#ff9896",
        "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
        "#e377c2", "#f7b6d9", "#7f7f7f", "#c7c7c7",
        "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"
      )
      colors <- rep(tab20_colors, length.out = length(sample_ids))

      # Before normalization
      plot_path_before <- paste0("distribution_", marker, "_before.png")
      png(plot_path_before, width = 800, height = 600, res = dpi)

      plot(0, 0, type = "n", xlim = intensity_ranges[[marker]], ylim = c(0, 0.3),
           xlab = paste0("log(", marker, ")"),
           ylab = "Density",
           main = paste(marker, "- Before Normalization"),
           cex.lab = 1.0,   # 10pt
           cex.main = 1.2,  # 12pt
           font.main = 2)   # bold

      for (i in seq_along(sample_ids)) {
        if (marker %in% names(feature_input[[i]])) {
          log_int <- log_transform_intensities(feature_input[[i]][[marker]])
          dens <- density(log_int)
          lines(dens$x, dens$y, col = colors[i], lwd = 2)
        }
      }
      legend("topright", legend = sample_ids, col = colors, lwd = 2, cex = 0.6, bty = "n")
      dev.off()

      # After normalization
      plot_path_after <- paste0("distribution_", marker, "_after.png")
      png(plot_path_after, width = 800, height = 600, res = dpi)

      plot(0, 0, type = "n", xlim = intensity_ranges[[marker]], ylim = c(0, 0.3),
           xlab = paste0("log(", marker, ")"),
           ylab = "Density",
           main = paste(marker, "- After Normalization"),
           cex.lab = 1.0,   # 10pt
           cex.main = 1.2,  # 12pt
           font.main = 2)   # bold

      for (i in seq_along(sample_ids)) {
        norm_col <- paste0(marker, "_normalized")
        if (norm_col %in% names(feature_input[[i]])) {
          log_int <- log_transform_intensities(feature_input[[i]][[norm_col]])
          dens <- density(log_int)
          lines(dens$x, dens$y, col = colors[i], lwd = 2)
        }
      }
      legend("topright", legend = sample_ids, col = colors, lwd = 2, cex = 0.6, bty = "n")
      dev.off()

      cat("    Saved distribution plots\n")
    }
  }

  # Save normalized data if requested
  if (save_normalized_features) {
    for (i in seq_along(sample_ids)) {
      sample_id <- sample_ids[i]
      output_file <- file.path(norm_dir, paste0(sample_id, "_normalized.rds"))
      saveRDS(feature_input[[i]], output_file)
      cat("  Saved normalized data for", sample_id, "\n")
    }
  }

  cat("Normalization complete!\n")

  return(invisible(feature_input))
}

#' Verify normalization results
#' @param original_data Original data
#' @param normalized_data Normalized data
#' @param markers Markers to check
#' @param sample_ids Sample IDs
#' @return Summary statistics
verify_normalization <- function(original_data, normalized_data, markers, sample_ids) {
  cat("\n=== Normalization Verification ===\n")

  results <- list()

  for (marker in markers) {
    cat("\nMarker:", marker, "\n")

    # Calculate statistics before normalization
    orig_means <- numeric(length(sample_ids))
    orig_sds <- numeric(length(sample_ids))

    for (i in seq_along(sample_ids)) {
      if (marker %in% names(original_data[[i]])) {
        log_int <- log_transform_intensities(original_data[[i]][[marker]])
        orig_means[i] <- mean(log_int)
        orig_sds[i] <- sd(log_int)
      }
    }

    cat("  Before normalization:\n")
    cat("    Mean range:", round(range(orig_means), 3), "\n")
    cat("    SD range:", round(range(orig_sds), 3), "\n")
    cat("    CV of means:", round(sd(orig_means) / mean(orig_means), 3), "\n")

    # Calculate statistics after normalization
    norm_means <- numeric(length(sample_ids))
    norm_sds <- numeric(length(sample_ids))

    norm_col <- paste0(marker, "_normalized")
    for (i in seq_along(sample_ids)) {
      if (norm_col %in% names(normalized_data[[i]])) {
        log_int <- log_transform_intensities(normalized_data[[i]][[norm_col]])
        norm_means[i] <- mean(log_int)
        norm_sds[i] <- sd(log_int)
      }
    }

    cat("  After normalization:\n")
    cat("    Mean range:", round(range(norm_means), 3), "\n")
    cat("    SD range:", round(range(norm_sds), 3), "\n")
    cat("    CV of means:", round(sd(norm_means) / mean(norm_means), 3), "\n")

    # Calculate improvement
    improvement <- (sd(orig_means) - sd(norm_means)) / sd(orig_means) * 100
    cat("  Improvement in consistency:", round(improvement, 1), "%\n")

    results[[marker]] <- list(
      original_means = orig_means,
      original_sds = orig_sds,
      normalized_means = norm_means,
      normalized_sds = norm_sds,
      improvement = improvement
    )
  }

  return(invisible(results))
}