# uniform_registration.R
# Registration functions for UniFORM normalization pipeline
# R implementation of Python UniFORM.registration module

#' Compute correlation-based shift between histograms
#' @param ref_hist Reference histogram counts
#' @param target_hist Target histogram counts
#' @param max_shift Maximum shift to test (in bins)
#' @return Optimal shift value
compute_correlation_shifts <- function(ref_hist, target_hist, max_shift = 50) {
  n_bins <- length(ref_hist)
  correlations <- numeric(2 * max_shift + 1)
  shifts <- seq(-max_shift, max_shift)

  for (i in seq_along(shifts)) {
    shift <- shifts[i]

    if (shift < 0) {
      # Shift left
      ref_portion <- ref_hist[1:(n_bins + shift)]
      target_portion <- target_hist[(-shift + 1):n_bins]
    } else if (shift > 0) {
      # Shift right
      ref_portion <- ref_hist[(shift + 1):n_bins]
      target_portion <- target_hist[1:(n_bins - shift)]
    } else {
      # No shift
      ref_portion <- ref_hist
      target_portion <- target_hist
    }

    # Calculate correlation
    if (length(ref_portion) > 0 && sd(ref_portion) > 0 && sd(target_portion) > 0) {
      correlations[i] <- cor(ref_portion, target_portion)
    } else {
      correlations[i] <- 0
    }
  }

  # Find shift with maximum correlation
  best_idx <- which.max(correlations)
  best_shift <- shifts[best_idx]

  return(best_shift)
}

#' Compute landmark-based shifts
#' @param landmarks List of landmark positions for each sample
#' @param reference_idx Index of reference sample
#' @param bin_width Width of each histogram bin
#' @return Vector of shifts for each sample
compute_landmark_shifts <- function(landmarks, reference_idx = 1, bin_width = 1) {
  n_samples <- length(landmarks)
  shifts <- numeric(n_samples)

  ref_landmark <- landmarks[[reference_idx]]
  if (is.na(ref_landmark) || is.null(ref_landmark)) {
    return(shifts)
  }

  for (i in 1:n_samples) {
    if (i == reference_idx) {
      shifts[i] <- 0
    } else {
      if (!is.na(landmarks[[i]]) && !is.null(landmarks[[i]])) {
        # Shift in bins
        shifts[i] <- round((ref_landmark - landmarks[[i]]) / bin_width)
      } else {
        shifts[i] <- 0
      }
    }
  }

  return(shifts)
}

#' Correct shifted histograms
#' @param histograms List of histogram data
#' @param shifts Vector of shift values
#' @return Corrected histograms
correct_shifted_histograms <- function(histograms, shifts) {
  corrected <- list()
  sample_names <- names(histograms)

  for (i in seq_along(sample_names)) {
    sample <- sample_names[i]
    shift <- shifts[i]
    hist_data <- histograms[[sample]]

    if (shift == 0) {
      corrected[[sample]] <- hist_data
    } else {
      # Apply shift to both bin edges and centers
      bin_width <- if (!is.null(hist_data$bin_edges) && length(hist_data$bin_edges) > 1) {
        diff(hist_data$bin_edges[1:2])
      } else if (!is.null(hist_data$bins) && length(hist_data$bins) > 1) {
        diff(hist_data$bins[1:2])
      } else {
        1
      }

      corrected_data <- list(
        counts = hist_data$counts
      )

      if (!is.null(hist_data$bin_edges)) {
        corrected_data$bin_edges <- hist_data$bin_edges + shift * bin_width
      }
      if (!is.null(hist_data$bin_centers)) {
        corrected_data$bin_centers <- hist_data$bin_centers + shift * bin_width
      }
      if (!is.null(hist_data$bins)) {
        corrected_data$bins <- hist_data$bins + shift * bin_width
      }

      corrected[[sample]] <- corrected_data
    }
  }

  return(corrected)
}

#' Plot histogram distributions
#' @param histograms List of histogram data
#' @param title Plot title
#' @param output_path Path to save plot
#' @param reference_id Reference sample ID to highlight
plot_histogram_distributions <- function(histograms, title = "", output_path = NULL, reference_id = NULL) {
  if (!is.null(output_path)) {
    png(output_path, width = 800, height = 600)
  }

  # Find global range - use bin_edges if available, otherwise bins
  all_x <- unlist(lapply(histograms, function(h) {
    if (!is.null(h$bin_edges)) h$bin_edges else h$bins
  }))
  all_counts <- unlist(lapply(histograms, function(h) h$counts))

  # Find max count for y-axis
  max_count <- max(all_counts)

  plot(0, 0, type = "n",
       xlim = range(all_x),
       ylim = c(0, max_count * 1.1),
       xlab = "Grid Points",
       ylab = "Pixel Counts",
       main = title,
       cex.lab = 1.2,  # 12pt equivalent
       font.lab = 2,   # bold
       font.main = 2)  # bold

  # Use tab20-like colors
  tab20_colors <- c(
    "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78",
    "#2ca02c", "#98df8a", "#d62728", "#ff9896",
    "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
    "#e377c2", "#f7b6d9", "#7f7f7f", "#c7c7c7",
    "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"
  )
  colors <- rep(tab20_colors, length.out = length(histograms))

  sample_names <- names(histograms)
  for (i in seq_along(histograms)) {
    hist_data <- histograms[[i]]
    # Use bin_edges if available, otherwise bins
    x_coords <- if (!is.null(hist_data$bin_edges)) hist_data$bin_edges else hist_data$bins
    # Highlight reference sample
    if (!is.null(reference_id) && sample_names[i] == reference_id) {
      lines(x_coords, hist_data$counts, col = "black", lwd = 3)
    } else {
      lines(x_coords, hist_data$counts, col = colors[i], lwd = 2, type = "l")
    }
  }

  # Custom legend with reference highlighted
  if (!is.null(reference_id)) {
    ref_idx <- which(sample_names == reference_id)
    legend_cols <- colors
    legend_cols[ref_idx] <- "black"
    legend_lwds <- rep(2, length(sample_names))
    legend_lwds[ref_idx] <- 3
    legend_labels <- sample_names
    legend_labels[ref_idx] <- paste0(reference_id, " (reference)")
    legend("topright", legend = legend_labels, col = legend_cols, lwd = legend_lwds, cex = 0.8, bty = "n")
  } else {
    legend("topright", legend = sample_names, col = colors, lwd = 2, cex = 0.8, bty = "n")
  }

  if (!is.null(output_path)) {
    dev.off()
  }
}

#' Automatic registration of histograms
#' @param histogram_data Histogram data from process_sample_distributions
#' @param all_markers Vector of all marker names
#' @param selected_markers Markers to process
#' @param sample_ids Sample IDs
#' @param reference_samples Optional reference samples for each marker
#' @param landmark_map Optional manual landmark positions
#' @param num_bins Number of histogram bins
#' @param dpi DPI for plots
#' @param x_limits X-axis limits for plots
#' @return List with shifts_map, chosen_references, implied_landmarks_map
automatic_registration <- function(
  histogram_data,
  all_markers,
  selected_markers = NULL,
  sample_ids = NULL,
  reference_samples = NULL,
  landmark_map = NULL,
  num_bins = 1024,
  dpi = 100,
  x_limits = NULL
) {

  if (is.null(selected_markers)) {
    selected_markers <- all_markers
  }

  if (is.null(sample_ids)) {
    # Get sample IDs from first marker
    first_marker <- names(histogram_data)[1]
    sample_ids <- names(histogram_data[[first_marker]])
  }

  shifts_map <- list()
  chosen_references <- list()
  implied_landmarks_map <- list()

  cat("Performing automatic registration...\n")

  for (marker in selected_markers) {
    if (!(marker %in% names(histogram_data))) {
      cat("  Marker", marker, "not found in histogram data\n")
      next
    }

    cat("  Processing marker:", marker, "\n")

    marker_hists <- histogram_data[[marker]]

    # Choose reference sample
    if (!is.null(reference_samples) && marker %in% names(reference_samples)) {
      ref_sample <- reference_samples[[marker]]
    } else {
      # Use first sample as reference by default
      ref_sample <- sample_ids[1]
    }

    chosen_references[[marker]] <- ref_sample
    ref_idx <- which(sample_ids == ref_sample)

    # Check for manual landmarks
    if (!is.null(landmark_map) && marker %in% names(landmark_map)) {
      landmarks <- landmark_map[[marker]]
      if (!is.null(landmarks) && length(landmarks) == length(sample_ids)) {
        # Use landmark-based registration
        cat("    Using landmark-based registration\n")
        bin_width <- diff(marker_hists[[1]]$bins[1:2])
        shifts <- compute_landmark_shifts(landmarks, ref_idx, bin_width)
        implied_landmarks_map[[marker]] <- landmarks
      } else {
        # Use correlation-based registration
        cat("    Using correlation-based registration\n")
        shifts <- numeric(length(sample_ids))
        ref_hist <- marker_hists[[ref_sample]]$counts

        for (i in seq_along(sample_ids)) {
          if (i == ref_idx) {
            shifts[i] <- 0
          } else {
            target_hist <- marker_hists[[sample_ids[i]]]$counts
            shifts[i] <- compute_correlation_shifts(ref_hist, target_hist)
          }
        }

        # Compute implied landmarks from shifts
        implied_landmarks <- numeric(length(sample_ids))
        for (i in seq_along(sample_ids)) {
          # Find mode of histogram as landmark
          hist_data <- marker_hists[[sample_ids[i]]]
          mode_idx <- which.max(hist_data$counts)
          implied_landmarks[i] <- mode_idx - shifts[i]
        }
        implied_landmarks_map[[marker]] <- implied_landmarks
      }
    } else {
      # Use correlation-based registration
      cat("    Using correlation-based registration\n")
      shifts <- numeric(length(sample_ids))
      ref_hist <- marker_hists[[ref_sample]]$counts

      for (i in seq_along(sample_ids)) {
        if (i == ref_idx) {
          shifts[i] <- 0
        } else {
          target_hist <- marker_hists[[sample_ids[i]]]$counts
          shifts[i] <- compute_correlation_shifts(ref_hist, target_hist)
        }
      }

      # Compute implied landmarks from shifts
      implied_landmarks <- numeric(length(sample_ids))
      for (i in seq_along(sample_ids)) {
        # Find mode of histogram as landmark
        hist_data <- marker_hists[[sample_ids[i]]]
        mode_idx <- which.max(hist_data$counts)
        implied_landmarks[i] <- mode_idx - shifts[i]
      }
      implied_landmarks_map[[marker]] <- implied_landmarks
    }

    # Store shifts
    names(shifts) <- sample_ids
    shifts_map[[marker]] <- shifts

    # Apply shifts and visualize
    corrected_hists <- correct_shifted_histograms(marker_hists, shifts)

    # Create before/after plots
    plot_path <- paste0("registration_", marker, ".png")
    png(plot_path, width = 1600, height = 600, res = dpi)
    par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

    # Before registration
    plot_histogram_distributions(marker_hists,
                                 paste(marker, "Original Distributions"),
                                 reference_id = ref_sample)

    # After registration
    plot_histogram_distributions(corrected_hists,
                                 paste(marker, "Normalized Distributions"),
                                 reference_id = ref_sample)

    dev.off()
    cat("    Saved registration plot to:", plot_path, "\n")
  }

  cat("Registration complete!\n")

  return(list(
    shifts_map = shifts_map,
    chosen_references = chosen_references,
    implied_landmarks_map = implied_landmarks_map
  ))
}