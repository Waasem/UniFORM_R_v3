# uniform_preprocessing.R
# Preprocessing functions for UniFORM normalization pipeline
# R implementation of Python UniFORM.preprocessing module

#' Log transform intensity values
#' @param intensities Numeric vector of raw intensity values
#' @param min_value Minimum threshold for values to be retained
#' @return Log-transformed intensity values
log_transform_intensities <- function(intensities, min_value = 1.0) {
  # Filter out values below threshold
  valid_mask <- intensities >= min_value
  filtered <- intensities[valid_mask]

  # Take natural log
  log_values <- log(filtered)

  # Replace any -Inf with 0
  log_values[is.infinite(log_values) & log_values < 0] <- 0.0

  return(log_values)
}

#' Fit two-component Gaussian Mixture Model
#' @param intensities 1D array of intensity values (log-transformed)
#' @param num_bins Number of evenly spaced points for evaluation
#' @param min_value Minimum value of the evaluation grid
#' @param max_value Maximum value of the evaluation grid
#' @param verbose Print messages if TRUE
#' @return List with threshold, eval_points, pdf_comp0, pdf_comp1
fit_two_component_gmm <- function(intensities, num_bins, min_value, max_value, verbose = FALSE) {

  # Remove NA values
  intensities <- intensities[!is.na(intensities)]

  if (length(intensities) < 10) {
    warning("Not enough data points for GMM fitting")
    eval_points <- seq(min_value, max_value, length.out = num_bins)
    return(list(
      threshold = NA,
      eval_points = eval_points,
      pdf_comp0 = rep(0, num_bins),
      pdf_comp1 = rep(0, num_bins),
      means = c(NA, NA),
      sds = c(NA, NA),
      weights = c(NA, NA)
    ))
  }

  # Try to fit GMM using simple EM algorithm
  tryCatch({
    # Initial guess using quantiles
    q25 <- quantile(intensities, 0.25)
    q75 <- quantile(intensities, 0.75)

    # Initialize two components
    mu1 <- q25
    mu2 <- q75
    sigma1 <- sd(intensities[intensities <= median(intensities)])
    sigma2 <- sd(intensities[intensities > median(intensities)])
    pi1 <- 0.5
    pi2 <- 0.5

    # EM algorithm (simplified)
    for (iter in 1:50) {
      # E-step: compute responsibilities
      dens1 <- dnorm(intensities, mu1, sigma1)
      dens2 <- dnorm(intensities, mu2, sigma2)

      resp1 <- (pi1 * dens1) / (pi1 * dens1 + pi2 * dens2)
      resp2 <- 1 - resp1

      # M-step: update parameters
      n1 <- sum(resp1)
      n2 <- sum(resp2)

      pi1 <- n1 / length(intensities)
      pi2 <- n2 / length(intensities)

      mu1 <- sum(resp1 * intensities) / n1
      mu2 <- sum(resp2 * intensities) / n2

      sigma1 <- sqrt(sum(resp1 * (intensities - mu1)^2) / n1)
      sigma2 <- sqrt(sum(resp2 * (intensities - mu2)^2) / n2)
    }

    # Ensure component 1 has lower mean
    if (mu1 > mu2) {
      temp <- mu1; mu1 <- mu2; mu2 <- temp
      temp <- sigma1; sigma1 <- sigma2; sigma2 <- temp
      temp <- pi1; pi1 <- pi2; pi2 <- temp
    }

    # Compute threshold as intersection point
    threshold <- (mu1 + mu2) / 2

    # Evaluate PDFs on grid
    eval_points <- seq(min_value, max_value, length.out = num_bins)
    pdf_comp0 <- pi1 * dnorm(eval_points, mu1, sigma1)
    pdf_comp1 <- pi2 * dnorm(eval_points, mu2, sigma2)

    return(list(
      threshold = threshold,
      eval_points = eval_points,
      pdf_comp0 = pdf_comp0,
      pdf_comp1 = pdf_comp1,
      means = c(mu1, mu2),
      sds = c(sigma1, sigma2),
      weights = c(pi1, pi2)
    ))

  }, error = function(e) {
    if (verbose) cat("GMM fitting failed:", e$message, "\n")
    eval_points <- seq(min_value, max_value, length.out = num_bins)
    return(list(
      threshold = median(intensities),
      eval_points = eval_points,
      pdf_comp0 = rep(0, num_bins),
      pdf_comp1 = rep(0, num_bins),
      means = c(NA, NA),
      sds = c(NA, NA),
      weights = c(NA, NA)
    ))
  })
}

#' Process sample distributions
#' @param feature_input List of data frames or AnnData object
#' @param sample_ids Vector of sample IDs
#' @param all_markers Vector of all marker names
#' @param markers_to_plot Markers to process (default: all_markers)
#' @param use_normalized Use normalized data if TRUE
#' @param num_bins Number of histogram bins
#' @param plots_per_row Plots per row for visualization
#' @param dpi DPI for plots
#' @param xlims X-axis limits
#' @param ylims Y-axis limits
#' @param output_figure_path Path to save figure
#' @param verbose Print progress messages
#' @return List with intensity_ranges, histograms, gmm_models
process_sample_distributions <- function(
  feature_input,
  sample_ids = NULL,
  all_markers = NULL,
  markers_to_plot = NULL,
  use_normalized = FALSE,
  num_bins = 1024,
  plots_per_row = 4,
  dpi = 100,
  xlims = NULL,
  ylims = NULL,
  output_figure_path = NULL,
  verbose = TRUE
) {

  # Determine input type and extract data
  if (is.list(feature_input) && !is.data.frame(feature_input)) {
    # List of data frames
    if (is.null(sample_ids)) {
      sample_ids <- paste0("Sample_", seq_along(feature_input))
    }
    if (is.null(all_markers)) {
      # Get markers from first sample
      all_markers <- setdiff(names(feature_input[[1]]), c("cell_id", "sample_id", "x", "y"))
    }
    data_type <- "list"
  } else {
    stop("Unsupported input type")
  }

  if (is.null(markers_to_plot)) {
    markers_to_plot <- all_markers
  }

  # Initialize output structures
  intensity_ranges <- list()
  histograms <- list()
  gmm_models <- list()

  for (marker in markers_to_plot) {
    if (verbose) cat("Processing marker:", marker, "\n")

    # Collect all intensities for this marker
    all_intensities <- c()

    if (data_type == "list") {
      for (i in seq_along(feature_input)) {
        if (marker %in% names(feature_input[[i]])) {
          if (use_normalized && paste0(marker, "_normalized") %in% names(feature_input[[i]])) {
            intensities <- feature_input[[i]][[paste0(marker, "_normalized")]]
          } else {
            intensities <- feature_input[[i]][[marker]]
          }
          log_int <- log_transform_intensities(intensities)
          all_intensities <- c(all_intensities, log_int)
        }
      }
    }

    if (length(all_intensities) == 0) {
      if (verbose) cat("  No data found for marker", marker, "\n")
      next
    }

    # Determine global range
    global_min <- min(all_intensities, na.rm = TRUE)
    global_max <- max(all_intensities, na.rm = TRUE)

    intensity_ranges[[marker]] <- c(global_min, global_max)

    # Create histogram bins
    hist_breaks <- seq(global_min, global_max, length.out = num_bins + 1)
    bin_centers <- (hist_breaks[-1] + hist_breaks[-length(hist_breaks)]) / 2
    bin_edges <- hist_breaks[-length(hist_breaks)]  # edges[:-1] in Python

    # Process each sample
    histograms[[marker]] <- list()
    gmm_models[[marker]] <- list()

    for (i in seq_along(sample_ids)) {
      sample_id <- sample_ids[i]

      if (data_type == "list") {
        if (marker %in% names(feature_input[[i]])) {
          if (use_normalized && paste0(marker, "_normalized") %in% names(feature_input[[i]])) {
            intensities <- feature_input[[i]][[paste0(marker, "_normalized")]]
          } else {
            intensities <- feature_input[[i]][[marker]]
          }
          log_int <- log_transform_intensities(intensities)

          # Calculate histogram
          hist_counts <- hist(log_int, breaks = hist_breaks, plot = FALSE)$counts
          histograms[[marker]][[sample_id]] <- list(
            bin_edges = bin_edges,  # Store edges for plotting
            bin_centers = bin_centers,  # Keep centers for compatibility
            counts = hist_counts
          )

          # Fit GMM
          gmm <- fit_two_component_gmm(log_int, num_bins, global_min, global_max, verbose = FALSE)
          gmm_models[[marker]][[sample_id]] <- gmm
        }
      }
    }
  }

  # Create visualization if requested
  if (!is.null(output_figure_path) && length(markers_to_plot) > 0) {
    if (verbose) cat("Creating distribution plots...\n")

    # Calculate grid dimensions
    n_markers <- length(markers_to_plot)
    n_rows <- ceiling(n_markers / plots_per_row)

    # Create plot - match Python figsize=(20, n_rows * 4)
    # Convert to pixels: 20 inches * 100 dpi = 2000 pixels width
    png(output_figure_path, width = 20 * dpi, height = n_rows * 4 * dpi, res = dpi)

    # Set up layout with proper margins for legend
    layout_matrix <- matrix(1:(n_rows * plots_per_row), nrow = n_rows, byrow = TRUE)
    layout(layout_matrix, widths = rep(1, plots_per_row))
    par(mar = c(5, 5, 3, 8), oma = c(0, 0, 0, 4))

    # Use tab20 color palette like Python
    # Approximate tab20 colors
    tab20_colors <- c(
      "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78",
      "#2ca02c", "#98df8a", "#d62728", "#ff9896",
      "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
      "#e377c2", "#f7b6d9", "#7f7f7f", "#c7c7c7",
      "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"
    )
    colors <- rep(tab20_colors, length.out = length(sample_ids))

    plot_idx <- 1
    for (marker in markers_to_plot) {
      if (marker %in% names(histograms)) {
        # Find y-axis range based on maximum frequency
        max_freq <- 0
        for (sample_id in sample_ids) {
          if (sample_id %in% names(histograms[[marker]])) {
            max_freq <- max(max_freq, max(histograms[[marker]][[sample_id]]$counts))
          }
        }

        # Plot histograms for all samples
        plot(0, 0, type = "n",
             xlim = if (!is.null(xlims)) xlims else intensity_ranges[[marker]],
             ylim = if (!is.null(ylims)) ylims else c(0, max_freq * 1.1),
             xlab = "Log(Cell Mean Intensity)",
             ylab = "Frequency",
             main = marker,
             cex.lab = 1.0,  # 10pt equivalent
             cex.main = 1.2, # 12pt equivalent
             font.lab = 2,   # bold
             font.main = 2)  # bold

        # Remove top and right spines
        box(bty = "l", lwd = 2)

        for (i in seq_along(sample_ids)) {
          sample_id <- sample_ids[i]
          if (sample_id %in% names(histograms[[marker]])) {
            hist_data <- histograms[[marker]][[sample_id]]
            # Plot raw counts (Frequency) using bin edges like Python
            x_coords <- if (!is.null(hist_data$bin_edges)) hist_data$bin_edges else hist_data$bins
            lines(x_coords, hist_data$counts,
                  col = colors[i], lwd = 2, type = "l")
          }
        }

        # Add legend only for first plot
        if (plot_idx == 1) {
          legend("topright", legend = sample_ids, col = colors,
                 lwd = 2, cex = 0.8, bty = "n")
        }

        plot_idx <- plot_idx + 1
      }
    }

    # Hide extra plots if needed
    while (plot_idx <= n_rows * plots_per_row) {
      plot.new()
      plot_idx <- plot_idx + 1
    }

    dev.off()
    if (verbose) cat("Saved distribution plots to:", output_figure_path, "\n")
  }

  return(list(
    intensity_ranges = intensity_ranges,
    histograms = histograms,
    gmm_models = gmm_models
  ))
}