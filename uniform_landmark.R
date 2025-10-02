# uniform_landmark.R
# Landmark refinement functions for UniFORM normalization pipeline
# R implementation of Python UniFORM.landmark module

#' Generate HTML for interactive landmark picking
#' @param histogram_data Histogram data
#' @param gmm_models GMM models
#' @param markers Markers to include
#' @param sample_ids Sample IDs
#' @param implied_landmarks Initial landmark positions
#' @param output_file HTML output file path
create_landmark_html <- function(
  histogram_data,
  gmm_models,
  markers,
  sample_ids,
  implied_landmarks = NULL,
  output_file = "landmark_picking.html"
) {

  html_content <- '
<!DOCTYPE html>
<html>
<head>
    <title>UniFORM Landmark Picker</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .marker-section { margin-bottom: 30px; border: 1px solid #ccc; padding: 15px; }
        .plot-container { width: 100%; height: 400px; margin: 10px 0; }
        .landmark-inputs { margin: 10px 0; }
        input[type="number"] { width: 80px; margin: 0 5px; }
        .sample-label { display: inline-block; width: 100px; }
        button { padding: 10px 20px; margin: 10px 5px; background: #4CAF50; color: white; border: none; cursor: pointer; }
        button:hover { background: #45a049; }
    </style>
</head>
<body>
    <h1>UniFORM Landmark Refinement</h1>
    <p>Click on the plots to set landmark positions, or enter values manually.</p>
'

  # Add plots for each marker
  for (marker in markers) {
    if (!(marker %in% names(histogram_data))) next

    html_content <- paste0(html_content, '
    <div class="marker-section">
        <h2>', marker, '</h2>
        <div id="plot_', marker, '" class="plot-container"></div>
        <div class="landmark-inputs">
            <h3>Landmark Positions:</h3>')

    # Add input fields for each sample
    for (i in seq_along(sample_ids)) {
      sample_id <- sample_ids[i]
      initial_val <- ""
      if (!is.null(implied_landmarks) && marker %in% names(implied_landmarks)) {
        if (i <= length(implied_landmarks[[marker]])) {
          initial_val <- as.character(round(implied_landmarks[[marker]][i]))
        }
      }

      html_content <- paste0(html_content, '
            <span class="sample-label">', sample_id, ':</span>
            <input type="number" id="landmark_', marker, '_', sample_id, '" value="', initial_val, '" />
            <br/>')
    }

    html_content <- paste0(html_content, '
        </div>
    </div>')
  }

  # Add save button and JavaScript
  html_content <- paste0(html_content, '
    <button onclick="saveLandmarks()">Save Landmarks to CSV</button>
    <button onclick="clearAll()">Clear All</button>

    <script>
    // Plot data for each marker
')

  # Generate JavaScript for plots
  for (marker in markers) {
    if (!(marker %in% names(histogram_data))) next

    html_content <- paste0(html_content, '
    (function() {
        var traces = [];
')

    # Add traces for each sample
    for (i in seq_along(sample_ids)) {
      sample_id <- sample_ids[i]
      if (sample_id %in% names(histogram_data[[marker]])) {
        hist_data <- histogram_data[[marker]][[sample_id]]
        # Use bin_edges if available, otherwise bins
        x_coords <- if (!is.null(hist_data$bin_edges)) hist_data$bin_edges else hist_data$bins

        html_content <- paste0(html_content, '
        traces.push({
            x: [', paste(x_coords, collapse = ","), '],
            y: [', paste(hist_data$counts, collapse = ","), '],
            type: "scatter",
            mode: "lines",
            name: "', sample_id, '",
            line: { width: 2 }
        });
')
      }
    }

    html_content <- paste0(html_content, '
        var layout = {
            title: "', marker, ' Distribution",
            xaxis: { title: "Log(Cell Mean Intensity)" },
            yaxis: { title: "Frequency" },
            hovermode: "x"
        };

        Plotly.newPlot("plot_', marker, '", traces, layout);

        // Add click handler
        document.getElementById("plot_', marker, '").on("plotly_click", function(data) {
            var x = data.points[0].x;
            var sample = data.points[0].data.name;
            document.getElementById("landmark_', marker, '_" + sample).value = Math.round(x * 100) / 100;
        });
    })();
')
  }

  # Add save function
  html_content <- paste0(html_content, '
    function saveLandmarks() {
        var csv = "Marker";
        var samples = [', paste0('"', sample_ids, '"', collapse = ","), '];

        // Header row
        for (var i = 0; i < samples.length; i++) {
            csv += "," + samples[i] + "_landmark";
        }
        csv += "\\n";

        // Data rows
        var markers = [', paste0('"', markers, '"', collapse = ","), '];
        for (var m = 0; m < markers.length; m++) {
            csv += markers[m];
            for (var s = 0; s < samples.length; s++) {
                var input = document.getElementById("landmark_" + markers[m] + "_" + samples[s]);
                csv += "," + (input ? input.value : "");
            }
            csv += "\\n";
        }

        // Download CSV
        var blob = new Blob([csv], { type: "text/csv" });
        var url = URL.createObjectURL(blob);
        var a = document.createElement("a");
        a.href = url;
        a.download = "landmark_annotations.csv";
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);

        alert("Landmarks saved to landmark_annotations.csv");
    }

    function clearAll() {
        if (confirm("Clear all landmark values?")) {
            var inputs = document.querySelectorAll("input[type=number]");
            for (var i = 0; i < inputs.length; i++) {
                inputs[i].value = "";
            }
        }
    }
    </script>
</body>
</html>')

  # Write HTML file
  writeLines(html_content, output_file)
  cat("Created landmark picking HTML:", output_file, "\n")
}

#' Landmark refinement workflow
#' @param histogram_data Histogram data from process_sample_distributions
#' @param gmm_models GMM models
#' @param markers Markers to process
#' @param sample_ids Sample IDs
#' @param data_source Original data source
#' @param num_bins Number of bins
#' @param dpi DPI for plots
#' @param x_limits X-axis limits
#' @param group_size Number of samples per group
#' @param implied_landmarks_map Initial landmarks from registration
#' @param verbose Print messages
#' @param output_directory Output directory for files
landmark_refinement <- function(
  histogram_data,
  gmm_models,
  markers = NULL,
  sample_ids = NULL,
  data_source = NULL,
  num_bins = 1024,
  dpi = 100,
  x_limits = NULL,
  group_size = 4,
  implied_landmarks_map = NULL,
  verbose = TRUE,
  output_directory = "landmark_picking"
) {

  # Create output directory
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }

  # Get markers and samples
  if (is.null(markers)) {
    markers <- names(histogram_data)
  }

  if (is.null(sample_ids)) {
    first_marker <- markers[1]
    sample_ids <- names(histogram_data[[first_marker]])
  }

  cat("Generating landmark refinement files...\n")

  # Create CSV template
  csv_file <- file.path(output_directory, "landmark_annotations.csv")
  csv_data <- data.frame(
    Marker = markers,
    stringsAsFactors = FALSE
  )

  # Add columns for each sample
  for (sample in sample_ids) {
    col_name <- paste0(sample, "_landmark")
    csv_data[[col_name]] <- NA

    # Add implied landmarks if available
    if (!is.null(implied_landmarks_map)) {
      for (i in seq_along(markers)) {
        marker <- markers[i]
        if (marker %in% names(implied_landmarks_map)) {
          landmarks <- implied_landmarks_map[[marker]]
          sample_idx <- which(sample_ids == sample)
          if (length(sample_idx) > 0 && sample_idx <= length(landmarks)) {
            csv_data[i, col_name] <- round(landmarks[sample_idx])
          }
        }
      }
    }
  }

  write.csv(csv_data, csv_file, row.names = FALSE)
  cat("  Created landmark CSV template:", csv_file, "\n")

  # Create filled version (for testing)
  csv_filled <- file.path(output_directory, "landmark_annotations_filled.csv")
  write.csv(csv_data, csv_filled, row.names = FALSE)

  # Generate HTML files for groups of samples
  n_groups <- ceiling(length(sample_ids) / group_size)

  for (g in 1:n_groups) {
    start_idx <- (g - 1) * group_size + 1
    end_idx <- min(g * group_size, length(sample_ids))
    group_samples <- sample_ids[start_idx:end_idx]

    html_file <- file.path(output_directory, paste0("group_", g, "_landmarks.html"))

    if (verbose) {
      cat("  Processing samples", paste(group_samples, collapse = " to "), "\n")
    }

    # Extract data for this group
    group_hists <- list()
    group_gmms <- list()

    for (marker in markers) {
      if (marker %in% names(histogram_data)) {
        group_hists[[marker]] <- list()
        group_gmms[[marker]] <- list()

        for (sample in group_samples) {
          if (sample %in% names(histogram_data[[marker]])) {
            group_hists[[marker]][[sample]] <- histogram_data[[marker]][[sample]]
          }
          if (!is.null(gmm_models) && marker %in% names(gmm_models)) {
            if (sample %in% names(gmm_models[[marker]])) {
              group_gmms[[marker]][[sample]] <- gmm_models[[marker]][[sample]]
            }
          }
        }

        if (verbose) {
          cat("    Marker:", marker, "\n")
        }
      }
    }

    # Create HTML file
    create_landmark_html(
      histogram_data = group_hists,
      gmm_models = group_gmms,
      markers = markers,
      sample_ids = group_samples,
      implied_landmarks = implied_landmarks_map,
      output_file = html_file
    )
  }

  cat("\nLandmark refinement files created in:", output_directory, "\n")
  cat("  1. Open the HTML files in a web browser\n")
  cat("  2. Click on plots or enter values to set landmarks\n")
  cat("  3. Save the landmarks using the 'Save Landmarks to CSV' button\n")
  cat("  4. Update 'landmark_annotations_filled.csv' with your adjustments\n")

  return(invisible(csv_data))
}

#' Read landmark annotations from CSV
#' @param csv_file Path to CSV file
#' @param sample_ids Sample IDs to extract
#' @return List of landmarks by marker
read_landmark_csv <- function(csv_file, sample_ids = NULL) {
  if (!file.exists(csv_file)) {
    stop("Landmark CSV file not found:", csv_file)
  }

  # Read CSV
  landmark_df <- read.csv(csv_file, stringsAsFactors = FALSE)

  # Extract sample names from columns
  landmark_cols <- grep("_landmark$", names(landmark_df), value = TRUE)
  csv_samples <- gsub("_landmark$", "", landmark_cols)

  if (is.null(sample_ids)) {
    sample_ids <- csv_samples
  }

  # Convert to list format
  landmark_map <- list()

  for (i in 1:nrow(landmark_df)) {
    marker <- landmark_df$Marker[i]
    landmarks <- numeric(length(sample_ids))

    for (j in seq_along(sample_ids)) {
      sample <- sample_ids[j]
      col_name <- paste0(sample, "_landmark")

      if (col_name %in% names(landmark_df)) {
        val <- landmark_df[i, col_name]
        if (!is.na(val) && val != "") {
          landmarks[j] <- as.numeric(val)
        } else {
          landmarks[j] <- NA
        }
      } else {
        landmarks[j] <- NA
      }
    }

    # Only add if at least one landmark is defined
    if (!all(is.na(landmarks))) {
      landmark_map[[marker]] <- landmarks
    } else {
      landmark_map[[marker]] <- NULL
    }
  }

  return(landmark_map)
}