#' Plot CpG changes around a central methylation heatmap
#'
#' @description Plot negative CpG changes to the left and positive changes to
#'   the right of a central methylation heatmap. Grey rectangles show the
#'   portions of estimate plus or minus one standard error on each side of
#'   zero. Vertical lines mark nonzero estimates (gray40 by default); black
#'   rectangles mark both edges of significant heatmap cells. This function
#'   copies the data preparation and smoothing of [plot_icr_cpg_diffs_rows()]
#'   but builds its plot independently, without calling that plotting function.
#' @param df_cpg Data frame containing `cpg_id` and the selected coefficient,
#'   standard error, methylation, and significance columns. Subset to one
#'   primary variable and model group before calling; duplicate CpG results
#'   within the ICR are rejected.
#' @param icr_id Single ICR identifier, such as `"ICR_10"`.
#' @param coeff_colname Name of the numeric coefficient column.
#' @param xlabel Optional x-axis label. Defaults to `coeff_colname` when NULL.
#' @param coeff_se_colname Numeric standard error column on the same scale as
#'   the coefficient. Standard errors must be nonnegative.
#' @param title_extra Text appended to the ICR identifier in the plot title.
#' @param spline_window_cpg_size Full smoothing window width in CpG-index units,
#'   an integer of at least 5. Defaults to 5; NULL disables CpG-mode smoothing.
#'   Larger values smooth more. The window is capped at the plotted CpG count.
#'   Uses normalized tricube-weighted means; the argument name is retained.
#' @param spline_window_bp_size Full genomic window width in base pairs, a
#'   positive finite number. Defaults to NULL. When supplied, overrides
#'   `spline_window_cpg_size`, including its validation. For example, 1000 uses
#'   sites within 500 bp on either side of each prediction position.
#' @param export_plot Logical; save the finished plot as a JPEG when TRUE.
#'   Defaults to FALSE. The ggplot object is returned in either case.
#' @param output_path Output directory, required when `export_plot = TRUE`.
#'   Missing directories are created. The filename is
#'   `diff_split_Est_<+ or ->ZF_<High, Medium, or Low>_<icr_id>.jpg`, using the same
#'   metadata as [plot_icr_dotplot()]. Unknown confidence is labeled Unknown.
#'   Images are 4 inches high and 4.5 inches wide. Existing files are overwritten.
#' @param methyl_colname Numeric beta-value column for the methylation heatmap.
#'   Defaults to `"contrast_fitted_low"`. Values are assigned to the nearest
#'   state among 0, 0.5, and 1, labeled `"0%"`, `"50%"`, and `"100%"`.
#'   Exact midpoint ties (0.25 and 0.75) use the lower state. Values outside
#'   the range 0 to 1 use the nearest endpoint. Missing or nonfinite values
#'   are blank.
#' @param methylation_colors Three fill colors for states 0%, 50%, and 100%,
#'   respectively. Defaults to `c("#FFFFFF", "#FFBAFC", "#E300DA")`.
#' @param df_cpg_sig_colname Numeric p-value or logical significance column
#'   for methylation-cell markers. Defaults to `"cpg_glm_adj_pval"`. Missing
#'   or nonfinite results are treated as nonsignificant.
#'   Logical TRUE is significant without applying a threshold.
#' @param cpg_sig_thresh P-value threshold, default 0.05. Numeric p-values must
#'   be strictly below this threshold. Ignored for logical significance flags.
#' @param border_color_methyl Border color for all methylation cells.
#'   Defaults to `"white"`; NA disables these borders. Line width is 0.15 mm.
#' @param border_color_sig Fill color of the rectangles at both edges of
#'   significant methylation cells. The argument name is retained for compatibility.
#'   Defaults to `"black"`; NA hides these rectangles.
#' @param delta_beta_fill_colors Two fill colors for the smoothed methylation
#'   changes: negative first, positive second. Defaults to `c("blue", "red")`.
#'   Accepts color names or hexadecimal colors. Fill alpha remains 0.4.
#'   Has no effect when smoothing is disabled.
#' @param zero_line_color Color of X = 0 reference lines. Defaults to NA,
#'   which hides them. Supply a color name or hexadecimal color to show them;
#'   line width is 1 mm.
#' @param estimate_line_color Color of the vertical per-CpG estimate markers
#'   for `coeff_colname` (for example, mean delta-beta). Defaults to `"gray40"`.
#'   Accepts a color name or hexadecimal color; NA hides these markers.
#' @details CpG membership and positions come from
#'   `tdhia::manifest_v1A2_design_scores`, using `icr_id` and `MAPINFO`.
#'   Only CpGs present in `df_cpg` are plotted. Rows with missing or nonfinite
#'   coefficients, standard errors, or positions are omitted with a warning.
#'   Indices number the remaining CpGs in genomic order from top to bottom;
#'   for large ICRs only a subset of index ticks is displayed.
#'   Smoothing mirrors values and coordinates about the first and last retained
#'   CpGs, without duplicating endpoint observations. Reflections repeat as
#'   needed to cover the smoothing radius. Padded observations are used only
#'   for smoothing; they do not change plotted CpGs, annotations, or counts.
#'   Both modes use a normalized kernel-weighted mean: sum(w * coefficient) /
#'   sum(w), with tricube weights w = (1 - abs(distance / radius)^3)^3 inside
#'   the window and zero weight at or beyond its boundary. CpG mode uses radius
#'   min(spline_window_cpg_size, n_cpgs) / 2 in CpG-index units. In bp mode,
#'   radius is spline_window_bp_size / 2 in genomic coordinates; prediction
#'   positions are interpolated from MAPINFO along the numbered CpG axis.
#'   Predictions include every retained CpG and a grid between them. A single
#'   contributing CpG returns its coefficient. Defined kernel means remain
#'   within the range of their contributing coefficients, avoiding regression
#'   overshoot. Locations with no positive kernel weights are connected by
#'   linear interpolation between neighboring defined means for display;
#'   these connections do not imply additional CpG support. Small ICRs,
#'   including a single valid CpG, are supported.
#'   Both window arguments set to NULL disable smoothing. SEs do not weight
#'   the smoother. Smoothed values are filled to zero using
#'   `delta_beta_fill_colors`, with alpha 0.4.
#'   The fill extends to the outer row boundaries, indices 0.5 and n_cpgs + 0.5,
#'   by holding the first and last smoothed values constant over each half-row.
#'   The two change panels have equal width and matched ranges: negative
#'   changes span -B to 0 on the left, and positive changes span 0 to B on the
#'   right. B is the largest absolute coefficient plus the mean standard error
#'   across the plotted CpGs, without extra axis padding.
#'   If all coefficients and SEs are zero, the half-range defaults to 1.
#'   SE rectangles and smoothed values outside the visible range are clipped
#'   at the panel boundary; their underlying values are retained.
#'   An SE interval crossing zero is split at zero and shown on both sides of
#'   the heatmap, regardless of the sign of the coefficient. Nonzero estimate
#'   lines appear only on their sign's side. Zero coefficients have no displaced
#'   estimate line; their SE intervals are still shown when nonzero.
#'   Smoothing is fitted once to all retained CpGs, then its negative and
#'   positive portions are drawn in their respective panels. No separate
#'   sign-specific model fits or significance calculations are performed.
#'   Each row extends from CpG index minus 0.5 to index plus 0.5. Estimate lines
#'   use the same bounds as the SE rectangles, so their heights adapt to the
#'   number of CpGs and panel height. Line width is fixed; line position marks
#'   the estimate. The y-axis has no expansion beyond the plotted geometry.
#'   Light grey SE rectangles and estimate lines are drawn beneath the
#'   smoothing layers and optional X = 0 reference lines. When enabled with
#'   `zero_line_color`, each change panel has its own 1 mm zero line at the
#'   edge adjoining the central heatmap.
#'   The methylation heatmap uses the same retained, ordered CpGs and row bounds
#'   as both change panels. Its width is 3.25% of their combined width, with no
#'   gap, axis ticks, or legend. Missing annotations do not remove CpGs.
#'   Cell fill denotes methylation state. Significant cells have a filled
#'   rectangle at each heatmap edge spanning the full row height. Each rectangle
#'   extends inward by 25% of the heatmap width, leaving the middle 50% visible.
#'   Rectangles have no outlines; their widths scale with the heatmap width.
#'   All cells use the same thin border, regardless of significance.
#' @returns A ggplot object, or NULL if no finite coefficients, standard errors,
#'   and positions remain. Its `data` contains the plotted rows, genomic
#'   positions, CpG indices, coefficients, standard errors, interval bounds,
#'   ordered `methylation_state`, and logical `significant`.
#'   Files are written only when `export_plot = TRUE`.
#' @export 
plot_icr_cpg_diffs_rows_split <- function(
    df_cpg, icr_id, coeff_colname = "cpg_glm_estimate", xlabel = NULL,
    coeff_se_colname = "cpg_glm_estimate_se", title_extra = "",
    spline_window_cpg_size = 5,
    spline_window_bp_size = NULL, export_plot = FALSE, output_path = NULL,
    methyl_colname = "contrast_fitted_low",
    methylation_colors = c("#FFFFFF", "#FFBAFC", "#E300DA"),
    df_cpg_sig_colname = "cpg_glm_adj_pval", cpg_sig_thresh = 0.05,
    border_color_methyl = "white", border_color_sig = "black",
    delta_beta_fill_colors = c("blue", "red"), zero_line_color = NA,
    estimate_line_color = "gray40") {  #("#1b9e77", "#d95f02"
  required <- c("cpg_id", coeff_colname, coeff_se_colname, methyl_colname,
                df_cpg_sig_colname)
  missing_cols <- setdiff(required, names(df_cpg))
  if (length(missing_cols)) {
    stop("Missing columns in df_cpg: ", paste(missing_cols, collapse = ", "))
  }
  if (length(icr_id) != 1L || is.na(icr_id)) {
    stop("icr_id must be a single nonmissing ICR identifier.")
  }
  if (length(methylation_colors) != 3L || anyNA(methylation_colors)) {
    stop("methylation_colors must contain three nonmissing fill colors.")
  }
  if (!is.character(delta_beta_fill_colors) ||
      length(delta_beta_fill_colors) != 2L || anyNA(delta_beta_fill_colors)) {
    stop("delta_beta_fill_colors must contain two nonmissing colors: ",
         "negative, then positive.")
  }
  if (length(border_color_methyl) != 1L || length(border_color_sig) != 1L) {
    stop("Specify one color for cell borders and one for significance marks, ",
         "or NA to disable either.")
  }
  sig_values <- df_cpg[[df_cpg_sig_colname]]
  if (!is.numeric(df_cpg[[methyl_colname]]) ||
      !(is.numeric(sig_values) || is.logical(sig_values))) {
    stop("Methylation must be numeric; significance must be numeric or logical.")
  }
  if (!is.logical(sig_values) &&
      (!is.numeric(cpg_sig_thresh) || length(cpg_sig_thresh) != 1L ||
       !is.finite(cpg_sig_thresh) || cpg_sig_thresh < 0 || cpg_sig_thresh > 1)) {
    stop("cpg_sig_thresh must be a finite number between 0 and 1.")
  }
  if (!is.null(spline_window_bp_size)) {
    if (!is.numeric(spline_window_bp_size) ||
        length(spline_window_bp_size) != 1L ||
        !is.finite(spline_window_bp_size) || spline_window_bp_size <= 0) {
      stop("spline_window_bp_size must be NULL or a positive finite number.")
    }
  } else if (!is.null(spline_window_cpg_size) &&
      (!is.numeric(spline_window_cpg_size) || length(spline_window_cpg_size) != 1L ||
       !is.finite(spline_window_cpg_size) || spline_window_cpg_size < 5 ||
       spline_window_cpg_size != floor(spline_window_cpg_size))) {
    stop("spline_window_cpg_size must be NULL or an integer of at least 5.")
  }

  # Use the package manifest for both membership and genomic ordering.
  manifest <- tdhia::manifest_v1A2_design_scores
  sites <- unique(manifest[which(manifest$icr_id == icr_id),
                           c("cpg_id", "MAPINFO")])
  if (anyDuplicated(sites$cpg_id)) {
    stop("The manifest contains conflicting positions for CpGs in this ICR.")
  }
  df_plot <- as.data.frame(df_cpg[df_cpg$cpg_id %in% sites$cpg_id, ])
  if (!nrow(df_plot)) stop("No CpG results found for ", icr_id, ".")
  if (anyDuplicated(df_plot$cpg_id)) {
    stop("Multiple results per CpG; subset df_cpg to one primary variable ",
         "and model group before plotting.")
  }

  df_plot$MAPINFO <- sites$MAPINFO[match(df_plot$cpg_id, sites$cpg_id)]
  df_plot$coefficient <- df_plot[[coeff_colname]]
  df_plot$coefficient_se <- df_plot[[coeff_se_colname]]
  if (!is.numeric(df_plot$coefficient) ||
      !is.numeric(df_plot$coefficient_se)) {
    stop("Coefficient and standard error columns must be numeric.")
  }
  if (any(df_plot$coefficient_se < 0, na.rm = TRUE)) {
    stop("Standard errors must be nonnegative.")
  }

  # Filter before numbering so visible CpG indices remain consecutive.
  keep <- is.finite(df_plot$coefficient) &
    is.finite(df_plot$coefficient_se) & is.finite(df_plot$MAPINFO)
  if (any(!keep)) warning("Omitting ", sum(!keep), " CpG(s) with missing ",
                          "or nonfinite coefficients, SEs, or positions.")
  df_plot <- df_plot[keep, , drop = FALSE]
  if (!nrow(df_plot)) return(NULL)
    
  df_plot <- df_plot[order(df_plot$MAPINFO, df_plot$cpg_id), , drop = FALSE]
  df_plot$cpg_index <- seq_len(nrow(df_plot))
  df_plot$lower <- df_plot$coefficient - df_plot$coefficient_se
  df_plot$upper <- df_plot$coefficient + df_plot$coefficient_se

  # Classify annotations after filtering and ordering to retain row alignment.
  methyl_values <- df_plot[[methyl_colname]]
  methyl_values[!is.finite(methyl_values)] <- NA_real_
  df_plot$methylation_state <- cut(
    methyl_values, breaks = c(-Inf, 0.25, 0.75, Inf),
    labels = c("0%", "50%", "100%"), ordered_result = TRUE)
  sig_values <- df_plot[[df_cpg_sig_colname]]
  df_plot$significant <- if (is.logical(sig_values)) {
    !is.na(sig_values) & sig_values
  } else {
    is.finite(sig_values) & sig_values < cpg_sig_thresh
  }

  # Keep ticks readable while retaining the full height of every CpG row.
  n_sites <- nrow(df_plot)
  ticks <- if (n_sites <= 20L) seq_len(n_sites) else
    sort(unique(c(1, pretty(c(1, n_sites), n = 8), n_sites)))
  ticks <- ticks[ticks >= 1 & ticks <= n_sites]
  bound <- max(abs(df_plot$coefficient)) + mean(df_plot$coefficient_se)
  if (bound == 0) bound <- 1
  strip_width <- 0.065 * bound
  strip_start <- 2 * bound
  panel_levels <- c("Negative", "Methylation", "Positive")
  negative_data <- positive_data <- df_plot
  negative_data$plot_panel <- factor("Negative", levels = panel_levels)
  positive_data$plot_panel <- factor("Positive", levels = panel_levels)
  change_data <- rbind(negative_data, positive_data)

  # Split intervals at zero so crossing SEs retain both visible portions.
  interval_data <- change_data[
    (change_data$plot_panel == "Negative" & change_data$lower < 0) |
      (change_data$plot_panel == "Positive" & change_data$upper > 0), ]
  interval_data$lower_plot <- ifelse(interval_data$plot_panel == "Positive",
                                     pmax(interval_data$lower, 0),
                                     interval_data$lower)
  interval_data$upper_plot <- ifelse(interval_data$plot_panel == "Negative",
                                     pmin(interval_data$upper, 0),
                                     interval_data$upper)
  estimate_data <- change_data[
    (change_data$plot_panel == "Negative" & change_data$coefficient < 0) |
      (change_data$plot_panel == "Positive" & change_data$coefficient > 0), ]
  panel_bounds <- data.frame(
    plot_panel = factor(rep(c("Negative", "Positive"), each = 2),
                        levels = panel_levels),
    x = c(-bound, 0, 0, bound), y = rep(c(0.5, n_sites + 0.5), 2))

  # Give each change panel its own zero edge, adjoining the center strip.
  gg <- ggplot2::ggplot(df_plot) +
    ggplot2::geom_rect(
      data = interval_data,
      ggplot2::aes(xmin = lower_plot, xmax = upper_plot,
                    ymin = cpg_index - 0.5, ymax = cpg_index + 0.5),
      inherit.aes = FALSE, fill = "grey80", alpha = 0.4, colour = NA) +
    ggplot2::geom_segment(
      data = estimate_data,
      ggplot2::aes(x = coefficient, xend = coefficient,
                    y = cpg_index - 0.5, yend = cpg_index + 0.5),
      inherit.aes = FALSE, colour = estimate_line_color, linewidth = 0.6,
      lineend = "butt", na.rm = TRUE) +
    # A separate x-domain identifies the center strip without changing labels.
    ggplot2::geom_blank(data = panel_bounds, ggplot2::aes(x = x, y = y),
                         inherit.aes = FALSE) +
    ggplot2::scale_x_continuous(
      limits = function(limits) {
        if (limits[1] < 0) c(-bound, 0) else if (limits[1] >= strip_start) {
          c(strip_start, strip_start + strip_width)
        } else c(0, bound)
      },
      breaks = function(limits) {
        if (limits[1] >= strip_start) numeric() else pretty(limits)
      },
      labels = function(values) {
        ifelse(values == 0, "0", format(values, trim = TRUE))
      }, oob = function(x, range) x) +
    ggplot2::scale_y_reverse(breaks = ticks, expand = c(0, 0)) +
    # Zoom without discarding SE rectangles that extend outside the panel.
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::labs(x = if (is.null(xlabel)) coeff_colname else xlabel,
                  y = "CpG index (genomic order)",
                  title = paste0(icr_id, " ", title_extra)) +
    ggplot2::facet_grid(cols = ggplot2::vars(plot_panel),
                        scales = "free_x", space = "free_x") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(panel.spacing.x = grid::unit(0, "pt"),
                    strip.background = ggplot2::element_blank(),
                    strip.text.x = ggplot2::element_blank())

  # Smooth with normalized kernel weights, then connect unsupported intervals.
  if (!is.null(spline_window_bp_size) || !is.null(spline_window_cpg_size)) {
    trend <- data.frame(cpg_index = sort(unique(c(
      seq(1, n_sites, length.out = max(200, n_sites)), df_plot$cpg_index))))
    # Reflect both values and positions; retain each source once at joins.
    smooth_data <- data.frame(
      position = if (is.null(spline_window_bp_size)) df_plot$cpg_index else
        df_plot$MAPINFO,
      coefficient = df_plot$coefficient, source_index = seq_len(n_sites))
    radius <- if (is.null(spline_window_bp_size)) {
      min(spline_window_cpg_size, n_sites) / 2
    } else spline_window_bp_size / 2
    left_edge <- min(smooth_data$position)
    extent <- max(smooth_data$position) - left_edge
    if (extent > 0) {
      n_mirrors <- ceiling(radius / extent)
      offset <- smooth_data$position - left_edge
      copies <- lapply(seq.int(-n_mirrors, n_mirrors), function(tile) {
        copy <- smooth_data
        copy$position <- left_edge + tile * extent +
          if (tile %% 2 == 0) offset else extent - offset
        copy
      })
      smooth_data <- unique(do.call(rbind, copies))
      smooth_data <- smooth_data[
        smooth_data$position >= left_edge - radius &
          smooth_data$position <= left_edge + extent + radius, ]
      smooth_data <- smooth_data[order(smooth_data$position), ]
    }
    positions <- if (is.null(spline_window_bp_size)) {
      trend$cpg_index
    } else if (n_sites == 1L) {
      rep(df_plot$MAPINFO, nrow(trend))
    } else stats::approx(df_plot$cpg_index, df_plot$MAPINFO,
                          xout = trend$cpg_index)$y
    # Positive normalized weights keep each mean within its contributors' range.
    trend$estimate <- vapply(positions, function(position) {
      distance <- abs(smooth_data$position - position) / radius
      inside <- distance < 1
      if (!any(inside)) return(NA_real_)
      weights <- (1 - distance[inside]^3)^3
      sum(weights * smooth_data$coefficient[inside]) / sum(weights)
    }, numeric(1))
    # Interpolate gaps for display, retaining all defined kernel estimates.
    supported <- is.finite(trend$estimate)
    if (any(!supported)) {
      trend$estimate <- stats::approx(trend$cpg_index[supported],
        trend$estimate[supported], xout = trend$cpg_index)$y
    }
    # Hold endpoint values through the outer half-rows for full-height fill.
    edges <- data.frame(cpg_index = c(0.5, n_sites + 0.5),
                        estimate = trend$estimate[c(1, nrow(trend))])
    trend <- rbind(edges[1, ], trend, edges[2, ])
    # Insert zero crossings so positive and negative fills meet at zero.
    crossing <- which(head(trend$estimate, -1) *
                        tail(trend$estimate, -1) < 0)
    if (length(crossing)) {
      left <- trend[crossing, ]
      right <- trend[crossing + 1L, ]
      zeros <- data.frame(cpg_index = left$cpg_index - left$estimate *
        (right$cpg_index - left$cpg_index) /
        (right$estimate - left$estimate), estimate = 0)
      trend <- rbind(trend, zeros)
      trend <- trend[order(trend$cpg_index), ]
    }
    # Connected estimates use one continuous fill group.
    trend$segment <- 1L
    positive_trend <- negative_trend <- trend
    positive_trend$plot_panel <- factor(rep("Positive", nrow(trend)),
                                        levels = panel_levels)
    negative_trend$plot_panel <- factor(rep("Negative", nrow(trend)),
                                        levels = panel_levels)
    gg <- gg + ggplot2::geom_ribbon(
      data = positive_trend,
      ggplot2::aes(y = cpg_index, xmin = 0, group = segment,
                    xmax = pmax(estimate, 0)),
      inherit.aes = FALSE, orientation = "y", fill = delta_beta_fill_colors[2],
      alpha = 0.4, colour = NA) +
      ggplot2::geom_ribbon(
        data = negative_trend,
        ggplot2::aes(y = cpg_index, group = segment,
                      xmin = pmin(estimate, 0), xmax = 0),
        inherit.aes = FALSE, orientation = "y", fill = delta_beta_fill_colors[1],
        alpha = 0.4, colour = NA)
  }
  # Mark both heatmap edges while keeping the center of each cell unobstructed.
  strip_data <- data.frame(
    plot_panel = factor("Methylation", levels = panel_levels),
    cpg_index = df_plot$cpg_index,
    fill_color = methylation_colors[as.integer(df_plot$methylation_state)])
  strip_data$fill_color[is.na(strip_data$fill_color)] <- "white"
  sig_data <- strip_data[rep(which(df_plot$significant), each = 2), ,
                         drop = FALSE]
  sig_data$xmin <- rep(strip_start + c(0, 0.75) * strip_width,
                       length.out = nrow(sig_data))
  sig_data$xmax <- rep(strip_start + c(0.25, 1) * strip_width,
                       length.out = nrow(sig_data))
  gg <- gg + ggplot2::geom_vline(
    data = data.frame(
      plot_panel = factor(c("Negative", "Positive"), levels = panel_levels),
      intercept = 0),
    ggplot2::aes(xintercept = intercept), colour = zero_line_color,
    linewidth = 1, na.rm = TRUE) +
    ggplot2::geom_rect(
      data = strip_data,
      ggplot2::aes(xmin = strip_start, xmax = strip_start + strip_width,
                    ymin = cpg_index - 0.5,
                    ymax = cpg_index + 0.5, fill = fill_color),
      inherit.aes = FALSE, colour = border_color_methyl, linewidth = 0.15) +
    # Use relative rectangle widths so image dimensions do not change coverage.
    ggplot2::geom_rect(
      data = sig_data,
      ggplot2::aes(xmin = xmin, xmax = xmax,
                    ymin = cpg_index - 0.5, ymax = cpg_index + 0.5),
      inherit.aes = FALSE, fill = border_color_sig, colour = NA,
      na.rm = TRUE) +
    ggplot2::scale_fill_identity()
  .export_icr_cpg_diffs_rows_split(gg, icr_id, export_plot, output_path)
  gg
}

# Export the row plot using the ICR metadata and established JPEG filename.
.export_icr_cpg_diffs_rows_split <- function(gg, icr_id, export_plot, output_path) {
  if (!is.logical(export_plot) || length(export_plot) != 1L ||
      is.na(export_plot)) stop("export_plot must be TRUE or FALSE.")
  if (!export_plot) return(invisible(NULL))
  if (!is.character(output_path) || length(output_path) != 1L ||
      is.na(output_path) || !nzchar(output_path)) {
    stop("output_path must be an output directory when export_plot is TRUE.")
  }
  icr_metadata <- tdhia::add_metadata_to_imp_sites(
    icr_id, imp_type = "icr", db_flag = FALSE)
  zinc_finger_str <- c("-", "+")[as.integer(icr_metadata$is_icr_zinc[1]) + 1L]
  icf_conf_str <- c("1" = "High", "2" = "Medium", "3" = "Low")[
    as.character(icr_metadata$icr_conf[1])]
  if (is.na(icf_conf_str)) icf_conf_str <- "Unknown"
  filename <- file.path(output_path, sprintf("diff_split_Est_%sZF_%s_%s.jpg",
                                            zinc_finger_str, icf_conf_str, icr_id))
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  cowplot::save_plot(filename = filename, plot = gg,
                     base_height = 4, base_width = 4.5)
  invisible(filename)
}

