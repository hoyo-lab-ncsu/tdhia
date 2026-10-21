#' Plot fitted low and high values at CpGs within an ICR
#'
#' @description Plot two supplied values per CpG, with low values in blue and
#'   high values in red. CpGs run from left to right in genomic order, and
#'   values run along the vertical axis. One asterisk marks each CpG below the
#'   supplied p-value threshold. No smoothing curves are drawn.
#'   Dots match [plot_icr_dotplot()]: size 1.5, alpha 0.5, and the default
#'   solid-circle shape and fill behavior.
#' @param df_cpg Data frame containing `cpg_id`, `is_icr_zinc`,
#'   `Nearest.Transcript`, and the selected value and p-value columns. Subset
#'   to one primary variable and model group before calling; duplicate CpG
#'   results within the ICR are rejected. Supply the QC-passing CpGs before
#'   restricting the displayed window so the title reports the full ICR count.
#' @param icr_id Single ICR identifier, such as `"ICR_10"`.
#' @param low_colname Numeric column plotted as blue dots. Defaults to
#'   `"contrast_fitted_low"`. Values are plotted directly, without averaging.
#' @param high_colname Numeric column plotted as red dots. Defaults to
#'   `"contrast_fitted_high"`. Low and high identify series; values are not
#'   reordered if the blue value exceeds the red value at a CpG.
#' @param cpg_sig_colname Numeric p-value column used to mark CpGs. Defaults to
#'   `"cpg_glm_adj_pval"`. Missing or nonfinite p-values receive no asterisk.
#' @param cpg_sig_thresh P-value threshold between 0 and 1. Defaults to 0.05.
#'   A CpG is marked when its supplied p-value is strictly below this threshold.
#'   No statistical tests or multiple-testing adjustments are performed.
#' @param title_extra Text appended to the ICR count, ZF status, and transcript
#'   title. Include any desired separator, such as `", p-val=0.001"`.
#' @param xlabel Optional vertical value-axis label. The argument name is
#'   retained so existing calls still label the fitted values. Defaults to
#'   `"Model fitted value"` when NULL. Expressions are supported.
#' @param max_sig_hwindow CpG window in genomic index order. NULL (default) or
#'   a single Inf plots all sites. Two nonzero integer endpoints select an
#'   inclusive index range, e.g. `c(10, 25)`. Negative endpoints count from the
#'   end, with -1 denoting the last site; `c(-10, -1)` selects the last 10 sites.
#'   Endpoints are clipped to the available sites. A single nonnegative integer
#'   adds that many neighboring CpGs on each side of the span from the first to
#'   last significant site. With no significant sites, all sites are plotted
#'   with a warning. Indices refer to the complete finite CpG pairs before
#'   windowing, not to rows of the unsorted input data.
#' @param export_plot Logical; save the finished plot as a JPEG when TRUE.
#'   Defaults to FALSE. The ggplot object is returned in either case.
#' @param output_path Output directory, required when `export_plot = TRUE`.
#'   Missing directories are created. The filename is
#'   `dot_Est_<+ or ->ZF_<High, Medium, or Low>_<icr_id>.jpg`, using the ICR
#'   metadata. Unknown confidence is labeled Unknown. Images are 2.5 inches high
#'   and 3 inches wide. Existing files are overwritten.
#' @param legend_position Legend placement passed to `theme(legend.position)`.
#'   Defaults to `"none"`.
#' @details This function prepares its data and builds its plot independently.
#'   CpG membership and positions come from
#'   `tdhia::manifest_v1A2_design_scores`, using `icr_id` and `MAPINFO`.
#'   The title follows `ICR_ID (N CpGs) +ZF/-ZF: transcripts`. N counts all
#'   input rows belonging to the requested ICR before finite-value filtering
#'   and windowing. ZF status comes from the first nonmissing `is_icr_zinc`
#'   value in those rows. Distinct nonmissing, nonempty `Nearest.Transcript`
#'   values are joined with commas. Missing annotations are labeled Unknown.
#'   Only CpGs present in `df_cpg` are plotted. Rows with missing or nonfinite
#'   low values, high values, or positions are omitted with a warning, so every
#'   retained CpG has both dots. Missing p-values do not remove CpGs.
#'   Windowing occurs after filtering and genomic ordering, without recomputing
#'   p-values or their significance. Original CpG indices are retained after
#'   windowing. When the final plot contains 23 or fewer CpGs, every site is
#'   labeled by its `cpg_id`, with rotated labels. Larger selections
#'   show a subset of numeric index ticks.
#'   The value axis includes both series within the selected window, with
#'   padding for asterisks above the values. It is not forced symmetric about
#'   zero because the plotted values are levels rather than differences.
#'   The former `spline_window_cpg_size` and `spline_window_bp_size` arguments
#'   have been removed; remove them from existing calls.
#' @returns A ggplot object, or NULL if no finite CpG pairs remain. Its `data`
#'   contains the retained input rows plus
#'   `MAPINFO`, `cpg_index`, `low_value`, `high_value`, and `significant`.
#'   Files are written only when `export_plot = TRUE`.
#' @export
plot_icr_cpg_dot <- function(
    df_cpg, icr_id, low_colname = "contrast_fitted_low",
    high_colname = "contrast_fitted_high", cpg_sig_colname = "cpg_glm_adj_pval",
    cpg_sig_thresh = 0.05, title_extra = "", xlabel = NULL,
    max_sig_hwindow = NULL, export_plot = FALSE, output_path = NULL, legend_position = "none") {
  
  required <- c("cpg_id", "is_icr_zinc", "Nearest.Transcript", low_colname,
                high_colname, cpg_sig_colname)
  missing_cols <- setdiff(required, names(df_cpg))
  if (length(missing_cols)) {
    stop("Missing columns in df_cpg: ", paste(missing_cols, collapse = ", "))
  }
  if (length(icr_id) != 1L || is.na(icr_id)) {
    stop("icr_id must be a single nonmissing ICR identifier.")
  }
  if (!is.numeric(cpg_sig_thresh) || length(cpg_sig_thresh) != 1L ||
      !is.finite(cpg_sig_thresh) || cpg_sig_thresh < 0 || cpg_sig_thresh > 1) {
    stop("cpg_sig_thresh must be a finite number between 0 and 1.")
  }
  if (!is.null(max_sig_hwindow)) {
    if (!is.numeric(max_sig_hwindow) ||
        !length(max_sig_hwindow) %in% c(1L, 2L) || anyNA(max_sig_hwindow)) {
      stop("max_sig_hwindow must be NULL, a padding size, or two index endpoints.")
    }
    if (length(max_sig_hwindow) == 1L) {
      if (max_sig_hwindow < 0 || max_sig_hwindow != floor(max_sig_hwindow)) {
        stop("Window padding must be a nonnegative integer or Inf.")
      }
    } else if (any(!is.finite(max_sig_hwindow)) || any(max_sig_hwindow == 0) ||
               any(max_sig_hwindow != floor(max_sig_hwindow))) {
      stop("Window endpoints must be finite nonzero integers.")
    }
  }

  # Align each pair of fitted values with its manifest position before sorting.
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

  # Keep the full ICR count and metadata before filtering values or zooming.
  total_cpgs <- nrow(df_plot)
  zinc_values <- df_plot$is_icr_zinc[!is.na(df_plot$is_icr_zinc)]
  zinc_label <- if (!length(zinc_values)) "Unknown ZF" else
    if (zinc_values[1]) "+ZF" else "-ZF"
  transcripts <- unique(trimws(as.character(df_plot$Nearest.Transcript)))
  transcripts <- transcripts[!is.na(transcripts) & nzchar(transcripts)]
  transcript_label <- if (length(transcripts)) {
    paste(transcripts, collapse = ", ")
  } else {"Unknown"}
  title <- sprintf("%s (%d CpGs) %s: %s", icr_id, total_cpgs, zinc_label,
                    transcript_label)

  df_plot$MAPINFO <- sites$MAPINFO[match(df_plot$cpg_id, sites$cpg_id)]
  df_plot$low_value <- df_plot[[low_colname]]
  df_plot$high_value <- df_plot[[high_colname]]
  if (!is.numeric(df_plot$low_value) || !is.numeric(df_plot$high_value) ||
      !is.numeric(df_plot[[cpg_sig_colname]])) {
    stop("Low, high, and p-value columns must be numeric.")
  }

  # Retain complete pairs; unavailable p-values only suppress their asterisks.
  keep <- is.finite(df_plot$low_value) & is.finite(df_plot$high_value) &
    is.finite(df_plot$MAPINFO)
  if (any(!keep)) warning("Omitting ", sum(!keep), " CpG(s) with missing ",
                          "or nonfinite low values, high values, or positions.")
  df_plot <- df_plot[keep, , drop = FALSE]
  if (!nrow(df_plot)) {
    return(NULL);
    save(list = ls(all.names = TRUE), file = "plot_icr_cpg_dot.RData")# load(file = "plot_icr_cpg_dot.RData")
    stop("No finite CpG results remain for ", icr_id, "."); 
  }
    
  df_plot <- df_plot[order(df_plot$MAPINFO, df_plot$cpg_id), , drop = FALSE]
  df_plot$cpg_index <- seq_len(nrow(df_plot))
  p_value <- df_plot[[cpg_sig_colname]]
  df_plot$significant <- is.finite(p_value) & p_value < cpg_sig_thresh

  # Select a genomic window while preserving each site's original index.
  total_sites <- nrow(df_plot)
  if (length(max_sig_hwindow) == 2L) {
    endpoints <- max_sig_hwindow
    endpoints[endpoints < 0] <- total_sites + endpoints[endpoints < 0] + 1
    if (endpoints[1] > endpoints[2]) {
      stop("The first window endpoint must not exceed the last.")
    }
    window <- c(max(1, endpoints[1]), min(total_sites, endpoints[2]))
    if (window[1] > window[2]) stop("max_sig_hwindow selects no CpGs.")
    df_plot <- df_plot[seq.int(window[1], window[2]), , drop = FALSE]
  } else if (length(max_sig_hwindow) == 1L && is.finite(max_sig_hwindow)) {
    sig_indices <- which(df_plot$significant)
    if (!length(sig_indices)) {
      warning("No significant CpGs; plotting all CpGs.")
    } else {
      window <- c(max(1, min(sig_indices) - max_sig_hwindow),
                   min(total_sites, max(sig_indices) + max_sig_hwindow))
      df_plot <- df_plot[seq.int(window[1], window[2]), , drop = FALSE]
    }
  }

  # Label small windows by CpG name; keep numeric indices for larger windows.
  named_sites <- nrow(df_plot) <= 27L
  index_range <- range(df_plot$cpg_index)
  ticks <- if (named_sites) df_plot$cpg_index else
    sort(unique(c(index_range, pretty(index_range, n = 8))))
  ticks <- ticks[ticks >= index_range[1] & ticks <= index_range[2]]
  tick_labels <- if (named_sites) as.character(df_plot$cpg_id) else
    as.character(ticks)

  # Reserve a separate margin for one significance mark per CpG.
  value_range <- range(c(df_plot$low_value, df_plot$high_value))
  padding <- if (diff(value_range) > 0) diff(value_range) * 0.1 else
    max(abs(value_range), 1) * 0.05
  star_position <- value_range[2] + padding
  gg <- ggplot2::ggplot(df_plot, ggplot2::aes(x = cpg_index)) +
    ggplot2::scale_x_continuous(limits = index_range + c(-0.5, 0.5),
                                breaks = ticks, labels = tick_labels,
                                expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = c(Low = "blue", High = "red"),
                                 breaks = c("Low", "High")) +
    ggplot2::coord_cartesian(ylim = c(value_range[1] - padding,
                                     value_range[2] + 2 * padding),
                             expand = FALSE) +
    ggplot2::labs(x = if (named_sites) "CpG site (genomic order)" else
                    "CpG index (genomic order)",
                  y = if (is.null(xlabel)) "Model fitted value" else xlabel,
                  colour = NULL,
                  title = paste0(title, title_extra)) +
    ggplot2::theme_classic(base_size = 7) +
    ggplot2::theme(legend.position = legend_position,
                    plot.title = ggplot2::element_text(size = 6))
  if (named_sites) {
    gg <- gg + ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 45, vjust = 1, hjust = 1))
  }
  gg <- gg + ggplot2::geom_point(ggplot2::aes(y = low_value, colour = "Low"),
                                 size = 1.5, alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(y = high_value, colour = "High"),
                        size = 1.5, alpha = 0.5) +
    ggplot2::geom_text(data = df_plot[df_plot$significant, , drop = FALSE],
                       y = star_position, label = "*", colour = "black", size = 4)

  # Export fitted values separately from the row-difference and raw-beta plots.
  if (!is.logical(export_plot) || length(export_plot) != 1L ||
      is.na(export_plot)) stop("export_plot must be TRUE or FALSE.")
  if (export_plot) {
    if (!is.character(output_path) || length(output_path) != 1L ||
        is.na(output_path) || !nzchar(output_path)) {
      stop("output_path must be an output directory when export_plot is TRUE.")
    }
    metadata <- add_metadata_to_imp_sites(icr_id, imp_type = "icr")
    zinc <- c("-", "+")[as.integer(metadata$is_icr_zinc[1]) + 1L]
    confidence <- c("1" = "High", "2" = "Medium", "3" = "Low")[
      as.character(metadata$icr_conf[1])]
    if (is.na(confidence)) confidence <- "Unknown"
    filename <- file.path(output_path, sprintf("dot_Est_%sZF_%s_%s.jpg",
                                              zinc, confidence, icr_id))
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
    cowplot::save_plot(filename = filename, plot = gg,
                       base_height = 2.5, base_width = 3)
  }
  gg
}
