#' Plot CpG estimates and standard errors as full-height rows
#'
#' @description Plot one coefficient per CpG in genomic order using full-height
#'   rows. Grey rectangles span estimate plus or minus one standard error.
#'   Black vertical lines mark the estimates. This function prepares its data,
#'   fits its optional local trend, and builds its plot independently.
#' @param df_cpg Data frame containing `cpg_id` and the selected coefficient
#'   and standard error columns. Subset to one primary variable and model group
#'   before calling; duplicate CpG results within the ICR are rejected.
#' @param icr_id Single ICR identifier, such as `"ICR_10"`.
#' @param coeff_colname Name of the numeric coefficient column.
#' @param xlabel Optional x-axis label. Defaults to `coeff_colname` when NULL.
#' @param coeff_se_colname Numeric standard error column on the same scale as
#'   the coefficient. Standard errors must be nonnegative.
#' @param title_extra Text appended to the ICR identifier in the plot title.
#' @param spline_window_cpg_size Approximate number of neighboring CpGs used by
#'   LOESS, an integer of at least 5. Defaults to 5; NULL disables smoothing.
#'   Larger values smooth more. The window is capped at the plotted CpG count.
#' @param spline_window_bp_size Full genomic window width in base pairs, a
#'   positive finite number. Defaults to NULL. When supplied, overrides
#'   `spline_window_cpg_size`, including its validation. For example, 1000 uses
#'   sites within 500 bp on either side of each prediction position.
#' @param export_plot Logical; save the finished plot as a JPEG when TRUE.
#'   Defaults to FALSE. The ggplot object is returned in either case.
#' @param output_path Output directory, required when `export_plot = TRUE`.
#'   Missing directories are created. The filename is
#'   `diff_Est_<+ or ->ZF_<High, Medium, or Low>_<icr_id>.jpg`, using the same
#'   metadata as [plot_icr_dotplot()]. Unknown confidence is labeled Unknown.
#'   Images are 4 inches high and 4.5 inches wide. Existing files are overwritten.
#' @details CpG membership and positions come from
#'   `tdhia::manifest_v1A2_design_scores`, using `icr_id` and `MAPINFO`.
#'   Only CpGs present in `df_cpg` are plotted. Rows with missing or nonfinite
#'   coefficients, standard errors, or positions are omitted with a warning.
#'   Indices number the remaining CpGs in genomic order from top to bottom;
#'   for large ICRs only a subset of index ticks is displayed.
#'   In CpG mode, LOESS fits coefficient against CpG index, using degree 1 and
#'   span = min(spline_window_cpg_size / n_cpgs, 1). Despite the argument name,
#'   this is a local regression, not a spline. With fewer than five valid CpGs,
#'   CpG-mode smoothing is skipped with a warning.
#'   In bp mode, local linear regression uses tricube distance weights within
#'   half the specified window width; sites at the boundary have zero weight.
#'   Predictions requiring fewer than two distinct positively weighted genomic
#'   positions are left as gaps, without widening the window. Prediction
#'   positions are interpolated from MAPINFO along the numbered CpG axis.
#'   Both window arguments set to NULL disable smoothing. SEs do not weight
#'   either smoother. Positive smoothed values have transparent red fill to
#'   zero, and negative values have transparent blue fill to zero.
#'   The x-axis is symmetric about zero, with a half-range equal to the largest
#'   absolute coefficient plus the mean standard error across the plotted CpGs,
#'   without extra axis padding.
#'   If all coefficients and SEs are zero, the half-range defaults to 1.
#'   SE rectangles and smoothed values outside the visible range are clipped
#'   at the panel boundary; their underlying values are retained.
#'   Each row extends from CpG index minus 0.5 to index plus 0.5. Estimate lines
#'   use the same bounds as the SE rectangles, so their heights adapt to the
#'   number of CpGs and panel height. Line width is fixed; line position marks
#'   the estimate. The y-axis has no expansion beyond the plotted geometry.
#'   Light grey SE rectangles and estimate lines are drawn beneath the
#'   smoothing layers and the zero reference line.
#' @returns A ggplot object. Its `data` contains the plotted rows, genomic
#'   positions, CpG indices, coefficients, standard errors, and interval bounds.
#'   Files are written only when `export_plot = TRUE`.
#' @export
plot_icr_cpg_diffs_rows <- function(
    df_cpg, icr_id, coeff_colname = "cpg_glm_estimate", xlabel = NULL,
    coeff_se_colname = "cpg_glm_estimate_se", title_extra = "",
    spline_window_cpg_size = 5,
    spline_window_bp_size = NULL, export_plot = FALSE, output_path = NULL) {
  required <- c("cpg_id", coeff_colname, coeff_se_colname)
  missing_cols <- setdiff(required, names(df_cpg))
  if (length(missing_cols)) {
    stop("Missing columns in df_cpg: ", paste(missing_cols, collapse = ", "))
  }
  if (length(icr_id) != 1L || is.na(icr_id)) {
    stop("icr_id must be a single nonmissing ICR identifier.")
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
  if (!nrow(df_plot)) {
    return(NULL);
    save(list = ls(all.names = TRUE), file = "plot_icr_cpg_diffs_rows.RData")# load(file = "plot_icr_cpg_diffs_rows.RData")
    stop("No finite CpG results remain for ", icr_id, "."); return(NULL)
  }
    
  df_plot <- df_plot[order(df_plot$MAPINFO, df_plot$cpg_id), , drop = FALSE]
  df_plot$cpg_index <- seq_len(nrow(df_plot))
  df_plot$lower <- df_plot$coefficient - df_plot$coefficient_se
  df_plot$upper <- df_plot$coefficient + df_plot$coefficient_se

  # Keep ticks readable while retaining the full height of every CpG row.
  n_sites <- nrow(df_plot)
  ticks <- if (n_sites <= 20L) seq_len(n_sites) else
    sort(unique(c(1, pretty(c(1, n_sites), n = 8), n_sites)))
  ticks <- ticks[ticks >= 1 & ticks <= n_sites]
  bound <- max(abs(df_plot$coefficient)) + mean(df_plot$coefficient_se)
  if (bound == 0) bound <- 1

  # Draw row geometry first so smoothing fills and the zero line stay on top.
  gg <- ggplot2::ggplot(df_plot, ggplot2::aes(x = coefficient, y = cpg_index)) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = lower, xmax = upper,
                    ymin = cpg_index - 0.5, ymax = cpg_index + 0.5),
      fill = "grey80", alpha = 0.4, colour = NA) +
    ggplot2::geom_segment(
      ggplot2::aes(x = coefficient, xend = coefficient,
                    y = cpg_index - 0.5, yend = cpg_index + 0.5),
      colour = "black", linewidth = 0.6, lineend = "butt") +
    ggplot2::scale_x_continuous(limits = function(limits) {
      extent <- max(abs(limits))
      if (extent == 0) extent <- 1
      c(-extent, extent)
    }) +
    ggplot2::scale_y_reverse(breaks = ticks, expand = c(0, 0)) +
    # Zoom without discarding SE rectangles that extend outside the panel.
    ggplot2::coord_cartesian(xlim = c(-bound, bound), expand = FALSE) +
    ggplot2::labs(x = if (is.null(xlabel)) coeff_colname else xlabel,
                  y = "CpG index (genomic order)",
                  title = paste0(icr_id, " ", title_extra)) +
    ggplot2::theme_classic(base_size = 10)

  # Fit the local trend directly from this function's ordered CpG data.
  if (!is.null(spline_window_bp_size) || !is.null(spline_window_cpg_size)) {
    if (is.null(spline_window_bp_size) && n_sites < 5L) {
      warning("Skipping LOESS: fewer than five valid CpGs.")
    } else {
      trend <- data.frame(cpg_index = seq(1, n_sites,
                                          length.out = max(200, n_sites)))
      if (is.null(spline_window_bp_size)) {
        fit <- stats::loess(coefficient ~ cpg_index, data = df_plot,
                            span = min(spline_window_cpg_size / n_sites, 1),
                            degree = 1, control = stats::loess.control(
                              surface = "direct"))
        trend$estimate <- as.numeric(stats::predict(fit, newdata = trend))
      } else {
        positions <- if (n_sites == 1L) rep(df_plot$MAPINFO, nrow(trend)) else
          stats::approx(df_plot$cpg_index, df_plot$MAPINFO,
                        xout = trend$cpg_index)$y
        # Center each local fit to avoid large genomic coordinates.
        trend$estimate <- vapply(positions, function(position) {
          distance <- (df_plot$MAPINFO - position) /
            (spline_window_bp_size / 2)
          inside <- abs(distance) < 1
          if (length(unique(df_plot$MAPINFO[inside])) < 2L) return(NA_real_)
          local_x <- distance[inside]
          weights <- (1 - abs(local_x)^3)^3
          fit <- stats::lm.wfit(cbind(1, local_x),
                               df_plot$coefficient[inside], w = weights)
          if (fit$rank < 2L) NA_real_ else unname(fit$coefficients[1])
        }, numeric(1))
      }
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
      # Leave gaps where a local fit cannot be estimated from nearby CpGs.
      trend$segment <- cumsum(!is.finite(trend$estimate))
      trend <- trend[is.finite(trend$estimate), ]
      gg <- gg + ggplot2::geom_ribbon(
        data = trend, ggplot2::aes(y = cpg_index, xmin = 0, group = segment,
                                    xmax = pmax(estimate, 0)),
        inherit.aes = FALSE, orientation = "y", fill = "red",
        alpha = 0.4, colour = NA) +
        ggplot2::geom_ribbon(
          data = trend, ggplot2::aes(y = cpg_index, group = segment,
                                      xmin = pmin(estimate, 0), xmax = 0),
          inherit.aes = FALSE, orientation = "y", fill = "blue",
          alpha = 0.4, colour = NA) +
        ggplot2::geom_path(
          data = trend, ggplot2::aes(x = estimate, y = cpg_index, group = segment),
          inherit.aes = FALSE, colour = "transparent", linewidth = 0.7)
    }
  }
  gg <- gg + ggplot2::geom_vline(xintercept = 0, colour = "black", linewidth = 1)
  .export_icr_cpg_diffs(gg, icr_id, export_plot, output_path)
  gg
}

# Export the row plot using the ICR metadata and established JPEG filename.
.export_icr_cpg_diffs <- function(gg, icr_id, export_plot, output_path) {
  if (!is.logical(export_plot) || length(export_plot) != 1L ||
      is.na(export_plot)) stop("export_plot must be TRUE or FALSE.")
  if (!export_plot) return(invisible(NULL))
  if (!is.character(output_path) || length(output_path) != 1L ||
      is.na(output_path) || !nzchar(output_path)) {
    stop("output_path must be an output directory when export_plot is TRUE.")
  }
  icr_metadata <- add_metadata_to_imp_sites(icr_id, imp_type = "icr")
  zinc_finger_str <- c("-", "+")[as.integer(icr_metadata$is_icr_zinc[1]) + 1L]
  icf_conf_str <- c("1" = "High", "2" = "Medium", "3" = "Low")[
    as.character(icr_metadata$icr_conf[1])]
  if (is.na(icf_conf_str)) icf_conf_str <- "Unknown"
  filename <- file.path(output_path, sprintf("diff_Est_%sZF_%s_%s.jpg",
                                            zinc_finger_str, icf_conf_str, icr_id))
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  cowplot::save_plot(filename = filename, plot = gg,
                     base_height = 4, base_width = 4.5)
  invisible(filename)
}


#' plot_icr_dotplot
#' @description produces a simple dot pot of beta values of CpG sites within a 
#' specified ICR, seperating patients into 2 groups.
#' 
#' @param mat_cpg_beta a matrix of beta values, cpg sites (rows) x patients (col)
#' @param sig_cpgs vector of cpg site IDs that are significant (can be across 
#' all ICRs and not just ICR being plotted).
#' @param df_patient_groups a dataframe with the sample_colname 
#' (default: patient_id) and a "group" column that specifies whether each patient 
#' is in the low (1, blue) or high (2, red) group.
#' @param icr_id the ID of the specific ICR to be plotted of the form ICR_#, i.e., ICR_10.
#' @param xlab_txt xlabel for plot 
#' @param plot_height_width vector of two numbers in inches, for plot width and 
#' height, as defined in cowplot::save_plot().
#' @param max_sig_hwindow dictates the window of CPG sites to be included,
#'   1) NULL: plot all cpg sites (default).
#'   2) Vector of 2 numerics: the lower and upper index of CpG sites to be 
#'   included. If positive, the index is from the start. If negative, index is 
#'   from the end. 
#'   3) Single numeric: how much padding is added around significant CpG sites. Recommended: 20.
#' @param output_path full file system path so save plot to.
#' @param db_flag when true, saves environment variables to disk to a file with 
#' same name as function.
#' @param filter_na_group boolean, when TRUE (default), removes NA group from plots.
#' @param legend.position specify position of legend as specified in ggplot theme
#'  (default = "none").
#' @param ytext text for y axis.
#' @param sample_colname name of column that refers to sample_id found with 
#' df_patient_groups, and also matches the column names in mat_cpg_beta.
#' @param overwrite_plot boolean, when true, overwrite plot to disk.
#' @returns a ggplot function handle to the plot.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
plot_icr_dotplot <- function(mat_cpg_beta, sig_cpgs = NA, df_patient_groups, icr_id, xlab_txt = "", 
                                 plot_height_width = c(5,3), output_path, max_sig_hwindow = NULL, db_flag = FALSE,
                                 filter_na_group = T, legend.position = "none", ytext = "Mean Beta Value", 
                                 sample_colname = "patient_id", overwrite_plot = F) {
  # Create output folder
  dir.create(path = output_path, recursive = TRUE, showWarnings = FALSE)
  if(db_flag) save(list = ls(all.names = TRUE), file = "plot_icr_methylation.RData")
  # load(file = "plot_icr_methylation.RData")
  
  # If sample_colname does not exist in df_patient_groups, use rownames
  if ( !(sample_colname %in% colnames(df_patient_groups))) {
    df_patient_groups  <- df_patient_groups %>% rownames_to_column(sample_colname)
  }
  
  # Get list of CpGs for specified icrs
  df = dplyr::left_join(x = data.frame(cpg_id = rownames(mat_cpg_beta)),
            y = tdhia::manifest_v1A2_design_scores %>% dplyr::select("cpg_id", "icr_id","MAPINFO") %>% dplyr::distinct(),
            by = dplyr::join_by("cpg_id"), keep = FALSE, na_matches = "never", 
            relationship = "one-to-one")
  
  # Subset the cpg_beta matrix
  sub_mat_cpg_beta = mat_cpg_beta[df$icr_id == icr_id,]
  
  
  total_cpgs <- nrow(sub_mat_cpg_beta)
  
  
  # Reorder rows to genomic location
  sub_mat_cpg_beta <- sub_mat_cpg_beta %>% dplyr::arrange(df[df$icr_id == icr_id,]$MAPINFO)
  # ordered_cpg_ids <- df[df$icr_id == icr_id,] %>% arrange(MAPINFO) %>% pull(cpg_id)
  
  # Subset the cpg sites to those around significant cpgs/ specified by user
  if (length(max_sig_hwindow)==2) {
    # If indices are negative, assume index is from end, convert to positive from start
    if (max_sig_hwindow[1] < 0) max_sig_hwindow[1] = nrow(sub_mat_cpg_beta) + max_sig_hwindow[1]
    if (max_sig_hwindow[2] < 0) max_sig_hwindow[2] = nrow(sub_mat_cpg_beta) + max_sig_hwindow[2]
    
    
    sub_mat_cpg_beta <- sub_mat_cpg_beta[max_sig_hwindow[1]:max_sig_hwindow[2],]
  } else if (!is.null(max_sig_hwindow) && (length(max_sig_hwindow)>0) &&
             (nrow(sub_mat_cpg_beta) > 2*max_sig_hwindow+1)) {
    
    sig_inds <- which(rownames(sub_mat_cpg_beta) %in% sig_cpgs)
    
    min_ind <- max(c(min(sig_inds) - max_sig_hwindow, 1))
    max_ind <- min(c(max(sig_inds) + max_sig_hwindow, nrow(sub_mat_cpg_beta)))
    
    sub_mat_cpg_beta <- sub_mat_cpg_beta[min_ind:max_ind,]
    
  }

  # Ordered list of cpg sites (for factor)
  # rownames(sub_mat_cpg_beta)
  # Convert data to long format
  df_long_cpg_beta <- sub_mat_cpg_beta %>% tibble::rownames_to_column("cpg_id") %>% 
    tidyr::pivot_longer(cols = -"cpg_id", names_to = sample_colname)
  
  df_long_cpg_beta <- dplyr::left_join(x = df_long_cpg_beta, y= df_patient_groups, 
                                by = dplyr::join_by({{sample_colname}}),
            keep = FALSE, na_matches = "never", relationship = "many-to-one")
  
  df_long_cpg_beta$cpg_id <- factor(df_long_cpg_beta$cpg_id, levels = rownames(sub_mat_cpg_beta), ordered=TRUE)
  df_long_cpg_beta$group <- factor(df_long_cpg_beta$group, levels = unique(df_long_cpg_beta$group) %>% sort(), ordered=TRUE)
  
  df_summary <- df_long_cpg_beta %>% dplyr::group_by(cpg_id, group) %>% 
    dplyr::summarize(beta_mean = mean(value, na.rm = T), beta_sd = stats::sd(value, na.rm = T),
              beta_sd = stats::sd(value, na.rm = T), beta_sem = stats::sd(value, na.rm = T)/sqrt(dplyr::n()))
  df_summary$cpg_id_rank = as.numeric(df_summary$cpg_id)
  df_summary$xmin = df_summary$cpg_id_rank -0.5
  df_summary$xmax = df_summary$cpg_id_rank +0.5
  df_summary$back_fill = df_summary$cpg_id_rank %% 2 == 0
  if (filter_na_group) df_summary <- df_summary %>% dplyr::filter(!is.na(group))
  
  # Get ICR metadata (closest genes, zinc finger)
  icr_metadata <- add_metadata_to_imp_sites(icr_id, imp_type = "icr")
  zinc_finger_str = c("-","+")[(as.numeric(icr_metadata$is_icr_zinc)+1)]
  icf_conf_str = c("High","Medium", "Low")[icr_metadata$icr_conf]
  
  
  exact_xlim = range(df_summary$beta_mean)
  padded_xlim  = c(exact_xlim[1] - 0.1*diff(exact_xlim),exact_xlim[2] + 0.1*diff(exact_xlim))
  
 
  # Plot methylation across cpg sites in ICR
  gg <- ggplot(data = df_summary, aes(x = cpg_id, y = beta_mean)) +
    geom_point(aes(color = group), size = 1.5, alpha = 0.5) + 
    geom_text(aes(label = ifelse(df_summary$cpg_id %in% sig_cpgs, "*", "")), 
              y = padded_xlim[2], size = 4, vjust=1) + 
    scale_color_manual(values=c("blue", "red")) + #, labels= c("Low","High"), guide = "none") + 
    scale_fill_manual(values=c("grey90", "white"), guide = "none") +
    # scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = (c(min(df_summary$cpg_id_rank), max(df_summary$cpg_id_rank))), ylim = padded_xlim) + #ylim=c(0,1)
      xlab(xlab_txt) + ylab(ytext) + 
    # geom_ribbon(aes(fill = ))+
    # ggtitle(sprintf("%s (%sZF, %s): %s", icr_id, zinc_finger_str,icf_conf_str, 
    #           icr_metadata$Nearest.Transcript)) +
    ggtitle(sprintf("%s (%d CpGs) %sZF: %s", icr_id, total_cpgs, zinc_finger_str,icr_metadata$Nearest.Transcript)) +
    theme_classic(base_size = 7) + theme(axis.text.x = element_text(
      angle = 45,vjust = 1, hjust = 1), plot.title = element_text(size = 6),
      legend.position = legend.position) 
  
  
  # # Connect group means and fill each curve down to zero without stacking.
  # gg2 <- ggplot(data = df_summary,
  #               aes(x = cpg_id_rank, y = beta_mean, group = group)) +
  #   geom_ribbon(aes(ymin = 0, ymax = beta_mean, fill = group),
  #               alpha = 0.2, colour = NA) +
  #   geom_line(aes(color = group), linewidth = 0.5) +
  #   scale_color_manual(values = c("blue", "red")) +
  #   scale_fill_manual(values = c("blue", "red")) +
  #   # Site ranks follow the same genomic order as the dot plot.
  #   scale_x_continuous(breaks = seq_len(nrow(sub_mat_cpg_beta))) +
  #   coord_cartesian(xlim = range(df_summary$cpg_id_rank),
  #                   ylim = padded_xlim) +
  #   xlab("CpG Sites") + ylab(ytext) +
  #   ggtitle(sprintf("%s (%d CpGs) %sZF: %s", icr_id, total_cpgs,
  #                   zinc_finger_str, icr_metadata$Nearest.Transcript)) +
  #   theme_classic(base_size = 7) +
  #   theme(plot.title = element_text(size = 6),
  #         legend.position = legend.position)

  # Export
  plot_path <-
  if (!is.na(output_path) && overwrite_plot) {
    # Print summaries of data to command line as well
    print(gg)
    # print(table(df_patient_groups$group))
    cowplot::save_plot(filename =  paste0(output_path, "/", sprintf(
      "dot_Beta_%sZF_%s_%s", zinc_finger_str, icf_conf_str, icr_id), ".jpg"), 
      plot = gg,base_height = 2, base_width = 2.5)
    # cowplot::save_plot(filename =  paste0(output_path, "/", sprintf(
    #   "line_Beta_%sZF_%s_%s", zinc_finger_str, icf_conf_str, icr_id), ".jpg"), 
    #   plot = gg2,base_height = 2, base_width = 2.5)
  }
  
  # cat(sprintf("%s, cpg_subset: %s\n", icr_id, paste0(max_sig_hwindow, collapse = ": ")))
  
  # Export
  return(list(dot_plot = gg, df_summary = df_summary, cpg_beta_plotted = sub_mat_cpg_beta))
}







#' plot_icr_diffbar
#' @description produces a simple dot pot of beta values of CpG sites within a 
#' specified ICR, seperating patients into 2 groups.
#' 
#' Note: cpg sites may appear out of order, but that is because they are sorted 
#' by their genomic position (MAPINFO) and not label
#' 
#' @param mat_cpg_beta a matrix of beta values, cpg sites (rows) x patients (col)
#' @param sig_cpgs vector of cpg site IDs that are significant (can be across 
#' all ICRs and not just ICR being plotted).
#' @param df_patient_groups a dataframe with the  
#' 1) a sample column  name (default: patient_id, specified by sample_colname) 
#' 2) diff_group: column that specifies whether each patient 
#' is in the low (1) or high (2) group.
#' 3) subset_group: column that specifies different subsets of patient population
#'  (combined will be plotted also).
#' @param icr_id the ID of the specific ICR to be plotted of the form ICR_#, i.e., ICR_10.
#' @param xlab_txt xlabel for plot 
#' @param plot_height_width vector of two numbers in inches, for plot width and 
#' height, as defined in cowplot::save_plot().
#' @param max_sig_hwindow dictates the window of CPG sites to be included,
#'   1) NULL: plot all cpg sites (default).
#'   2) Vector of 2 numerics: the lower and upper index of CpG sites to be 
#'   included. If positive, the index is from the start. If negative, index is 
#'   from the end. 
#'   3) Single numeric: how much padding is added around significant CpG sites.0
#' @param output_path full file system path so save plot to.
#' @param db_flag when true, saves environment variables to disk to a file with 
#' same name as function.
#' @param filter_na_group boolean, when TRUE (default), removes NA group from plots.
#' @param legend.position specify position of legend as specified in ggplot theme
#'  (default = "none").
#' @param sample_colname name of column that refers to sample_id found with 
#' df_patient_groups, and also matc the column names in mat_cpg_beta.
#' @param overwrite_plot boolean, when true, plot is overwritten to disk.
#' @returns a ggplot function handle to the plot.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
plot_icr_diffbar <- function(mat_cpg_beta, sig_cpgs = NA, df_patient_groups, icr_id, xlab_txt = "", 
                             plot_height_width = c(3,2), output_path, max_sig_hwindow = NULL, db_flag = FALSE,
                             filter_na_group = T, legend.position = "none", 
                             sample_colname = "patient_id", overwrite_plot = T) {
  # Create output folder
  dir.create(path = output_path, recursive = TRUE, showWarnings = FALSE)
  if(db_flag) save(list = ls(all.names = TRUE), file = "plot_icr_diffbar.RData")
  # load(file = "plot_icr_diffbar.RData")
  
  # Get list of CpGs for specified icrs
  df = dplyr::left_join(x = data.frame(cpg_id = rownames(mat_cpg_beta)),
                 y = tdhia::manifest_v1A2_design_scores %>% dplyr::select("cpg_id", "icr_id","MAPINFO") %>% dplyr::distinct(),
                 by = dplyr::join_by("cpg_id"), keep = FALSE, na_matches = "never", 
                 relationship = "one-to-one")
  
  # Subset the cpg_beta matrix
  sub_mat_cpg_beta = mat_cpg_beta[df$icr_id == icr_id,]
  
  total_cpgs <- nrow(sub_mat_cpg_beta)
  
  # Reorder rows to genomic location, assuming same chromosomes
  sub_mat_cpg_beta <- sub_mat_cpg_beta %>% dplyr::arrange(df[df$icr_id == icr_id,]$MAPINFO)
  # ordered_cpg_ids <- df[df$icr_id == icr_id,] %>% arrange(MAPINFO) %>% pull(cpg_id)
  
  # Subset the cpg sites to those around significant cpgs/ specified by user
  if (length(max_sig_hwindow)==2) {
    # If indices are negative, assume index is from end, convert to positive from start
    if (max_sig_hwindow[1] < 0) max_sig_hwindow[1] = nrow(sub_mat_cpg_beta) + max_sig_hwindow[1]
    if (max_sig_hwindow[2] < 0) max_sig_hwindow[2] = nrow(sub_mat_cpg_beta) + max_sig_hwindow[2]
    
    
    sub_mat_cpg_beta <- sub_mat_cpg_beta[max_sig_hwindow[1]:max_sig_hwindow[2],]
  } else if (!is.null(max_sig_hwindow) && (length(max_sig_hwindow)>0) &&
             (nrow(sub_mat_cpg_beta) > 2*max_sig_hwindow+1)) {
    
    sig_inds <- which(rownames(sub_mat_cpg_beta) %in% sig_cpgs)
    
    min_ind <- max(c(min(sig_inds) - max_sig_hwindow, 1))
    max_ind <- min(c(max(sig_inds) + max_sig_hwindow, nrow(sub_mat_cpg_beta)))
    
    sub_mat_cpg_beta <- sub_mat_cpg_beta[min_ind:max_ind,]
    
  }

  # Ordered list of cpg sites (for factor)
  # rownames(sub_mat_cpg_beta)
  # Convert data to long format
  df_long_cpg_beta <- sub_mat_cpg_beta %>% tibble::rownames_to_column("cpg_id") %>% 
    tidyr::pivot_longer(cols = -"cpg_id", names_to = sample_colname)
  
  df_long_cpg_beta <- dplyr::left_join(x = df_long_cpg_beta, y= df_patient_groups, by = dplyr::join_by({{sample_colname}}),
                                keep = FALSE, na_matches = "never", relationship = "many-to-one")
  
  df_long_cpg_beta$cpg_id <- factor(df_long_cpg_beta$cpg_id, 
                                    levels = rev(rownames(sub_mat_cpg_beta)), ordered=TRUE)
  df_long_cpg_beta$diff_group <- factor(df_long_cpg_beta$diff_group,
                                        levels = unique(df_long_cpg_beta$diff_group) %>% sort(), ordered=TRUE)
  df_long_cpg_beta$subset_group <- factor(df_long_cpg_beta$subset_group,
                                          levels = unique(df_long_cpg_beta$subset_group) %>% sort(), ordered=TRUE)
  
  
  # Calculate beta difference across all patients
  df_summary_all <- df_long_cpg_beta %>% dplyr::group_by(cpg_id) %>% 
    dplyr::summarize(
      subset_group = "All",
      beta_mean_diff = mean(value[diff_group==2], na.rm = T) - mean(value[diff_group==1], na.rm = T),
      n1 = sum(diff_group==1), n2 = sum(diff_group==2),
      beta_sd_diff = sqrt(stats::sd(value[diff_group==1], na.rm = T)^2 +
                            stats::sd(value[diff_group==2], na.rm = T)^2),
      beta_sem_diff = sqrt(stats::sd(value[diff_group==1], na.rm = T)^2/n1 +
                             stats::sd(value[diff_group==2], na.rm = T)^2/n2))
  df_summary_all$cpg_id_rank = as.numeric(df_summary_all$cpg_id)
  df_summary_all$xmin = df_summary_all$cpg_id_rank -0.5
  df_summary_all$xmax = df_summary_all$cpg_id_rank +0.5
  df_summary_all$back_fill = df_summary_all$cpg_id_rank %% 2 == 0
  df_summary_all$cpg_sig = df_summary_all$cpg_id %in% sig_cpgs
  
  # Get difference in beta for each of the sample subgroups
  df_summary_sub <- df_long_cpg_beta %>% dplyr::group_by(cpg_id, subset_group) %>% 
    dplyr::summarize(
      beta_mean_diff = mean(value[diff_group==2], na.rm = T) - mean(value[diff_group==1], na.rm = T),
      n1 = sum(diff_group==1), n2 = sum(diff_group==2),
      beta_sd_diff = sqrt(stats::sd(value[diff_group==1], na.rm = T)^2 +
                            stats::sd(value[diff_group==2], na.rm = T)^2),
      beta_sem_diff = sqrt(stats::sd(value[diff_group==1], na.rm = T)^2/n1 +
                             stats::sd(value[diff_group==2], na.rm = T)^2/n2))
  
  temp <- df_summary_all %>% dplyr::select(c("cpg_id_rank", "xmin", "xmax", "back_fill", "cpg_sig"))
  df_summary_sub <- cbind(df_summary_sub, rbind(temp,temp))
  # Bind summary stat of all group and subset groups
  df_summary <- rbind(df_summary_all, df_summary_sub)
  # Enforce factor level order (reverse because y-axis inverted in plotting)
  df_summary$subset_group <- factor(
    df_summary$subset_group, levels = rev(c("All", levels(df_long_cpg_beta$subset_group))), ordered = TRUE)

  
  # Get ICR metadata (closest genes, zinc finger)
  icr_metadata <- add_metadata_to_imp_sites(icr_id, imp_type = "icr")
  zinc_finger_str = c("-","+")[(as.numeric(icr_metadata$is_icr_zinc)+1)]
  icf_conf_str = c("High","Medium", "Low")[icr_metadata$icr_conf]

  exact_xlim = range(df_summary$beta_mean_diff)
  padded_xlim  = c(exact_xlim[1] - 0.5*diff(exact_xlim),exact_xlim[2] + 0.5*diff(exact_xlim))
  

  # Plot methylation across cpg sites in ICR
  gg <- ggplot(data = df_summary, aes(y = cpg_id, x = beta_mean_diff, )) +
    geom_rect(data = df_summary_all,  aes(ymin = xmin, ymax = xmax, xmin = -Inf, xmax = +Inf), 
               fill = ifelse(df_summary_all$back_fill, "grey92", "white"),
               color = ifelse(df_summary_all$cpg_sig, "black", NA), linewidth = 0.25) +
    geom_col(aes(fill = subset_group),alpha = 1, position = "dodge", width = 1) + 
    scale_fill_manual(values = c("All" = "black", "1" = "#f03b20", "2" = "#67a9cf"))+
    coord_cartesian(ylim = c(min(df_summary$cpg_id_rank)-0.5,
                             max(df_summary$cpg_id_rank)+0.5), xlim = padded_xlim, expand = c(0,0)) +
    geom_vline(xintercept = 0, color = "black")+
    ylab(xlab_txt) + xlab("Mean Beta Value") + 
    ggtitle(sprintf("%s (%d CpGs): %s", icr_id, total_cpgs, icr_metadata$Nearest.Transcript)) +
    theme_classic(base_size = 7) + 
    theme(axis.text.x = element_text(vjust = 0.5, hjust = 1),
          axis.text.y = element_text(colour = ifelse( df_summary$cpg_sig, "black", "grey50"),
                                     face = "bold"),
      plot.title = element_text(size = 7), legend.position = legend.position)
  gg
  
  
  plot_path <- paste0(output_path, "/", 
                      sprintf("Beta_%sZF_%s_%s", zinc_finger_str, icf_conf_str, icr_id), ".jpg")
  if (overwrite_plot | !file.exists(plot_path)) {
    # Print summaries of data to command line as well
    print(gg)
    cowplot::save_plot(filename = plot_path, plot = gg, base_height = plot_height_width[1], base_width = plot_height_width[1])
  }

  cat(sprintf("%s, cpg_subset: %s", icr_id, max_sig_hwindow))
  # Export
  return(list(plot = gg, df_summary = df_summary, cpg_beta_plotted = sub_mat_cpg_beta))
}

#' Count shared and exclusive significant ICRs between response variables
#'
#' @param df_icr ICR results containing `icr_id`, `primary_var`, `model_group`,
#'   and the significance column selected by `icr_sig_colname`.
#' @param primary_var_ref Single reference response variable name, matched
#'   against `primary_var` in `df_icr`.
#' @param primary_vars Character vector of response variable names to compare
#'   with the reference. Comparisons retain the supplied order.
#' @param icr_sig_colname Name of a logical, numeric binary, or numeric p-value
#'   column used for selection. Logical TRUE and numeric 1 mean significant.
#'   Numeric columns whose observed values are all 0 or 1 are treated as binary.
#'   Otherwise, supply p-values (adjusted p-values for adjusted significance).
#' @param icr_sig_pval_threh P-value threshold, default 0.05. P-values must be
#'   strictly below this threshold. Ignored for binary significance columns.
#' @param model_group Single model group to select from `df_icr`.
#' @param output_path Optional output directory. NULL (default) disables CSV
#'   export. Otherwise, missing directories are created and both tables are
#'   saved as `<primary_var_ref>_shared_icr_summary.csv` and
#'   `<primary_var_ref>_shared_icrs_long.csv`, without row names. Existing files
#'   with these names are overwritten.
#' @details Selection uses the supplied significance results without
#'   recalculating them. Binary detection uses the full supplied column,
#'   ignoring missing values, before selecting the model group. Missing results,
#'   nonfinite p-values, and missing or empty ICR identifiers are excluded.
#'   Each significant ICR is counted once per response variable,
#'   even if it appears in multiple rows. Requested responses must be present
#'   in the selected model group; responses with no significant ICRs are
#'   represented by empty sets.
#' @returns A named list with two data frames:
#'   \describe{
#'     \item{summary}{One row per entry in `primary_vars`, with columns
#'       `primary_a` (the reference), `primary_b` (the comparison), `n_a` (ICRs
#'       significant only for A), `n_b` (ICRs significant only for B), and
#'       `n_shared` (ICRs significant for both).}
#'     \item{shared_icrs}{Long table with columns `primary_a`, `primary_b`, and
#'       `shared_icr_id`, followed by metadata from
#'       [add_metadata_to_imp_sites()] with `imp_type = "icr"`: `is_icr_zinc`,
#'       `ZFP57`, `ZFP445`, `icr_conf`, `Genomic.Coordinates`,
#'       `Nearest.Transcript`, and `Distance.to.Nearest.Transcript`.
#'       Each response-pair/ICR combination appears once.
#'       ICRs shared across several comparisons appear once per pair. Pairs
#'       without shared ICRs have no rows. An empty result retains the key and
#'       metadata columns. ICRs without annotations are retained, using the
#'       metadata helper's missing-value and membership defaults. Comparisons
#'       retain their first supplied order, and ICRs follow their order in the
#'       reference results.}
#'   }
#'   Access the counts with `$summary`; the function previously returned the
#'   counts data frame directly. Both tables are returned regardless of whether
#'   CSV export is enabled.
#' @examples
#' df_icr <- data.frame(
#'   primary_var = c("response_a", "response_a", "response_b", "response_b"),
#'   model_group = "all", icr_id = c("ICR_1", "ICR_2", "ICR_2", "ICR_3"),
#'   adjusted_pval = c(0.01, 0.02, 0.03, 0.04))
#' overlap <- shared_icrs_reponse_vars(
#'   df_icr, "response_a", "response_b", "adjusted_pval", 0.05, "all")
#' overlap$summary
#' overlap$shared_icrs
#' df_icr$significant <- df_icr$adjusted_pval < 0.05
#' shared_icrs_reponse_vars(
#'   df_icr, "response_a", "response_b", "significant", model_group = "all")
#' @export
shared_icrs_reponse_vars <- function(
    df_icr, primary_var_ref, primary_vars, icr_sig_colname,
    icr_sig_pval_threh = 0.05, model_group, output_path = NULL) {
  required <- c("icr_id", "primary_var", "model_group", icr_sig_colname)
  missing_cols <- setdiff(required, names(df_icr))
  if (length(missing_cols)) {
    stop("Missing columns in df_icr: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.character(primary_var_ref) || length(primary_var_ref) != 1L ||
      !is.character(primary_vars) ||
      anyNA(c(primary_var_ref, primary_vars))) {
    stop("Supply one reference name and a vector of nonmissing response names.")
  }
  if (length(model_group) != 1L || is.na(model_group)) {
    stop("model_group must select one nonmissing model group.")
  }
  sig_values <- df_icr[[icr_sig_colname]]
  if (!is.numeric(sig_values) && !is.logical(sig_values)) {
    stop("icr_sig_colname must select a logical or numeric column.")
  }
  observed_values <- sig_values[!is.na(sig_values)]
  is_binary <- is.logical(sig_values) ||
    (length(observed_values) > 0L && all(observed_values %in% c(0, 1)))
  if (!is_binary && (!is.numeric(icr_sig_pval_threh) ||
      length(icr_sig_pval_threh) != 1L || !is.finite(icr_sig_pval_threh) ||
      icr_sig_pval_threh < 0 || icr_sig_pval_threh > 1)) {
    stop("icr_sig_pval_threh must be a single number between 0 and 1.")
  }

  # Restrict the model before identifying each response's significant ICRs.
  df_group <- df_icr[which(df_icr$model_group == model_group), , drop = FALSE]
  if (!nrow(df_group)) stop("No results for model_group: ", model_group)
  missing_vars <- setdiff(c(primary_var_ref, primary_vars), df_group$primary_var)
  if (length(missing_vars)) {
    stop("Response variables absent from selected model group: ",
         paste(missing_vars, collapse = ", "))
  }
  # Binary flags bypass the threshold, including its validation.
  sig_values <- df_group[[icr_sig_colname]]
  significant <- if (is_binary) {
    !is.na(sig_values) & sig_values == 1
  } else {
    is.finite(sig_values) & sig_values < icr_sig_pval_threh
  }
  keep <- significant &
    !is.na(df_group$icr_id) & nzchar(as.character(df_group$icr_id))
  df_sig <- df_group[which(keep), , drop = FALSE]
  reference_icrs <- unique(as.character(df_sig$icr_id[
    which(df_sig$primary_var == primary_var_ref)]))

  # Preallocate to preserve column types, including for an empty comparison list.
  result <- data.frame(
    primary_a = rep(primary_var_ref, length(primary_vars)),
    primary_b = primary_vars, n_a = integer(length(primary_vars)),
    n_b = integer(length(primary_vars)), n_shared = integer(length(primary_vars)),
    stringsAsFactors = FALSE)
  shared_rows <- vector("list", length(primary_vars))
  for (i in seq_along(primary_vars)) {
    comparison_icrs <- unique(as.character(df_sig$icr_id[
      which(df_sig$primary_var == primary_vars[i])]))
    shared_icrs <- intersect(reference_icrs, comparison_icrs)
    result$n_a[i] <- length(setdiff(reference_icrs, comparison_icrs))
    result$n_b[i] <- length(setdiff(comparison_icrs, reference_icrs))
    result$n_shared[i] <- length(shared_icrs)
    shared_rows[[i]] <- data.frame(
      primary_a = rep(primary_var_ref, length(shared_icrs)),
      primary_b = rep(primary_vars[i], length(shared_icrs)),
      shared_icr_id = shared_icrs, stringsAsFactors = FALSE)
  }

  # Keep one row per pair/ICR, including when comparison names are repeated.
  shared_table <- if (length(shared_rows)) {
    unique(do.call(rbind, shared_rows))
  } else {
    data.frame(primary_a = character(), primary_b = character(),
               shared_icr_id = character(), stringsAsFactors = FALSE)
  }
  rownames(shared_table) <- NULL

  # Annotate each distinct ICR once, then match it across response pairs.
  metadata <- tdhia::add_metadata_to_imp_sites(
    unique(shared_table$shared_icr_id), imp_type = "icr", db_flag = FALSE)
  shared_table <- dplyr::left_join(
    shared_table, metadata, by = c("shared_icr_id" = "icr_id"),
    relationship = "many-to-one")

  # Export the returned tables, retaining column headers for empty results.
  if (!is.null(output_path)) {
    if (!is.character(output_path) || length(output_path) != 1L ||
        is.na(output_path) || !nzchar(output_path)) {
      stop("output_path must be NULL or a single nonempty directory path.")
    }
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(result, file.path(output_path,
      paste0(primary_var_ref, "_shared_icr_summary.csv")), row.names = FALSE)
    utils::write.csv(shared_table, file.path(output_path,
      paste0(primary_var_ref, "_shared_icrs_long.csv")), row.names = FALSE)
  }
  return(list(summary = result, shared_icrs = shared_table))
}
