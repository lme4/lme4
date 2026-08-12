## Shared, parameterized plotting code for the paramsurvey comparison
## plots. Original `08_summary_plots.R` had a fixed 6-method/5-dataset
## config baked in; this library factors the actual plot-building logic
## out into make_summary_plots(), which takes the METHOD and DATASET
## subsets (and an output filename prefix) as arguments, so different
## callers can compare different subsets of a potentially larger set of
## methods/datasets without editing shared code or overwriting each
## other's PNGs.
##
## Master registries below (`.method_registry`, `.method_palette`,
## `.dataset_label_registry`) are the single place to add a new
## method/dataset -- extend them there, never reassign an existing
## method's color, so plots comparing different subsets stay visually
## consistent (same method = same color everywhere).

suppressMessages({
  library(ggplot2)
  library(patchwork)
  library(ggh4x)
})

.method_registry <- c(
  glmmTMB          = "glmmTMB",
  jointphi         = "joint-phi (R)",
  pirlsdigamma     = "PIRLS/digamma (R)",
  pirlsmoment      = "PIRLS/moment (R)",
  lme4current      = "PIRLS/moment (C++)",
  lme4old          = "PIRLS/fixed-phi (CRAN)",
  juliaMixedModels = "MixedModels.jl (pa/dispersion-again)"
)

## Okabe-Ito, excluding black only (7 colours for 7 registered methods).
## Indexed by method name (not position) so a subset always gets the same
## colour a method has everywhere else.
.method_palette <- setNames(
  c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
  names(.method_registry)
)

.dataset_label_registry <- c(
  epil2_simple  = "epil2 (simple)",
  epil2_complex = "epil2 (complex)",
  epil2_phigt1  = "epil2 (phi>1)",
  report4bb     = "whale_crate",
  schizophrenia = "schizophrenia"
)

build_long <- function(sim, results_list, method_levels, method_labels) {
  beta_names <- names(sim$pretty$beta)
  param_cols <- c("sd1", "sd2", "corr", "phi", beta_names)
  out <- do.call(rbind, lapply(names(results_list), function(m) {
    df <- results_list[[m]]
    do.call(rbind, lapply(param_cols, function(p) {
      data.frame(method = m, i = df$i, parameter = p, value = df[[p]])
    }))
  }))
  out$parameter <- factor(out$parameter, levels = param_cols)
  out$method <- factor(out$method, levels = method_levels, labels = method_labels[method_levels])
  out
}

build_truth <- function(sim) {
  sd <- sim$pretty$sd
  parts <- list(sd1 = unname(sd[1]))
  if (length(sd) > 1) parts$sd2 <- unname(sd[2])
  parts$corr <- sim$pretty$corr
  parts$phi <- sim$pretty$phi
  for (nm in names(sim$pretty$beta)) parts[[nm]] <- unname(sim$pretty$beta[[nm]])
  data.frame(parameter = names(parts), true_value = unlist(parts))
}

make_plot <- function(sim, results_list, title, method_levels, method_labels, palette,
                       show_ylab = TRUE) {
  long <- build_long(sim, results_list, method_levels, method_labels)
  long <- long[!is.na(long$value), ]
  truth <- build_truth(sim)
  truth <- truth[!is.na(truth$true_value) & truth$parameter %in% unique(as.character(long$parameter)), ]
  truth$parameter <- factor(truth$parameter, levels = levels(long$parameter))

  ggplot(long, aes(x = factor(1), y = value, color = method, fill = method)) +
    geom_violin(position = position_dodge(width = 0.8), alpha = 0.25, color = NA, width = 0.9) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.15,
                 alpha = 0.6, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.8), alpha = 0.2, size = 1.4, show.legend = FALSE) +
    geom_hline(data = truth, aes(yintercept = true_value), linetype = "dashed", color = "black",
               inherit.aes = FALSE) +
    facet_wrap(~parameter, ncol = 1, scales = "free_y", strip.position = "right") +
    scale_color_manual(values = palette, drop = FALSE) +
    scale_fill_manual(values = palette, drop = FALSE) +
    labs(title = title, x = NULL, y = if (show_ylab) "estimate" else NULL,
         color = "method", fill = "method") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(), strip.text.y.right = element_text(angle = 0))
}

summarize_long <- function(long) {
  agg <- aggregate(value ~ method + parameter, data = long,
                    FUN = function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x))))
  data.frame(method = agg$method, parameter = agg$parameter,
             mean = agg$value[, "mean"], se = agg$value[, "se"])
}

make_plot_pointrange <- function(sim, results_list, title, method_levels, method_labels, palette,
                                  show_ylab = TRUE) {
  long <- build_long(sim, results_list, method_levels, method_labels)
  long <- long[!is.na(long$value), ]
  summ <- summarize_long(long)
  truth <- build_truth(sim)
  truth <- truth[!is.na(truth$true_value) & truth$parameter %in% unique(as.character(summ$parameter)), ]
  truth$parameter <- factor(truth$parameter, levels = levels(long$parameter))

  ggplot(summ, aes(x = factor(1), y = mean, ymin = mean - 2 * se, ymax = mean + 2 * se,
                    color = method)) +
    geom_pointrange(position = position_dodge(width = 0.8), size = 0.4) +
    geom_hline(data = truth, aes(yintercept = true_value), linetype = "dashed", color = "black",
               inherit.aes = FALSE) +
    facet_wrap(~parameter, ncol = 1, scales = "free_y", strip.position = "right") +
    scale_color_manual(values = palette, drop = FALSE) +
    labs(title = title, x = NULL, y = if (show_ylab) "estimate (mean ± 2 SE)" else NULL,
         color = "method") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(), strip.text.y.right = element_text(angle = 0))
}

load_results <- function(wd, example, methods) {
  sim <- readRDS(file.path(wd, paste0(example, "_simdata.rds")))
  results_list <- setNames(lapply(methods, function(m) {
    readRDS(file.path(wd, paste0(example, "_results_", m, ".rds")))
  }), methods)
  list(sim = sim, results = results_list)
}

build_metric_long <- function(name, results_list, extract) {
  do.call(rbind, lapply(names(results_list), function(m) {
    df <- results_list[[m]]
    data.frame(dataset = name, method = m, i = df$i, value = extract(df))
  }))
}

build_negll_diff_long <- function(name, results_list, ref_method) {
  ref <- results_list[[ref_method]]
  ref_negll <- setNames(ref$negll, ref$i)
  methods <- setdiff(names(results_list), ref_method)
  do.call(rbind, lapply(methods, function(m) {
    df <- results_list[[m]]
    data.frame(dataset = name, method = m, i = df$i,
               value = df$negll - ref_negll[as.character(df$i)])
  }))
}

make_metric_plot <- function(d, xlab, palette, logx = FALSE, vline0 = FALSE, annot = NULL) {
  p <- ggplot(d, aes(y = method, x = value, color = method, fill = method)) +
    geom_violin(alpha = 0.25, color = NA, width = 0.9) +
    geom_boxplot(width = 0.15, alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
    geom_point(alpha = 0.2, size = 1.4, show.legend = FALSE) +
    facet_wrap(~dataset, nrow = 1, scales = "free_x") +
    scale_color_manual(values = palette, drop = FALSE) +
    scale_fill_manual(values = palette, drop = FALSE) +
    labs(x = xlab, y = NULL) +
    theme_bw(base_size = 10)
  if (logx) p <- p + scale_x_log10()
  if (vline0) p <- p + geom_vline(xintercept = 0, linetype = "dashed", color = "black")
  if (!is.null(annot) && nrow(annot) > 0) {
    p <- p + geom_text(data = annot, aes(x = x, y = method, label = label),
                        inherit.aes = FALSE, vjust = -0.6, size = 5, fontface = "bold")
  }
  p
}

#' Build and save the full set of paramsurvey comparison plots for a
#' chosen subset of methods and datasets.
#'
#' @param examples character vector of dataset names (matching
#'   `<example>_simdata.rds`/`<example>_results_<method>.rds` file
#'   naming), in the order they should appear left-to-right.
#' @param methods character vector of method names (must be registered
#'   in `.method_registry` above), in the order/colour-group they should
#'   appear.
#' @param out_prefix filename prefix for all four PNGs produced, so
#'   re-runs with a different method/dataset subset never overwrite an
#'   earlier comparison's output. E.g. out_prefix="julia_compare_" ->
#'   "julia_compare_param_summary_distrib.png", etc.
#' @param negll_ref_method reference method for the paired
#'   Delta(-2*logLik) panel (must be one of `methods`); excluded from
#'   that panel's own method axis (it would trivially be all zeros).
#' @param negll_diff_outliers optional data.frame(example, method,
#'   threshold) of Delta(-2*logLik) points to drop from that panel (per
#'   dataset/method, values > threshold), so a single extreme replicate
#'   doesn't compress the rest of that facet's scale. Each dropped
#'   method/facet gets a "*" annotated just above its remaining points.
#' @param dataset_labels named character vector, dataset name -> display
#'   title; defaults to `.dataset_label_registry`, falling back to the
#'   raw name for anything not registered there.
make_summary_plots <- function(examples, methods,
                                wd = here::here("misc/Gamma_GLMM/paramsurvey"),
                                out_prefix = "",
                                negll_ref_method = "glmmTMB",
                                negll_diff_exclude = character(0),
                                negll_diff_outliers = NULL,
                                dataset_labels = NULL) {
  stopifnot(all(methods %in% names(.method_registry)))
  stopifnot(negll_ref_method %in% methods)
  method_labels <- .method_registry[methods]
  ## palette must be named by the *display labels*, since that's what the
  ## `method` factor's levels are set to everywhere below -- ggplot's
  ## scale_*_manual() matches `values` names against factor levels, not
  ## against the internal method keys.
  palette <- setNames(.method_palette[methods], method_labels)

  if (is.null(dataset_labels)) {
    dataset_labels <- ifelse(examples %in% names(.dataset_label_registry),
                              .dataset_label_registry[examples], examples)
    names(dataset_labels) <- examples
  }

  data_list <- setNames(lapply(examples, load_results, wd = wd, methods = methods), examples)

  ## ---- full per-replicate distribution (violin + boxplot + points) ----
  plots <- lapply(seq_along(examples), function(k) {
    ex <- examples[k]
    make_plot(data_list[[ex]]$sim, data_list[[ex]]$results, dataset_labels[[ex]],
              methods, method_labels, palette, show_ylab = (k == 1))
  })
  combined <- Reduce(`|`, plots) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  outfile <- file.path(wd, paste0(out_prefix, "param_summary_distrib.png"))
  ggsave(outfile, combined, width = 4.2 * length(examples), height = 9, dpi = 130, limitsize = FALSE)
  cat("saved to", outfile, "\n")

  ## ---- mean +/- 2SE version ----
  plots_pr <- lapply(seq_along(examples), function(k) {
    ex <- examples[k]
    make_plot_pointrange(data_list[[ex]]$sim, data_list[[ex]]$results, dataset_labels[[ex]],
                          methods, method_labels, palette, show_ylab = (k == 1))
  })
  combined_pr <- Reduce(`|`, plots_pr) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  outfile_pr <- file.path(wd, paste0(out_prefix, "param_summary_stderr.png"))
  ggsave(outfile_pr, combined_pr, width = 3.2 * length(examples), height = 9, dpi = 130, limitsize = FALSE)
  cat("saved to", outfile_pr, "\n")

  ## ---- elapsed time + paired Delta(-2*logLik) vs negll_ref_method ----
  long_time <- do.call(rbind, lapply(examples, function(ex) {
    build_metric_long(dataset_labels[[ex]], data_list[[ex]]$results, function(df) df$time_sec)
  }))
  long_time <- long_time[!is.na(long_time$value), ]
  long_time$dataset <- factor(long_time$dataset, levels = unname(dataset_labels[examples]))
  long_time$method <- factor(long_time$method, levels = methods, labels = method_labels)

  long_negll_diff <- do.call(rbind, lapply(examples, function(ex) {
    build_negll_diff_long(dataset_labels[[ex]], data_list[[ex]]$results, negll_ref_method)
  }))
  long_negll_diff <- long_negll_diff[!is.na(long_negll_diff$value), ]
  methods_diff <- setdiff(methods, c(negll_ref_method, negll_diff_exclude))
  long_negll_diff <- long_negll_diff[long_negll_diff$method %in% methods_diff, ]

  ## drop flagged outliers (per dataset/method, value > threshold) before
  ## the free_x per-facet scaling is fixed, tracking a "*" annotation
  ## positioned at each dropped method's remaining max.
  outlier_annot <- NULL
  if (!is.null(negll_diff_outliers) && nrow(negll_diff_outliers) > 0) {
    for (k in seq_len(nrow(negll_diff_outliers))) {
      ex_lab <- dataset_labels[[negll_diff_outliers$example[k]]]
      m <- negll_diff_outliers$method[k]
      thresh <- negll_diff_outliers$threshold[k]
      sel <- long_negll_diff$dataset == ex_lab & long_negll_diff$method == m
      excl <- sel & long_negll_diff$value > thresh
      if (any(excl)) {
        remaining_max <- max(long_negll_diff$value[sel & !excl])
        outlier_annot <- rbind(outlier_annot, data.frame(
          dataset = ex_lab, method = method_labels[[m]], x = remaining_max, label = "*"
        ))
        long_negll_diff <- long_negll_diff[!excl, ]
      }
    }
  }

  long_negll_diff$dataset <- factor(long_negll_diff$dataset, levels = unname(dataset_labels[examples]))
  long_negll_diff$method <- factor(long_negll_diff$method, levels = methods_diff,
                                    labels = method_labels[methods_diff])
  if (!is.null(outlier_annot)) {
    outlier_annot$dataset <- factor(outlier_annot$dataset, levels = levels(long_negll_diff$dataset))
    outlier_annot$method <- factor(outlier_annot$method, levels = levels(long_negll_diff$method))
  }

  p_time <- make_metric_plot(long_time, "elapsed time (s)", palette, logx = TRUE) +
    theme(legend.position = "bottom")
  p_negll_diff <- make_metric_plot(long_negll_diff,
                                    paste0("Δ (-2*logLik) vs ", method_labels[[negll_ref_method]],
                                           " (paired by replicate)"),
                                    palette[method_labels[methods_diff]], vline0 = TRUE,
                                    annot = outlier_annot) +
    theme(legend.position = "none")
  p_time_negll <- p_time / p_negll_diff

  outfile2 <- file.path(wd, paste0(out_prefix, "time-negll_summary.png"))
  ggsave(outfile2, p_time_negll, width = 3 * length(examples), height = 9, dpi = 130, limitsize = FALSE)
  cat("saved to", outfile2, "\n")

  invisible(list(distrib = outfile, stderr = outfile_pr, time_negll = outfile2))
}
