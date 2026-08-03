## Rough summary plots: for each of the four examples (epil2_simple,
## epil2_complex, report4bb, schizophrenia), one ggplot with parameter
## value on the y-axis, fitting method mapped to colour (dodged
## horizontally so the six methods' B=500 replicate estimates are visually
## separated), and one facet row per parameter (scales="free_y", since
## sd/corr/phi/beta are all on different scales). A dashed horizontal line
## marks the true ("pretty") value used to simulate the data, per facet.
## The four per-dataset plots are combined side by side with patchwork
## (not faceted together, since each dataset has a different parameter
## set).

suppressMessages({
  library(ggplot2)
  library(patchwork)
  library(ggh4x)
})

wd <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/param_survey"

method_levels <- c("glmmTMB", "jointphi", "pirlsdigamma", "pirlsmoment", "lme4current", "lme4old")
method_labels <- c(glmmTMB = "glmmTMB", jointphi = "joint-phi (R)",
                    pirlsdigamma = "PIRLS/digamma (R)", pirlsmoment = "PIRLS/moment (R)",
                    lme4current = "PIRLS/moment (C++)", lme4old = "PIRLS/fixed-phi (CRAN)")

## Okabe-Ito palette, excluding black and yellow (6 colours for 6 methods)
okabe_ito6 <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
names(okabe_ito6) <- method_labels[method_levels]

build_long <- function(sim, results_list) {
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

make_plot <- function(sim, results_list, title) {
  long <- build_long(sim, results_list)
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
    scale_color_manual(values = okabe_ito6) +
    scale_fill_manual(values = okabe_ito6) +
    labs(title = title, x = NULL, y = "estimate", color = "method", fill = "method") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(), strip.text.y.right = element_text(angle = 0))
}

load_results <- function(example) {
  sim <- readRDS(file.path(wd, paste0(example, "_simdata.rds")))
  results_list <- setNames(lapply(method_levels, function(m) {
    readRDS(file.path(wd, paste0(example, "_results_", m, ".rds")))
  }), method_levels)
  list(sim = sim, results = results_list)
}

d_simple <- load_results("epil2_simple")
d_complex <- load_results("epil2_complex")
d_report4bb <- load_results("report4bb")
d_schizo <- load_results("schizophrenia")

p1 <- make_plot(d_simple$sim, d_simple$results, "epil2 (simple)")
p2 <- make_plot(d_complex$sim, d_complex$results, "epil2 (complex)")
p3 <- make_plot(d_report4bb$sim, d_report4bb$results, "Report4BB")
p4 <- make_plot(d_schizo$sim, d_schizo$results, "schizophrenia")

combined <- (p1 | p2 | p3 | p4) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

outfile <- file.path(wd, "param_summary_distrib.png")
ggsave(outfile, combined, width = 17, height = 9, dpi = 130)
cat("saved to", outfile, "\n")

## ---- same parameter summary, but mean +/- 2 SE (geom_pointrange) ----
## instead of the full per-replicate distribution.
summarize_long <- function(long) {
  agg <- aggregate(value ~ method + parameter, data = long,
                    FUN = function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x))))
  data.frame(method = agg$method, parameter = agg$parameter,
             mean = agg$value[, "mean"], se = agg$value[, "se"])
}

make_plot_pointrange <- function(sim, results_list, title) {
  long <- build_long(sim, results_list)
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
    scale_color_manual(values = okabe_ito6) +
    labs(title = title, x = NULL, y = "estimate (mean ± 2 SE)", color = "method") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(), strip.text.y.right = element_text(angle = 0))
}

p1b <- make_plot_pointrange(d_simple$sim, d_simple$results, "epil2 (simple)")
p2b <- make_plot_pointrange(d_complex$sim, d_complex$results, "epil2 (complex)")
p3b <- make_plot_pointrange(d_report4bb$sim, d_report4bb$results, "Report4BB")
p4b <- make_plot_pointrange(d_schizo$sim, d_schizo$results, "schizophrenia")

combined_pr <- (p1b | p2b | p3b | p4b) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

outfile_pr <- file.path(wd, "param_summary_stderr.png")
ggsave(outfile_pr, combined_pr, width = 17, height = 9, dpi = 130)
cat("saved to", outfile_pr, "\n")

## ---- same again, excluding lme4 2.0-6 (the old/buggy method) ----
## lme4 2.0-6 is so far off on several parameters that it stretches each
## panel's free_y scale and makes the five "current" methods' much smaller
## differences hard to see; dropping it lets those differences show up.
drop_old <- function(results_list) results_list[setdiff(names(results_list), "lme4old")]

p1c <- make_plot_pointrange(d_simple$sim, drop_old(d_simple$results), "epil2 (simple)")
p2c <- make_plot_pointrange(d_complex$sim, drop_old(d_complex$results), "epil2 (complex)")
p3c <- make_plot_pointrange(d_report4bb$sim, drop_old(d_report4bb$results), "Report4BB")
p4c <- make_plot_pointrange(d_schizo$sim, drop_old(d_schizo$results), "schizophrenia")

combined_pr_new <- (p1c | p2c | p3c | p4c) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

outfile_pr_new <- file.path(wd, "param_summary_stderr_newonly.png")
ggsave(outfile_pr_new, combined_pr_new, width = 17, height = 9, dpi = 130)
cat("saved to", outfile_pr_new, "\n")

## ---- second plot: elapsed time and -2*logLik, all datasets/methods ----
## Two separate plots (elapsed time, Δ(-2*logLik)) combined vertically
## with patchwork, each facet_wrap'd by dataset with its own free_y scale
## per panel. A single facet_grid(metric ~ dataset) with a shared numeric
## axis per row doesn't work here: even after log-scaling time and
## recentring negll, the ranges differ too much across datasets (e.g.
## report4bb fits take ~0.05-0.7s vs ~100s for epil2-complex PIRLS fits)
## for one axis per row to stay readable -- splitting into two independent
## plots, each free per dataset panel, avoids that.
dataset_levels <- c("epil2 (simple)", "epil2 (complex)", "Report4BB", "schizophrenia")

build_metric_long <- function(name, results_list, extract) {
  do.call(rbind, lapply(names(results_list), function(m) {
    df <- results_list[[m]]
    data.frame(dataset = name, method = m, i = df$i, value = extract(df))
  }))
}

datasets <- list("epil2 (simple)" = d_simple, "epil2 (complex)" = d_complex,
                  "Report4BB" = d_report4bb, "schizophrenia" = d_schizo)

long_time <- do.call(rbind, lapply(names(datasets), function(nm) {
  build_metric_long(nm, datasets[[nm]]$results, function(df) df$time_sec)
}))

## -2*logLik values live on wildly different absolute scales across
## datasets (epil2 ~O(1000), Report4BB ~O(10-50)). Recentre to "excess
## over the best fit found for this dataset, across all methods/
## replicates" -- 0 = as good as the best fit found -- so each dataset's
## panel is comparably interpretable even before the free_y scale kicks in.
long_negll <- do.call(rbind, lapply(names(datasets), function(nm) {
  d <- build_metric_long(nm, datasets[[nm]]$results, function(df) df$negll)
  d$value <- d$value - min(d$value, na.rm = TRUE)
  d
}))

prep_long <- function(d) {
  d <- d[!is.na(d$value), ]
  d$dataset <- factor(d$dataset, levels = dataset_levels)
  d$method <- factor(d$method, levels = method_levels, labels = method_labels[method_levels])
  d
}
long_time <- prep_long(long_time)
long_negll <- prep_long(long_negll)

## Paired, per-replicate Δ(-2*logLik) vs glmmTMB (the reference method):
## for each simulated dataset i, method_negll[i] - glmmTMB_negll[i].
## Unlike the marginal Δ(-2*logLik) above -- which is dominated by
## between-replicate sampling variation in how hard/easy each simulated
## dataset is to fit -- pairing within replicate cancels that shared
## variation out and isolates the much smaller, systematic difference
## between each method and glmmTMB on the *same* data. PIRLS/fixed-phi
## (CRAN)/lme4old is excluded here (not just from the reference role): its
## gap is so much larger than the others' that including it would stretch
## each panel's free_x scale right back to hiding the differences this
## panel exists to show, and it isn't a fair phi-estimator comparison
## anyway (see README).
build_negll_diff_long <- function(name, results_list) {
  ref <- results_list[["glmmTMB"]]
  ref_negll <- setNames(ref$negll, ref$i)
  methods <- setdiff(names(results_list), c("glmmTMB", "lme4old"))
  do.call(rbind, lapply(methods, function(m) {
    df <- results_list[[m]]
    data.frame(dataset = name, method = m, i = df$i,
               value = df$negll - ref_negll[as.character(df$i)])
  }))
}

long_negll_diff <- do.call(rbind, lapply(names(datasets), function(nm) {
  build_negll_diff_long(nm, datasets[[nm]]$results)
}))

method_levels_diff <- setdiff(method_levels, c("glmmTMB", "lme4old"))

prep_long_diff <- function(d) {
  d <- d[!is.na(d$value), ]
  d$dataset <- factor(d$dataset, levels = dataset_levels)
  d$method <- factor(d$method, levels = method_levels_diff, labels = method_labels[method_levels_diff])
  d
}
long_negll_diff <- prep_long_diff(long_negll_diff)

## One replicate (i=3) is a joint-phi convergence failure specific to
## epil2 (simple) -- phi estimate >2 vs ~0.2 for every other replicate/
## method (see earlier discussion) -- with Δ(-2*logLik) ~326, which would
## otherwise stretch that panel's free_x scale and hide everything else.
## Drop it and mark its position with a "*" at the panel's (now-reduced)
## max x-value, lined up with joint-phi's row.
outlier_mask <- long_negll_diff$dataset == "epil2 (simple)" & long_negll_diff$value > 100
long_negll_diff <- long_negll_diff[!outlier_mask, ]

outlier_ann <- data.frame(
  dataset = factor("epil2 (simple)", levels = dataset_levels),
  method = factor("joint-phi (R)", levels = levels(long_negll_diff$method)),
  value = max(long_negll_diff$value[long_negll_diff$dataset == "epil2 (simple)"], na.rm = TRUE))

## methods on the y-axis (mapped, not coord_flip), value on the x-axis;
## facet_wrap (not facet_grid) gives each dataset panel its own free_x
## scale, so flipping the axes doesn't reintroduce the shared-numeric-axis
## problem the facet_grid version had.
make_metric_plot <- function(d, xlab, logx = FALSE, vline0 = FALSE) {
  p <- ggplot(d, aes(y = method, x = value, color = method, fill = method)) +
    geom_violin(alpha = 0.25, color = NA, width = 0.9) +
    geom_boxplot(width = 0.15, alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
    geom_point(alpha = 0.2, size = 1.4, show.legend = FALSE) +
    facet_wrap(~dataset, nrow = 1, scales = "free_x") +
    scale_color_manual(values = okabe_ito6) +
    scale_fill_manual(values = okabe_ito6) +
    labs(x = xlab, y = NULL) +
    theme_bw(base_size = 10)
  if (logx) p <- p + scale_x_log10()
  if (vline0) p <- p + geom_vline(xintercept = 0, linetype = "dashed", color = "black")
  p
}

p_time <- make_metric_plot(long_time, "elapsed time (s)", logx = TRUE)
p_negll <- make_metric_plot(long_negll, "Δ (-2*logLik)")
p_negll_diff <- make_metric_plot(long_negll_diff, "Δ (-2*logLik) vs glmmTMB (paired by replicate)",
                                  vline0 = TRUE) +
  geom_text(data = outlier_ann, aes(x = value, y = method, label = "*"),
            inherit.aes = FALSE, size = 6, fontface = "bold", color = "black", hjust = 0)

## p_negll_diff's data spans only 4 of the 6 methods (glmmTMB and lme4old
## excluded); even with identical scale colours/labels across all three
## plots, patchwork's guides="collect" doesn't reliably merge a legend
## built from 4 present levels with ones built from 6, so instead keep
## exactly one legend (from the top panel, which has all 6 methods) and
## suppress the other two explicitly. (Don't use patchwork's `&` to set
## legend.position here -- it applies to, and would override, all three
## subplots' settings.)
p_time <- p_time + theme(legend.position = "bottom")
p_negll <- p_negll + theme(legend.position = "none")
p_negll_diff <- p_negll_diff + theme(legend.position = "none")
p_time_negll <- p_time / p_negll / p_negll_diff

outfile2 <- file.path(wd, "time-negll_summary.png")
ggsave(outfile2, p_time_negll, width = 12, height = 13, dpi = 130)
cat("saved to", outfile2, "\n")
