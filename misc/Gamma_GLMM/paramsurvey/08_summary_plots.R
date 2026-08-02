## Rough summary plots: for each of the three examples (epil2_simple,
## epil2_complex, report4bb), one ggplot with parameter value on the
## y-axis, fitting method mapped to colour (jitter-dodged horizontally so
## the six methods' B=10 replicate estimates are visually separated), and
## one facet row per parameter (scales="free_y", since sd/corr/phi/beta
## are all on different scales). A dashed horizontal line marks the true
## ("pretty") value used to simulate the data, per facet. The three
## per-dataset plots are combined side by side with patchwork (not
## faceted together, since each dataset has a different parameter set).

suppressMessages({
  library(ggplot2)
  library(patchwork)
  library(ggh4x)
})

wd <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/param_survey"

method_levels <- c("glmmTMB", "jointphi", "pirlsdigamma", "pirlsmoment", "lme4current", "lme4old")
method_labels <- c(glmmTMB = "glmmTMB", jointphi = "joint-phi",
                    pirlsdigamma = "PIRLS/digamma", pirlsmoment = "PIRLS/moment",
                    lme4current = "lme4 current", lme4old = "lme4 2.0-6")

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

p1 <- make_plot(d_simple$sim, d_simple$results, "epil2 (simple)")
p2 <- make_plot(d_complex$sim, d_complex$results, "epil2 (complex)")
p3 <- make_plot(d_report4bb$sim, d_report4bb$results, "Report4BB")

combined <- (p1 | p2 | p3) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

outfile <- file.path(wd, "summary_plot.png")
ggsave(outfile, combined, width = 13, height = 9, dpi = 130)
cat("saved to", outfile, "\n")

## ---- second plot: elapsed time and -2*logLik, all datasets/methods ----
## facet_grid this time (unlike the per-parameter plot above) since time
## and negll are defined uniformly for every dataset/method -- rows =
## metric (time, -2*logLik), columns = dataset, x = method (jittered),
## colour = method, both axes free per panel given the huge scale
## differences (e.g. report4bb fits take ~0.05-0.7s, epil2-complex PIRLS
## fits take ~100s; epil2 negll is ~O(1000), report4bb negll is ~O(10-50)).
dataset_levels <- c("epil2 (simple)", "epil2 (complex)", "Report4BB")

build_time_negll_long <- function(name, results_list) {
  out <- do.call(rbind, lapply(names(results_list), function(m) {
    df <- results_list[[m]]
    rbind(
      data.frame(dataset = name, method = m, i = df$i,
                 metric = "elapsed time (s)", value = df$time_sec),
      data.frame(dataset = name, method = m, i = df$i,
                 metric = "Δ (-2*logLik)", value = df$negll)
    )
  }))
  ## -2*logLik values live on wildly different absolute scales across
  ## datasets (epil2 ~O(1000), Report4BB ~O(10-50)), which makes them
  ## unreadable sharing a facet_grid row's y-axis with other datasets.
  ## Recentre to "excess over the best fit found for this dataset, across
  ## all methods/replicates" -- comparable scale across datasets, and
  ## more directly informative (0 = as good as the best fit found).
  is_ll <- out$metric == "Δ (-2*logLik)"
  out$value[is_ll] <- out$value[is_ll] - min(out$value[is_ll], na.rm = TRUE)
  out
}

long2 <- rbind(
  build_time_negll_long("epil2 (simple)", d_simple$results),
  build_time_negll_long("epil2 (complex)", d_complex$results),
  build_time_negll_long("Report4BB", d_report4bb$results)
)
long2 <- long2[!is.na(long2$value), ]
long2$dataset <- factor(long2$dataset, levels = dataset_levels)
long2$metric <- factor(long2$metric, levels = c("elapsed time (s)", "Δ (-2*logLik)"))
long2$method <- factor(long2$method, levels = method_levels, labels = method_labels[method_levels])

## ggh4x::facetted_pos_scales() lets a facet_grid row get its own scale
## (here: log10 for elapsed time, linear for Δ(-2*logLik)) while keeping
## everything in one faceted plot -- avoids splitting into two separate
## patchwork panels. With scales="free_y", facet_grid creates one y-scale
## *group per row* (shared across columns), not one per panel, so the
## list here has 2 elements (one per metric row), not 6.
p_time_negll <- ggplot(long2, aes(x = method, y = value, color = method, fill = method)) +
  geom_violin(alpha = 0.25, color = NA, width = 0.9) +
  geom_boxplot(width = 0.15, alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
  geom_point(alpha = 0.2, size = 1.4, show.legend = FALSE) +
  facet_grid(metric ~ dataset, scales = "free_y") +
  facetted_pos_scales(y = list(scale_y_log10(), scale_y_continuous())) +
  scale_color_manual(values = okabe_ito6) +
  scale_fill_manual(values = okabe_ito6) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "none")

outfile2 <- file.path(wd, "summary_plot_time_negll.png")
ggsave(outfile2, p_time_negll, width = 11, height = 6, dpi = 130)
cat("saved to", outfile2, "\n")
