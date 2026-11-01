## Usage: Rscript plot_dispsim.R
## Relative bias and relative RMSE (+/- 2 SE) of beta0, beta1, the RE sd and
## phi, by method and scenario (true phi x cluster size), from all
## dispsim_g*.rds files. Writes dispsim_bias.png and dispsim_rmse.png.
## Reporting-only phi estimators are shown for the glmer_n fit only (they
## hardly depend on the glmer arm). GLMMadaptive_AGQ1 is left out: its
## biases (beta0 about -0.5) would set the axis scales.

library(ggplot2)
theme_set(theme_bw())
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

fit_methods <- c("glmer_n" = "glmer dev/n", "glmer_n-p" = "glmer dev/(n-p)",
                 "glmer_n-q" = "glmer dev/(n-q)", "glmmTMB_ML" = "glmmTMB ML",
                 "glmmTMB_REML" = "glmmTMB REML",
                 "GLMMadaptive_AGQ11" = "GLMMadaptive AGQ(11)")
alt_methods <- c(phi_digamma = "glmer + digamma", phi_digamma_q = "glmer + digamma, n-q",
                 phi_pearson_p = "glmer + Pearson, n-p",
                 phi_pearson_q = "glmer + Pearson, n-q")
params <- c(beta0 = "beta0 (intercept)", beta1 = "beta1 (slope)",
            sd_re = "RE sd", phi = "phi")

## relative bias and RMSE with SEs; SE(RMSE) by the delta method,
## SE(MSE) / (2 RMSE)
stat <- function(est, true) {
    e <- (est - true) / true
    n <- length(e)
    rmse <- sqrt(mean(e^2))
    data.frame(bias = mean(e), se_bias = sd(e) / sqrt(n),
               rmse = rmse, se_rmse = sd(e^2) / sqrt(n) / (2 * rmse), n = n)
}

one <- function(f) {
    res <- readRDS(f)
    truth <- attr(res, "truth"); setup <- attr(res, "setup")
    res <- res[res$error == "" & res$method %in% names(fit_methods), ]
    out <- list()
    for (d in split(res, res$method)) {
        for (p in names(params))
            out[[length(out) + 1]] <- cbind(param = p, method = fit_methods[[d$method[1]]],
                                            type = "fit", stat(d[[p]], truth[[p]]))
    }
    gl <- res[res$method == "glmer_n", ]
    for (a in intersect(names(alt_methods), names(gl)))
        out[[length(out) + 1]] <- cbind(param = "phi", method = alt_methods[[a]],
                                        type = "reporting only",
                                        stat(gl[[a]], truth[["phi"]]))
    cbind(phi_true = 1 / setup$shape, nper = setup$nper, do.call(rbind, out))
}

files <- list.files(pattern = "^dispsim_g.*\\.rds$")
all <- do.call(rbind, lapply(files, one))
all <- within(all, {
    method <- factor(method, levels = rev(c(fit_methods, alt_methods)))
    param <- factor(params[param], levels = params)
    phi_lab <- factor(sprintf("true phi = %g", phi_true))
    cluster_size <- factor(nper)
})

pd <- position_dodge(width = 0.6, orientation = "y")
plot_metric <- function(what, xlab) {
    se <- all[[paste0("se_", what)]]
    ggplot(all, aes(x = .data[[what]], y = method, colour = cluster_size, shape = type)) +
        geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
        geom_pointrange(aes(xmin = .data[[what]] - 2 * se, xmax = .data[[what]] + 2 * se),
                        position = pd, size = 0.25) +
        facet_grid(phi_lab ~ param, scales = "free_x") +
        scale_x_continuous(n.breaks = 4) +
        scale_colour_manual(values = okabe_ito, name = "cluster size") +
        scale_shape_manual(values = c(fit = 16, "reporting only" = 1), name = NULL) +
        labs(x = xlab, y = NULL) +
        theme(panel.spacing = grid::unit(14, "pt"), legend.position = "bottom")
}

ggsave("dispsim_bias.png",
       plot_metric("bias", "relative bias, (estimate - true)/true (+/- 2 SE)"),
       width = 10, height = 7, dpi = 150)
ggsave("dispsim_rmse.png",
       plot_metric("rmse", "relative RMSE (+/- 2 SE)"),
       width = 10, height = 7, dpi = 150)
