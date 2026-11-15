## Usage: Rscript combine_dispsim.R
## Combine all dispsim_*.rds results into dispsim_combined_summary.csv (long
## format: setup x method/estimator x parameter) and print relative bias and
## relative RMSE tables for phi and the RE sd.

files <- list.files(pattern = "^dispsim_g.*\\.rds$")
alt <- c("phi_digamma", "phi_digamma_q", "phi_pearson_p", "phi_pearson_q")

one <- function(f) {
    res <- readRDS(f)
    truth <- attr(res, "truth"); setup <- attr(res, "setup")
    res <- res[res$error == "", ]
    q <- setup$ngrp + 1       # rank of [X Z] for y ~ x + (1|g)
    n <- setup$ngrp * setup$nper
    stat <- function(est, p, method) {
        e <- est - truth[[p]]
        data.frame(ngrp = setup$ngrp, nper = setup$nper, phi_true = 1 / setup$shape,
                   q_over_n = q / n, param = p, method = method,
                   relbias = mean(e) / truth[[p]],
                   se_relbias = sd(e) / sqrt(length(e)) / truth[[p]],
                   relrmse = sqrt(mean(e^2)) / truth[[p]], n = length(e))
    }
    out <- list()
    for (d in split(res, res$method)) {
        for (p in names(truth)) out[[length(out) + 1]] <- stat(d[[p]], p, d$method[1])
        if (grepl("^glmer", d$method[1]))
            for (a in intersect(alt, names(d)))
                out[[length(out) + 1]] <- stat(d[[a]], "phi",
                                               paste0(d$method[1], ":", sub("phi_", "", a)))
    }
    do.call(rbind, out)
}
all <- do.call(rbind, lapply(files, one))
all <- all[order(all$nper, all$phi_true, all$param, all$method), ]
write.csv(all, "dispsim_combined_summary.csv", row.names = FALSE)

wide <- function(p, what) {
    d <- all[all$param == p, ]
    d$cell <- sprintf("phi=%g,nper=%d", d$phi_true, d$nper)
    tab <- tapply(d[[what]], list(d$method, d$cell), function(x) round(x, 3))
    tab[, order(as.numeric(sub(".*nper=", "", colnames(tab))),
                as.numeric(sub("phi=([^,]*),.*", "\\1", colnames(tab))))]
}
cat("phi: relative bias\n"); print(wide("phi", "relbias"))
cat("\nphi: relative RMSE\n"); print(wide("phi", "relrmse"))
cat("\nRE sd: relative bias\n"); print(wide("sd_re", "relbias"))
cat("\nRE sd: relative RMSE\n"); print(wide("sd_re", "relrmse"))
cat("\nbeta0: bias (relative to 1)\n"); print(wide("beta0", "relbias"))
