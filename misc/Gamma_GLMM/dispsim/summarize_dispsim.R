## Usage: Rscript summarize_dispsim.R dispsim_<setup>.rds
## Bias, RMSE and failure counts per method and parameter.

args <- commandArgs(trailingOnly = TRUE)
res <- readRDS(args[1])
truth <- attr(res, "truth")
ok <- res$error == ""
cat(sprintf("%s\n", args[1]))
print(attr(res, "setup")[c("nsim", "ngrp", "nper", "shape", "sd_re", "wall")])
cat("\nerrors / warnings / singular / mean time (s):\n")
print(do.call(rbind, lapply(split(res, res$method), function(d)
    data.frame(errors = sum(d$error != ""), warnings = sum(d$warn != ""),
               singular = sum(d$singular, na.rm = TRUE),
               time = round(mean(d$time), 3)))))
stats <- do.call(rbind, lapply(names(truth), function(p) {
    do.call(rbind, lapply(split(res[ok, ], res$method[ok]), function(d) {
        e <- d[[p]] - truth[[p]]
        data.frame(param = p, method = d$method[1], true = truth[[p]],
                   mean = mean(d[[p]]), bias = mean(e),
                   se_bias = sd(e) / sqrt(length(e)),
                   rmse = sqrt(mean(e^2)), n = length(e))
    }))
}))
## reporting-only phi alternatives, computed post hoc from each glmer fit
alt <- intersect(c("phi_digamma", "phi_digamma_q", "phi_pearson_p", "phi_pearson_q"),
                names(res))
gl <- res[ok & grepl("^glmer", res$method), ]
stats <- rbind(stats, do.call(rbind, lapply(alt, function(a) {
    do.call(rbind, lapply(split(gl, gl$method), function(d) {
        e <- d[[a]] - truth[["phi"]]
        data.frame(param = "phi", method = paste0(d$method[1], ":", sub("phi_", "", a)),
                   true = truth[["phi"]], mean = mean(d[[a]]), bias = mean(e),
                   se_bias = sd(e) / sqrt(length(e)),
                   rmse = sqrt(mean(e^2)), n = length(e))
    }))
})))
rownames(stats) <- NULL
cat("\n")
print(stats, digits = 3)
