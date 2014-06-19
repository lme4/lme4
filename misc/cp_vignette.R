paper.dir <- "~/R/lme4_misc/papers/lmer_JSS"
## assume we start in misc
vignette.dir <- "../vignettes"
setwd(vignette.dir)
system(sprintf("sed -e \"s/\\\\documentclass{jss}/\\documentclass[nojss]{jss}/\" -e \"s/%%VIGNETTE_NOTE/Submitted to \\\\\\emph{Journal of Statistical Software}/\" %s/lmer.Rnw >%s/lmer.Rnw",paper.dir,"."))
knitr::knit2pdf("lmer.Rnw")
tools::compactPDF("lmer.pdf",gs_quality="ebook")
file.copy("lmer.pdf","../inst/doc",overwrite=TRUE)
unlink("lmer.pdf")
