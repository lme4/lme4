paper.dir <- "~/R/D/bitbucket/lme4_misc/papers/lmer_JSS"
paper.dir <- "~/R/lme4_misc/papers/lmer_JSS"

## assume we start in lme4/misc, but lme4/vignettes works, too:
vignette.dir <- "../vignettes"
lmerSrc <- file.path(paper.dir, "lmer.Rnw")
if(!file.exists(lmerSrc)) {
    stop("paper.dir directory wrong? File ", lmerSrc," does not exist")
} else if(inherits(tryCatch(owd <- setwd(vignette.dir), error= function(e)e), "error")) {
    stop("Cannot setwd(",vignette.dir,"). Must run inside  lme4/misc/")
} else {
    system(sprintf("sed -e \"s/\\\\documentclass{jss}/\\documentclass[nojss]{jss}/\" -e \"s/%%VIGNETTE_NOTE/Submitted to \\\\\\emph{Journal of Statistical Software}/\" %s > lmer.Rnw",
                   lmerSrc))
    knitr::knit2pdf("lmer.Rnw")
    ## for arXiv:
    if(!file.exists("jss.cls"))
        file.symlink(file.path(R.home(), "share/texmf/tex/latex", "jss.cls"), ".")
    ## unfortunately, R's tar() is not usable for this ... {aargh!}
    ## tar("lmer_4_arxiv.tar",
    ##     files = c("lmer.tex", "lmer.bbl", "figure", "jss.cls"))
    system("tar cfzh lmer_4_arxiv.tar.gz lmer.tex lmer.bbl figure jss.cls figure")
    untar("lmer_4_arxiv.tar.gz", list=TRUE)
    ##
    tools::compactPDF("lmer.pdf",gs_quality="ebook")
    file.copy("lmer.pdf","../inst/doc", overwrite=TRUE)
    ls.lmer <- function() list.files(pattern="^lmer[.]...$")
    fls <- local({L <- ls.lmer(); L[!(L %in% paste0("lmer.", c("bib", "Rnw")))]})
    cat("--> tar for arXiv:  ", list.files(pattern="^lmer_4_arxiv"),"\n")
    cat("removing ",fls,": .."); unlink(fls); cat("\n")
    cat(" --> lmer.* in vignettes:\n"); print(ls.lmer())
    setwd(owd)# back
}
