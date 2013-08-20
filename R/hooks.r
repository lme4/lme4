.onUnload <- function(libpath) {
    gc()
    if (is.loaded("lmer_Deviance",PACKAGE="lme4"))
        library.dynam.unload("lme4", libpath)
}
