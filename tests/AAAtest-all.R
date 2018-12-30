if(require("testthat", quietly = TRUE)) {
    pkg   <- "lme4"
    require(pkg, character.only=TRUE, quietly=TRUE)
    if(getRversion() < "3.5.0") { withAutoprint <- identity ; prt <- print } else { prt <- identity }
    if(Sys.getenv("USER") %in% c("maechler", "bbolker")) withAutoprint({
        ## for developers' sake:
        lP <- .libPaths() # ---- .libPaths() : ----
        prt(lP)
        ## ---- Entries in .libPaths()[1] : ----
        prt(list.files(lP[1], include.dirs=TRUE))
        prt(sessionInfo())
        prt(packageDescription("Matrix"))
        ## 'lme4' from packageDescription "file" :
        prt(attr(packageDescription("lme4"), "file"))
    })
    test_check(pkg)
    ##======== ^^^
    print(warnings()) # TODO? catch most of these by expect_warning(..)
} else {
    cat( "package 'testthat' not available, cannot run unit tests\n" )
}
