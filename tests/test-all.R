if(require("testthat", quietly = TRUE)) {
    pkg   <- "lme4"
    require(pkg, character.only=TRUE, quietly=TRUE)
    test_package(pkg)
} else {
    print( "package 'testthat' not available, cannot run unit tests" )
}
