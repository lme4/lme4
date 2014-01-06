if(require("testthat", quietly = TRUE)) {
    pkg   <- "lme4"
    require(pkg, character.only=TRUE, quietly=TRUE)
    test_package(pkg)
    print(warnings()) # TODO? catch most of these by expect_warning(..)
} else {
    cat( "package 'testthat' not available, cannot run unit tests\n" )
}
