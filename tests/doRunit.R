#### doRUnit.R --- Run RUnit tests
####------------------------------------------------------------------------

### MM: Vastly changed:  This should also be "runnable" for *installed*
##              package which has no ./tests/
## ----> put the bulk of the code e.g. in  ../inst/unitTests/runTests.R :


if(require("RUnit", quietly = TRUE)) {
    pkg <- "lme4Eigen"

    require( pkg, character.only=TRUE)

    path <- system.file("unitTests", package = pkg)

    stopifnot(file.exists(path), file.info(path.expand(path))$isdir)

    source(file.path(path, "runTests.R"), echo = TRUE)
} else {
    print( "package RUnit not available, cannot run unit tests" )
}

