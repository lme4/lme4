cat ../../tests/testthat/$1 | \
    sed -e 's/test_that(/## test_that(/' | \
    sed -e 's/context(/## context(/'  |  \
    sed -e 's/})/## })/'  | \
    sed -e 's/library("testthat")//' |
    sed -e 's/require("testthat")//' \
> $1


