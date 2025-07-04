    # Dummy functions to allow formula parsing of special covariance structures.
    # These functions do not perform any calculations; they simply act as
    # tags for the `reformulas::splitForm` function to identify them as "special"
    # terms. They accept a formula argument and return it immediately.

    us <- function(x) {
      x
    }

    cs <- function(x) {
      x
    }

    ar1 <- function(x) {
      x
    }

    # We define our own `diag` to prevent a name collision with the base R
    # `diag()` function during formula evaluation. Our version will be found
    # first and will correctly handle the formula syntax.
    diag <- function(x) {
      x
    }
