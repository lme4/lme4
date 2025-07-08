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

 dcov <- function(x) {
      x
    }
 
