library(lme4Eigen)

allEQ <- function(x,y, tolerance = 4e-4, ...)
    all.equal.numeric(x,y, tolerance=tolerance, ...)

(nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~
              0 + Asym + xmid + scal + (0 + Asym|Tree),
              Orange, 
              start = c(Asym = 200, xmid = 725, scal = 350),
              verbose = 1L))
fixef(nm1)

## 'Theoph' Data modeling

Th.start <- c(lKe = -2.5, lKa = 0.5, lCl = -3)
system.time(nm2 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
                         0 + lKe + lKa + lCl + 
                         (0 + lKe+lKa+lCl|Subject), verb = 1,
                         Theoph, start = Th.start))  # ~ 5.7s {dual-opteron 2814, on 64b, no optim.}
fixef(nm2)

system.time(nm3 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
                          0 + lKe + lKa + lCl + (0 + lKe|Subject) +
                          (0 + lKa|Subject) + (0 + lCl|Subject),
                          Theoph, start = Th.start,
                          verbose = 1L)) # ~ 3.2s
fixef(nm3)

## dropping   lKe  from random effects:
system.time(nm4 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
                         0 + lKe + lKa + lCl + (0+lKa+lCl|Subject),
                         Theoph, start = Th.start, verbose = 1L))
fixef(nm4)
sigma(nm4)

system.time(nm5 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
                         0 +lKe + lKa + lCl + (0 + lKa|Subject) +
                         (0 + lCl|Subject), verbose = 1L,
                         Theoph, start = Th.start))
fixef(nm5)

#e3 <- expand(nm3)
#stopifnot(identical(sapply(e3, class),
#                    c(sigma = "numeric", P = "pMatrix",
#                      T = "dtCMatrix", S = "ddiMatrix"))
#          , allEQ(e3$sigma, c(sigmaML = 0.70777))
#          , all(e3$P@perm == outer(12*(0:2), 1:12, "+"))
#          , identical(as(e3$T, "diagonalMatrix"), Diagonal(3*12))
#          , allEQ(e3$S@x, rep(c(0, 0.92746, 0.23667), each=12))
#          )

## e2 <- expand(nm2) # -> gave error!
## stopifnot(identical(sapply(e2, class),
##                     c(sigma = "numeric", P = "pMatrix",
##                       T = "dtCMatrix", S = "ddiMatrix"))
## #          , allEQ(e2$sigma, c(sigmaML = 0.70777))
##           , all(e2$P@perm == outer(12*(0:2), 1:12, "+"))
##           , all(diag(e2$T == 1))
##           , nnzero(e2$T) == 36 + 24 + 12
## #          , allEQ(unique(e2$T@x),
## #                  c(1, 0.0310, 0.0909, -0.00187), tol = .009)
## #          , allEQ(e2$S@x, rep(c(0.0000000, 0.9274296, 0.2366269), each=12))
##           )

#showProc.time()


