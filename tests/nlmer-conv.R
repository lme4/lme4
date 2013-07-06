### nlmer() convergence testing / monitoring / ...
##  -------------------
(testLevel <- if (nzchar(s <- Sys.getenv("LME4_TEST_LEVEL"))) as.numeric(s) else 1)

### The output of tests here are *not* 'diff'ed  (<==> no *.Rout.save file)
library(lme4)

## 'Theoph' Data modeling

if (testLevel > 1) {
    Th.start <- c(lKe=-2.5, lKa=0.5, lCl=-3)

    (nm2 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~ lKe+lKa+lCl|Subject,
                  Theoph, start = Th.start))
    (nm3 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
                  (lKe|Subject)+(lKa|Subject)+(lCl|Subject),
                  Theoph, start = Th.start))
    ## dropping   lKe  from random effects:
    (nm4 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~lKa+lCl|Subject,
                  Theoph, start = Th.start, tolPwrss=1e-8))
    (nm5 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~(lKa|Subject)+(lCl|Subject),
                  Theoph, start = Th.start))
}
