### nlmer() convergence testing / monitoring / ...
##  -------------------
### The output of tests here are *not* 'diff'ed  (<==> no *.Rout.save file)
library(lme4a)

## 'Theoph' Data modeling

Th.start <- c(lKe = -2.5, lKa = 0.5, lCl = -3)

(nm2 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
              0+lKe+lKa+lCl+(0+lKe+lKa+lCl|Subject),
              Theoph, start = Th.start, verb = TRUE))
(nm3 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
              0+lKe+lKa+lCl+(0+lKe|Subject)+(0+lKa|Subject)
              +(0+lCl|Subject),
              Theoph, start = Th.start, verbose = 1))
## dropping   lKe  from random effects:
(nm4 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
              0+lKe+lKa+lCl+(0+lKa+lCl|Subject),
              Theoph, start = Th.start, verbose = 1))
(nm5 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
              0+lKe+lKa+lCl+(0+lKa|Subject)+(0+lCl|Subject),
              Theoph, start = Th.start, verbose = 1))
