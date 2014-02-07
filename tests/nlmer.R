library(lme4)

allEQ <- function(x,y, tolerance = 4e-4, ...)
    all.equal.numeric(x,y, tolerance=tolerance, ...)

(nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ (Asym|Tree),
              Orange, start = c(Asym = 200, xmid = 725, scal = 350)))
fixef(nm1)

if (lme4:::testLevel() > 2) {
    ## 'Theoph' Data modeling
    Th.start <- c(lKe = -2.5, lKa = 0.5, lCl = -3)

system.time(nm2 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
                         (lKe+lKa+lCl|Subject),
                         Theoph, start = Th.start, tolPwrss=1e-8))
print(nm2, corr=FALSE)

system.time(nm3 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
                         (lKe|Subject) + (lKa|Subject) + (lCl|Subject),
                         Theoph, start = Th.start))
print(nm3, corr=FALSE)

## dropping   lKe  from random effects:
system.time(nm4 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~ (lKa+lCl|Subject),
                         Theoph, start = Th.start, tolPwrss=1e-8))
print(nm4, corr=FALSE)

system.time(nm5 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
                         (lKa|Subject) + (lCl|Subject),
                         Theoph,
                         start = Th.start, tolPwrss=1e-8))
print(nm5, corr=FALSE)

    if (require("PKPDmodels")) {
        oral1cptSdlkalVlCl <-
            PKmod("oral", "sd", list(ka ~ exp(lka), Cl ~ exp(lCl), V ~ exp(lV)))
        if (FALSE) {
            ## FIXME: Error in get(nm, envir = nlenv) : object 'k' not found
            ##  probably with environments/call stack etc.?
            ## 'pnames' is c("lV","lka","k") -- not ("lV","lka","lCl")
            ##  nlmer -> nlformula -> MkRespMod
            ##  pnames are OK in nlformula, but in MkRespMod we try to recover
            ##   them from the column names of the gradient attribute of the
            ##   model evaluated in nlenv -- which are wrong.
            system.time(nm2a <- nlmer(conc ~ oral1cptSdlkalVlCl(Dose, Time, lV, lka, lCl) ~
                                      (lV+lka+lCl|Subject),
                                      Theoph, start = c(lV=-1, lka=-0.5, lCl=-3), tolPwrss=1e-8))
            print(nm2a, corr=FALSE)
        }
    }
}  ## testLevel > 2
