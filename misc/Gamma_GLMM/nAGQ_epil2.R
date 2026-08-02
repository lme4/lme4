library(lme4)
library(glmmTMB)
suppressMessages(library(dplyr))
data(epil2, package="glmmTMB")
epil2_red <- filter(epil2, y>0)

fam <- Gamma(link="log")
form <- y ~ trt + (1|subject)

nAGQvals <- c(1,5,15,25)
res <- lapply(nAGQvals, function(n) {
  t0 <- Sys.time()
  m <- glmer(form, data=epil2_red, family=fam, nAGQ=n)
  list(nAGQ=n, time=as.numeric(Sys.time()-t0),
       sd_RE=as.numeric(attr(VarCorr(m)$subject,"stddev")),
       fixef=fixef(m), sigma=sigma(m), logLik=as.numeric(logLik(m)))
})

for (r in res) {
  cat(sprintf("nAGQ=%2d  sd(RE)=%.5f  sigma(disp)=%.5f  (Intercept)=%.5f  trt=%.5f  logLik=%.4f  time=%.1fs\n",
              r$nAGQ, r$sd_RE, r$sigma, r$fixef[1], r$fixef[2], r$logLik, r$time))
}

m_glmmTMB <- glmmTMB(form, data=epil2_red, family=fam)
cat(sprintf("glmmTMB  sd(RE)=%.5f  sigma(disp)=%.5f  (Intercept)=%.5f  trt=%.5f  logLik=%.4f\n",
            as.numeric(attr(VarCorr(m_glmmTMB)$cond$subject,"stddev")),
            sigma(m_glmmTMB), fixef(m_glmmTMB)$cond[1], fixef(m_glmmTMB)$cond[2], as.numeric(logLik(m_glmmTMB))))
