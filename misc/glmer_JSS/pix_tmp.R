library(lme4)
load("pix_tmp.rda")
library(ggplot2)
theme_set(theme_bw())
revguide <- guide_legend(reverse = TRUE)
library(RColorBrewer)

avg_size_2sd <- 2*sd(cbpp2$avg_size)
scfun <- function(v) {
  combCI[combCI$var == "avg_size", v] <-
    combCI[combCI$var == "avg_size", v]*avg_size_2sd
  assign("combCI", combCI, parent.frame())
}
scfun("lwr")
scfun("upr")
scfun("est")
combCI <- transform(combCI,
            vartype = factor(ifelse(grepl("^sd", var), "random effects", "fixed effects"))
)                    

ggplot(combCI, aes(var, est, colour=type, shape=model)) +
  geom_pointrange(aes(ymin=lwr, ymax=upr),
                  position=position_dodge(width=0.6)) +
  labs(x="", y="Estimate (log-odds of seropositivity)") +
  geom_hline(yintercept=0, lty=2) +
  coord_flip() +
  scale_colour_brewer(palette="Dark2", guide = revguide) +
    scale_shape(guide=revguide, labels = mnames) +
    facet_wrap(~ vartype, ncol = 1, scale = "free")
