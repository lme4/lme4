library("testthat")
library("lme4")

context("glmer.nb")
test_that("basic", {
   set.seed(101)
   dd <- expand.grid(f1 = factor(1:3),
                     f2 = LETTERS[1:2], g=1:9, rep=1:15,
                     KEEP.OUT.ATTRS=FALSE)
   mu <- 5*(-4 + with(dd, as.integer(f1) + 4*as.numeric(f2)))
   dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)
   require("MASS")
   m.glm <- glm.nb(y ~ f1*f2, data=dd)
   m.nb <- glmer.nb(y ~ f1*f2 + (1|g), data=dd)
   expect_is(m.nb,"glmerMod")
   ## 'family' properly quoted/not expanded in call?
   expect_equal(gsub("[0-9]+\\.[0-9]+","NUM",deparse(m.nb@call$family)),
          "negative.binomial(theta = NUM)")
   expect_null(m.nb@call$verbose)  ## check: GH #321
   expect_equal(fixef(m.nb), coef (m.glm), tol=1e-5)    ## GH #319

   ## GH #285
   expect_error(glmer(Reaction > 250 ~ Days + (1|Subject),
                         data = sleepstudy, family=poisson),
                "must be numeric")

   m.nb2 <- glmer.nb(y ~ f1*f2 + (1|g), data=dd,
                     subset = g!=8)
   ## expect parameters, ngrps *not* to equal full model
   expect_equal(unname(fixef(m.nb2)),
                c(1.629240234, 0.76028840, 1.008629913, 1.6172507, 
                  -0.6814426, -0.66468330),tol=1e-5)
   expect_equal(unname(ngrps(m.nb2)),8)

   ## control handling ... this should suppress warnings ...
   old.opts <- options(warning=2)
   m.nb2 <- glmer.nb(round(Reaction) ~ Days + (1|Subject),
                     data = sleepstudy, subset = Subject != 370,
                     control=glmerControl(check.conv.grad="ignore"))
   expect_is(m.nb2,"glmerMod")
   options(old.opts)

   set.seed(101)
   dd <- expand.grid(f1 = factor(1:3),
                  f2 = LETTERS[1:2], g=1:9, rep=1:15,
                     KEEP.OUT.ATTRS=FALSE)
   dd$y <- rnbinom(nrow(dd),mu=3,size=1)
   mu <- 5*(-4 + with(dd, as.integer(f1) + 4*as.numeric(f2)))
   m.nb2 <- glmer.nb(y~f1+(1|g),
                       data=dd,
                       ## offset is broken:
                       ## offset=rep(0,nrow(dd)),
                       ## offset=rep(log(2),nrow(dd)),
                       ##  -> Downdated VtV is not positive definite
                       contrasts=list(f1=contr.sum))

   m.nb2 <- glmer.nb(y~f1+(1|g), dd)
   expect_equal(names(m.nb2@call),c("","formula","data","family"))

   if (FALSE) {
       m.nb2 <- glmer.nb(y~f1+(1|g),
                       data=dd,
                         offset=rep(0,nrow(dd)))

}
}
)
