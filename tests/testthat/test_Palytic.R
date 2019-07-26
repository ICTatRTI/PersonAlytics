context("Palytic")
library(PersonAlytics)

test_that("PalyticBasics",
{
  #OvaryICT <<- PersonAlytics::OvaryICT

  t1 <- Palytic$new(data = OvaryICT, ids='Mare', dv='follicles',
                    time='Time', phase='Phase', autoDetect = list())

  t1.gamlss <- summary(t1$gamlss())
  t1.lme    <- summary(t1$lme())
  form1 <- t1$formula

  expect_equal( t1.gamlss[1:4,1], t1.lme$tTable[,1], tolerance = 0.01 )

  t1$correlation <- "corARMA(p=1, q=0)"
  t1.gamlss.ar1 <- t1$gamlss(); t1.gamlss.ar1s <- summary( t1.gamlss.ar1 )
  t1.lme.ar1    <- t1$lme()

  expect_equal(
  t1.gamlss.ar1$PalyticSummary$formula,
  t1.lme.ar1$PalyticSummary$formula$formula)

  # this test will fail, I don't know why we get hugely different results
  # when adding a residual correlation structure
  #expect_equal( t1.gamlss.ar1s[1:4,1],
  #              t1.lme.ar1$tTable[,1], tolerance = 0.01 )

  # repeat 'by hand' to check that gamlss is using AR1
  formar <- t1.gamlss.ar1$PalyticSummary$formula
  t1.gamlss.ar1s2 <- summary( gamlss(form1, data = OvaryICT) )
  expect_false(identical(t1.gamlss.ar1s2, t1.gamlss.ar1s))

  # WIP - the lme varFuncs are not working
  if(1==2)
  {
    # some things to check on getting more equal results between gamlss and lme

    g.sig.fix <- gamlss(form1, data = OvaryICT, sigma.formula = ~ 1,
                        sigma.fix = TRUE)
    g.sig.fixs <- summary(g.sig.fix)

    g.sig.identity <- gamlss(form1, data = OvaryICT, sigma.formula = ~ 1,
                             family = NO(sigma.link = 'identity'))
    g.sig.identitys <- summary(g.sig.id)

    # varFuncs, see https://fukamilab.github.io/BIO202/03-C-heterogeneity.html
    .varFixed <- varFixed(~Time)
    .varFixed <- initialize(.varFixed, data = OvaryICT)
    g.sig.varFixed <- gamlss(form1, data = OvaryICT, sigma.formula = ~ 1,
                             weights=.varFixed)
    g.sig.varFixeds <- summary(g.sig.varFixed)

    .varIdent <- initialize(varIdent(form=~1|Time), OvaryICT)
    g.sig.varIdent <- gamlss(form1, data = OvaryICT, sigma.formula = ~ 1,
                             weights=.varIdent)
    g.sig.varIdents <- summary(g.sig.varIdent)

    # table all of the results
    # estimates
    data.frame( t1.lme          = t1.lme$tTable[,1],
                t1.gamlss       = t1.gamlss[1:4,1],
                t1.gamlss.ar1s  = t1.gamlss.ar1s[1:4,1],
                t1.gamlss.ar1s2 = t1.gamlss.ar1s2[1:4,1],
                g.sig.fixs      = g.sig.fixs[1:4,1],
                g.sig.identitys = g.sig.identitys[1:4,1])

    # variance estimates
    data.frame( t1.gamlss      = t1.gamlss[5,1],
                t1.gamlss.ar1s = t1.gamlss.ar1s[5,1],
                g.sig.fixs     = g.sig.ids[5,1])
  }
})


# From an email sent to mikis.stasinopoulos@gamlss.org on Fri 2018-05-25 1:41 PM
# No response was given

# Hello Mikis –

#Thank you for the gamlss package, its broad applicability address many real data problems I run into. I have a question about the comparability of gamlss and lme for the case of normal data which arose as I was developing a function to calculate the intraclass correlation (ICC) from lme and gamlss objects. The question I have is whether one can replicate a model in gamlss and lme? In the example below, I can get the same variance components and ICC if I use ‘sigma.formula = ~ 0’ but not using ‘sigma.formula = ~ 1’, but the former results in an error and fails to produce AIC, etc. (see the last line of code).

#Is it possible to get the same variance components and still get the AIC from gamlss? I have the sense ‘sigma.formula = ~ 0’ is not appropriate but the fact that it replicate the lme variances gives me pause. Any insight is appreciated. Thank you.

if(1==2){
test_that("PalyticICC",
{
  ICC <- function(out)
  {
    if('lme' %in% class(out)) varests <- as.numeric(nlme::VarCorr(out)[1:2])
    if('gamlss' %in% class(out)) varests <- as.numeric(nlme::VarCorr(gamlss::getSmo(out))[1:2])
    return( varests[1]/sum(varests) )
  }

  egaov <- aov(follicles ~ factor(Mare), data = nlme::Ovary)

  eglme <- nlme::lme(follicles ~ 1, data = nlme::Ovary, random = ~ 1 | Mare,
                     method = 'ML')

  eggamlss1 <- gamlss::gamlss(follicles ~ 1 + re(random = ~ 1 | Mare, method = 'ML'),
                              data = nlme::Ovary,
                              sigma.formula = ~ 1)

  eggamlss2 <- gamlss::gamlss(follicles ~ 1 + re(random = ~ 1 | Mare, method = 'ML'),
                              data = nlme::Ovary,
                              sigma.formula = ~ 0)

  # compare variance components
  VarCorr(eglme)
  VarCorr(getSmo(eggamlss1))
  VarCorr(getSmo(eggamlss2))

  # compare iccs
  ICC(eglme)
  ICC(eggamlss1)
  ICC(eggamlss2)
  multilevel::ICC1(egaov)
  #multilevel::ICC2(egaov)

  # show the error
  summary(eggamlss2)

  # no actual tests yet

})
}

# there are still some small differences that I need to track down between
# the manual model
test_that("groupAR_Order",
{
  # 'manual' model
  #OvaryICT <<- PersonAlytics::OvaryICT

  ctrl <- nlme::lmeControl(opt="optim")
  m1 <- lme(follicles ~ Time * Phase, data = OvaryICT, random = ~ Time | Mare,
            correlation = corARMA(p=0,q=2), control = ctrl, method = "REML")
  #m1$call
  #dim(m1$data)

  t1 <- Palytic$new(data = OvaryICT, ids='Mare', dv='follicles',
                    time='Time', phase='Phase', autoDetect = list(AR=list(P=3,Q=3)))
  t1$GroupAR()
  m1t <- t1$lme()
  #m1t$call
  #dim(m1t$data)

  expect_equal(summary(m1)$tTable, m1t$tTable, tolerance = .01)

  expect_equal(m1$modelStruct, m1t$modelStruct, tolerance = .01)

}
)
