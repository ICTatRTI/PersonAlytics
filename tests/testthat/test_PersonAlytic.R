context("PersonAlytic")
library(PersonAlytics)
library(nlme)


test_that("ManyV1",
{
  OvaryICT <- PersonAlytics::OvaryICT

  t0 <- PersonAlytic(output  = NULL        ,
                     data    = OvaryICT    ,
                     ids     = "Mare"      ,
                     dvs     = "follicles" ,
                     phase   = "Phase"     ,
                     time    = "Time"      ,
                     package = "nlme"      ,
                     autoSelect = list(TO=list(polyMax=3))   )

  testthat::expect_equal(t0$tTable[,1],
               c( `(Intercept)`     =   10.487701 ,
                  Time              =  -13.174889 ,
                  Phase             =  -91.115823 ,
                  `I(Time^2)`       =   31.066552 ,
                  `I(Time^3)`       =    7.217541 ,
                  `Time:Phase`      =  354.693415 ,
                  `Phase:I(Time^2)` = -416.363710 ,
                  `Phase:I(Time^3)` =  128.882955
                  )
  )

  mod.lme <- lme(t0$PalyticSummary$fixed, data = OvaryICT,
             random = t0$PalyticSummary$random)

  testthat::expect_equal(summary(mod.lme)$tTable, t0$tTable)

  # individual models
  t1 <- PersonAlytic(output          = NULL        ,
                     data            = OvaryICT    ,
                     ids             = "Mare"      ,
                     dvs             = "follicles" ,
                     phase           = "Phase"     ,
                     time            = "Time"      ,
                     package         = "arma"      ,
                     individual_mods = TRUE        )


  # verify a single model run against batch model run in t1 (`individual_mods = TRUE`)
  mare1 <- PersonAlytic(output   = NULL             ,
                        data     = OvaryICT         ,
                        ids      = "Mare"           ,
                        dvs      = "follicles"      ,
                        phase    = "Phase"          ,
                        time     = "Time"           ,
                        package  = "arma"           ,
                        subgroup = OvaryICT$Mare==1
  )

  testthat::expect_equal( c(t(mare1$tTable[1:5,])),
             unname(unlist(t1[t1$Mare==1,44:63])) )

  # verify against a 'manual' auto.arima run
  m1   <- OvaryICT$Mare==1
  xr   <- all.vars(mare1$PalyticSummary$fixed)
  xreg <- model.matrix(mare1$PalyticSummary$fixed,  OvaryICT[m1,xr])[,-1]

  m1ar <-
  forecast::auto.arima(y = OvaryICT$follicles[m1],
                       xreg = data.matrix(xreg))

  testthat::expect_equal(m1ar$coef, mare1$arima$coef)

  # verify against saved values
  testthat::expect_equal(m1ar$coef,   c(ma1          =  0.4407915 ,
                              intercept    =  6.2146980 ,
                              Time         = -5.9190285 ,
                              Phase        = -2.0250462 ,
                              `Time:Phase` = 10.7527872 ))

  # clean up
  file.remove( dir(getwd(), glob2rx("*.txt")) )
  file.remove( dir(getwd(), glob2rx("*.csv")) )

})
