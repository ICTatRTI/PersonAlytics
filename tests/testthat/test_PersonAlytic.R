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
                     interactions = list(c("Time", "Phase"),
                                         c("Phase", "I(Time^2)"),
                                         c("Phase", "I(Time^3)")),
                     autoSelect = list(TO=list(polyMax=3))   )

  testthat::expect_equal(t0$tTable[,1],
               c( `(Intercept)`     =   10.487701 ,
                  Time              =  -13.174886 ,
                  `I(Time^2)`       =   31.066515 ,
                  `I(Time^3)`       =    7.217603 ,
                  Phase             =  -91.115799 ,
                  `Time:Phase`      =  354.693314 ,
                  `I(Time^2):Phase` = -416.363545 ,
                  `I(Time^3):Phase` =  128.882840
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
                     interactions    = list(c("Time", "Phase")),
                     individual_mods = TRUE        )


  # verify a single model run against batch model run in t1 (`individual_mods = TRUE`)
  mare1 <- PersonAlytic(output   = NULL             ,
                        data     = OvaryICT         ,
                        ids      = "Mare"           ,
                        dvs      = "follicles"      ,
                        phase    = "Phase"          ,
                        time     = "Time"           ,
                        package  = "arma"           ,
                        interactions    = list(c("Time", "Phase")),
                        subgroup = OvaryICT$Mare==1
  )

  l <- which(names(t1)=="ma1.Estimate")
  u <- which(names(t1)=="Time.Phase.Pr...z..")
  testthat::expect_equal( c(t(mare1$tTable[1:5,])),
             unname(unlist(t1[t1$Mare==1,l:u])) )



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

  # target ivs + psuite
  t2 <- PersonAlytic(output          = "targets"   ,
                     data            = OvaryICT    ,
                     ids             = "Mare"      ,
                     dvs             = "follicles" ,
                     phase           = "Phase"     ,
                     time            = "Time"      ,
                     package         = "arma"      ,
                     target_ivs      = paste("Target", 1:2, sep=""),
                     nbest           = 2
                     )
  # no real test, but it needs to run to check directory changes

  # just clean up
  file.remove( dir(getwd(), glob2rx("*.txt")) )
  file.remove( dir(getwd(), glob2rx("*.csv"), recursive = TRUE) )
  dirs <- list.dirs(getwd())
  dirs <- dirs[grepl("PersonAlytics output for", dirs)]
  unlink( dirs, recursive = TRUE, force = TRUE)

})
