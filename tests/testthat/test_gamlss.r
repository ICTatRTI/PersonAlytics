if(1==2)
{
  library(PersonAlytics)

  # set up the palytic object, turning off autoDetection for now
  t1 <- Palytic$new(data = OvaryICT, ids='Mare', dv='follicles',
                    time='Time', phase='Phase', autoDetect = list(Done=TRUE))
  t1$plot()

  # get the best distribution via palytic
  dev.new()
  t1$dist()

  # try gamlss.demo - the density plot doesn't match the theoretical density,
  # but that is limited to the count distributions
  dev.new()
  gamlss.demo::demoDist()


  # compare normal, Double Poisson, and Weibull2
  t1.NO   <- t1$gamlss(family = 'NO') # this is the default
  t1.WEI2 <- t1$gamlss(family = 'WEI2')
  t1.DPO  <- t1$gamlss(family = 'DPO')
  t1.NO$tTable
  t1.WEI2$tTable
  t1.DPO$tTable

  # Demo fully automatic DIST, TO, AR
  t2 <- Palytic$new(data = OvaryICT, ids='Mare', dv='follicles',
                    time='Time', phase='Phase',
                    autoDetect=list(AR=list(P=3, Q=3)     ,
                                    TO=list(polyMax=3)    ,
                                    DIST=list()))
  m2 <- t2$gamlss() # this will take some time to run
  m2$tTable

}









if(1==2)
{
  library(PersonAlytics)

  OvaryICT$follicles01 <- to01(OvaryICT$follicles)
  hist(OvaryICT$follicles01)

  set.seed(4269)
  OvaryICT$folliclesbinary <- as.numeric(OvaryICT$follicles +
                                           rnorm(nrow(OvaryICT), 0, sd(OvaryICT$follicles)) > median(OvaryICT$follicles, na.rm=TRUE))
  table(OvaryICT$folliclesbinary)

  start <- Sys.time()
  fbi <- PersonAlytic("LogisticTest", data = OvaryICT, ids = "Mare", dvs = "folliclesbinary",
                      time = "Time", phase = "Phase", family = BI(), package = "gamlss",
                      standardize = FALSE, alignPhase = 'none', sigma.formula = ~0,
                      detectAR = FALSE, detectTO = FALSE)
  fbi$PalyticSummary$formula
  fbi.sum <- summary(fbi)
  Sys.time() - start

  # validate againts lme4
  start <- Sys.time()
  library(lme4)
  fbi4 <- glmer(folliclesbinary ~ Time*Phase + (Time|Mare), data = OvaryICT,
                family = binomial)
  fbi4@call
  summary(fbi4)
  Sys.time() - start

  # results are close(ish) when sigma.formual =~0; gamlss has lower SE values
  fbi.sum
  summary(fbi4)$coefficients


  fbb <- PersonAlytic("LogisticTest", data = OvaryICT, ids = "Mare", dvs = "follicles01",
                      time = "Time", phase = "Phase", family = BEINF(), package = "gamlss"  )
  summary(fbb)
}




