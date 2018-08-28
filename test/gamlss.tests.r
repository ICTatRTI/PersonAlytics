if(1==2)
{
  library(PersonAlytics)

OvaryICT$follicles01 <- to01(OvaryICT$follicles)
hist(OvaryICT$follicles01)

OvaryICT$folliclesbinary <- as.numeric(OvaryICT$follicles +
                                         rnorm(nrow(OvaryICT), 0, sd(OvaryICT$follicles)) > median(OvaryICT$follicles, na.rm=TRUE))
table(OvaryICT$folliclesbinary)

fbi <- PersonAlytic("LogisticTest", data = OvaryICT, ids = "Mare", dvs = "folliclesbinary",
             time = "Time", phase = "Phase", family = BI(), package = "gamlss"  )
summary(fbi)

# validate againts lme4, non-intercept results are close
library(lme4)
fbi4 <- glmer(folliclesbinary ~ Time*Phase + (Time|Mare), data = OvaryICT,
              family = binomial)
summary(fbi4)


fbb <- PersonAlytic("LogisticTest", data = OvaryICT, ids = "Mare", dvs = "follicles01",
                    time = "Time", phase = "Phase", family = BEINF(), package = "gamlss"  )
summary(fbb)

}
