
library(PersonAlytics)

OvaryICT$follicles01 <- to01(OvaryICT$follicles)
hist(OvaryICT$follicles01)

set.seed(4269)
OvaryICT$folliclesbinary <- as.numeric(OvaryICT$follicles +
                                         rnorm(nrow(OvaryICT), 0, sd(OvaryICT$follicles)) > median(OvaryICT$follicles, na.rm=TRUE))
table(OvaryICT$folliclesbinary)

fbi <- PersonAlytic("LogisticTest", data = OvaryICT, ids = "Mare", dvs = "folliclesbinary",
             time = "Time", phase = "Phase", family = BI(), package = "gamlss",
             standardize = FALSE, alignPhase = FALSE, sigma.formula = ~0,
             detectAR = FALSE, detectTO = FALSE)
fbi$PalyticSummary$formula
fbi.sum <- summary(fbi)

# validate againts lme4
library(lme4)
fbi4 <- glmer(folliclesbinary ~ Time*Phase + (Time|Mare), data = OvaryICT,
              family = binomial)
fbi4@call
summary(fbi4)

# results are close(ish) when sigma.formual =~0; gamlss has lower SE values
fbi.sum
summary(fbi4)$coefficients


fbb <- PersonAlytic("LogisticTest", data = OvaryICT, ids = "Mare", dvs = "follicles01",
                    time = "Time", phase = "Phase", family = BEINF(), package = "gamlss"  )
summary(fbb)


