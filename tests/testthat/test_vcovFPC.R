context("vcovFPC")
library(PersonAlytics)
library(nlme)
#library(gamlss)
library(lme4)

test_that("Ovary",
          {
            mod.lme <- lme(follicles ~ Time, data = Ovary, random = ~ Time | Mare,
                           method = 'ML')
            mod.merMod <- lmer(follicles ~ Time + (Time | Mare), data = Ovary,
                               REML = FALSE)
            popsize2 <- seq(100,1000,by=100)
            fpc.se.merMod <- fpc.se.lme <- list()
            for(i in seq_along(popsize2))
            {
              fpc.se.merMod[[i]] <- sqrt(diag(vcovFPC(mod.merMod, popsize2 = popsize2[i])))
              fpc.se.lme[[i]] <- sqrt(diag(vcovFPC(mod.lme, popsize2 = popsize2[[i]])))
            }

            testthat::expect_equal(
              do.call(rbind, fpc.se.merMod),
              do.call(rbind, fpc.se.lme),
              tolerance = 9e-05
            )

            # use getLambdat to create Lambdat and show equivalence to lme4 (within rounding)
            Lambdat <- getLambdat(mod.lme)
            testthat::expect_equal(Lambdat, getME(mod.merMod, "Lambdat"), tolerance = 9e-05)

            # use getZt to create Zt and show equivalence to lme4
            Zt <- getZt(mod.lme)
            testthat::expect_equal(unname(Zt), unname(getME(mod.merMod, "Zt")))
          }

)

if(1==2)
{
# 20190726 update to lme4 causes lmer to fail to converge, troubleshoot this later

test_that("BodyWeight",
          {
            mod.lme <- lme(weight ~ Time, data = BodyWeight, random = ~ Time | Rat,
                           method = 'ML')
            mod.merMod <- lmer(weight ~ Time + (Time | Rat), data = BodyWeight,
                               REML = FALSE)
            popsize2 <- seq(100,1000,by=100)
            fpc.se.merMod <- fpc.se.lme <- list()
            for(i in seq_along(popsize2))
            {
              fpc.se.merMod[[i]] <- sqrt(diag(vcovFPC(mod.merMod, popsize2 = popsize2[i])))
              fpc.se.lme[[i]] <- sqrt(diag(vcovFPC(mod.lme, popsize2 = popsize2[[i]])))
            }
            testthat::expect_equal(
              do.call(rbind, fpc.se.merMod),
              do.call(rbind, fpc.se.lme),
              tolerance = 9e-06
            )

            # use getLambdat to create Lambdat and show equivalence to lme4 (within rounding)
            Lambdat <- getLambdat(mod.lme)
            testthat::expect_equal(Lambdat, getME(mod.merMod, "Lambdat"), tolerance = 9e-06)
            # use getZt to create Zt and show equivalence to lme4
            Zt <- getZt(mod.lme)
            testthat::expect_equal(unname(Zt), unname(getME(mod.merMod, "Zt")))
          }
)

}
