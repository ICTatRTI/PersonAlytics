library(PersonAlytics)

head(OvaryICT)
hist(OvaryICT$follicles)

# if autoSelect=list(DIST=list()), best model is WEI
t0.wei <- PersonAlytic(output  = 'Test0'     ,
                       data    = OvaryICT    ,
                       ids     = "Mare"      ,
                       dvs     = "follicles" ,
                       phase   = "Phase"     ,
                       time    = "Time"      ,
                       package = "gamlss"     ,
                       family  = WEI(),
                       autoSelect  = list())
summary(t0.wei)

# but...NBI does fit better, check autoselect and fitdist options
t0.nbi <- PersonAlytic(output  = 'Test0'     ,
                       data    = OvaryICT    ,
                       ids     = "Mare"      ,
                       dvs     = "follicles" ,
                       phase   = "Phase"     ,
                       time    = "Time"      ,
                       package = "gamlss"     ,
                       family  = NBI(),
                       autoSelect  = list())
summary(t0.nbi)

# compare normal - even this fits better than weibull, needs qc
t0.no <- PersonAlytic(output  = 'Test0'     ,
                       data    = OvaryICT    ,
                       ids     = "Mare"      ,
                       dvs     = "follicles" ,
                       phase   = "Phase"     ,
                       time    = "Time"      ,
                       package = "gamlss"     ,
                       family  = NO(),
                       autoSelect  = list())
summary(t0.no)

# this fails, still WEI best,
t0.nbi <- PersonAlytic(output  = 'Test0'     ,
                       data    = OvaryICT    ,
                       ids     = "Mare"      ,
                       dvs     = "follicles" ,
                       phase   = "Phase"     ,
                       time    = "Time"      ,
                       package = "gamlss"   ,
                       autoSelect  = list(DIST=list(type="counts")))
summary(t0.nbi)
