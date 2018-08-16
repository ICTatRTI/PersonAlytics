# easy run defaults for testing
pah <- function()
{
  data=PersonAlyticsLite::OvaryICT
  ids="Mare"
  dvs=list("follicles")
  phase="Phase"
  time="TimeSin"
  ivs=NULL
  ivsl=NULL
  interactions=NULL
  time_power=1
  correlation=NULL
  family=gamlss.dist::NO()
  subgroup=NULL
  standardize=TRUE
  package='gamlss'
  ind.mods=TRUE
  grp.mod=FALSE
  PalyticObj=NULL
  detectAR=TRUE
  detectTO=TRUE
  maxOrder=3
  dv=NULL
  fixed=NULL
  random=NULL
  formula=NULL
  charSub=NULL
  t1=PalyticHTP$new(data=PersonAlyticsLite::OvaryICT,
                ids="Mare",
                dv=dvs[[1]],
                time="TimeSin",
                phase="Phase",
                ivs=NULL,
                interactions=NULL,
                time_power=1,
                correlation=NULL,
                family=gamlss.dist::NO()
  )
}

pahdefault <- function()
{
  file=NULL
  data=PersonAlyticsLite::OvaryICT
  ids="Mare"
  dvs="follicles"
  time="TimeSin"
  phase=NULL
  ivs=NULL
  ivsl=NULL
  interactions=NULL
  time_power=1
  correlation=NULL
  family=gamlss.dist::NO()
  subgroup=NULL
  standardize=TRUE
  package='nlme'
  ind.mods=TRUE
  grp.mod=FALSE
  PalyticObj=NULL
  detectAR=TRUE
  detectTO=TRUE
  maxOrder=3
  charSub=NULL
  sigma.formula=~1
  debugforeach = FALSE
  p.method = "BY"
  alpha = .05
}