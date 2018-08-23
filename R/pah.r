# easy run defaults for testing
pah <- function()
{
  data=OvaryICT
  ids="Mare"
  dvs=list("follicles")
  phase="Phase"
  time="TimeSin"
  ivs=NULL
  target_ivs=NULL
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
  t1=Palytic$new(data=OvaryICT,
                ids="Mare",
                dv=dvs[[1]],
                time="Time",
                phase="Phase",
                ivs=NULL,
                interactions=NULL,
                time_power=1,
                correlation=NULL,
                family=gamlss.dist::BEINF()
  )
}

pahdefault <- function()
{
  file=NULL
  data=OvaryICT
  phase=NULL
  ivs=NULL
  target_ivs=NULL
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
  charSub=NULL
  sigma.formula=~1
  debugforeach = FALSE
  p.method = "BY"
  alpha = .05
}
