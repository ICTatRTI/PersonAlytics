library(PersonAlytics)

# multiple target ivs
t2 <- PersonAlytic(output='Test1'                                     ,
                   data=OvaryICT                                      ,
                   ids="Mare"                                         ,
                   dvs="follicles"                                    ,
                   phase="Phase"                                      ,
                   time="Time"                                        ,
                   target_ivs=as.list( paste('Target', 1:6, sep='') ) ,
                   package="nlme"                                     ,
                   individual_mods=TRUE                               )
