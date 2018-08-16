# this is supposedly taboo but I've been unable to work around it b/c we
# cannot :: the %dopar% operator
require(foreach)
# parralelization
funcs <- c("mean") # c(".eds") -- not importing from PersonAlytic correctly
cl    <- snow::makeCluster(parallel::detectCores(), type="SOCK")
snow::clusterExport(cl, funcs)
doSNOW::registerDoSNOW(cl)
pkgs  <- c("gamlss", "nlme")
