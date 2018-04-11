phaseLocations <- function(ids, phases, times)
{
  id.min     <- aggregate(times, list(ids, phases), min)
  phase.info <- aggregate(id.min$x, list(id.min$Group.2), median)
  phase.info <- phase.info[order(phase.info$x),]
  phase.nms  <- phase.info$Group.1
  ### center the phase locations
  phase.mt <- c(phase.info$x, max(times, na.rm=TRUE))
  phase.mm <- cbind(phase.mt[1:(length(phase.mt)-1)], phase.mt[2:length(phase.mt)])
  phase.m  <- apply(phase.mm, 1, mean)
  ### test results
  if(length(phase.m)!=length(phase.nms)) stop(paste('in phase.loc(): there are', length(phase.m),
                                                    'phase locations but', length(phase.nms),
                                                    'phase names.'))
  if((1+length(phase.m))!=length(phase.mt)) stop(paste('in phase.loc(): there are', length(phase.mt),
                                                       'phase transition points but there should be', 1+length(phase.m),
                                                       'points.'))
  ### output
  return(list(phase_loc=phase.m, phase_tic=phase.mt, phase_nms=phase.nms))
}
