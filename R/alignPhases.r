### align the data at the transition between the first and second phase
alignPhases <- function(dat, id, phase, time, do.plot=FALSE)
{
  time.old <- dat[,time]

  ### only align if phase order is the same for everyone
  phase.order <- aggregate(dat[,time], by=list(dat[,id], dat[,phase]), min, na.rm=TRUE)
  phase.order <- aggregate(phase.order$x, by=list(phase.order$Group.1), order)
  phase.is.ordered <- all(phase.order$x[1,1]==phase.order$x[,1])

  u.id <- unique(dat[,id])

  ### order phases if they are unordered
  if(!phase.is.ordered)
  {
    phase.Levels <- levels(as.factor(dat[,phase]))
    ### this automatically takes the levels alphabetically except for the
    #   selected phase.Level; user input will be needed in the future
    phase.Levels <- c(phase.Level, phase.Levels[phase.Levels!=phase.Level])
    min.times <- max.times <- vector('list', 0)

    ### rescale intraphase times to start at 0
    for(i in 1:length(phase.Levels))
    {
      temp      <- dat[dat[,phase]==phase.Levels[i], c(id,time)]
      min.times[[i]] <- aggregate(temp[,time], list(temp[,id]), min, na.rm=TRUE)
      max.times[[i]] <- aggregate(temp[,time], list(temp[,id]), max, na.rm=TRUE)
      for(j in 1:length(u.id))
      {
        temp.time <- temp[temp[,id]==u.id[j],time]
        temp.time <- temp.time - min.times[[i]]$x[min.times[[i]]$Group.1==u.id[j]]
        temp[temp[,id]==u.id[j], time] <- temp.time
        rm(temp.time)
      }
      dat[dat[,phase]==phase.Levels[i], time] <- temp[,time]
    }

    ### reaccumulate time for the phase order in phase.Levels
    for(i in 2:length(phase.Levels))
    {
      ### get the interphase gap times
      max.pmi <- cbind(min.times[[i-1]]$Group.1, apply(cbind(min.times[[i-1]]$x, min.times[[i]]$x),1,max))
      min.pmi <- cbind(max.times[[i-1]]$Group.1, apply(cbind(max.times[[i-1]]$x, max.times[[i]]$x),1,min))
      dif.pmi <- data.frame(id=min.times[[i-1]]$Group.1, diff=max.pmi[,2] - min.pmi[,2])

      ### accumulate
      temp_im1  <- dat[dat[,phase]==phase.Levels[i-1], c(id,time)]
      max.times <- aggregate(temp_im1[,time], list(temp_im1[,id]), max, na.rm=TRUE)
      temp      <- dat[dat[,phase]==phase.Levels[i], c(id,time)]
      min.times <- aggregate(temp[,time], list(temp[,id]), min, na.rm=TRUE)
      for(j in 1:length(u.id))
      {
        temp.time <- temp[temp[,id]==u.id[j],]
        temp.time <- temp.time + (max.times$x[max.times$Group.1==u.id[j]] + dif.pmi$diff[dif.pmi$id==u.id[j]])
        temp[temp[,id]==u.id[j], time] <- temp.time[,time]
        rm(temp.time)
      }
      dat[dat[,phase]==phase.Levels[i], time] <- temp[,time]
    }
    ### re-sort the data
    dat <- dat[order(dat[,id], dat[,time]),]
  }

  ### alignment to the first phase transition
  if(do.plot)
  {
    par(mfrow=c(2,2))
    hist(dat[,time])
    time.old <- dat[,time]
  }
  for(i in u.id)
  {
    w       <- dat[,id]==i
    phase.i <- factor(dat[w,phase])
    if(is.numeric(phase.i)) phase.t <- table(phase.i!=min(phase.i))
    if(is.factor(phase.i))  phase.t <- table(phase.i)
    time.i  <- dat[w, time]
    time.p  <- time.i - time.i[ phase.t[1]+1 ]
    #time.p2  <- -(phase.t[1]):(phase.t[2]-1)
    dat[w,time] <- time.p
    #cbind(time.i, time.p, time.p2, phase.i)
  }
  if(do.plot)
  {
    plot(time.old, dat[,time])
    hist(dat[,time])
    par(mfrow=c(1,1))
  }


  return(dat)
}
