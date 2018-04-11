## ------------------------------------------------------------------------
# change data to standard data.frame, and create a phase variable
Ovary <- as.data.frame(nlme::Ovary)
Ovary$Mare <- factor(Ovary$Mare, ordered = FALSE)
Ovary$Phase <- as.numeric(Ovary$Time > .5)
Ovary$TimeSin <- sin(2*pi*Ovary$Time) # the model won't converge without this!
#par(mfrow=c(2,2))
#hist(Ovary$Time, xlim=c(-1,1.5))
#plot(Ovary$Time, Ovary$TimeSin)
#hist(Ovary$TimeSin, xlim=c(-1,1.5))
#par(mfrow=c(1,1))
# figure 5.10 page 240
#ggplot2::ggplot(Ovary, aes(x = Time, y = follicles)) + geom_line() + geom_point() +
#  facet_grid( ~ Mare)
## rescaled time
#ggplot2::ggplot(Ovary, aes(x = TimeSin, y = follicles)) + geom_line() + geom_point() +
#  facet_grid( ~ Mare)
frm <- as.formula('follicles ~ Phase + TimeSin + Phase*TimeSin')

