% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{decompFormula}
\alias{decompFormula}
\title{decompFormula - function to deconstruct a return formula}
\usage{
decompFormula(formula = NULL)
}
\arguments{
\item{See}{\code{\link{Palytic}}}
}
\description{
decompFormula - function to deconstruct a return formula
}
\examples{

formula1 <- formula(y~x+z+a*b+re(random=time + I(time^2)|id,
            correlation=corARMA(p=1,q=1), method="ML") + k + l)

formula2 <- follicles ~ Time * Phase + re(random = ~1 | Mare, method = "REML",
            correlation = corARMA(p = 1, q = 0))

formula3 <- PASAT ~ (re(random = ~1 | ID, method = "REML",
            correlation = nlme::corARMA(p = 2, q = 2))) +
            Time2 + Tx + Time2:Tx

formula4 <- PASAT ~ (re(random = ~1 | ID, method = "REML", correlation = NULL)) +
            Time2 + Tx + X19.26_888.2232n + Time2:Tx

formula5 <- PASAT ~ Time2 + Tx + (re(random = ~1 | ID, method = "REML",
            correlation = nlme::corARMA(p = 1, q = 1))) +
             X19.26_888.2232n + Time2:Tx

formula6 <- DSST ~ Time2 + Tx + I(Time2^2) + I(Time2^3) +
            (re(random = ~Time2 + I(Time2^2) + I(Time2^3) | ID,
            method = "REML", correlation = nlme::corARMA(p = 3, q = 3))) +
            X3.81_442.2293m.z + Batch + Session2 + Time2:Tx + Tx:I(Time2^2) +
            Tx:I(Time2^3)

formula7 <- 1 ~ Time * Phase + re(random = ~1 | Mare, method = "REML",
            correlation = NULL, family = )

 PersonAlytics:::decompFormula(formula1)
 PersonAlytics:::decompFormula(formula2)
 PersonAlytics:::decompFormula(formula3)
 PersonAlytics:::decompFormula(formula4)
 PersonAlytics:::decompFormula(formula5)
 PersonAlytics:::decompFormula(formula6)
 PersonAlytics:::decompFormula(formula7)
}
\author{
Stephen Tueller \email{stueller@rti.org}
}
\keyword{internal}
