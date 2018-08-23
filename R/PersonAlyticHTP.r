#' High through-put options for Palytic objects.
#'
#' @export
#' @import snow
#' @import doSNOW
#' @import foreach
#' @import plyr
#'
#' @description A user interface for creating a Palytic object
#' (which inherits from \code{\link{Palytic}} and invoking
#' high through-put options including automated autocorrelation detection, automated
#' detection of a polynomial for time, and high throughput analyses of all individuals
#' in the data set, multiple dependent variables, and target independent variables.
#' Type I error corrections are implemented via [pending].
#'
#' @param file The file name (or full path with `/` instead of `\`) where output should
#' be saved. If left \code{NULL}, the date and time will prefix `PersonAlyticHTP_Output.csv`.
#' @param data See \code{\link{PersonAlytic}}.
#' @param ids See \code{\link{PersonAlytic}}.
#' @param dvs A list of one or more character dependent variable names in \code{data}.
#' The linear mixed effects model \code{dvs[d] ~ phase + time + phase*time + target_ivs[c] + ivs}
#' with random effects \code{~ time | ids[i]} will be fit to the data using
#' \code{\link{gamlss}}. The iterators \code{[d]}, \code{[d]}, and \code{[d]}
#' indicates the model will be fit for each combination of dependent variable in \code{dvs},
#' independent variable in \code{target_ivs}, and each unique ID in \code{ids}
#' (overridden using \code{ind.mods=FALSE}) controlling for indepented variables in
#' \code{ivs}. For more options submit a \code{PalyticObj}.
#' @param time See \code{\link{PersonAlytic}}.
#' @param phase See \code{\link{PersonAlytic}}.
#' @param ivs See \code{\link{PersonAlytic}}, noting that the variables in \code{ivs}
#' cannot also be in \code{target_ivs}.
#' @param target_ivs Independent variables that are iterated over one at a time. Effects for
#' these variables are labeled as 'target predictor' in the output.
#' @param interactions See \code{\link{PersonAlytic}}.
#' @param time_power See \code{\link{PersonAlytic}}. Ignored if \code{detectAR=TRUE}
#' @param correlation See \code{\link{PersonAlytic}}. Ignored if \code{detectTO=TRUE}
#' @param family See \code{\link{PersonAlytic}}.
#' @param subgroup See \code{\link{PersonAlytic}}.
#' @param standardize See \code{\link{PersonAlytic}}. The default is \code{TRUE}
#' which makes parameter estimate magnitudes comparable across individuals, outcomes in
#' \code{dvs}, and covariates in \code{target_ivs}.
#' @param package See \code{\link{PersonAlytic}}.
#' @param ind.mods Logical, defaults to \code{TRUE}. Should individual models be
#' fit for each ID?
#' @param PalyticObj See \code{\link{Palytic}}. If \code{PalyticObj} is submitted
#' then only \code{dvs}, \code{target_ivs}, and \code{ind.mods} will be
#' used. This allows users access to additional options including generalized linear
#' mixed effects models via the \code{family} option, user specified \code{correlation}
#' structures (in non-\code{NULL} this will override the automated correlation structure
#' search), and user specified models via \code{formula}.
#' @param detectAR Logical, defaults to \code{TRUE}. Should the autoregressive structure
#' be automatically detected? If the \code{time} variable is equally spaced, this is
#' done using the function \code{\link{forecast}}. If the \code{time} variable is not
#' equally spaced, this is done using likelihood ratio tests with mixed effects models
#' using the specified \code{package} under maximum likelihood. This is done separately
#' for each case in \code{ids} if \code{ind.mods=TRUE}.
#' @param detectTO Logical, defaults to \code{TRUE}. Should the \code{time_power} value
#' be automatically detected? Values from 1 to \code{maxOrder} will be tested using
#' likelihood ratio tests with mixed effects models using the specified \code{package}
#' under maximum likelihood. This is done separately for each case in \code{ids}
#' if \code{ind.mods=TRUE}.
#' @param  maxOrder Numeric, defaults to 3. What is the highest order of \code{time}
#' that sholud be tested for both fixed and randmo effects, e.g., \code{time +
#' I(time^2)+...+I(time^maxOrder)}.
#' @param charSub list of paired character strings for character substitution.
#' If the names of the target predictors
#' in \code{target_ivs} had to be edited to make valid variable names, this parameter allows
#' users put the illegal characters back in. For example, if the original variable name
#' was "17.00_832.2375m/z", a letter would need to prefix the variable name and the
#' "/" would need to be replaced with another character, e.g., "X17.00_832.2375m.z".
#' To get the row names of the output back to original varibale name, use
#' \code{charSub=list(c("X", ""), c("m.z", "m/z"))}. Note that inputs to charSub
#' must be in double quotes and are case sensitive. All duplicates will be substituted.
#' For example, if the variable name was "X1X23.x" and \code{charSub=list(c("X", ""))},
#' the resulting row label for this variable would be "123.x".
#' @param sigma.formula A formula for the variance under \code{\link{gamlss}}. Static.
#' It will not change dynamically over iterations nor will it be updated by \code{time_power}
#' or \code{detectTO}. If model fitting using this option fails, another attempt will be
#' made after reseting it to its defaul, i.e., \code{~1}.
#'
#' @examples
#' # group model
#' t0 <- PersonAlyticHTP(data=OvaryICT,
#'                  ids="Mare",
#'                  dvs="follicles",
#'                  phase="Phase",
#'                  time="TimeSin",
#'                  package='nlme',
#'                  ind.mods=FALSE)
#' # individual models (using defaults)
#' t1 <- PersonAlyticHTP(data=OvaryICT,
#'                  ids="Mare",
#'                  dvs="follicles",
#'                  phase="Phase",
#'                  time="TimeSin",
#'                  package='nlme')
#'
#' summary(t0)
#' summary(t1)

# \dontrun{
# # if you wish to delete the automatically created csv file, run
# #NOT IMPLEMENTED YET
# }

PersonAlyticHTP <- function(file=NULL                ,
                            data=NULL                ,
                            ids                      ,
                            dvs                      ,
                            time                     ,
                            phase=NULL               ,
                            ivs=NULL                 ,
                            target_ivs=NULL          ,
                            interactions=NULL        ,
                            time_power=1             ,
                            correlation=NULL         ,
                            family=gamlss.dist::NO() ,
                            subgroup=NULL            ,
                            standardize=TRUE         ,
                            package='nlme'           ,
                            ind.mods=TRUE            ,
                            PalyticObj=NULL          ,
                            detectAR=TRUE            ,
                            detectTO=TRUE            ,
                            maxOrder=3               ,
                            charSub=NULL             ,
                            sigma.formula=~1         ,
                            debugforeach = FALSE     ,
                            p.method = "BY"          ,
                            alpha = .05               )
{

}
