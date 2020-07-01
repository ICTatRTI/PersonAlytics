######################################################################################################################################
######################################################################################################################################
## This R script is trying to implement the MLM with Cholesky Transformation (MLM-CT) by Jahug & Wood (2017) that 
## was written in SAS into R code. The goal is to incorporate heterogenous level-1 error covariance structure 
## (e.g., residual variance sigma, correlation parameter rho, etc.). The method first involves fitting individual 
## time series data with appropriate structure to account for autocorrelation. Essentially, a transformation
## on each single time series regression equation is performed so that one can assume no serial correlation
## after the transformation (ID error covariance structure). Using the transformed data we can easily fit 
## a MLM in the 'lme' function with all persons having homogeneous level-1 error structure, e.g., ID error covariance matrix
## To obtain the transformation matrix A for each individual, one can use a Cholesky composition of the covariance matrix sigma.
## Although the estimation of Sigma can get complicated if a complex autocorrelation structure is imposed.

######################################################################################################################################
## One drawback of the method is the step to compose all indivudal time series models and compute sigma for the 
## transformation matricies A's. In theory the individual time series model can have different autocorrelation struture
## if necessary, however, this is a bit unrealistic to do it automatically, particularly if there are many people with 
## supposedly very different dynamic patterns (not just the value of parameters).  

######################################################################################################################################
## Another note on this is that even the simulation study in this article found that negnecting the heterogeneous error structure
## is NOT really causing problems to EITHER the parameter estimation of the fixed effect, NOR the variance or covariance of
## of random effect!! The correction with heterogeneous error structure only moderately corrected the standard error of 
## fixed intercept and SE of fixed effect of level-2 covariates WHEN the average autocorrelation is really high, and such correction 
## method would require sufficient # of observation (N>50)!! A large # of timepoints is also required to compute 
## the Sigma and A matrix for the first step above, another drawback of the MLM-CT method.  This altogether to me seemed a trivial 
## issue that might need tedious or (time-wise) expansive adjustment to do. What do you think??

######################################################################################################################################
#######################################################################################################################################

# Data simulation from the R file named 'group specific correlation structures in R.r'
# imported notes
# try getting indidual correlation structures
# https://stackoverflow.com/questions/11819720/converting-repeated-measures-mixed-model-formula-from-sas-to-r
## create reproducible example fake panel data set: <-- from stackoverflow -->

set.seed(94); subject.id = abs(round(rnorm(10)*10000,0)) # 10 subjects
set.seed(99); sds = rnorm(10,15,5);means = 1:10*runif(10,7,13);trends = runif(10,0.5,2.5)
this = NULL; set.seed(98)
for(i in 1:10) { this = c(this,rnorm(6, mean = means[i], sd = sds[i])*trends[i]*1:6)} # 6 timepoints
set.seed(97)
that = sort(rep(rnorm(10,mean = 20, sd = 3),6))
df1 = data.frame(day = rep(1:6,10), GROUP = c(rep('TEST',30),rep('CONTROL',30)), # GROUP factor, 2 values
                 Y = this, # response *60 values
                 X1 = that, # individual-level covariates *10 values
                 person = sort(rep(subject.id,6))) # *10 values
write.csv(df1, 'df1.csv', row.names=FALSE)

######################################################################################################################################
# PART I  What is capable in package 'nmle'
## Before I jumped into other methods, I want to test how 'nmle' handles heteregeneity both in the 
## random variance matrix G and in the within-person error variance. The functions to play are 'correlation' and 'weight'.  
## In order to have a G matrix that is group-specific and/or a R matrix that is person-specific (exploratory 
## dynamic structure), one needs to specific the group variable in correlation or weight functions, respectively.

model.01 = lme(fixed = Y ~ GROUP + X1,  
               random = ~ 1 | person,  # this is the same to random=list(person = pdSymm(form = ~ 1)),
               data = df1,
               method = 'REML')
summary(model.01)

VarCorr(model.01)

# with compound symmetric error matrix (homogeneous for all persons)
model.02 = lme(fixed = Y ~ GROUP + X1, random = ~ 1 | person,
               correlation=corCompSymm(form=~day|person), na.action = na.exclude,
               data = df1, method='REML')
summary(model.02) # fit almost the same but AIC/BIC are a bit worse, not an appropriate error correlation
VarCorr(model.02)

anova(model.01, model.02) 

# with AR(1) error matrix (homogeneous for all persons)
model.03 = lme(fixed = Y ~ GROUP + X1, random = ~ 1 | person,
               correlation=corAR1(form=~day|person), na.action = na.exclude,
               data = df1, method='REML')
summary(model.03) 
VarCorr(model.03)
anova(model.01, model.03) # fit is much better, use this error correlation structure

all.equal( summary(model.02)$tTable, summary(model.01)$tTable ) # TRUE
all.equal( VarCorr(model.02), VarCorr(model.01) ) # slight difference "4 string mismatches"

# 2. Heterogeneity for between variance G (by GROUP) yet homogeneous within person variance R 
model.11 = lme(fixed = Y ~ GROUP + X1, random = list(person = pdDiag(form = ~ factor(GROUP))),
               data = df1, method='REML')
summary(model.11) # fit is similar to model.02
VarCorr(model.11)

anova(model.01, model.11)  # fit almost the same but AIC/BIC are a bit worse
# heteregeneous random variance is not necessary for Identity error matrix

# Heterogeneity for between variance G (by treatment GROUP) and homogeneous within person variance R (compound symmetry structure)
model.12 = lme(fixed = Y ~ GROUP + X1, random = list(person = pdDiag(form = ~ factor(GROUP))),
               correlation=corAR1(form=~day|person), na.action = na.exclude,
               data = df1, method='REML')
summary(model.12) 
VarCorr(model.12)

anova(model.02, model.12) # fit is much better, heteregeneous random variance 
# when using the correct error correlation structure (AR-1)

anova(model.11, model.12)
# again compound symmetry structure is NOT appropriate for the within-person residual variance

# Heterogeneity for between variance G (by individual-level continuous covariate X1 or Person) 
# and homogeneous within person variance R
model.13 = lme(fixed = Y ~ GROUP + X1, random = list(person = pdDiag(form = ~ factor(person))),
               data = df1, method='REML')
# Error in lme.formula(fixed = Y ~ 1, random = list(person = pdDiag(form = ~factor(X1))),  : 
#                       fewer observations than random effects in all level 1 groups
# continuous group covariate not estimable

# 3. Everything the same, heterogeneity in within-person residual variance by GROUP (Level-1 random effects).
# using the weights = formula in the same way for the dyadic models, in order to get the desired structure
model.21 = lme(fixed = Y ~ GROUP + X1, random = list(person = pdDiag(form = ~ factor(GROUP))),
               #correlation=corCompSymm(form=~day|person), na.action = na.exclude,
               weights = varIdent(form = ~ 1 | factor(GROUP)),
               data = df1, method='REML')
summary(model.21) # fit is similar to model.02
VarCorr(model.21)

## For the within-person variances we need to do some multiplication
# The Residual StdDev is
summary(model.21)$sigma
# and the weights are 
coef(model.21$modelStruct$varStruct, unconstrained=FALSE)
# So the estimated Residual variance for group = 0, σ2eg=0 is 
(summary(model.21)$sigma*1.0000)^2
# and the estimated variance for group = 1, σ2eg=1
(summary(model.21)$sigma*coef(model.21$modelStruct$varStruct, uncons=FALSE))^2

# To formally test if they heterogenous variances are needed we would do a log-likelihood test, e.g., using the anova() function.
anova(model.11, model.21)

# Heterogeneity for between variance G (by GROUP) AS WELL AS within person variance R by GROUP (compound symmetry structure)
model.22 = lme(fixed = Y ~ GROUP + X1, random = list(person = pdDiag(form = ~ factor(GROUP))),
               correlation=corCompSymm(form=~day|person), na.action = na.exclude,
               weights = varIdent(form = ~ 1 | factor(GROUP)),
               data = df1, method='REML')
summary(model.22) 
VarCorr(model.22)

anova(model.12, model.22)
# Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# model.12     1  7 807.4750 821.7764 -396.7375                        
# model.22     2  8 777.8113 794.1557 -380.9057 1 vs 2 31.66369  <.0001
## heterogeneity in within-person residual variance IS NEEDED!!
anova(model.21, model.22)
# compound symmetry structure is NOT appropriate for the within-person residual variance. DROP

# TEST Heterogeneous within person variance R by Person or Individual-level covariate X1 (same model)
model.23 = lme(fixed = Y ~ GROUP + X1, random = list(person = pdDiag(form = ~ factor(GROUP))),
               # correlation=corCompSymm(form=~day|person), na.action = na.exclude,
               weights = varIdent(form = ~ 1 | person),
               data = df1, method='REML')
summary(model.23) 
VarCorr(model.23)
anova(model.21, model.23) # Significant improvement based on LRT, but worse AIC/BIC upon model.21

# The Residual StdDev is
summary(model.23)$sigma # the first person is the reference group
# and the weights are 
coef(model.23$modelStruct$varStruct, unconstrained=FALSE)
# So the estimated Residual variance for group = 0, σ2eg=0 is 
(summary(model.23)$sigma*1.0000)^2
# and the estimated variance for group = 1, σ2eg=1
(summary(model.23)$sigma*coef(model.23$modelStruct$varStruct, uncons=FALSE))^2


# TEST Heterogeneous within person variance R by BOTH GROUP & Person 
#model.24 = lme(fixed = Y ~ GROUP + X1, random = list(person = pdDiag(form = ~ factor(GROUP))),
#               # correlation=corCompSymm(form=~day|person), na.action = na.exclude,
#               weights = varIdent(form = ~ 1+ 1|GROUP + 1|person), # NOT WORKING, this is for random function ONLY RECOGNIZE PERSON
#               data = df1, method='REML')
#summary(model.24) 
#VarCorr(model.24)
#anova(model.23, model.24) 

model.24 = lme(fixed = Y ~ GROUP + X1, random = list(person = pdDiag(form = ~ factor(GROUP))),
               # correlation=corCompSymm(form=~day|person), na.action = na.exclude,
               weights = varIdent(form = ~ 1 | GROUP * person), 
               data = df1, method='REML')
summary(model.24) 
VarCorr(model.24)
anova(model.23, model.24)   
# same model, because this is a fully crossed design person either fall under treatment or control

all.equal( summary(model.23)$tTable, summary(model.24)$tTable ) # TRUE
all.equal( VarCorr(model.23), VarCorr(model.24) ) # TRUE

# However, the autocorrelation parameter is the same for everything still,
# can i specify person-specific compound coefficient rho??
#model.25 = lme(fixed = Y ~ GROUP + X1, random = list(person = pdDiag(form = ~ factor(GROUP))),
#               correlation=corCompSymm(form= ~ 1|day), na.action = na.exclude,
#               weights = varIdent(form = ~ 1 | person), 
#               data = df1, method='REML')
##Error in lme.formula(fixed = Y ~ GROUP + X1, random = list(person = pdDiag(form = ~factor(GROUP))),  : 
##                       incompatible formulas for groups in 'random' and 'correlation'
# no higher order group variable in correlation function than the random function

model.25 = lme(fixed = Y ~ GROUP + X1, random = ~ 1 | person,
               correlation=corAR1(form= ~ day| person ), na.action = na.exclude,
               weights = varIdent(form = ~ 1 | person ), 
               data = df1, method='REML')

summary(model.25) 
VarCorr(model.25)
random=list(person = pdSymm(form = ~ 1))
df1$GROUP = as.numeric(df1$GROUP)
model.26 = lme(fixed = Y ~ GROUP + X1, random = list(person = pdDiag(form = ~ factor(GROUP))),
               na.action = na.exclude, correlation=croCAR1( form= ~ day | GROUP),
               weights = varIdent(form = ~ 1 | person ), 
               data = df1, method='REML')

#summary(model.26) 
model.26$modelStruct$varStruct
model.26$modelStruct$corStruct
VarCorr(model.26)
anova(model.25, model.26)   

all.equal( summary(model.25)$tTable, summary(model.26)$tTable ) # TRUE
all.equal( VarCorr(model.25), VarCorr(model.26) ) # TRUE


######################################################################################################################################
## PART II Implement the MLM-CT Procedure Using A Simple Model (unconditional) on the Simulated Data
df1$person = as.factor(df1$person)
arima.fit = list() # list to store individual arima model result object
y_hat = list() # list to store individual fitted response values 
int = vector() # vector to store individual fitted intercept
rho1 = vector() # vector to store individual autocorrelation parameter
sigma2_w = vector() # vector to store individual residual variance
trans.Y = vector() # vector to store the transformed response variable
trans.Int = vector() # vector to store the transformed intercepts
A = list() # list to store individual transformation matrix A's

## Using a for loop to repeat a individual arima (1,0,0) model for all parcitipants, respectively.
for(i in 1:length(levels(df1$person))){
  arima.fit[[i]] = arima(df1$Y[df1$person == levels(df1$person)[i]],  order = c(1,0,0))
  y_hat[[i]] = fitted(arima.fit[[i]])
  rho1[i] = arima.fit[[i]]$coef[1]
  int[i] = arima.fit[[i]]$coef[2]
  sigma2_w[i] = arima.fit[[i]]$sigma2
  # Calculate first-order autocovariance matrix Sigma
  struc.matrix = matrix(c(1, rho1[i], rho1[i]^2, rho1[i]^3, rho1[i]^4, rho1[i]^5, 
                          rho1[i], 1, rho1[i], rho1[i]^2, rho1[i]^3, rho1[i]^4, 
                          rho1[i]^2, rho1[i], 1, rho1[i], rho1[i]^2, rho1[i]^3, 
                          rho1[i]^3, rho1[i]^2, rho1[i], 1, rho1[i], rho1[i]^2, 
                          rho1[i]^4, rho1[i]^3, rho1[i]^2, rho1[i], 1, rho1[i], 
                          rho1[i]^5, rho1[i]^4, rho1[i]^3, rho1[i]^2, rho1[i], 1), 6, 6)
  autocov_Sigma = sigma2_w[i]/(1-rho1[i]^2) * struc.matrix
  # Calculate A matrix (inverse of the Cholesky factor of the autocovariance matrix)
  autocov_Sigma_inv = solve(autocov_Sigma)
  v_inv = sigma2_w[i] * autocov_Sigma_inv
  L.chol <- chol(v_inv) 
  # Verify: t(L.chol) %*% L.chol == v_inv
  A[[i]] = solve(t(L.chol))
  # create individual set of transformed vectors 
  obs = sum(df1$person == levels(df1$person)[i])
  index = ((i-1)*obs+1):((i-1)*obs+obs)
  trans.Y[index] = A[[i]] %*% df1$Y[df1$person == levels(df1$person)[i]]
  trans.Int[index] = A[[i]] %*% rep(int[i], obs)
}

df1$trans.Y = trans.Y 
df1$trans.Int = trans.Int

# Fit MLM assumming the ID independent residual structure
# using trans.Y, every individual should have error residual that is approximately a white noise series
# it is now okay to assume the ID or other independent resitual structure across people.

model.trans.01 = lme(fixed = trans.Y ~ trans.Int, random = ~ trans.Int | person,
                     correlation=corAR1(form=~day|person),na.action = na.exclude,
                     data = df1, method='REML')

# This contrast the MLM model using original data/variable that assumes the homogeneous AR(1) error correlation for all individuals
#model.orig.01 = lme(fixed = Y ~ 1, random = ~ 1 | person,
#                    correlation=corAR1(form=~day|person), na.action = na.exclude,
#                    data = df1, method='REML')

summary(model.trans.01) 
summary(model.orig.01)
VarCorr(model.trans.01)
VarCorr(model.orig.01)
anova(model.trans.01, model.orig.01)
# cannot test LRT on models using different response variables
