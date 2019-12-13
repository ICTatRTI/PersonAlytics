---
title: "PersonAlytics&copy; User's Guide"
author: "Stephen Tueller, Ty Ridenour, Derek Ramirez"
date: "2019-12-13"
output: 
  rmarkdown::html_vignette:
    number_sections: true
    toc: true
urlcolor: blue
vignette: >
  %\VignetteIndexEntry{PersonAlytics User's Guide}
  %\VignetteEngine{knitr::knitr}
  %\usepackage[UTF-8]{inputenc}
---



<style>
p.comment {
background-color: #DBDBDB;
padding: 10px;
border: 1px solid black;
margin-left: 25px;
border-radius: 5px;
font-style: italic;
}
</style>



# Introduction

The purpose of this user's guide is to illustrate how to use `PersonAlytics` to analyze 

1. Single case time series data, also known as single subject or N-of-1 studies. Interrupted time series design designs can be used to introduce one (or more) interventions which results in two (or more) phases (e.g., pre-intervention and post-intervention). 

2. Small sample intensive longitudinal design data.

3. Ideographic clinical trial data (ICT), which a single-case or small sample longitudinal data set where all participants experience two or more treatment phases (e.g., baseline and follow-up). 

These three types of data are analyzed using longitudinal mixed effects models, also known as hierarchical linear models, multilevel models, latent growth curve models, or mixed method trajectory analysis. There are already several R packages such as `nlme` and `gamlss` that can be used to implement these models. 

The purpose of the PersonAlytics package is to automate the following aspect of analyzing these type of data: 

1. Model selection for detecting the shape of the trajectory over time. Using a good fitting trajectory shape yields parameters that correspond to the shape and improve interpretability. Current options include polynomial growth (e.g., linear, quadratic, cubic, etc.) and piecewise growth (e.g., linear growth within each phase). 

2. Model selection for detecting an appropriate time series model for the residuals. Getting a good fitting residual correlation model improves standard error estimation.

3. Automated, paralellized, high throughput analyses of mixed effects models in situations where the combination of predictors, outcomes, and/or individuals (i.e., in the case individual level models are desired as may be the case in personalized medicine) yields a number of analyses to unweildly implement manually. This is especially true if the model selection process for the trajectory shape and the residual correlation model is conducted for each combination of predictors, outcomes, and/or individuals. Included in PersonAlytics are options for the user to specify Type I error rate or False Discovery Rate (FDR) corrections. 

<p class="comment">
**Parallelization**. High throughput analyses are achieved through parallelization. Jobs are split among two or more processors on a computer and are run in parallel to each other. Results are recombined at the end of the process. 
</p>

<p class="comment">
**High Throughput Example 1: Migraine Triggers**. 346 migraine patients were followed for 90 days. They recorded information on 71 potential migraine and non-migraine headache triggers such as alcohol, weather, and exercise. Individual models were required to determine the five person-specific triggers with the largest effect size. Ever though 90 time points is large, it was deemed to small to estimate all 71 trigger effects simultaneously, so trigger effects were estimated one at a time. The analysis required 346 patients X 2 outcomes X 71 triggers = 49,132 PersonAlytic runs.
</p>

<p class="comment"> 
**High Throughput Example 2: THC Metabolomics**. 17 patients participated in a two-phase design. The baseline phase had 2 hours of observation. The intervention was 25mg of THC and the intervention phase had 6 hours of observation. A total of 20 time points generated blood metabolomic data and outcomes including sleepiness, reaction time, memory, and behavior. The analysis required 18,023 chemical compounds (the metabolites) X 8 outcomes = 144,184 PersonAlytic runs.  
</p>


# The `PersonAlytic` framework 

Ty to write content here

# Installing `PersonAlytics`

First the user must install `R`, available at `https://cran.r-project.org/`. It is also suggested that the user use a modern code editor such as Rstudio (`https://www.rstudio.com/`). It is assumed that the reader is familiar with basic `R` use. If needed, an internet search will provide tutorials to help the user become familiar with R. 

Installing `PersonAlytics` is done by using the `install_github` function of the `devtools` package. First install `devtools` using


```r
install.packages('devtools')
```

Then use


```r
devtools::install_github("https://github.com/ICTatRTI/PersonAlytics")
```

The `install_github` function may give the user the option update other R packages needed by `PersonAlytics` that the user already have. Unless the user are experienced with R package versions and updates, it is recommended that the user ignore this step and press enter without selecting an option to update packages. If the user do decide to update packages and run into an error, the user may need to restart `R`, manually update the package in the error message, and then rerun the `install_github` command. 

Once `install_github` starts running, it will install all of the other `R` packages required by `PersonAlytics` if the user do not already have them. 

# Starting `PersonAlytics`

Now we can load `PersonAlytics`:


```r
library(PersonAlytics)
```

# Basic `PersonAlytics` Use

We will illustrate the basic options of the `PersonAlytic()` function using the Ovary data from the `nlme` package which have been modified to represent an ICT. The `PersonAlytic()` function is the primary user interface for the `PersonAlytics` package. the user can access the documentation for `PersonAlytic()` by typing


```r
?PersonAlytic
```

## The OvaryICT Data

The first six rows of the OvaryICT data are shown in Table x where it can be seen that the original Ovary data set has been augmented with Phase variables which are an essential component of an ICT but which are not required by the PersonAlytics package. Note that the data are in 'long format'. In the example below, Mare 1's first six time points are represented in separate rows. The time points are days prior to ovulation (for negative values), at ovulation (for values of 0 or 1) or between ovulations (for positive values other than 0 and 1). There are also six randomly generated predictors named Target1 to Target6 which will be used to illustrate other `PersonAlytics` features. 


|   | Mare|  Time| follicles| Phase|Phase2 |Target1 |Target2 |Target3 | Target4| Target5| Target6|
|:--|----:|-----:|---------:|-----:|:------|:-------|:-------|:-------|-------:|-------:|-------:|
|83 |    1| -0.14|         9|     0|1      |2       |1       |4       |   -0.53|    0.18|    0.74|
|84 |    1| -0.09|         9|     0|1      |2       |1       |2       |   -1.37|    1.34|    0.63|
|85 |    1| -0.05|         7|     0|1      |4       |4       |4       |   -2.21|    1.20|    1.24|
|86 |    1|  0.00|         6|     0|1      |1       |3       |4       |    1.82|    0.87|    0.23|
|87 |    1|  0.05|         7|     0|1      |3       |1       |4       |   -0.65|   -0.12|   -0.31|
|88 |    1|  0.09|         6|     0|1      |1       |2       |3       |   -0.28|    0.34|    1.50|

## Required `PersonAlytics` Parameters

The four required parameters for `PersonAlytic()` are

1. `data`: the name of the user's data set. This data needs to have been read into R prior to using `PersonAlytic()`. A web search can be used to learn how to read in multiple types of data into R, including (but not limited to) csv, xlsx, sas, stata, spss, and most database formats. The data must be structured in 'long format' where time points are repeated within individual as was illustrated above for the `OvaryICT` data. In the example below, we use the `OvaryICT` data set.

2. `ids`: the name of the identification variable for individuals. This must be a quoted variable name matching the user's ID variable in the data set that the user provided to the parameter `data`. 

3. `dvs`: the name of the dependent variable. This must be a quoted variable name matching the user's dependent variable in the data set that the user provided to the parameter `data`. If the user have multiple dependent variables, `PersonAlytics` will iterate over them one at a time. Multiple dependence variables are specified as a character list, e.g., `dvs=list('dv1', 'dv2', etc.)`. This is a high throughput option.

4. `time`: the name of the time variable. This must be a quoted variable name matching the user's dependent variable in the data set that the user provided to the parameter `data`.

Here is an illustration using the four basic parameters. This example also includes setting `autoDetect=NULL`, which turns off automated model comparisons for selecting the trajectory shape and residual correlation structure until this is discussed in its own section. 


```r
eg_required <- PersonAlytic(data=OvaryICT,
                ids="Mare",
                dvs="follicles",
                time="Time",
                autoDetect=NULL)
```

It is important to assign the results of a call to `PersonAlytic` to an R object. In the example above, the object is named `eg_required`. If the user neglect to assign the results to an R object, minimal output will be printed to the screen and the remaining results will be lost. If high throughput is initialized by including more than one dependent variable (discussed above), more than one target independent variables (see parameter `target_ivs` below), or if individual models are requested (see parameter `individual_mods` below), the resulting object is a `data.frame` concatenating results across of the options submitted to the high throughput run. Reading the output is detailed in the next section.

# Reading the Output

The argument `output` is character string that will be used to name a file for saving output. If left `NULL`, the default is 'PersonAlytic_Output'. Do not give a file extension, these will be added automatically depending on whether a single analysis was run versus when high throughput options are invoked.

## Single Analysis Output

Here we fit the same model using `nlme` and `gamlss` and illustrate options for viewing and manipulating the output. If the user are familiar with the `lme` package, this is an implementation of the linear (normal) mixed effects model. Using the `package` parameter, the `gamlss` package can alternatively be used in conjunction with the `family` parameter for outcomes with non-normal distributions. The `family` parameter is described in more detail later, the default distribution assumed when using `package="gamlss"` is the normal distribution.  

First, we run the model using `package="nlme"`. Since this is the default, we need to explicitly specify the `package` parameter. The output will be saved the the R object `eg_nlme` for further use inside the R console. In addition, `output="nlme_example"` will create a file in the working directory named 'nlme_example.txt'. If the user are unsure what the user's current working directory is, type `getwd()` into the R console. A full path with '/' instead of '\' can also be used. For example, `output='C:/MyResults'`. For convenience, we will also turn off the `autoDetect` parameter by setting it to an empty list. The `autoDetect` parameter is discussed in detail below. 


```r
eg_nlme <- PersonAlytic(output="nlme_example",
                data=OvaryICT,
                ids="Mare",
                dvs="follicles",
                phase="Phase",
                time="Time",
                autoDetect=NULL)
```

Then we run the same model using `package="gamlss"`. The output will be saved the the R object `eg_gamlss` for further use inside the R console. In addition, `output="gamlss_example"` will create a file in the working directory named 'gamlss_example.txt'.


```r
eg_gamlss <- PersonAlytic(output="gamlss_example",
                data=OvaryICT,
                ids="Mare",
                dvs="follicles",
                phase="Phase",
                time="Time",
                package="gamlss",
                autoDetect=NULL)
```

In these examples, the R object `eg_nlme` is a class `lme` object, and the R object `eg_gamlss` is a class `gamlss` object.


```r
class(eg_nlme)
```

```
## [1] "lme"
```

```r
class(eg_gamlss)
```

```
## [1] "gamlss" "gam"    "glm"    "lm"
```

Both `lme` and `gamlss` objects have a `summary()` method which prints detailed results to the R console. For the `nlme` example we get the following.


```r
summary(eg_nlme)
```

```
Linear mixed-effects model fit by REML
 Data: tempData 
       AIC      BIC    logLik
  1670.449 1700.185 -827.2244

Random effects:
 Formula: ~Time | Mare
 Structure: General positive-definite, Log-Cholesky parametrization
            StdDev   Corr  
(Intercept) 2.751709 (Intr)
Time        3.871676 -0.202
Residual    3.317972       

Fixed effects: follicles ~ Time * Phase 
               Value Std.Error  DF   t-value p-value
(Intercept) 10.66163 0.9016267 294 11.824883  0.0000
Time        -0.86898 1.7667760 294 -0.491845  0.6232
Phase       10.97209 1.3126329 294  8.358845  0.0000
Time:Phase  -8.64385 1.9742483 294 -4.378299  0.0000
 Correlation: 
           (Intr) Time   Phase 
Time       -0.318              
Phase      -0.105  0.135       
Time:Phase  0.175 -0.504 -0.817

Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max 
-2.49764683 -0.65159342 -0.01330142  0.61976543  2.61092061 

Number of Observations: 308
Number of Groups: 11 
```

The fixed effects results in the section titled 'Fixed Effects' is what gets saved to the file 'nlme_example.txt'. For the `gamlss` example we get the following: 


```r
summary(eg_gamlss)
```

```
******************************************************************
Family:  c("NO", "Normal") 

Call:  
gamlss::gamlss(formula = self$formula, sigma.formula = sigma.formula,  
    family = currentFamily, data = tempData, control = ctrl) 

Fitting method: RS() 

------------------------------------------------------------------
Mu link function:  identity
Mu Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)  10.6898     0.3404  31.402  < 2e-16 ***
Time         -0.8637     1.2787  -0.675      0.5    
Phase        10.9725     1.2654   8.671 3.28e-16 ***
Time:Phase   -8.6439     1.9033  -4.541 8.25e-06 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

------------------------------------------------------------------
Sigma link function:  log
Sigma Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)  1.16333    0.04029   28.87   <2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

------------------------------------------------------------------
NOTE: Additive smoothing terms exist in the formulas: 
 i) Std. Error for smoothers are for the linear effect only. 
ii) Std. Error for the linear terms may not be reliable. 
------------------------------------------------------------------
No. of observations in the fit:  308 
Degrees of Freedom for the fit:  23.39381
      Residual Deg. of Freedom:  284.6062 
                      at cycle:  2 
 
Global Deviance:     1590.679 
            AIC:     1637.467 
            SBC:     1724.728 
******************************************************************
```

The fixed effects results in the section titled 'Mu Coefficients' and the variance coefficients in the section titled 'Sigma Coefficients' is what gets saved to the file 'gamlss_example.txt'. 

It is beyond the scope of this document to detail the differences between the `nlme` and `gamlss` approaches, but we will highlight one important difference. Notice that the parameter estimates in the 'Value' column of the 'Fixed Effects' section of the `eg_nlme` output and the 'Estimate' column of the 'Mu Coefficients' section of the `eg_gamlss` output are very similar. However, the standard errors are of the `gamlss` model are lower. This is because the `gamlss` approach models not only the mean, but also the variance. If there is any heteroscedasticity of variance, the `gamlss` results will consequently have lower standard errors than the `nlme` results. In the `gamlss` output, the variance parameter(s) are in the 'Sigma Coefficients' section.   

## High throughput output

Here we illustrate a high throughput example by first creating two additional outcomes and fit a basic growth model to each outcome. Note that creating variables that are the square and the root of the outcome is not advised unless the user has substantive reasons to do so, this is simply for illustration purposes. 


```r
OvaryICT$follicles2 <- OvaryICT$follicles^2
OvaryICT$folliclesr <- sqrt(OvaryICT$follicles2)
eg_htp <- PersonAlytic(output="htp_example",
                data=OvaryICT,
                ids="Mare",
                dvs=list("follicles", "follicles2", "folliclesr"),
                phase="Phase",
                time="Time",
                autoDetect=NULL)
```

```
## 
## 
## Model fitting starting...
```

```
##   |                                                                         |======================                                           |  33%  |                                                                         |===========================================                      |  67%  |                                                                         |=================================================================| 100%
```

```
## 
## 
## Model fitting took:
##  Time difference of 4.909943 secs.
```

The output will be saved the the R object `eg_htp` for further use inside the R console and `output="htp_example"` will create a file in the working directory named 'htp_example.csv' that is the same thing as the R object `eg_htp`. This csv file can be opened in any spreadsheet program. The rows are all possible combinations of `dvs`, `target_ivs`, and `ids` (if `individual_mods=TRUE`). The column variables are list below and come in five sets:

1. Variables identify the combination of user inputs that lead to a given row's analysis. Most column names correspond to their respective parameter name in `PersonAlytic`. Other variables include information on the version of 'PersonAlytics', date, time, and the directory in which the models were run.

    * 'Mare': This is is the name of the id variable passed to the `ids` parameter. In the current example, the value is 'All Cases'. If indivdiual models are requested by setting `individual_mods=TRUE`, this column will give the id for each case in the data set.
    
    * 'ids': This is the name of the `ids` variable which is the column name of column 1. 
    
    * 'dv': The name of the dependent variable.
    
    * 'time': The name of the time variable.
    
    * 'ivs': The names of the indepedent variables (describe below).
    
    * 'target_iv': The name of a target independent variable (described below).
    
    * 'interactions': The names of the interaction terms (described below). 
    
    * 'time_power': The shape of the trajectories over time (described below).
    
    * 'alignPhase': How time was realigned by phase, if any (described below). 
    
    * 'correlation': The residual correlation structure (described below). If the value is `NULL`, this  corresponds to no within-group correlations (see `?lme`). 
    
    * 'family': The distribution for the dependent variable (described below).
    
    * 'standardize': Which variables were standardized (described below).
    
    * 'method': The estimation method.
    
    * 'package': Which R package was used to fit the models.
    
    * 'PersonAylitics': The version of the `PersonAylitics` used for the analysis.
    
    * 'Date_Time': The date and time the model was run. 
    
    * 'estimator': What model estimator was used (described below).
    
    * 'analyzed_N': The number of observations analyzed. This is the sum of time points (per case if multiple cases were analyzed). 
    
    * 'call': The model formula created from the user inputs.
    
    * 'wasLRTrun': If 'target_ivs' were provided, a likelihood ratio test (LRT) for models with and without the target independent variable will be attempted, and if succesfull, 'wasLRTrun' will be 'TRUE'. 
    
    * 'targ_ivs_lrt_pvalue': If the LRT was run, the p-value is recorded here. 
    
    * 'fixed': The fixed effects portion of the 'call'.
    
    * 'random': The random effects portion of the 'call'.
    
    * 'formula': This is similar to the 'call' variable which gives the intended formula. If the model will not converge using the intended formula, simplifications of the formula are attempted in the following order: 
    
        + No correlation structure
        
        + No correlation structure and no random slopes
        
        + If the model is piecewise, drop the phase by time interaction
        
    * 'correlation0': The correlation portion of the 'call'.
    
    * 'directory': The directory where model output is saved.
    
    * 'date': The date the output was saved (which may be different from the 'Date_Time'  model was run for long runs).

2. Variables describing whether the analysis converged and variables to help diagnose a failed run (e.g., zero variance in an outcome or independent variable). 

    * 'N_participants': The number of unique cases in the `ids` variable.
    
    * 'N_time_points': The number of unique time points. 
    
    * 'Nobs': The number of individuals across all unique time points. If the time points are all unique, this number will be the same as 'N_time_points'. If cases share time points, 'Nobs' will be smaller than 'N_time_points'. 
    
    * 'dvVar': The variance of the dependent variable. If this is zero, the model will not converge.
    
    * 'timeVar': A check whether the time variable is monotonically increasing. If it is not, the `PersonAlytic` function will stop with an error.
    
    * 'ivVar': The variance of the independent variables (if any are included in the model).
    
    * 'target_ivVar': The variance of the target independent variable (if one is included in the model).
    
    * 'converge': Model convergence status.

3. Descriptive statistics in pairs with the first column describing the statistic with the prefix `statName` and the second column in each pair with the prefix `statValue` giving the statistic's value. 
4. Model results with a parameter estimates, standard error, t-value, degrees of freedom, and p-value. 

5. If the finite population correction (FPC) is specified (see below for details), the model results repeated with FPCs for the standard errors (and consequently, the p-values). The example below doesn't include the FPC, but if it did the parameter estimates, standard errors, t-values, degrees of freedom, and p-values would be repeated with the suffix `fpc`. 


```
 [1] "Mare"                   "ids"                   
 [3] "dv"                     "time"                  
 [5] "phase"                  "ivs"                   
 [7] "target_iv"              "interactions"          
 [9] "time_power"             "alignPhase"            
[11] "correlation"            "family"                
[13] "standardize"            "method"                
[15] "package"                "Personalytics"         
[17] "Date_Time"              "N_participants"        
[19] "N_time_points"          "Nobs"                  
[21] "dvVar"                  "timeVar"               
[23] "ivVar"                  "target_ivVar"          
[25] "converge"               "estimator"             
[27] "analyzed_N"             "call"                  
[29] "wasLRTrun"              "targ_ivs_lrt_pvalue"   
[31] "fixed"                  "random"                
[33] "formula"                "correlation0"          
[35] "directory"              "date"                  
[37] "statName1"              "statValue1"            
[39] "statName2"              "statValue2"            
[41] "statName3"              "statValue3"            
[43] "X.Intercept..Value"     "X.Intercept..Std.Error"
[45] "X.Intercept..DF"        "X.Intercept..t.value"  
[47] "X.Intercept..p.value"   "Time.Value"            
[49] "Time.Std.Error"         "Time.DF"               
[51] "Time.t.value"           "Time.p.value"          
[53] "Phase.Value"            "Phase.Std.Error"       
[55] "Phase.DF"               "Phase.t.value"         
[57] "Phase.p.value"          "Time.Phase.Value"      
[59] "Time.Phase.Std.Error"   "Time.Phase.DF"         
[61] "Time.Phase.t.value"     "Time.Phase.p.value"    
```

# Autodetection of the Residual Autocorrelation Structure, Time Order, and Distribution

The autodetection options automate the tedious process of conducting ML model comparisons to determine any of three options described in this section. The value `NULL` (`autoDetect=NULL`), or an empty list (`autoDetect=NULL`) turns all options off. Leaving any option out of the list will turn that option off (e.g., `autoDetect=list(TO=list(polyMax=3))` will only implement autodetection of the time order).

## Residual Correlation Structure `AR`

The autoregressive moving-average (ARMA) order of the residual correlation structure. Determining a good fitting correlation structure will provide more accurate standard errors. The default, `AR=list(P=3, Q=3)` will search all combinations of `p=c(0,1,2,3)` and `q=c(0,1,2,3)`, as well as the default with no within-group correlations (`correlation=NULL`), for the best fitting ARMA correlation structure. This is done with a fit index instead of ML likelihood ratio tests (LRT) because not all ARMA models are nested (a requirement of the LRT). The fit index that will be used is set using the `whichIC` parameter with options `BIC` or `AIC`. If $n=1$ or `individual_mods=TRUE`, correlation model selection is implemented using the `auto.arima` function of the `forecast` package. See `?auto.arima`.

Here is an example of using the `AR` option of the `autoDetect` parameter without the `TO` or `DIST` options:















































