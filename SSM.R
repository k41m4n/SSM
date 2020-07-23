 
# R code for 
#  ****************************************************************
#  *  An Introduction to State Space Time Series Analysis (2007). *
#  *  Jacques J.F. Commandeur and Siem Jan Koopman.			          *
#  *  Oxford: Oxford University Press.							              *
#  ****************************************************************

#R LIBRARIES, FOLDERS and FUNCTIONS####

#Installing, if needed, packages and loading the installed pacakages from the library
if(!(require(normtest))){install.packages('normtest')}
library(normtest)
if(!(require(KFAS))){install.packages('KFAS')}
library(KFAS)
if(!(require(rstudioapi))){install.packages('rstudioapi')}
library(rstudioapi)
if(!(require(knitr))){install.packages('knitr')}
library(knitr)
if(!(require(eurostat))){install.packages('eurostat')}
library(eurostat)
if(!(require(dplyr))){install.packages('dplyr')}
library(dplyr)
if(!(require(forecast))){install.packages('forecast')}
library(forecast)

#Cleaning workspace
rm(list=ls())
#Setting directory for files with data; they should be in the same directory as the files of source code
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())


#Function for Q-statistic
#Q-statistic is a general omnibus test that can be used to check whether 
#the combined first k autocorrelations significantly deviate from 0,
#meaning the null hypothesis of independence must be rejected
#We assume that the residuals are independent if the test statistic does not exceed the critical value.
qStatistic <- function(predResid, k, w) {
  #Standardised residuals as predResid should be submitted into this function!
  #k - first k autocorrelations to be used in test; w - number of the disturbance variances
  #(see Commandeur and Koopman, p.90-96)
  value <- Box.test(predResid, lag = k, type = "Ljung")$statistic #Q-statistic based on the statistic calulated in the Ljung-Box test for lags from 1 to k
  criticalValue <- qchisq(0.95, k-w+1) #Critical value corresponding to the upper 5% in the chi-square-distribution with k-w+1 degrees of freedom
  list(#List of values provided by the function
    k = k, #First k autocorrelations
    value = unname(value), #Value of the test statistic
    criticalValue = criticalValue #Critical value
  )
}

#Function for r-statistic
#r-statistic checks independence of one-step-ahead prediction residuals 
#It provides values of the autocorrelations at lags 1 and l, together with the 95% confidence limits. 
#We assume that the residuals are independent if the values of the autocorrelations 
#do not exceed the critical values for the 95% confidence limits 
#(see Commandeur and Koopman, p.90-96)
rStatistic <- function(predResid, d, l) {
  #Standardised residuals as predResid should be submitted into this function!
  #d - diffuse initial value of the state, l - autocorrelation at lag l to be provided by the function
  n <- (length(predResid)-d) #Length of the series after subtracting d, the number of diffuse initial elements of the state
  acfValues <- acf(predResid[-(1:d)], plot = FALSE)$acf[-1]# List of the values of the autocorrelations for the series without the first d values
  criticalValue <- 2 / sqrt(n) # +/- critical value for the 95% confidence limits
  list( #List of values provided by the function
    l = c(1,l), #lags 1 and l 
    value1 = acfValues[1], #Value of the autocorrelation at lag 1
    value2 = acfValues[l], #Value of the autocorrelation at lag l
    criticalValue = criticalValue # +/- critical value for 95% confidence limits
  )
}

#Function for H-statistic
#H-statistic checks homoscedasticity of one-step-ahead prediction residuals 
#This is done by testing the null hypothesis of the equal variances of the residuals 
#in the first third part of the series and the last third part of the series.
#The ratio between these two variances is tested against an F-distribution with (h,h) degrees of freedom 
#applying the usual 5% rule for rejection of the null hypothesis of equal variances, for a two-tailed test.  
#We must find critical values corresponding to the upper and lower 2.5% in the two tails of the F-distribution;
#If, however, the tested statistic is larger than or equal to 1, it is enough to check 
#whether it is lower than the critical value #corresponding to the upper 2.5% in the F-distribution; 
#On the other hand, if the statistic is lower than 1 we have to test 
#if its reciprocal value (1/ratio) is lower than the above-mentioned critical value.
#We assume that the residuals are homoscedastic if the test statistic does not exceed the critical value.
#(see Commandeur and Koopman, p.90-96)
hStatistic <- function(predResid, d) {
  #Standardised residuals as predResid should be submitted into this function!
  #d - number of diffuse initial values in the state,
  n <- length(predResid) # Number of observations/residuals
  h <- round((n-d)/3, digits = 0) #One third of the series: nearest integer to (n-d)/3; also degrees of freedom for the test
  ratio <- sum(predResid[(n-h+1):n]^2) / sum(predResid[(d+1):(d+h)]^2) #Ratio between the variance of the residuals in the last third part of the series and the variance of residuals in the first third part of the series
  value <- ifelse(ratio >= 1, ratio, 1/ratio) # Value of the test statistic; if the ratio is smaller than 1 then the reciprocal value is used for testing (1/ratio)
  criticalValue <- qf(0.975, h, h) # Critical value corresponding to the upper 2.5% in the F-distribution with (h,h) degrees of freedom
  list( #List of values provided by the function
    h = h, #Degrees of freedom
    ratio = ratio, #Ratio between the two variances
    value = value, #Value of the test statistic
    criticalValue = criticalValue #Critical value 
  )
}

#Function for N-statistic
#H-statistic checks normality of one-step-ahead prediction residuals 
#This is done by testing the null hypothesis of normality
#We assume that the residuals are normally distributed if the test statistic does not exceed the critical value at 5% level
#(see Commandeur and Koopman, p.90-96)
nStatistic <- function(predResid, d) {
  #Standardised residuals as predResid should be submitted into this function!
  #d - number of diffuse initial values in the state
  value <- jb.norm.test(predResid[-(1:d)])$statistic #N-statistic based on the statistic calculated in the Jarque and Bera or Shenton and Bowman test;
  criticalValue <- qchisq(0.95,2) #Critical value corresponding to the upper 5% in the chi-square-distribution with 2 degrees of freedom
  list(#List of values provided by the function
    value = unname(value), #Value of the test statistic
    criticalValue = criticalValue #Critical value
  )
}

#Function to create a table with statistics
dTable <- function(qStatistic, rStatistic, hStatistic, nStatistic, title){

cat(title)
cat("\n")

diagnosticTemplateTable <- c(
  
  "-----------------------------------------------------------------------------",    
  "                    statistic    value   critical value   asumption satisfied",    
  "-----------------------------------------------------------------------------",    
  "independence           Q(%2d)   %7.3f            %5.2f        %1s",  # Q-statistic, 4 args    
  "                        r(%1d)   %7.3f           +-%4.2f        %1s", # r-statistics,      4 args    
  "                       r(%2d)   %7.3f           +-%4.2f        %1s", # r,      4 args    
  "homoscedasticity     %-3s(%2d)   %7.3f            %5.2f        %1s",  # Homo,     5 args    
  "normality                  N   %7.3f            %5.2f        %1s",    # N,        3 args    
  "-----------------------------------------------------------------------------"  
) 

cat(    sprintf(      paste(diagnosticTemplateTable, collapse = "\n"),       
                      # Q-statistic, 4 args      
                      qStatistic$k,          
                      qStatistic$value,      
                      qStatistic$criticalValue,      
                      ifelse(qStatistic$value < qStatistic$criticalValue, "+", "-"),     
                      # r-statistic, 4 args      
                      rStatistic$l[1],       
                      rStatistic$value1,      
                      rStatistic$criticalValue,      
                      ifelse(abs(rStatistic$value1) < rStatistic$criticalValue, "+", "-"),      
                      # r-statistic, 4 args   
                      rStatistic$l[2],   
                      rStatistic$value2,
                      rStatistic$criticalValue,      
                      ifelse(abs(rStatistic$value2) < rStatistic$criticalValue, "+", "-"),      
                      # H-statistic, 5 args      
                      ifelse(hStatistic$ratio > 1, "  H", "1/H"),      
                      hStatistic$h,       
                      hStatistic$value,       
                      hStatistic$criticalValue,       
                      ifelse(hStatistic$value < hStatistic$criticalValue, "+", "-"),      
                      # N, 3 args      
                      nStatistic$value,        
                      nStatistic$criticalValue,      
                      ifelse(nStatistic$value < nStatistic$criticalValue, "+", "-")    )  )  

}

#Function to find best initial values for optim ver. 1 
initValOpt <- function(w_ = w , model_ = model, updatefn_ = ownupdatefn, method = "Nelder-Mead", maxLoop = 100){
  results  <- matrix(NA, maxLoop, 2) %>% 
    data.frame() %>%
    `colnames<-`(c("Initial.value", "Log.likelihood"))
  #set.seed(123)
  cat("Loop: ")
  for (j in 1:maxLoop){
    cat(paste(j, " "))
    x <- runif(1, min = 0.00001, max = 2) %>% round(3)
    fit <- fitSSM(inits = log(rep(x, w_)), model = model_, updatefn = updatefn_, method = method)
    maxLik <- logLik(fit$model, method = method)/n
    results[j, ] <- c(x, maxLik)
  }     
  cat("\n")
  results %>% arrange(desc(Log.likelihood)) %>% arrange(Initial.value) %>% print()
  return(results[1,1])
}

#Function to find best initial values for optim ver. 2
initValOpt2 <- function(formula = "log(rep(x, 3))", model_ = model, updatefn_ = ownupdatefn, method = "Nelder-Mead", maxLoop = 100){
  results  <- matrix(NA, maxLoop, 2) %>% 
    data.frame() %>%
    `colnames<-`(c("Initial.value", "Log.likelihood"))
  #set.seed(123)
  cat("Loop: ")
  for (j in 1:maxLoop){
    cat(paste(j, ""))
    x <- runif(1, min = 0.00001, max = 2) %>% round(3)
    fit <- fitSSM(inits = eval(parse(text = formula)), model = model_, updatefn = updatefn_, method = method)
    maxLik <- logLik(fit$model, method = method)/n
    results[j, ] <- c(x, maxLik)
  }     
  cat("\n")
  results %>% arrange(desc(Log.likelihood)) %>% arrange(Initial.value) %>% print()
  return(results[1,1])
}


#CHAPTER 1: Introduction####

data <- log(read.table("UKdriversKSI.txt"))
colnames(data) <- "logUKdriversKSI"
time <- 1:nrow(data)
fit <- lm(data$logUKdriversKSI~time)
(coef <- fit$coefficients) # Coefficients of regrerssion
f.stat <- summary(fit)$fstatistic
(f.stat.val <- f.stat[1]) # F-test value
(f.stat.p <- pf(f.stat[1], f.stat[2], f.stat[3], lower.tail = F)) # p-value for F-test
(error.var <- summary(fit)$sigma^2) # Error variance

#Figure 1.1. Scatter plot of the log of the number of UK drivers KSI
#against time (in months), including regression line
plot(data$logUKdriversKSI,  col = "darkgrey", xlab = "",ylab = "log UK drivers KSI",
    pch = 3,cex = 0.5, cex.lab = 0.8,cex.axis = 0.9,xlim = c(0,200))
abline(coefs , col = "blue", lwd  = 2, lty = 2)
title(main = "Figure 1.1. Scatter plot of the log of the number of UK drivers KSI
against time (in months), including regression line", cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI against time (in months)",
    "regression line"), cex = 0.5,lty = c(0, 2), col = c("darkgrey","blue"),
    pch = c(3,NA), bty = "y",horiz = T)

#Figure 1.2. Log of the number of UK drivers KSI plotted as a time series
plot(ts(data$logUKdriversKSI),ylab = "",xlab = "",xlim = c(0,200), col = "darkgrey")
title(main = "Figure 1.2. Log of the number of UK drivers KSI plotted as a time series", 
    cex.main = 0.8)
legend("topright",leg = "log UK drivers KSI",cex = 0.5,lty = 1, col = "darkgrey",horiz = T)

#Figure1.4. Correlogram of random time series
random.series <- rnorm(nrow(data))
Acf(c(random.series), 15, main = "", ylab = "")
title(main = "Figure1.4. Correlogram of random time series", 
      cex.main = 0.8)
legend("topright",leg = "ACF - random residuals",cex = 0.5,lty = 1, col = "black",horiz = T)

#Figure 1.5. Correlogram of calssical regression residuals
residuals <- residuals(fit)
Acf(c(residuals), 15,main = "", ylab = "")
title(main="Figure 1.5. Correlogram of calssical regression residuals", 
      cex.main=0.8)
legend("topright",leg = "ACF - regression residuals",cex = 0.5,lty = 1, col = "black",horiz = T)


 #CHAPTER 2: The local level model####

#2.1 Deterministic level####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% ts(start = 1969,frequency = 12)

#Defining model
model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree = 1, Q = list(matrix(0))), H = matrix(NA))

ownupdatefn <- function(pars,model){
model$H[,,1] <- exp(pars[1])
model
}

d <- q <- 1 #Number of diffuse initial values in the state 
w <- 1 #Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15 #First k autocorrelations to be used in Q-statistic
n <- 192 #Number of observations

#Fitting model and getting output

fit <- fitSSM(model, inits = 0.001, updatefn = ownupdatefn, method = "BFGS")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states (level)
level <- coef(outKFS$model)

#Initial value of level
(initLevel <- coef(outKFS$model)[1])

#Figure 2.1. Deterministic level
plot(dataUKdriversKSI , xlab = "", ylab = "", lty = 1)
lines(level, lty = 3)
title(main = "Figure 2.1. Deterministic level", cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI", "deterministic level"), 
       cex = 0.5, lty = c(1, 3), horiz = T)

#Auxiliary irregular residuals (non-standardised)
irregResid <- residuals(outKFS, "pearson") 

#Figure 2.2. Irregular component for deterministic level model
plot(irregResid  , xlab = "", ylab = "", lty = 2)
abline(h = 0, lty = 1)
title(main = "Figure 2.2. Irregular component for deterministic level model", cex.main = 0.8)
legend("topright",leg = "irregular",cex = 0.5, lty = 2, horiz = T)

#Diagnostic for one-step-ahead prediction residuals (standardised)
predResid <- rstandard(outKFS) 
qStat <- qStatistic(predResid, k, w)
rStat <- rStatistic(predResid, d, l)
hStat <- hStatistic(predResid, d)
nStat <- nStatistic(predResid, d)

#Table 2.1. Diagnostic tests for deterministic level model and log UK drivers KSI
title = "Table 2.1. Diagnostic tests for deterministic level model and log UK drivers \nKSI"
dTable(qStat, rStat, hStat, nStat, title)


#2.2 Stochastic level####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% ts(start = 1969,frequency = 12)

#Defining model
model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree = 1, Q = list(matrix(NA))), H = matrix(NA))

ownupdatefn <- function(pars,model){
  model$H[,, 1] <- exp(pars[1])
  model$Q[,, 1] <- exp(pars[2])
  model
}

d <- q <- 1 #Number of diffuse initial values in the state 
w <- 2#Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15#First k autocorrelations to be used in Q-statistic
n <- 192 #Number of observations
 
#Fitting model and getting ouput
fit <- fitSSM(model, inits = log(c(0.001, 0.001)), updatefn = ownupdatefn, method = "BFGS")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model, marginal = F)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states
level <- coef(outKFS$model)

#Initial value of level
(initLevel <- level[1])

#Figure 2.3. Stochastic level
plot(dataUKdriversKSI , xlab = "", ylab = "", lty = 1)
lines(level, lty = 3)
title(main = "Figure 2.3. Stochastic level", cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI", "stochastic level"), 
       cex = 0.5, lty = c(1, 3), horiz = T)

#Auxiliary irregular residuals (non-standardised)
irregResid <- residuals(outKFS, "pearson") 

#Figure 2.4. Irregular component for local level model
plot(irregResid  , xlab = "", ylab = "", lty = 2)
abline(h = 0, lty = 1)
title(main = "Figure 2.4. Irregular component for local level model", cex.main = 0.8)
legend("topright",leg = "irregular",cex = 0.5, lty = 2, horiz = T)

#Diagnostic for one-step-ahead prediction residuals (standardised)
predResid <- rstandard(outKFS) 
qStat <- qStatistic(predResid, k, w)
rStat <- rStatistic(predResid, d, l)
hStat <- hStatistic(predResid, d)
nStat <- nStatistic(predResid, d)

#Table 2.2. Diagnostic tests for local level model and log UK drivers KSI
title = "Table 2.2. Diagnostic tests for local level model and log UK drivers KSI"
dTable(qStat, rStat, hStat, nStat, title)


#2.3 The local level model and Norwegian fatalities####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataNOfatalities <- log(read.table("NorwayFinland.txt")[,2]) %>% ts(start = 1970, frequency = 1)

#Defining model
model <- SSModel(dataNOfatalities ~ SSMtrend(degree = 1, Q = list(matrix(NA))), H = matrix(NA))
ownupdatefn <- function(pars,model){
  model$H[,,1] <- exp(pars[1])
  model$Q[,,1] <- exp(pars[2])
  model
}

d <- q <- 1 #Number of diffuse initial values in the state 
w <- 2#Number of estimated hyperparameters (i.e. disturbance variances)
l <- 4 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 10#First k autocorrelations to be used in Q-statistic
n <- 34 #Number of observations

#Fitting model and getting output
fit <- fitSSM(model, inits = log(c(0.001,0.001)), updatefn = ownupdatefn, method = "BFGS")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states (level)
level <- coef(outKFS$model)

#Initial value of level
(initLevel <- level[1])

#Figure 2.5. Stochastic level for Norwegian fatalities
plot(dataNOfatalities, xlab = "", ylab = "", lty = 1)
lines(level, lty = 3)
title(main = "Figure 2.5. Stochastic level for Norwegian fatalities", cex.main = 0.8)
legend("topright",leg = c("log fatalities in Norway", "stochastic level"), 
       cex = 0.5, lty = c(1, 3), horiz = T)

#Auxiliary irregular residuals (non-standardised)
irregResid <- residuals(outKFS, "pearson") 

#Figure 2.6. Irregular component for Norwegian fatalities
plot(irregResid  , xlab = "", ylab = "", lty = 2)
abline(h = 0, lty = 1)
title(main = "Figure 2.6. Irregular component for Norwegian fatalities", cex.main = 0.8)
legend("topleft",leg = "irregular",cex = 0.5, lty = 2, horiz = T)

#Diagnostic for one-step-ahead prediction residuals (standardised)
predResid <- rstandard(outKFS) 
qStat <- qStatistic(predResid, k, w)
rStat <- rStatistic(predResid, d, l)
hStat <- hStatistic(predResid, d)
nStat <- nStatistic(predResid, d)

#Table 2.3. Diagnostic tests for local level model and log Norwegian fatalities
title = "Table 2.3. Diagnostic tests for local level model and log Norwegian \nfatalities"
dTable(qStat, rStat, hStat, nStat, title)


#Chapter 3: The local linear trend model####

#3.1 Deterministic level and slope####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% ts(start = 1969, frequency = 12)

#Defining model
model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree = 2, Q = list(matrix(0), matrix(0))), H = matrix(NA))

d <- q <- 2 #Number of diffuse initial values in the state 
w <- 1#Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15#First k autocorrelations to be used in Q-statistic
n <- 192 #Number of observations

#Fitting model and getting output
fit <- fitSSM(model, inits = log(0.001) ,method = "BFGS")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
outKFS$model$Q

#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states
smoothEstStat <- outKFS$alphahat

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

#Diagnostic for one-step-ahead prediction residuals (standardised)
predResid <- rstandard(outKFS) 
qStat <- qStatistic(predResid, k, w)
rStat <- rStatistic(predResid, d, l)
hStat <- hStatistic(predResid, d)
nStat <- nStatistic(predResid, d)

#Table 3.1 Diagnostic tests for deterministic linear trend model and log UK drivers KSI
title = "Table 3.1 Diagnostic tests for deterministic linear trend model and log UK \ndrivers KSI"
dTable(qStat, rStat, hStat, nStat, title)


 

#3.2 Stochastic level and slope####
#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% ts(start = 1969, frequency = 12)

#Defining model
model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree = 2, Q = list(matrix(NA), matrix(NA))),  H = matrix(NA))
ownupdatefn <- function(pars, model){
  model$H[,,1] <- exp(pars[1])
  diag(model$Q[,,1]) <- exp(pars[2:3])
  model
}

d <- q <- 2 #Number of diffuse initial values in the state 
w <- 3#Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15#First k autocorrelations to be used in Q-statistic
n <- 192 #Number of observations

#Fitting model
fit <- fitSSM(model, inits = log(c(0.001, 0001, 0001)) ,method = "BFGS")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states
smoothEstStat <- outKFS$alphahat

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

#Figure 3.1. Trend of stochastic linear trend model
plot(dataUKdriversKSI , xlab = "", ylab = "", lty = 1)
lines(smoothEstStat[, "level"], lty = 3)
title(main = "Figure 3.1. Trend of stochastic linear trend model", cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI", "stochastic level and slope"), 
       cex = 0.5, lty = c(1, 3), horiz = T)

#Figure 3.2. Slope of stochastic linear trend model
plot(smoothEstStat[, "slope"], xlab = "", ylab = "", lty = 1)
title(main = "Figure 3.2. Slope of stochastic linear trend model", cex.main = 0.8)
legend("topleft",leg = "stochastic slope", 
       cex = 0.5, lty = 1, horiz = T)

#Auxiliary irregular residuals (non-standardised)
irregResid <- residuals(outKFS, "pearson") 

#Figure 3.3. Irregular component of stochastic trend model
plot(irregResid  , xlab = "", ylab = "", lty = 2)
abline(h = 0, lty = 1)
title(main = "Figure 3.3. Irregular component of stochastic trend model", cex.main = 0.8)
legend("topright",leg = "irregular",cex = 0.5, lty = 2, horiz = T)

#Diagnostic for one-step-ahead prediction residuals (standardised)
predResid <- rstandard(outKFS) 
qStat <- qStatistic(predResid, k, w)
rStat <- rStatistic(predResid, d, l)
hStat <- hStatistic(predResid, d)
nStat <- nStatistic(predResid, d)

#Table 3.2 Diagnostic tests for the local linear trend model applied 
#to the log of the UK drivers KSI
title = "Table 3.2 Diagnostic tests for the local linear trend model applied to \nthe log of the UK drivers KSI"
dTable(qStat, rStat, hStat, nStat, title)


#3.3 Stochastic level and deterministic slope####
#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% ts(start = 1969, frequency = 12)

#Defining model
model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree = 2, Q = list(matrix(NA), matrix(0))),  H = matrix(NA))
ownupdatefn <- function(pars, model){
  model$H[,,1] <- exp(pars[1])
  diag(model$Q[,,1]) <- c(exp(pars[2]), 0)
  model
}

d <- q <- 2 #Number of diffuse initial values in the state 
w <- 2#Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15#First k autocorrelations to be used in Q-statistic
n <- 192 #Number of observations

#Fitting model and getting output
fit <- fitSSM(model, inits = log(c(0.001, 0.001)), updatefn = ownupdatefn, method = "L-BFGS-B")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states
smoothEstStat <- outKFS$alphahat

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

#Figure 3.4. Trend of stochastic level and deterministic slope model
plot(dataUKdriversKSI, xlab = "", ylab = "", lty = 1)
lines(smoothEstStat[, "level"] , lty = 3)
title(main = "Figure 3.4. Trend of stochastic level and deterministic slope model", cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI", "stochastic level and deterministic slope"), 
       cex = 0.5, lty = c(1, 3), horiz = T)

#Diagnostic for one-step-ahead prediction residuals (standardised)
#predResid <- rstandard(outKFS) 
#qStat <- qStatistic(predResid, k, w)
#rStat <- rStatistic(predResid, d, l)
#hStat <- hStatistic(predResid, d)
#nStat <- nStatistic(predResid, d)
#dTable(qStat, rStat, hStat, nStat, title = "")


#3.4 The local linear trend model and Finnish fatalities####

#A) Both the level and the slope component vary####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataFIfatalities <- log(read.table("NorwayFinland.txt")[,3]) %>% ts(start = 1970, frequency = 1)

#Defining model
model <- SSModel(dataFIfatalities ~ SSMtrend(degree = 2, Q = list(matrix(NA), matrix(NA))),  H = matrix(NA))
ownupdatefn <- function(pars, model){
  model$H[,,1] <- exp(pars[1])
  diag(model$Q[,,1]) <- exp(pars[2:3])
  model
}

d <- q <- 2 #Number of diffuse initial values in the state 
w <- 3#Number of estimated hyperparameters (i.e. disturbance variances)
l <- 4 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 10#First k autocorrelations to be used in Q-statistic
n <- 34 #Number of observations

#Fitting model and getting output
fit <- fitSSM(model, inits = log(c(0.001, 0.001, 0.001)), updatefn = ownupdatefn, method = "L-BFGS-B")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states
smoothEstStat <- outKFS$alphahat

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

#One-step-ahead prediction residuals (standardised)
#predResid <- rstandard(outKFS) 


#B) Deterministic level and stochastic slope####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataFIfatalities <- log(read.table("NorwayFinland.txt")[,3]) %>% ts(start = 1970, frequency = 1)

#Defining model
model <- SSModel(dataFIfatalities ~ SSMtrend(degree = 2, Q = list(matrix(0), matrix(NA))),  H = matrix(NA))
ownupdatefn <- function(pars, model){
  model$H[,,1] <- exp(pars[1])
  diag(model$Q[,,1]) <- c(0, exp(pars[2]))
  model
}

d <- q <- 2 #Number of diffuse initial values in the state 
w <- 2#Number of estimated hyperparameters (i.e. disturbance variances)
l <- 4 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 10#First k autocorrelations to be used in Q-statistic
n <- 34 #Number of observations

#Fitting model and getting output
fit <- fitSSM(model, inits = log(c(0.001, 0.001)), updatefn = ownupdatefn, method = "L-BFGS-B")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states
smoothEstStat <- outKFS$alphahat

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

#Figure 3.5. Trend of deterministic level and stochastic slope model
  # for Finnish fatalities (top) and stochastic slope component (bottom) linear trend model
par(mfrow = c(2, 1), mar = c(1.5, 4, 4, 4))
plot(dataFIfatalities, xlab = "", ylab = "", lty = 1)
lines(smoothEstStat[, "level"], lty = 3)
title(main = "Figure 3.5. Trend of deterministic level and stochastic slope model 
      for Finnish fatalities (top) and stochastic slope component (bottom) linear trend model", 
      cex.main = 0.8)
legend("topright",leg = c("log fatalities Finland", "deterministic level, stochastic slope"), 
       cex = 0.5, lty = c(1, 3), horiz = T)
par(mar = c(4, 4, 1.5, 4))
plot(smoothEstStat[, "slope"], xlab = "", ylab = "", lty = 3)
abline(h = 0, lty = 1)
legend("topright",leg = "stochastic slope", cex = 0.5, lty = 1, horiz = T)
par(mfrow=c(1, 1))

#Auxiliary irregular residuals (non-standardised)
irregResid <- residuals(outKFS, "pearson") 

#Figure 3.6. Irregular component for Finish fatalities
par(mfrow = c(1, 1), mar = c(4, 4, 4, 4))
plot(irregResid  , xlab = "", ylab = "", lty = 2)
abline(h = 0, lty = 1)
title(main = "Figure 3.6. Irregular component for Finish fatalities", cex.main = 0.8)
legend("topright",leg = "irregular",cex = 0.5, lty = 2, horiz = T)

#Diagnostic for one-step-ahead prediction residuals (standardised)
predResid <- rstandard(outKFS) 
qStat <- qStatistic(predResid, k, w)
rStat <- rStatistic(predResid, d, l)
hStat <- hStatistic(predResid, d)
nStat <- nStatistic(predResid, d)

#Table 3.3. Diagnostic tests for deterministic level and stochastic slope
#model, and log Finnish fatalities
title = "Table 3.3. Diagnostic tests for deterministic level and stochastic slope \nmodel, and log Finnish fatalities"
dTable(qStat, rStat, hStat, nStat, title)


#CHAPTER 4: Local level model with seasonal####

#4.1 Deterministic level and seasonal####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% ts(start = 1969, frequency=12)

#Figure 4.1. Log of the number of UK drivers KSI with time lines for years
plot(dataUKdriversKSI, ylab = "",xlab = "")
abline(v=seq(1969, 1985, 1), lty = 3)
title(main = "Figure 4.1. Log of the number of UK drivers KSI with time lines for years", 
      cex.main = 0.8)
legend("topright", leg = "log UK drivers KSI",cex = 0.5, lty = 1, horiz = T)

#Defining model
model <- SSModel(dataUKdriversKSI ~ SSMtrend(1, Q=0) + SSMseasonal(12, sea.type='dummy', Q = 0),  H=NA)

d <- q <- 12 #Number of diffuse initial values in the state 
w <- 1#Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15#First k autocorrelations to be used in Q-statistic
n <- 192 #Number of observations

#Fitting model and getting output
fit <- fitSSM(inits = c(0.001), model = model, method = "BFGS")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
outKFS$model$Q

#Maximum likelihood 
(maxLik <- logLik(fit$model, method = "BFGS")/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states (level and seasonal components)
smoothEstStat <- coef(outKFS$model)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

#Figure 4.2. Combined deterministic level and seasonal
plot(dataUKdriversKSI, xlab = "", ylab = "", lty = 1)
lines(smoothEstStat[, "level"] + smoothEstStat[, "sea_dummy1"], lty = 3)
title(main = "Figure 4.2. Combined deterministic level and seasonal", cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI", "deterministic level + seasonal"), 
       cex = 0.5, lty = c(1, 3), horiz = T)

#Figure 4.3. Deterministic level
plot(dataUKdriversKSI, xlab = "", ylab = "", lty = 1)
lines(smoothEstStat[, "level"], lty = 3)
title(main = "Figure 4.3. Deterministic level", cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI", "deterministic level"), 
       cex = 0.5, lty = c(1, 3), horiz = T)

#Figure 4.4. Deterministic seasonal
plot(smoothEstStat[, "sea_dummy1"], xlab = "", ylab = "", lty = 1)
abline(h = 0, lty = 3)
title(main = "Figure 4.4. Deterministic seasonal", cex.main = 0.8)
legend("topleft",leg = "deterministic seasonal", 
       cex = 0.5, lty = 1, horiz = T)

#Auxiliary irregular residuals (non-standardised)
irregResid <- residuals(outKFS, "pearson") 

#Figure 4.5. Irregular component for deterministic level and seasonal model
plot(irregResid  , xlab = "", ylab = "", lty = 2)
abline(h = 0, lty = 1)
title(main = "Figure 4.5. Irregular component for deterministic level and seasonal model", cex.main = 0.8)
legend("topright",leg = "irregular",cex = 0.5, lty = 2, horiz = T)

#Diagnostic for one-step-ahead prediction residuals (standardised)
predResid <- rstandard(outKFS) 
qStat <- qStatistic(predResid, k, w)
rStat <- rStatistic(predResid, d, l)
hStat <- hStatistic(predResid, d)
nStat <- nStatistic(predResid, d)

#Table 4.1. Diagnostic tests for deterministic level and seasonal 
#model and log UK drivers KSI
title = "Table 4.1. Diagnostic tests for deterministic level and seasonal model \nand log UK drivers KSI"
dTable(qStat, rStat, hStat, nStat, title)


#4.2 Stochastic level and seasonal####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% ts(start = 1969, frequency = 12)

#Defining model
model <- SSModel(dataUKdriversKSI ~ SSMtrend(1, Q = NA) + SSMseasonal(12, sea.type = 'dummy', Q = NA),  H = NA)
ownupdatefn <- function(pars,model){
model$H[,,1] <- exp(pars[1])
diag(model$Q[,,1]) <- exp(c(pars[2], pars[3]))
model
}

d <- q <- 12 #Number of diffuse initial values in the state 
w <- 3 #Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15#First k autocorrelations to be used in Q-statistic
n <- 192 #Number of observations

#Fitting model
#x <- initValOpt() #Finding best initial values for optim
fit <- fitSSM(inits = log(rep(0.887, w)), model = model, updatefn = ownupdatefn, method = "L-BFGS-B")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model, method = method)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states (level and seasonal components)
smoothEstStat <- coef(outKFS$model)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

#Figure 4.6. Stochastic level
plot(dataUKdriversKSI, xlab = "", ylab = "", lty = 1)
lines(smoothEstStat[, "level"], lty = 3)
title(main = "Figure 4.6. Stochastic level", cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI", "stochastic level"), 
       cex = 0.5, lty = c(1, 3), horiz = T)

#Figure 4.7. Stochastic seasonal
plot(smoothEstStat[, "sea_dummy1"], xlab = "", ylab = "", lty = 1)
abline(h = 0, lty = 3)
title(main = "Figure 4.4. Stochastic seasonal", cex.main = 0.8)
legend("topleft",leg = "stochastic seasonal", 
       cex = 0.5, lty = 1, horiz = T)

#Figure 4.8. Stochastic seasonal for the year 1969
plot(window(smoothEstStat[, "sea_dummy1"], start = c( 1969, 1), end = c( 1969, 12)), xlab = "", ylab = "", lty = 1, xaxt = "n")
axis(1, seq(1969, 1970-1/12, length.out = 12), c("1969-Jan", "", "", "1969-Apr", "", "", "1969-July", "", "", "1969-Oct", "", ""))
abline(h = 0, lty = 3)
title(main = "Figure 4.8. Stochastic seasonal for the year 1969", cex.main = 0.8)
legend("topleft",leg = "stochastic seasonal", 
       cex = 0.5, lty = 1, horiz = T)

#Auxiliary irregular residuals (non-standardised)
irregResid <- residuals(outKFS, "pearson") 

#Figure 4.9. Irregular component for stochastic level and seasonal model
plot(irregResid  , xlab = "", ylab = "", lty = 2)
abline(h = 0, lty = 1)
title(main = "Figure 4.9. Irregular component for stochastic level and seasonal model",
      cex.main = 0.8)
legend("topright",leg = "irregular",cex = 0.5, lty = 2, horiz = T)


#Diagnostic for one-step-ahead prediction residuals (standardised)
predResid <- rstandard(outKFS) 
qStat <- qStatistic(predResid, k, w)
rStat <- rStatistic(predResid, d, l)
hStat <- hStatistic(predResid, d)
nStat <- nStatistic(predResid, d)

#Table 4.2. Diagnostic tests for stochastic level and seasonal 
#model and log UK drivers KSI
title = "Table 4.2. Diagnostic tests for stochastic level and seasonal model \nand log UK drivers KSI"
dTable(qStat, rStat, hStat, nStat, title)


#4.3 Stochastic level and deterministic seasonal####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% ts(start = 1969, frequency=12)

#Defining model
model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree=1, Q=list(matrix(NA))) + SSMseasonal(12, sea.type='dummy', Q=matrix(0)),  H=matrix(NA))

ownupdatefn <- function(pars,model,...){
  model$H[,,1] <- exp(pars[1])
  diag(model$Q[,,1]) <- c(exp(pars[2]), 0)
  model
}

d <- q <- 12 #Number of diffuse initial values in the state 
w <- 2 #Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15#First k autocorrelations to be used in Q-statistic
n <- 192 #Number of observations

#Fitting model
x <- initValOpt() #Finding best initial values for optim
fit <- fitSSM(inits = log(rep(0.936, w)), model = model, updatefn = ownupdatefn, method = "Nelder-Mead")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model, method = method)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states (level and seasonal components)
smoothEstStat <- coef(outKFS$model)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

#Auxiliary irregular residuals (non-standardised)
#irregResid <- residuals(outKFS, "pearson") 

#One-step-ahead prediction residuals (standardised)
#predResid <- rstandard(outKFS) 


#4.4 The local level and seasonal model and UK inflation####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKinflation <- read.table("UKinflation.txt") %>% ts(start=1950, frequency=4)

#Fitting model
model <- SSModel(dataUKinflation ~ SSMtrend(degree=1, Q=list(matrix(NA))) + SSMseasonal(period=4, sea.type='dummy', Q=matrix(NA)),  H=matrix(NA))

ownupdatefn <- function(pars,model){
  model$H[,,1] <- exp(pars[1])
  diag(model$Q[,,1]) <- exp(c(pars[2], pars[3]))
  model
}

d <- q <- 4 #Number of diffuse initial values in the state 
w <- 3#Number of estimated hyperparameters (i.e. disturbance variances)
l <- 4 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 10#First k autocorrelations to be used in Q-statisticlogLik <- logLik( ) dlmLL(dataUKdriversKSI, mod)
n <- 208 #Number of observations

#Fitting model
#x <- initValOpt() #Finding best initial values for optim
fit <- fitSSM(inits = log(rep(1.207, w)), model = model, updatefn = ownupdatefn, method = "L-BFGS-B")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model, method = method)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states (level and seasonal components)
smoothEstStat <- coef(outKFS$model)

#Last values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[208,])

#Auxiliary irregular residuals (non-standardised)
irregResid <- residuals(outKFS, "pearson") 

#Figure 4.10 Stochastic level, seasonal and irregular in UK inflation series
par(mfrow = c(3, 1), mar = c(2, 2, 2, 2))
plot(dataUKinflation, xlab = "", ylab = "", lty = 1)
lines(smoothEstStat[, "level"], lty = 3)
title(main = "Figure 4.10 Stochastic level, seasonal and irregular in UK inflation series", 
      cex.main = 1)
legend("topleft",leg = c("quarterly price changes in UK", "stochastic level"), 
       cex = 0.8, lty = c(1, 3), horiz = T)

plot(smoothEstStat[, "sea_dummy1"], xlab = "", ylab = "", lty = 1, ylim = c(-0.006, 0.008))
abline(h = 0, lty = 3)
legend("topleft",leg = "stochastic seasonal", 
       cex = 0.8, lty = 1, horiz = T)

plot(irregResid  , xlab = "", ylab = "", lty = 2)
abline(h = 0, lty = 3)
legend("topleft",leg = "irregular",cex = 0.8, lty = 2, horiz = T)
par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1))

#Diagnostic for one-step-ahead prediction residuals (standardised)
predResid <- rstandard(outKFS) 
qStat <- qStatistic(predResid, k, w)
rStat <- rStatistic(predResid, d, l)
hStat <- hStatistic(predResid, d)
nStat <- nStatistic(predResid, d)

#Table 4.3. Diagnostic tests for local level and seasonal 
#model and UK inflation series
title = "Table 4.3. Diagnostic tests for local level and seasonal model \nand UK inflation series"
dTable(qStat, rStat, hStat, nStat, title)


#CHAPTER 5: The local level model with explanatory variable####

#5.1 Deterministic level and explanatory variable####

#A) Time as explanatory variable####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% ts(start = 1969,frequency=12)
time <- as.numeric(seq(from=1, to=192))

#Defining model
model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree=1, Q=list(matrix(0))) + SSMregression(~ time, Q=matrix(0)), H = matrix(NA))

ownupdatefn <- function(pars,model){
  model$H[,,1] <- exp(pars[1])
  diag(model$Q[,,1]) <- c(0, 0)
  model
}

d <- q <- 2 #Number of diffuse initial values in the state 
w <- 1#Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15#First k autocorrelations to be used in Q-statistic
n <- 192 #Number of observations

#Fitting model
#x <- initValOpt() #Finding best initial values for optim
fit <- fitSSM(model, inits = log(rep(0.286, w)), updatefn = ownupdatefn, method = "L-BFGS-B")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
 
#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states 
smoothEstStat <- coef(outKFS$model)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

#One-step-ahead prediction residuals (standardised)
#predResid <- rstandard(outKFS) 

#Auxiliary irregular residuals (non-standardised)
#irregResid <- residuals(outKFS, "pearson") 


#B) Petrol prices as explanatory variable####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% ts(start = 1969,frequency=12)
petrolPrices <- read.table("logUKpetrolprice.txt")[,1]

#Defining model
model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree=1, Q=list(matrix(0))) + SSMregression(~ petrolPrices, Q=matrix(0)), H = matrix(NA))

ownupdatefn <- function(pars,model){
  model$H[,, 1] <- exp(pars[1])
  diag(model$Q[,, 1]) <- c(0, 0)
  model
}

d <- q <- 2 #Number of diffuse initial values in the state 
w <- 1#Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15#First k autocorrelations to be used in Q-statistic
n <- 192 #Number of observations

#Fitting model
#x <- initValOpt() #Finding best initial values for optim
fit <- fitSSM(model, inits = log(rep(1.1, w)), updatefn = ownupdatefn, method = "L-BFGS-B")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states 
smoothEstStat <- coef(outKFS$model)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1, 1])

#Elasticity: percentage change in the number of UK drivers KSI due to petrol price
#1% increase in the petrol price results in a x% change in the number of drivers KSI
(elast <- smoothEstStat[1, 1] %>% unname)

#Figure 5.1. Deterministic level and explanatory variable 'log petrol price'
plot(dataUKdriversKSI, xlab = "", ylab = "", lty = 1)
lines(smoothEstStat[, "level"] + smoothEstStat[1, "petrolPrices"] * petrolPrices, lty = 3)
title(main = "Figure 5.1. Deterministic level and explanatory variable 'log petrol price'", 
      cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI", "deterministic level + beta*log(PETROL PRICE)"), 
       cex = 0.5, lty = c(1, 3), horiz = T)

#Figure 5.2. Conventional classical regression representation of deterministic level 
#and explanatory variable 'log petrol price'
plot(petrolPrices, dataUKdriversKSI,  xlab = "", ylab = "", pch = 3, cex = 0.5, cex.lab = 0.8, cex.axis = 0.9)
abline(a = smoothEstStat[1, "level"], b = smoothEstStat[1, "petrolPrices"], lty = 3)
title(main = "Figure 5.2. Conventional classical regression representation of deterministic level and explanatory variable 'log petrol price'", 
      cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI against log PETROL PRICE", "deterministic level + beta*log(PETROL PRICE)"), 
       cex = 0.5, lty = c(1, 3), horiz = T)

#Auxiliary irregular residuals (non-standardised)
irregResid <- residuals(outKFS, "pearson") 

#Figure 5.3. Irregular  component for deterministic level model 
#with deterministic explanatory variable 'log petrol price'
plot(irregResid  , xlab = "", ylab = "", lty = 2)
abline(h = 0, lty = 1)
title(main = "Figure 5.3. Irregular  component for deterministic level model with \n deterministic explanatory variable 'log petrol price'",
      cex.main = 0.8)
legend("topleft",leg = "irregular",cex = 0.5, lty = 2, horiz = T)

#One-step-ahead prediction residuals (standardised)
#predResid <- rstandard(outKFS) 

#5.2 Stochastic level and explanatory variable####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% ts(start = 1969,frequency=12)
petrolPrices <- read.table("logUKpetrolprice.txt")[,1]

#Defining model
model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree=1, Q=list(matrix(0))) + SSMregression(~ petrolPrices, Q=matrix(0)), H = matrix(NA))

ownupdatefn <- function(pars,model){
  model$H[,,1] <- exp(pars[1])
  diag(model$Q[,,1]) <- c(0, exp(pars[2]))
  model
}

d <- q <- 2 #Number of diffuse initial values in the state 
w <- 2#Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15#First k autocorrelations to be used in Q-statistic
n <- 192 #Number of observations

#Fitting model
#x <- initValOpt() #Finding best initial values for optim
fit <- fitSSM(model, inits = log(rep(0.971, w)), updatefn = ownupdatefn, method = "L-BFGS-B")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states 
smoothEstStat <- coef(outKFS$model)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

#Elasticity: percentage change in the number of UK drivers KSI due to petrol price
#1% increase in the petrol price results in a x% change in the number of drivers KSI
(elast <- smoothEstStat[1, 1] %>% unname)

#Figure 5.4. Stochastic level and deterministic explanatory variable 'log petrol price'
plot(dataUKdriversKSI, xlab = "", ylab = "", lty = 1)
lines(smoothEstStat[, "level"] + smoothEstStat[1, "petrolPrices"] * petrolPrices, lty = 3)
title(main = "Figure 5.4. Stochastic level and deterministic explanatory variable 'log petrol price'", 
      cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI", "stochastic level + beta*log(PETROL PRICE)"), 
       cex = 0.5, lty = c(1, 3), horiz = T)

#Auxiliary irregular residuals (non-standardised)
irregResid <- residuals(outKFS, "pearson") 

#Figure 5.5. Irregular component for stochastic level model 
#with deterministic explanatory variable 'log petrol price'
plot(irregResid, xlab = "", ylab = "", lty = 2)
abline(h = 0, lty = 1)
title(main = "Figure 5.5. Irregular  component for stochastic level model with \n deterministic explanatory variable 'log petrol price'",
      cex.main = 0.8)
legend("topleft",leg = "irregular",cex = 0.5, lty = 2, horiz = T)

#One-step-ahead prediction residuals (standardised)
#predResid <- rstandard(outKFS) 


#Chapter 6: The local level model with intervention variable

#6.1 Deterministic level and intervention variable####
#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% ts(start = 1969,frequency=12)
seatbeltLaw <- as.numeric(rep(c(0, 1), times=c(169, 23)))  #Intervention variable

#Defining model
model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree = 1, Q = list(matrix(0))) + SSMregression(~ seatbeltLaw, Q = matrix(0), P1 = 100, P1inf = 0), H = matrix(NA))
#Non-exact difuse initialisation is used for the regression coefficient state of the intervention variable. 
#The standardized residuals are not defined for the diffuse phase, which does not end before there is a change in the intervention variable, 
#and thus we only have standardized residual estimates available after the intervention. It is suggested to switch to non-exact diffuse prior 
#for the intervention variable: SSMregression(~ seatbeltLaw, Q=matrix(0), P1 = 100, P1inf = 0)

ownupdatefn <- function(pars,model){
  model$H[,, 1] <- exp(pars[1])
  model
}

d <- 1 #Number of the elements of the initial state vector with exact difuse initialization
q <- 2 #Number of the elements of the initial state vector with  difuse initialization (exact or non-exact)
w <- 1#Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15#First k autocorrelations to be used in Q-statistic
n <- 192 #Number of observations

#Fitting model
#x <- initValOpt() #Finding best initial values for optim
fit <- fitSSM(model, inits = log(rep(0.595, w)), updatefn = ownupdatefn, method = "L-BFGS-B")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states 
smoothEstStat <- coef(outKFS$model)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

#Elasticity: percentage change in the number of UK drivers KSI due to the seat belt law
(elast <- 100*(exp(smoothEstStat[1, 1])-1) %>% unname)

#Figure 6.1. Deterministic level and intervention variable
plot(dataUKdriversKSI, xlab = "", ylab = "", lty = 1)
lines(smoothEstStat[, "level"] + smoothEstStat[1, "seatbeltLaw"] * seatbeltLaw, lty = 3)
title(main = "Figure 6.1. Deterministic level and intervention variable", 
      cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI", "deterministic level + lambda*(SEATBELT LAW)"), 
       cex = 0.5, lty = c(1, 3), horiz = T)

#Figure 5.2. Conventional classical regression representation of deterministic level 
#and intervention variable
plot(seatbeltLaw, dataUKdriversKSI,  xlab = "", ylab = "", pch = 3, cex = 0.5, cex.lab = 0.8, cex.axis = 0.9)
abline(a = smoothEstStat[1, "level"], b = smoothEstStat[1, "seatbeltLaw"], lty = 3)
title(main = "Figure 5.2. Conventional classical regression representation of deterministic level and intervention variable", 
      cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI against log PETROL PRICE", "deterministic level + lambda*(SEATBELT LAW)"), 
       cex = 0.5, lty = c(1, 3), horiz = T)

#Auxiliary irregular residuals (non-standardised)
irregResid <- residuals(outKFS, "pearson") 

#Figure 5.3. Irregular  component for deterministic level model 
#with intervention variable
plot(irregResid  , xlab = "", ylab = "", lty = 2)
abline(h = 0, lty = 1)
title(main = "Figure 5.3. Irregular  component for deterministic level model with intervention variable",
      cex.main = 0.8)
legend("topleft",leg = "irregular",cex = 0.5, lty = 2, horiz = T)

#One-step-ahead prediction residuals (standardised)
#predResid <- rstandard(outKFS) 

#6.2 Stochastic level and intervention variable####
#Removing all objects except functions It does not work!!!
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% ts(start = 1969,frequency=12)
seatbeltLaw <- as.numeric(rep(c(0, 1), times=c(169, 23))) #Intervention variable

#Fitting model
model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree=1, Q=list(matrix(NA))) + SSMregression(~ seatbeltLaw, Q=matrix(0), P1=100, P1inf=0), H = matrix(NA))
#Non-exact difuse initialisation is used for the regression coefficient state of the intervention variable. 
#The standardized residuals are not defined for the diffuse phase, which does not end before there is a change in the intervention variable, 
#and thus we only have standardized residual estimates available after the intervention. It is suggested to switch to non-exact diffuse prior 
#for the intervention variable: SSMregression(~ seatbeltLaw, Q=matrix(0), P1 = 100, P1inf = 0)

ownupdatefn <- function(pars,model){
  model$H[,,1] <- exp(pars[1])
  diag(model$Q[,,1]) <- c(0, exp(pars[2]))
  model
}

d <- 1 #Number of the elements of the initial state vector with exact difuse initialization
q <- 2 #Number of the elements of the initial state vector with  difuse initialization (exact or non-exact)
w <- 2#Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15#First k autocorrelations to be used in Q-statistic
n <- 192 #Number of observations

#Fitting model
#x <- initValOpt() #Finding best initial values for optim
fit <- fitSSM(model, inits = log(rep(0.904, w)), updatefn = ownupdatefn, method = "L-BFGS-B")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states 
smoothEstStat <- coef(outKFS$model)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

#Elasticity: percentage change in the number of UK drivers KSI due to the seat belt law
(elast <- 100*(exp(smoothEstStat[1, 1])-1) %>% unname)

#Figure 6.4. Stochastic level and intervention variable
plot(dataUKdriversKSI, xlab = "", ylab = "", lty = 1)
lines(smoothEstStat[, "level"] + smoothEstStat[1, "seatbeltLaw"] * seatbeltLaw, lty = 3)
title(main = "Figure 6.4. Stochastic level and intervention variable", 
      cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI", "stochastic level + lambda*(SEATBELT LAW)"), 
       cex = 0.5, lty = c(1, 3), horiz = T)

#Auxiliary irregular residuals (non-standardised)
irregResid <- residuals(outKFS, "pearson") 

#Figure 5.3. Irregular  component for stochastic level model 
#with intervention variable
plot(irregResid  , xlab = "", ylab = "", lty = 2)
abline(h = 0, lty = 1)
title(main = "Figure 5.3. Irregular  component for stochastic level model with intervention variable",
      cex.main = 0.8)
legend("topleft",leg = "irregular",cex = 0.5, lty = 2, horiz = T)

#One-step-ahead prediction residuals (standardised)
#predResid <- rstandard(outKFS) 


#CHAPTER 7: The UK seatbelt and inflation models####

#7.1 Deterministic level and seasonal####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% ts(start = 1969, frequency=12)
petrolPrices <- read.table("logUKpetrolprice.txt")[,1] #Explanatory variable
seatbeltLaw <- as.numeric(rep(c(0, 1), times=c(169, 23)))  #Intervention variable

#Defining model
model <- SSModel(dataUKdriversKSI ~ SSMregression(~ petrolPrices, Q = matrix(0)) + SSMregression(~ seatbeltLaw, Q=matrix(0), P1=100, P1inf=0) + SSMtrend(1, Q=0) + SSMseasonal(12, sea.type='dummy', Q = 0),  H=NA)
ownupdatefn <- function(pars,model){
  model$H[,, 1] <- exp(pars[1])
  model
}

d <- 13 #Number of the elements of the initial state vector with exact difuse initialization
q <- 14 #Number of the elements of the initial state vector with  difuse initialization (exact or non-exact)
w <- 1 #Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15 #First k autocorrelations to be used in Q-statisticlogLik <- logLik( ) dlmLL(dataUKdriversKSI, mod)
n <- 192 #Number of observations

#Fitting model
#x <- initValOpt() #Finding best initial values for optim
fit <- fitSSM(model, inits = log(rep(0.303, w)), updatefn = ownupdatefn, method = "L-BFGS-B")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states 
smoothEstStat <- coef(outKFS$model)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

#Elasticity: percentage change in the number of UK drivers KSI due to the seat belt law
(elastSeatbeltLaw <- 100*(exp(smoothEstStat[1, "seatbeltLaw"])-1) %>% unname)

#Elasticity: percentage change in the number of UK drivers KSI due to petrol price
#1% increase in the petrol price results in a x% change in the number of drivers KSI
(elastPetrolPrices <- smoothEstStat[1, "petrolPrices"] %>% unname)

#Figure 7.1. Deterministic level plus variable log petrol price and seat belt law
plot(dataUKdriversKSI, xlab = "", ylab = "", lty = 1)
lines(smoothEstStat[, "level"] + smoothEstStat[1, "petrolPrices"] * petrolPrices + 
        smoothEstStat[1, "seatbeltLaw"] * seatbeltLaw, lty = 3)
title(main = "Figure 7.1. Deterministic level plus variable log petrol price and seat belt law", 
      cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI", "deterministic level + beta*log(PETROL PRICE) + lambda*(SEATBELT LAW)"), 
       cex = 0.5, lty = c(1, 3), horiz = T)

#Diagnostic for one-step-ahead prediction residuals (standardised)
predResid <- rstandard(outKFS) 
qStat <- qStatistic(predResid, k, w)
rStat <- rStatistic(predResid, d, l)
hStat <- hStatistic(predResid, d)
nStat <- nStatistic(predResid, d)

#Table 7.1. Diagnostic tests for deterministic model applied to the UK drivers KSI series
title = "Table 7.1. Diagnostic tests for deterministic model applied to the UK drivers KSI series"
dTable(qStat, rStat, hStat, nStat, title)

#Auxiliary irregular residuals (non-standardised)
irregResid <- residuals(outKFS, "pearson") 

#Figure 7.5. Correlogram of irregular component of completely deterministic level and seasonal model (for point 7.3)
Acf(as.numeric(irregResid), 15, main = "", ylab = "")
title(main="Figure 7.5. Correlogram of irregular component of completely deterministic level and seasonal model", 
      cex.main=0.8)
legend("topright",leg = "ACF - deterministic level and seasonal model residuals",cex = 0.5,lty = 1, col = "black",horiz = T)

#Regression estimates with standard errors to calculate t-ratio (for point 7.3)
outKFS 


#7.2 Stochastic level and seasonal####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% ts(start = 1969, frequency=12)
petrolPrices <- read.table("logUKpetrolprice.txt")[,1] #Explanatory variable
seatbeltLaw <- as.numeric(rep(c(0, 1), times=c(169, 23)))  #Intervention variable

#Defining model
model <- SSModel(dataUKdriversKSI ~ SSMregression(~ petrolPrices, Q = matrix(0)) + SSMregression(~ seatbeltLaw, Q = matrix(0), P1=100, P1inf = 0) + SSMtrend(1, Q = NA) + SSMseasonal(12, sea.type = 'dummy', Q = NA),  H = NA)
ownupdatefn <- function(pars, model){
  model$H[,, 1] <- exp(pars[1])
  diag(model$Q[,, 1]) <- c(0, 0, exp(pars[2]), exp(pars[3]))
  model
}

d <- 13 #Number of the elements of the initial state vector with exact difuse initialization
q <- 14 #Number of the elements of the initial state vector with  difuse initialization (exact or non-exact)
w <- 3#Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15#First k autocorrelations to be used in Q-statistic
n <- 192 #Number of observations

#Fitting model
#x <- initValOpt() #Finding best initial values for optim
fit <- fitSSM(model, inits = log(rep(0.846, w)), updatefn = ownupdatefn, method = "L-BFGS-B")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states 
smoothEstStat <- coef(outKFS$model)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

#Elasticity: percentage change in the number of UK drivers KSI due to the seat belt law
(elastSeatbeltLaw <- 100*(exp(smoothEstStat[1, "seatbeltLaw"])-1) %>% unname)

#Elasticity: percentage change in the number of UK drivers KSI due to petrol price
#1% increase in the petrol price results in a x% change in the number of drivers KSI
(elastPetrolPrices <- smoothEstStat[1, "petrolPrices"] %>% unname)

#Figure 7.2. Stochastic level plus variable log petrol price and seat belt law
plot(dataUKdriversKSI, xlab = "", ylab = "", lty = 1)
lines(smoothEstStat[, "level"] + smoothEstStat[1, "petrolPrices"] * petrolPrices + 
        smoothEstStat[1, "seatbeltLaw"] * seatbeltLaw, lty = 3)
title(main = "Figure 7.1. Stochastic level plus variable log petrol price and seat belt law", 
      cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI", "stochastic level + beta*log(PETROL PRICE) + lambda*(SEATBELT LAW)"), 
       cex = 0.5, lty = c(1, 3), horiz = T)

#Figure 7.3. Stochastic seasonal
plot(smoothEstStat[, "sea_dummy1"], xlab = "", ylab = "", lty = 1)
abline(h = 0, lty = 3)
title(main = "Figure 4.4. Stochastic seasonal", cex.main = 0.8)
legend("topleft",leg = "stochastic seasonal", 
       cex = 0.5, lty = 1, horiz = T)

#Auxiliary irregular residuals (non-standardised)
irregResid <- residuals(outKFS, "pearson") 

#Figure 7.4. Irregular component for stochastic level and seasonal model
plot(irregResid  , xlab = "", ylab = "", lty = 2)
abline(h = 0, lty = 1)
title(main = "Figure 7.4. Irregular component for stochastic level and seasonal model",
      cex.main = 0.8)
legend("topright",leg = "irregular",cex = 0.5, lty = 2, horiz = T)

#Diagnostic for one-step-ahead prediction residuals (standardised)
predResid <- rstandard(outKFS) 
qStat <- qStatistic(predResid, k, w)
rStat <- rStatistic(predResid, d, l)
hStat <- hStatistic(predResid, d)
nStat <- nStatistic(predResid, d)

#Table 7.2. Diagnostic tests for the stochastic model applied to the UK drivers KSI series
title = "Table 7.2. Diagnostic tests for the stochastic model applied to the UK drivers KSI series"
dTable(qStat, rStat, hStat, nStat, title)

#7.3 Stochastic level and deterministic seasonal####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% ts(start = 1969, frequency=12)
petrolPrices <- read.table("logUKpetrolprice.txt")[,1] #Explanatory variable
seatbeltLaw <- as.numeric(rep(c(0, 1), times=c(169, 23)))  #Intervention variable

#Defining model
model <- SSModel(dataUKdriversKSI ~ SSMregression(~ petrolPrices, Q = matrix(0)) + SSMregression(~ seatbeltLaw, Q=matrix(0), P1 = 100, P1inf = 0) + SSMtrend(1, Q = NA) + SSMseasonal(12, sea.type ='dummy', Q = 0),  H = NA)
ownupdatefn <- function(pars, model){
  model$H[,, 1] <- exp(pars[1])
  diag(model$Q[,, 1]) <- c(0, 0, exp(pars[2]), 0)
  model
}

d <- 13 #Number of the elements of the initial state vector with exact difuse initialization
q <- 14 #Number of the elements of the initial state vector with  difuse initialization (exact or non-exact)
w <- 2#Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15#First k autocorrelations to be used in Q-statisticlogLik <- logLik( ) dlmLL(dataUKdriversKSI, mod)
n <- 192 #Number of observations

#Fitting model
#x <- initValOpt() #Finding best initial values for optim
fit <- fitSSM(model, inits = log(rep(1.124, w)), updatefn = ownupdatefn, method = "L-BFGS-B")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states 
smoothEstStat <- coef(outKFS$model)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

#Elasticity: percentage change in the number of UK drivers KSI due to the seat belt law
(elastSeatbeltLaw <- 100*(exp(smoothEstStat[1, "seatbeltLaw"])-1) %>% unname)

#Elasticity: percentage change in the number of UK drivers KSI due to petrol price
#1% increase in the petrol price results in a x% change in the number of drivers KSI
(elastPetrolPrices <- smoothEstStat[1, "petrolPrices"] %>% unname)

#Auxiliary irregular residuals (non-standardised)
irregResid <- residuals(outKFS, "pearson") 

#Figure 7.6. Correlogram of irregular component of stochastic level and deterministic seasonal model
Acf(as.numeric(irregResid), 15, main = "", ylab = "")
title(main="Figure 7.6. Correlogram of irregular component of stochastic level and deterministic seasonal model", 
      cex.main=0.8)
legend("topright",leg = "stochastic level and deterministic seasonal model residuals",cex = 0.5,lty = 1, col = "black",horiz = T)

#Figure 7.5. Correlogram of irregular component of completely deterministic level and seasonal model 
#(see in point 7.1)

#One-step-ahead prediction residuals (standardised)
#predResid <- rstandard(outKFS) 

#Regression estimates with standard errors to calculate t-ratio 
outKFS 
#(see point 7.1 for deterministic model)


#7.4 The UK inflation model####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKinflation <- read.table("UKinflation.txt") %>% ts(start=1950, frequency=4)

n <- 208 #Number of observations
pulse1 <- ts(rep(0, n), frequency=4, start=c(1950, 1)) # Pulse intervention variable 1
pulse1[which(time(pulse1) == 1975.25)] <- 1  
pulse2 <- ts(rep(0, n), frequency=4, start=c(1950, 1)) # Pulse intervention variable 2
pulse2[which(time(pulse2) == 1979.50)] <- 1  

#Defining model
model <- SSModel(dataUKinflation ~ SSMregression(~ pulse1, Q = matrix(0), P1 = 100, P1inf = 0) + SSMregression(~ pulse2, Q = matrix(0), P1 = 100, P1inf = 0) + SSMtrend(degree = 1, Q = list(matrix(NA))) + SSMseasonal(period= 4, sea.type = 'dummy', Q = matrix(NA)),  H = matrix(NA))

ownupdatefn <- function(pars, model){
  model$H[,, 1] <- exp(pars[1])
  diag(model$Q[,, 1]) <- c(0, 0, exp(pars[2]), exp(pars[3]))
  model
}

d <- 4 # #Number of exact diffuse initial values in the state
q <- 6 #Number of diffuse initial values (exact and non-exact) in the state 
w <- 3 #Number of estimated hyperparameters (i.e. disturbance variances)
l <- 4 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 10 #First k autocorrelations to be used in Q-statisticlogLik <- logLik( ) dlmLL(dataUKdriversKSI, mod)

#Fitting model
#x <- initValOpt() #Finding best initial values for optim
x <- initValOpt(method = "Nelder-Mead") 

#fit <- fitSSM(inits = rep(0.273, w), model = model, method = "L-BFGS-B") 
#the Nelder-Mead algorithm that can be more robust than BFGS although it may converge more slowly.
fit <- fitSSM(inits = rep(1.442, w), model = model, method = "Nelder-Mead") 

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states 
smoothEstStat <- coef(outKFS$model)

#Initial values of the smoothed estimates of states
#(initSmoothEstStat <- smoothEstStat[1,])

#Auxiliary irregular residuals (non-standardised)
irregResid <- residuals(outKFS, "pearson") 

#Figure 7.7 Local level (including pulse interventions), local seasonal and irregular in UK inflation time series data
par(mfrow = c(3, 1), mar = c(2, 2, 2, 2))
plot(dataUKinflation, xlab = "", ylab = "", lty = 1)
lines(smoothEstStat[, "level"], lty = 3)
title(main = "Figure 7.7. Local level (including pulse interventions), local seasonal and irregular in UK inflation time series data", 
      cex.main = 1)
legend("topleft",leg = c("quarterly price changes in UK", "stochastic level + pulse intervention variables"), 
       cex = 0.8, lty = c(1, 3), horiz = T)

plot(smoothEstStat[, "sea_dummy1"], xlab = "", ylab = "", lty = 1, ylim = c(-0.006, 0.008))
abline(h = 0, lty = 3)
legend("topleft",leg = "stochastic seasonal", 
       cex = 0.8, lty = 1, horiz = T)

plot(irregResid  , xlab = "", ylab = "", lty = 2)
abline(h = 0, lty = 3)
legend("topleft",leg = "irregular",cex = 0.8, lty = 2, horiz = T)
par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1))

#Diagnostic for one-step-ahead prediction residuals (standardised)
predResid <- rstandard(outKFS) 
qStat <- qStatistic(predResid, k, w)
rStat <- rStatistic(predResid, d, l)
hStat <- hStatistic(predResid, d)
nStat <- nStatistic(predResid, d)

#Table 7.3. Diagnostic tests for the local level and seasonal 
#model incuding pulse intervention variables for UK inflation series
title = "Table 4.3. Diagnostic tests for local level and seasonal model \nincuding pulse intervention variables for UK inflation series"
dTable(qStat, rStat, hStat, nStat, title)


#8. General treatment of univariate state space models####

#8.3 Confidence intervals####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% ts(start = 1969, frequency=12)

#Defining model
model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree=1, Q=list(matrix(NA))) + SSMseasonal(12, sea.type='dummy', Q=matrix(0)),  H=matrix(NA))

ownupdatefn <- function(pars,model,...){
  model$H[,,1] <- exp(pars[1])
  diag(model$Q[,,1]) <- c(exp(pars[2]), 0)
  model
}

#d <- q <- 12 #Number of diffuse initial values in the state 
w <- 2 #Number of estimated hyperparameters (i.e. disturbance variances)
#l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
#k <- 15#First k autocorrelations to be used in Q-statistic
n <- 192 #Number of observations

#Fitting model
#x <- initValOpt()
fit <- fitSSM(inits = log(rep(0.009, w)), model = model, updatefn = ownupdatefn, method = "Nelder-Mead")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))



#Maximum likelihood 
#(maxLik <- logLik(fit$model, method = method)/n)

#Akaike information criterion (AIC)
#(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
#(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
#(Q <- fit$model$Q)

#Smoothed estimates of states (level and seasonal components)
smoothEstStat <- coef(outKFS$model)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

#Auxiliary irregular residuals (non-standardised)
#irregResid <- residuals(outKFS, "pearson") 

#One-step-ahead prediction residuals (standardised)
#predResid <- rstandard(outKFS) 


#Interesting! to be reviewed?
outKFS$d #The last time index of diffuse phase, i.e.  the non-diffuse phase began at timed+ 1. 

#Level estimation error variance
levEstErVar <- outKFS$V[1, 1, ] %>% ts(start = 1969, frequency=12)
outKFS$V[1, 1, ]

#Figure 8.1. Level estimation error variance for stochastic level and deterministic seasonal model applied to the log of UK drivers KSI
plot(levEstErVar, xlab = "", ylab = "", lty = 1)
title(main = "Figure 8.1. Level estimation error variance for stochastic level and deterministic seasonal model \n applied to the log of UK drivers KSI", 
      cex.main = 0.8)
legend("topright",leg = "level estimation error variance", 
       cex = 0.5, lty = 1, horiz = T)


#Predict out (naming to be reviewed)
outPredictLev <- predict(fit$model, states = "level", interval = "confidence", level = 0.90)
outPredictSeas <- predict(fit$model, states = "seasonal", interval = "confidence", level = 0.90)
outPredictSig <- predict(fit$model, states = "all", interval = "confidence", level = 0.90)

#Figure 8.2. Stochastic level and its 90% confidence interval for stochastic level and deterministic seasonal model applied to the log of UK drivers KSI
plot(dataUKdriversKSI, xlab = "", ylab = "", lty = 1)
lines(pred[, "fit"], lty = 3)
lines(pred[, "lwr"], lty = 3)
lines(pred[, "upr"], lty = 3)
title(main = "Figure 8.2. Stochastic level and its 90% confidence interval for stochastic level \nand deterministic seasonal model applied to the log of UK drivers KSI", 
      cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI", "stochastic level +/- 1.64SE"), 
       cex = 0.5, lty = c(1, 3), horiz = T)

#Figure 8.3. Deterministic seasonal and its 90% confidence interval for stochastic level and deterministic seasonal model applied to the log of UK drivers KSI
plot(window(outPredictSeas[, "fit"], start = 1981), xlab = "", ylab = "", lty = 3, ylim = c(-0.2, 0.3))
lines(window(outPredictSeas[, "lwr"], start = 1981), lty = 3)
lines(window(outPredictSeas[, "upr"], start = 1981), lty = 3)
title(main = "Figure 8.3. Deterministic seasonal and its 90% confidence interval for stochastic level \nand deterministic seasonal model applied to the log of UK drivers KSI", 
      cex.main = 0.8)
legend("topright",leg = "deterministic seasonal +/- 1.64SE", 
       cex = 0.5, lty = 3, horiz = T)

#Figure 8.4. Stochastic level plus deterministic seasonal and its 90% confidence interval for stochastic level and deterministic seasonal model applied to the log of UK drivers KSI
plot(window(outPredictSeas[, "fit"], start = 1981), xlab = "", ylab = "", lty = 1, ylim = c(-0.2, 0.3))
lines(window(outPredictSeas[, "lwr"], start = 1981), lty = 1)
lines(window(outPredictSeas[, "upr"], start = 1981), lty = 1)
title(main = "Figure 8.4. Stochastic level plus deterministic seasonal and its 90% confidence interval for stochastic level \nand deterministic seasonal model applied to the log of UK drivers KSI", 
      cex.main = 0.8)
legend("topright",leg = "signal +/- 1.64SE", 
       cex = 0.5, lty = 1, horiz = T)



#9. Multivariate time series analysis####

#9.4 An illustration of multivariate state space analysis####

#A) Model without rank restriction####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
seatbeltLaw <- rep(c(0, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
data <- read.table("UKfrontrearseatKSI.txt", header = TRUE)[-1] %>% 
  cbind(., seatbeltLaw) %>%
  ts(start = 1969,frequency=12) 

#Defining model
model <- SSModel(cbind(frontseatKSI, rearseatKSI) ~ kmDriven + petrolPrice + 
                   seatbeltLaw +
                   SSMtrend(degree = 1, Q = matrix(NA, 2, 2)) +
                   SSMseasonal(period = 12, sea.type = 'dummy'), H = matrix(NA, 2, 2), data = data)

#SSMregression(~ seatbeltLaw, Q = matrix(0), P1 = 100, P1inf = 0)

ownupdatefn <- function(pars, model){
  L1T <- diag(exp(pars[1:2])) # see explanationhttps://mc-stan.org/docs/2_18/reference-manual/covariance-matrices-1.html#fn16
  L1T[upper.tri(L1T)] <- pars[3]
  model["Q"] <- crossprod(L1T) #crossprod (X,Y) = t(X) %*% Y or crossprod (X) = t(X) %*% X
  L2T <- diag(exp(pars[4:5]))
  L2T[upper.tri(L2T)] <- pars[6]
  model["H"] <- crossprod(L2T)
  model
}

#d <- ?? # #Number of exact diffuse initial values in the state
#q <- ?? #Number of diffuse initial values (exact and non-exact) in the state 
w <- 6 #Number of estimated hyperparameters (i.e. disturbance variances)
#l <- ?? #Autocorrelation at lag l to be provided by r-statistic / ACF function
#k <- ?? #First k autocorrelations to be used in Q-statisticlogLik <- logLik( ) dlmLL(dataUKdriversKSI, mod)
n <- 2*192 #Number of observations

#Fitting model
#x <- initValOpt2(method = "Nelder-Mead", formula = "c(log(x), log(x), x, log(x), log(x), x)")
x <- 0.148
fit <- fitSSM(inits = c(log(x), log(x), x, log(x), log(x), x), model = model, updatefn = ownupdatefn, method = "Nelder-Mead")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
#All exact diffused explanatory variablea and intervention variable give good resutls
# Please not that difussion ends after the intervention variable setbeltLaw changes from 0 to 1
# This is gives very short time series of predictive residuals

#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
#(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)
 
#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states 
smoothEstStat <- coef(outKFS$model)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])


#predResid <- rstandard(outKFS) #One-step-ahead prediction residuals (standardised)
#irregResid <- rstandard(outKFS, "pearson") #Auxiliary irregular  residuals (standardised)

#Figure 9.1 Log of monthly numbers of front seat passengers (top) and rear seat passangers (bottom)
#killed or seriously injured in the UK in the period 1969-1984
plot(data[, "frontseatKSI"], xlab= "", ylab = "", lty=1, ylim=c(5.4, 7.2))
lines(data[, "rearseatKSI"], lty=5)
title(main = "Figure 9.1. Log of monthly numbers of fron seat passengers (top) and rear seat passengers (bottom)\n killed or seriously injured in the UK in the period 1969-1984", 
      cex.main = 0.8)
legend("topright",leg = c("log(front seat KSI)", "log(rear seat KSI)"), 
       lty = c(1, 5), cex = 0.55, horiz=T)
 
#Auxiliary level  residuals (standardised)
levResid <- residuals(outKFS, "state") 

#Figure 9.2 Level disturbances for rear seat (horizontal) versus front seat KSI (vertical)
#in a seamingly unrelated model
#Fit model before

plot(x=levResid[,"rearseatKSI"], y=levResid[,"frontseatKSI"], xlab = "Level disturbances for rear seat KSI", ylab = "Level disturbance for front seat KSI",
     pch = 3, cex = 0.5, cex.lab = 0.8, cex.axis = 0.9)
abline(lm(levResid[,"frontseatKSI"] ~ levResid[,"rearseatKSI"]))
abline(h = 0, lty = 3)
title(main = "Figure 9.2 Level disturbances for rear seat (horizontal) versus front seat KSI (vertical) \nin a seamingly unrelated model", 
      cex.main = 0.8) 
legend("topleft",leg = c("Level disturbances rear against front seat KSI", "Regression line"), 
       lty = c(NA, 1), pch=c(3, NA), cex = 0.55, horiz=T)

#Figure 9.3 Levels of treatment and control series in the seemingly unrelated model
#(levels of rear and front seat KSI)

par(mfrow=c(2,1), mar=c(2.1, 4.1, 2.1, 2.1))
plot(smoothEstStat[, "level.frontseatKSI"], xlab= "", ylab = "", lty=1)
title(main = "Figure 9.3 Levels of treatment and control series in the seemingly unrelated model", 
      cex.main = 0.8)
legend("topright",leg = "level front", lty = 1, cex = 0.55, horiz=T)
plot(smoothEstStat[, "level.rearseatKSI"], xlab= "", ylab = "",lty=1)
legend("topright",leg = "level rear", lty = 1, cex = 0.55, horiz=T)
par(mfrow=c(1,1), par(mar=c(5.1, 4.1, 4.1, 2.1)))

#Figure 9.4 Level of treatment against level of control series in the seemingly unrelated model
plot(x=smoothEstStat[, "level.rearseatKSI"], y=smoothEstStat[, "level.frontseatKSI"], xlab = "Level for rear seat KSI", ylab = "Level for front seat KSI",
     pch = 3, cex = 0.5, cex.lab = 0.8, cex.axis = 0.9)
abline(lm(smoothEstStat[, "level.frontseatKSI"] ~ smoothEstStat[, "level.rearseatKSI"]))
title(main = "Figure 9.4 Level of treatment against level of control series\nin the seemingly unrelated model", 
      cex.main = 0.8)
legend("topleft",leg = c("Level  rear against front seat KSI", "Regression line"), 
       lty = c(NA, 1), pch=c(3, NA), cex = 0.55, horiz=T)


#B) Model with rank restriction####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
seatbeltLaw <- rep(c(0, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
data <- read.table("UKfrontrearseatKSI.txt", header = TRUE)[-1] %>% 
  cbind(., seatbeltLaw) %>%
  ts(start = 1969,frequency=12) 

#Defining model
model <- SSModel(cbind(frontseatKSI, rearseatKSI) ~ kmDriven + petrolPrice + 
                   SSMregression(~ seatbeltLaw, index = 1) +
                   SSMtrend(degree = 1, Q = matrix(NA, 2, 2)) +
                   SSMseasonal(period = 12, sea.type = 'dummy'), H = matrix(NA, 2, 2), data = data)

#SSMregression(~ seatbeltLaw, Q = matrix(0), P1 = 100, P1inf = 0)

ownupdatefn <- function(pars, model){
  L1T <- diag(c(exp(pars[1]), 0)) # see explanationhttps://mc-stan.org/docs/2_18/reference-manual/covariance-matrices-1.html#fn16
  L1T[upper.tri(L1T)] <- pars[2]
  model["Q"] <- crossprod(L1T) #crossprod (X,Y) = t(X) %*% Y or crossprod (X) = t(X) %*% X
  L2T <- diag(exp(pars[3:4]))
  L2T[upper.tri(L2T)] <- pars[5]
  model["H"] <- crossprod(L2T)
  model
}

#d <- ?? # #Number of exact diffuse initial values in the state
#q <- ?? #Number of diffuse initial values (exact and non-exact) in the state 
w <- 5 #Number of estimated hyperparameters (i.e. disturbance variances)
#l <- ?? #Autocorrelation at lag l to be provided by r-statistic / ACF function
#k <- ?? #First k autocorrelations to be used in Q-statisticlogLik <- logLik( ) dlmLL(dataUKdriversKSI, mod)
n <- 2*192 #Number of observations

#Fitting model
#x <- initValOpt2(method = "Nelder-Mead", formula = "c(log(x), log(x), x, log(x), log(x), x)")
x <- 0.360
fit <- fitSSM(inits = c(log(x), log(x), x, log(x), log(x), x), model = model, updatefn = ownupdatefn, method = "Nelder-Mead")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
#All exact diffused explanatory variablea and intervention variable give good resutls
# Please not that difussion ends after the intervention variable setbeltLaw changes from 0 to 1
# This is gives very short time series of predictive residuals

#Maximum likelihood 
(maxLik <- logLik(fit$model)/n)

#Akaike information criterion (AIC)
#(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states 
smoothEstStat <- coef(outKFS$model)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])


#predResid <- rstandard(outKFS) #One-step-ahead prediction residuals (standardised)
#irregResid <- rstandard(outKFS, "pearson") #Auxiliary irregular  residuals (standardised)

#Auxiliary level  residuals (standardised)
levResid <- residuals(outKFS, "state") 

#Figure 9.5 Level disturbances for rear seat (horizontal) versus front seat KSI (vertical)
#in a seamingly unrelated model
plot(x=levResid[,"rearseatKSI"], y=levResid[,"frontseatKSI"], xlab = "Level disturbances for rear seat KSI", ylab = "Level disturbance for front seat KSI",
     pch = 3, cex = 0.5, cex.lab = 0.8, cex.axis = 0.9)
abline(lm(levResid[,"frontseatKSI"] ~ levResid[,"rearseatKSI"]))
abline(h = 0, lty = 3)
title(main = "Figure 9.5 Level disturbances for rear seat (horizontal) versus front seat KSI (vertical) \nin a seamingly unrelated model", 
      cex.main = 0.8)
legend("topleft",leg = c("Level disturbances rear against front seat KSI", "Regression line"), 
       lty = c(NA, 1), pch=c(3, NA), cex = 0.55, horiz=T)

#Figure 9.6 Level of treatment against level of control series in rank one model
plot(x = smoothEstStat[, "level.rearseatKSI"], y = smoothEstStat[, "level.frontseatKSI"], xlab = "Level for rear seat KSI", ylab = "Level for front seat KSI",
     pch = 3, cex = 0.5, cex.lab = 0.8, cex.axis = 0.9)
abline(lm(smoothEstStat[, "level.frontseatKSI"] ~ smoothEstStat[, "level.rearseatKSI"]))
title(main = "Figure 9.6 Level of treatment against level of control series\nin rank one model", 
      cex.main = 0.8)
legend("topleft",leg = c("Level  rear against front seat KSI", "Regression line"), 
       lty = c(NA, 1), pch=c(3, NA), cex = 0.55, horiz=T)

#Figure 9.7 Levels of treatment and control series in rank one model
par(mfrow=c(2,1), mar=c(2.1, 4.1, 2.1, 2.1))
plot(smoothEstStat[, "level.frontseatKSI"], xlab= "", ylab = "", lty=1)
title(main = "Figure 9.7 Levels of treatment and control series in rank one model", 
      cex.main = 0.8)
legend("topright",leg = "level front", lty = 1, cex = 0.55, horiz=T)
plot(smoothEstStat[, "level.rearseatKSI"], xlab= "", ylab = "",lty=1)
legend("topright",leg = "level rear", lty = 1, cex = 0.55, horiz=T)
par(mfrow=c(1,1), par(mar=c(5.1, 4.1, 4.1, 2.1)))

#Figure 9.8 Level of treatment series plus intervention and level of control series, rank one model
par(mfrow=c(2,1), mar=c(2.1, 4.1, 2.1, 2.1))
plot(smoothEstStat[, "level.frontseatKSI"] + data[, "seatbeltLaw"]*smoothEstStat[, "seatbeltLaw.frontseatKSI"], xlab= "", ylab = "", lty=1)
title(main = "Figure 9.8 Level of treatment series plus intervention and level of control series, \nrank one model", 
      cex.main = 0.8)
legend("topright",leg = "level front", lty = 1, cex = 0.55, horiz=T)
plot(smoothEstStat[, "level.rearseatKSI"], xlab= "", ylab = "",lty=1)
legend("topright",leg = "level rear", lty = 1, cex = 0.55, horiz=T)
par(mfrow=c(1,1), par(mar=c(5.1, 4.1, 4.1, 2.1)))

#Figure 9.9 Deterministic seasonal of treatment and control series, rank one model
par(mfrow=c(2,1), mar=c(2.1, 4.1, 2.1, 2.1))
plot(smoothEstStat[, "sea_dummy1.frontseatKSI"], xlab= "", ylab = "", lty=1)
title(main = "Figure 9.9 Deterministic seasonal of treatment and control series, rank one model", 
      cex.main = 0.8)
legend("topright",leg = "seasonal front", lty = 1, cex = 0.55, horiz=T)
plot(smoothEstStat[, "sea_dummy1.rearseatKSI"], xlab= "", ylab = "",lty=1)
legend("topright",leg = "seasonal rear", lty = 1, cex = 0.55, horiz=T)
par(mfrow=c(1,1), par(mar=c(5.1, 4.1, 4.1, 2.1)))


#BOOK END####  ########################################################



#Version 00 ####
rm(list = setdiff(ls(), lsf.str())) 

UKfrontrearseatKSI <- read.table("UKfrontrearseatKSI.txt", head=TRUE)  
data <- ts(UKfrontrearseatKSI[,2:3], start = 1969,frequency=12) #Two time series: front seat passengers killed or seriously injured and rear seat passengers killed or seriously injured

seatbeltLaw <- rep(c(0, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
X <- cbind(seatbeltLaw, UKfrontrearseatKSI[, c(5,4)])
colnames(X)[2:3] <- c("petrolPrice", "kmDriven" )
X <- ts(X, start = 1969,frequency=12)

model <- SSModel(data ~ kmDriven + petrolPrice + SSMregression (~ seatbeltLaw, index = 1) + 
                   SSMtrend(degree = 1, Q=matrix(NA, 2, 2)) +
                   SSMseasonal(period = 12, sea.type ='dummy'), H = matrix(NA, 2, 2), data = X)


ownupdatefn <- function(pars, model){
  H <- diag(exp(pars[1:2]))
  H[upper.tri(H)] <- pars[3]
  model["H"] <- crossprod(H)
  Q <- matrix(c(exp(pars[4]), pars[5]), 1, 2) 
  model["Q"] <- crossprod(Q) #crossprod (X,Y) = t(X) %*% Y or crossprod (X) = t(X) %*% X
  model
}

set.seed(123)
for (j in 1:50)
{
  #x <- runif(6, 0.0000001, 0.5)
  #fit <- fitSSM(inits = x, model = model, updatefn = ownupdatefn, method = "BFGS")
  x <- round(runif(1, 0.0001, 0.2),3)
  fit <- fitSSM(inits = c(log(x), log(x), x, log(x), x), model = model, updatefn = ownupdatefn, method = "BFGS")
  
  outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
  print(x)
  print(logLik(fit$model))
  
  print(fit$model$Q)
  levelResid <- residuals(outKFS, "state") #Auxiliary level  residuals (standardised)
  plot(levelResid[, 1], levelResid[, 2])
}

init = 0.028
fit <- fitSSM(inits = c(log(init), log(init), init, log(init), init), model = model, updatefn = ownupdatefn, method = "BFGS") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

(H <- fit$model$H)  

(Q <- fit$model$Q)

logLik(fit$model)

#Figure 9.5 Level disturbances for rear (horizontal) versus front seat KSI (vertical)
#rank one model

#Fit model before
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
levelResid <- residuals(outKFS, "state") #Auxiliary level  residuals (standardised)
colnames(levelResid) <- c("Front.seats", "Rear.seats")
plot(x=levelResid[,"Rear.seats"], y=levelResid[,"Front.seats"], xlab = "Level disturbances for rear seat KSI", ylab = "Level disturbance for front seat KSI")

abline(abline(lm(levelResid[,"Front.seats"] ~ levelResid[,"Rear.seats"])))

legend("topleft",leg = c("Level disturbances rear against front seat KSI", "Regression line"), 
       lty = c(NA, 1), pch=c(1, NA), cex = 0.55, horiz=T)


#Figure 9.6 Level of control (rear seats) against level of treatment series (front seat) in rank one model
plot(x=levels[, "level.Rear.seats.KSI"], y=levels[, "level.Front.seats.KSI"], xlab = "Level for rear seat KSI", ylab = "Level for front seat KSI")

abline(abline(lm(levels[, "level.Front.seats.KSI"] ~ levels[, "level.Rear.seats.KSI"])))

legend("topleft",leg = c("Level  rear against front seat KSI", "Regression line"), 
       lty = c(NA, 1), pch=c(1, NA), cex = 0.55, horiz=T)


plot(data$logUKdriversKSI, xlab="",ylab = "",pch=3, xlim=c(0,200))
abline(coefs)
legend("topright",leg = c("log UK drivers KSI against time", "regression line"), 
       lty = c(NA, 1), pch=c(3,NA), cex = 0.55, horiz=T)


#Figure 9.7 Levels of treatment and control series, rank one model
#(levels of front and rear seat KSI)
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
levels <- outKFS$alphahat[, c("level.Front.seats.KSI", "level.Rear.seats.KSI")]
treatLevelplusInterv <- outKFS$alphahat[, "level.Front.seats.KSI"] + X[, "seatbeltLaw"] * outKFS$alphahat[, "seatbeltLaw.Front.seats.KSI"]
controlLevel <- outKFS$alphahat[, "level.Rear.seats.KSI"]

par(mfrow=c(2,1))
plot(treatLevelplusInterv, xlab= "", ylab = "", lty=1)
legend("topright",leg = "level + intervention front", lty = 1, cex = 0.55, horiz=T)

plot(controlLevel, xlab= "", ylab = "",lty=1)
legend("topright",leg = "level rear", lty = 1, cex = 0.55, horiz=T)
par(mfrow=c(1,1))
  
#Figure 9.8 Level of treatment series plus intervention, and level of control series, rank one model
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
levels <- outKFS$alphahat[, c("level.Front.seats.KSI", "level.Rear.seats.KSI")]

par(mfrow=c(2,1))
plot(levels[, "level.Front.seats.KSI"], xlab= "", ylab = "", lty=1)
legend("topright",leg = "level front", lty = 1, cex = 0.55, horiz=T)

plot(levels[, "level.Rear.seats.KSI"], xlab= "", ylab = "",lty=1)
legend("topright",leg = "level rear", lty = 1, cex = 0.55, horiz=T)
par(mfrow=c(1,1))


#Figure 9.9 Deterministic seasonal of treatment ans control series, rank one model
#(levels of front and rear seat KSI)
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
seasonal <- outKFS$alphahat[, c("sea_dummy1.Front.seats.KSI", "sea_dummy1.Rear.seats.KSI")]
treatLevelplusInterv <- outKFS$alphahat[, "level.Front.seats.KSI"] + X[, "seatbeltLaw"] * outKFS$alphahat[, "seatbeltLaw.Front.seats.KSI"]
controlLevel <- outKFS$alphahat[, "level.Rear.seats.KSI"]

par(mfrow=c(2,1))
plot(seasonal[, "sea_dummy1.Front.seats.KSI"], xlab= "", ylab = "", lty=1)
legend("topright",leg = "seasonal front", lty = 1, cex = 0.55, horiz=T)

plot(seasonal[, "sea_dummy1.Rear.seats.KSI"], xlab= "", ylab = "",lty=1)
legend("topright",leg = "seasonal rear", lty = 1, cex = 0.55, horiz=T)
par(mfrow=c(1,1))

#Version 0 (It works ) ####

rm(list = setdiff(ls(), lsf.str())) 

UKfrontrearseatKSI <- read.table("UKfrontrearseatKSI.txt", head=TRUE)  
data <- ts(UKfrontrearseatKSI[,2:3], start = 1969,frequency=12) #Two time series: front seat passengers killed or seriously injured and rear seat passengers killed or seriously injured

seatbeltLaw <- rep(c(0, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
X <- cbind(seatbeltLaw, UKfrontrearseatKSI[, c(5,4)])
colnames(X)[2:3] <- c("petrolPrice", "kmDriven" )
X <- ts(X, start = 1969,frequency=12)

model <- SSModel(data ~ -1 + kmDriven + petrolPrice + SSMregression(~seatbeltLaw, P1inf = 0, P1=100, data=X, index=1) + 
                   SSMcustom(Z = diag(2), T = diag(2), R = matrix(1, 2, 1), Q = matrix(1), P1inf = diag(2)) +
                   SSMseasonal(period=12, sea.type='dummy'), H = matrix(NA, 2, 2), data=X)


ownupdatefn <- function(pars, model){
  Qsd <- c(exp(pars[1:2]))
  #Qcorr <- tanh(pars[3])
  Q <- Qsd %o% Qsd
  #Q[1,2] <- Q[2,1] <- Q[1,2] * Qcorr
  Q[1,2] <- Q[2,1] <- Q[1,2]
  model$Q[,,1] <- Q
  Hsd <- c(exp(pars[3:4]))
  Hcorr <- tanh(pars[5])
  H <- Hsd %o% Hsd
  H[1,2] <- H[2,1] <- H[1,2] * Hcorr
  model$H[,,1] <- H
  model
}

fit <- fitSSM(inits = rep(0.001 , 4), model = model,  method = "L-BFGS-B")

fit <- fitSSM(inits = rep(1e-10 , 5), model = model, updatefn = ownupdatefn, method = "L-BFGS-B") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

fit <- fitSSM(inits = rep(1e-10 , 5), model = model, updatefn = ownupdatefn, method = "L-BFGS-B") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))


#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)



#Version 1 (wit set-off)####
rm(list = setdiff(ls(), lsf.str())) 

UKfrontrearseatKSI <- read.table("UKfrontrearseatKSI.txt", head=TRUE)  
data <- ts(UKfrontrearseatKSI[,2:3], start = 1969,frequency=12) #Two time series: front seat passengers killed or seriously injured and rear seat passengers killed or seriously injured

seatbeltLaw <- rep(c(0.001, 1), times=c(169, 23))  # Introduction of the seat belt law  (intervention variable)
X <- cbind(seatbeltLaw, UKfrontrearseatKSI[, c(5,4)])
colnames(X)[2:3] <- c("petrolPrice", "kmDriven" )
X <- ts(X, start = 1969,frequency=12)

model <- SSModel(data ~ SSMregression(~ 1, index = 2, remove.intercept = FALSE) + kmDriven + petrolPrice + SSMregression (~ seatbeltLaw, index = 1)   + 
                   SSMtrend(degree = 1, type = "common", state_names = "common_level", Q=matrix(1)) +
                   SSMseasonal(period = 12, sea.type ='dummy'), H = matrix(NA, 2, 2), data = X) #SSMregression(~ 1, index = 2, remove.intercept = FALSE) +

ownupdatefn <- function(pars, model){
  #model$H[,,1] <- diag(exp(pars[1:2])); model$H[1,2,1] <- model$H[2,1,1]<- pars[3]
  model$H[,,1] <- exp(0.5*pars[1:2])
  model$H[1,2,1] <- model$H[2,1,1] <- tanh(pars[3])*prod(sqrt(exp(0.5*pars[1:2])))
  model$Z[c("Front.seats.KSI", "Rear.seats.KSI"), "common_level",] <- pars[4:5]

model
}


fit <- fitSSM(inits = rep(1e-10, 5), model = model, updatefn = ownupdatefn, method = "CG") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))


#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)


#Version 2 (demeaned)####
rm(list = setdiff(ls(), lsf.str())) 

UKfrontrearseatKSI <- read.table("UKfrontrearseatKSI.txt", head=TRUE)  
data <- scale(ts(UKfrontrearseatKSI[,2:3], start = 1969,frequency=12), scale=F) #Two time series: front seat passengers killed or seriously injured and rear seat passengers killed or seriously injured

seatbeltLaw <- rep(c(0, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
X <- cbind(seatbeltLaw, UKfrontrearseatKSI[, c(5,4)])
colnames(X)[2:3] <- c("petrolPrice", "kmDriven" )
X <- ts(X, start = 1969,frequency=12)

model <- SSModel(data ~ kmDriven + petrolPrice + SSMregression (~ seatbeltLaw, index = 1)   + 
                   SSMtrend(degree = 1, type = "common", state_names = "common_level", Q=matrix(1)) +
                   SSMseasonal(period = 12, sea.type ='dummy'), H = matrix(NA, 2, 2), data = X) #SSMregression(~ 1, index = 2, remove.intercept = FALSE) +

ownupdatefn <- function(pars, model){
  #model$H[,,1] <- diag(exp(pars[1:2])); model$H[1,2,1] <- model$H[2,1,1]<- pars[3]
  model$H[,,1] <- exp(0.5*pars[1:2])
  model$H[1,2,1] <- model$H[2,1,1] <- tanh(pars[3])*prod(sqrt(exp(0.5*pars[1:2])))
  model$Z[c("Front.seats.KSI", "Rear.seats.KSI"), "common_level",] <- pars[4:5]
  
  model
}


fit <- fitSSM(inits = rep(0.1, 5), model = model, updatefn = ownupdatefn, method = "L-BFGS-B") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))


#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)





#Version 3 (it does not work)####
rm(list = setdiff(ls(), lsf.str())) 

UKfrontrearseatKSI <- read.table("UKfrontrearseatKSI.txt", head=TRUE)  
data <- ts(UKfrontrearseatKSI[,2:3], start = 1969,frequency=12) #Two time series: front seat passengers killed or seriously injured and rear seat passengers killed or seriously injured

seatbeltLaw <- rep(c(0, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
X <- cbind(seatbeltLaw, UKfrontrearseatKSI[, c(5,4)])
colnames(X)[2:3] <- c("petrolPrice", "kmDriven" )
X <- ts(X, start = 1969,frequency=12)

model <- SSModel(data ~ -1 + kmDriven + petrolPrice + SSMregression (~ seatbeltLaw, index = 1)   + 
                   SSMcustom(Z = matrix(NA, 2,1),T = 1, R = 1, Q = matrix(1), P1inf = 1) +
                   SSMseasonal(period = 12, sea.type ='dummy'), H = matrix(NA, 2, 2), data = X) #SSMregression(~ 1, index = 2, remove.intercept = FALSE) +

ownupdatefn <- function(pars, model){
  model$Z[c("Front.seats.KSI", "Rear.seats.KSI"), "custom1",] <- matrix(params[1:2],2,1)
  model$H[,,1] <- diag(exp(pars[3:4])); model$H[1,2,1] <- model$H[2,1,1]<- pars[5]
  print(model$Z)
  model
}


fit <- fitSSM(inits = rep(0.1, 5), model = model, method = "BFGS") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))


#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)



#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)


#Version5 (diagonal identity)####
rm(list = setdiff(ls(), lsf.str())) 

UKfrontrearseatKSI <- read.table("UKfrontrearseatKSI.txt", head=TRUE)  
data <- ts(UKfrontrearseatKSI[,2:3], start = 1969,frequency=12) #Two time series: front seat passengers killed or seriously injured and rear seat passengers killed or seriously injured

seatbeltLaw <- rep(c(0, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
X <- cbind(seatbeltLaw, UKfrontrearseatKSI[, c(5,4)])
colnames(X)[2:3] <- c("petrolPrice", "kmDriven" )
X <- ts(X, start = 1969,frequency=12)

model <- SSModel(data ~ kmDriven + petrolPrice + SSMregression (~ seatbeltLaw, index = 1, P1=100, P1inf=0)   + 
                   SSMtrend(degree = 1, Q=matrix(0, 2, 2)) +
                   SSMseasonal(period = 12, sea.type ='dummy'), H = matrix(0, 2, 2), data = X)

ownupdatefn <- function(pars, model){
  model$H[,,1] <- matrix(c(exp(pars[1]), pars[2], pars[2], exp(pars[3])), 2,2)
  Q <- matrix(c(sqrt(pars[4]), 0, pars[5]*sqrt(pars[4]), 0), 2,2)
  model$Q[,,1] <- crossprod(Q)
  model
}

fit <- fitSSM(inits = rep(0.001, 5), model = model, updatefn = ownupdatefn, method = "L-BFGS-B") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

(H <- fit$model$H)  

(Q <- fit$model$Q)


 

#Model with a cumulator variables####
#see article: MESSY TIME SERIES: A  UNIFIED  APPROACH, on https://www.stat.berkeley.edu/~brill/Stat248/messyts.pdf

#Version 1####
rm(list = setdiff(ls(), lsf.str()))
if(!(require(eurostat))){install.packages('eurostat')}
library(eurostat)
if(!(require(KFAS))){install.packages('KFAS')}
library(KFAS)

#if(!(require(lubridate))){install.packages('lubridate')}
#library(lubridate)
#dataTemp2$time <- dataTemp2$time %m+% months(1)

#label_eurostat(x="namq_10_gdp")
label_eurostat_tables("namq_10_gdp", lang = "en")
namq_10_gdp <- get_eurostat("namq_10_gdp", time_format = "date_last", 
                            filters = list(geo = "EA19", na_item = "B1GQ", unit = "CLV15_MEUR", s_adj = "SCA"))

namq_10_gdp <- get_eurostat("namq_10_gdp", time_format = "date_last") 

dataTemp <- subset(namq_10_gdp, geo == "EA19" & na_item == "B1GQ" & unit == "CLV15_MEUR" & s_adj == "SCA" & time >= "1995-03-31", select=c("time", "values"))                         

dataGDP <- c(NA, as.vector(sapply(dataTemp$values, function(x) c(NA, NA, x))))# To covert to monthly data, two NA are added before the quarterly value; and NA is added at the begining to shift the data 
dataGDP <- ts(log(dataGDP), start=c(1995, 1), frequency=12)# Quarterly values instead of beining in the last month of each quarter are shifted forward by one month
                                                          # This is for a cumulator variable which is agregated in April, July, October and then January, April ...

model <- SSModel(dataGDP ~ SSMtrend(degree=1, Q=list(matrix(NA))) + SSMcustom(Z = 1, T = 1, R = 1, a1 = 0, P1inf = 1, P1 = 0, Q = matrix(NA), state_names = "cumulator"), H = matrix(0))
#Warning: if you change name of data change also this name in the part defining matrix T

ownupdatefn <- function(pars,model){
  #model$H[1,1,] <- exp(pars[1])
  
  diag(model$Q[,,1]) <- exp(pars[1:2])
  
  model$Z[1,1,1] <- 0
  
  T <- model$T #Creating time varying matrix T
  T[2,1,] <- 1
  T <- array(T, c(2,2, length(dataGDP))) #Make sure that this name of data is the same as in the definition of the model above
  T[2, 2, 1:dim(T)[3] %% 3 == 1]  <- 0 #Changing to 0 the value in the second row and second column of the matrix for time t4, t8, t12 ...
  model$T <- T                        # array dim=c(n.rows, n.columns, n.third.dimension)
  
  model$R[2,1,1] <- 1
  
  #attr(model, "tv") <- as.integer(rep(1,5))
  attr(model, "tv") <- as.integer(c(0, 0, 1, 0, 0))
  #Integer vector stating whether Z,H,T,R or Q is time-varying (indicated by 1 intvand 0 otherwise).  
  #If you manually change the dimensions of the matrices youmust change this attribute also.
  #Make sure that 1 is integer! 
  model
}

fit <- fitSSM(model, inits = log(c(0.001,0.001)), updatefn = ownupdatefn ,method = "BFGS")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
outKFS$model$H
outKFS$model$Q

#Version 2####
rm(list = setdiff(ls(), lsf.str()))
if(!(require(eurostat))){install.packages('eurostat')}
library(eurostat)
if(!(require(KFAS))){install.packages('KFAS')}
library(KFAS)


label_eurostat_tables("namq_10_gdp", lang = "en")
namq_10_gdp <- get_eurostat("namq_10_gdp", time_format = "date_last", 
                            filters = list(geo = "EA19", na_item = "B1GQ", unit = "CLV15_MEUR", s_adj = "SCA"))

namq_10_gdp <- get_eurostat("namq_10_gdp", time_format = "date_last") 

dataTemp <- subset(namq_10_gdp, geo == "EA19" & na_item == "B1GQ" & unit == "CLV15_MEUR" & s_adj == "NSA" & time >= "1995-03-31", select=c("time", "values"))                         

dataGDP <- c(NA, as.vector(sapply(dataTemp$values, function(x) c(NA, NA, x))))# To covert to monthly data, two NA are added before the quarterly value; and NA is added at the begining to shift the data 
dataGDP <- ts(log(dataGDP), start=c(1995, 1), frequency=12)# Quarterly values instead of beining in the last month of each quarter are shifted forward by one month
# This is for a cumulator variable which is agregated in April, July, October and then January, April ...

model <- SSModel(dataGDP ~ SSMtrend(degree=1, Q=list(matrix(NA)), a1 = 0, P1inf = 0, P1 = 50) + SSMseasonal(period = 12, Q=matrix(NA), a1 = rep(0, 11), P1inf = diag(rep(0, 11)), P1 = diag(rep(50, 11))) + SSMcustom(Z = 1, T = 1, R = 1, a1 = 0, P1inf = 0, P1 = 50, Q = matrix(NA), state_names = "cumulator"), H = matrix(0))
# It works with the model defined above; we have to adjust P1 for cumulator variable in SSMcustom between 0 and 100
#model <- SSModel(dataGDP ~ SSMtrend(degree=1, Q=list(matrix(NA))) + SSMseasonal(period = 12, Q=matrix(NA)) + SSMcustom(Z = 1, T = 1, R = 1, a1 = 0, P1inf = 0, P1 = 50, Q = matrix(NA), state_names = "cumulator"), H = matrix(0))

#Warning: if you change name of data change also this name in the part defining matrix T

ownupdatefn <- function(pars,model){
  
  diag(model$Q[, , 1]) <- exp(pars[1:3])#Third element of the diagonal matrix Q is the irregular disturbance for the cumulator variable
  
  model$T["cumulator", , 1] <-   model$Z[1, , 1] %*% model$T[, , 1]
  model$T <- replicate(length(dataGDP), model$T[,,1])#Make sure that this name of data is the same as in the definition of the model above
  model$T["cumulator", "cumulator", 1:dim(model$T)[3] %% 3 == 1]  <- 0 #Changing to 0 the value in the row and column of the cumulator variable of matriz T for time t4, t8, t12 ...
  model$T["cumulator", "cumulator", 1]  <- 0
  
  model$R["cumulator", , 1] <- model$Z[1, , 1] %*% model$R[, , 1]
  
  model$Z[1, colnames(model$Z)!="cumulator", 1] <- 0   # All the elements of matriz Z should be 0  except the last element that is 1 for the cumulator variable
  
  attr(model, "tv") <- as.integer(c(0, 0, 1, 0, 0)) #attr(model, "tv") <- as.integer(rep(1,5))
     #Integer vector stating whether Z,H,T,R or Q is time-varying (indicated by 1 intvand 0 otherwise).  
  #If you manually change the dimensions of the matrices youmust change this attribute also.
  #Make sure that 1 is integer!
  model
}

#is.SSModel(model)
fit <- fitSSM(model, inits = log(c(0.1, 0.1, 0.1)), updatefn = ownupdatefn, method = "BFGS") #"SANN"
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
outKFS$model$H
outKFS$model$Q



#Eurostat data####
if(!(require(knitr))){install.packages('knitr')}
library(knitr)
if(!(require(eurostat))){install.packages('eurostat')}
library(eurostat)

#Downloading a list of datasets with their titles and codes
datasets <- get_eurostat_toc()
View(datasets)

#label_eurostat_tables("ei_bsin_m_r2", lang = "en")

#Downloading the specific dataset/table
ei_bsin_m_r2 <- get_eurostat("ei_bsin_m_r2", time_format = "date_last")

#Listing codes with labels of all the dimensions/variables of the dataset
labels <- function (dataset) {
  kable(cbind(names(dataset)[!(names(dataset)) %in% "values"], label_eurostat_vars(names(dataset), lang = "en")), col.names=c("Code of dimension", "Label of dimension"))
  }

labels(ei_lmjv_q_r2) 

#Listing codes with lables for the factors of the chosen dimenions 
labels2 <- function (dataset, dimension) {
  datasetDimension <- as.name(paste(dataset,"$",dimension, sep=""))
  kable(cbind(levels(datasetDimension), label_eurostat(levels(datasetDimension), dic=dimension)), col.names=c("Code", "Label"))
  
}

labels2(ei_lmjv_q_r2, "nace_r2")


#data <-data.frame(ei_bssi_m_r2)

label1 <- "Industry - monthly data, EA19 -  Employment expectations over the next 3 months (BS-IEME)"
ei_bsin_m_r2 <- get_eurostat("ei_bsin_m_r2", time_format = "date_last")
data <- subset(ei_bsin_m_r2, indic == "BS-IEME" & s_adj == "SA" & geo == "EA19" & time >= "2010-01-01",  select=c("time", "values"))
View(data)
datats <- ts(data[,"values"], start=2010, frequency=12)
plot(datats)


label2 <- "Job vacancy rate, EA19, Industry (B-E)"
ei_lmjv_q_r2 <- get_eurostat("ei_lmjv_q_r2", time_format = "date_last")
data2 <- subset(ei_lmjv_q_r2, indic == "JOBRATE" & nace_r2 == "B-N" & sizeclas == "TOTAL" & s_adj == "SA" & geo == "EA19" & time >= "2010-01-01",  select=c("time", "values"))
View(data2)
data2ts <- ts(data2[,"values"], start=2010, frequency=4)
plot(data2ts)

par(mfrow=c(2,1))
plot(datats, main=label1)
plot(data2ts, main=label2)

########################################################

#Ver 2 with cumulator variable####

#Removing all objects except func tions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt"))
dataUKdriversKSI[193,] <- NA
dataUKdriversKSI[2:nrow(dataUKdriversKSI),] <- dataUKdriversKSI[1:nrow(dataUKdriversKSI)-1,] 
dataUKdriversKSI[1,] <- NA 
dataUKdriversKSI <- ts(dataUKdriversKSI, start = 1969, frequency=12)

#Fitting model

#model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree=1, Q=list(matrix(NA)), a1 = 0, P1inf = 0, P1 = 100) + SSMseasonal(period = 12, sea.type='dummy', Q=matrix(NA), a1 = rep(0, 11), P1inf = diag(rep(0, 11)), P1 = diag(rep(100, 11))) + SSMcustom(Z = 1, T = 1, R = 1, a1 = 0, P1inf = 0, P1 = 0, Q = matrix(NA), state_names = "cumulator"), H = matrix(0))

#Important: in SSMcustom P1inf = 0 and P1 = 0
model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree=1, Q=list(matrix(NA)), a1 = 0, P1inf = 1, P1 = 0) + SSMseasonal(period = 12, sea.type='dummy', Q=matrix(NA), a1 = rep(0, 11), P1inf = diag(rep(1, 11)), P1 = diag(rep(0, 11))) + SSMcustom(Z = 1, T = 1, R = 1, a1 = 0, P1inf = 0, P1 = 0, Q = matrix(NA), state_names = "cumulator"), H = matrix(0))
#Important: in SSMcustom P1inf = 0 and P1 = 0 to 50
ownupdatefn <- function(pars,model){
  
  diag(model$Q[, , 1]) <- exp(pars[1:3])#Third element of the diagonal matrix Q is the irregular disturbance for the cumulator variable
  
  model$T["cumulator", , 1] <-   model$Z[1, , 1] %*% model$T[, , 1]
  model$T <- replicate(length(dataUKdriversKSI), model$T[,,1])#Make sure that this name of data is the same as in the definition of the model above
  model$T["cumulator", "cumulator", ]  <- 0 
  #model$T["cumulator", "cumulator", 1:dim(model$T)[3] %% 3 == 1]  <- 0 #Changing to 0 the value in the row and column of the cumulator variable of matriz T for time t4, t8, t12 ...
  #model$T["cumulator", "cumulator", 1]  <- 1
  
  model$R["cumulator", , 1] <- model$Z[1, , 1] %*% model$R[, , 1]  
  
  model$Z[1, colnames(model$Z)!="cumulator", 1] <- 0 # All the elements of matriz Z should be 0  except the last element that is 1 for the cumulator variable
  
  attr(model, "tv") <- as.integer(c(0, 0, 1, 0, 0)) #attr(model, "tv") <- as.integer(rep(1,5))
  #Integer vector stating whether Z,H,T,R or Q is time-varying (indicated by 1 intvand 0 otherwise).  
  #If you manually change the dimensions of the matrices youmust change this attribute also.
  #Make sure that 1 is integer!
  model
}

fit <- fitSSM(inits = log(c(0.001, 0.001, 0.001)), model = model, updatefn = ownupdatefn, method = "BFGS") #"Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

#model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree=2, Q=list(matrix(NA), matrix(NA))), H = matrix(NA))
#fit <- fitSSM(model, inits = c(0.001,0.001, 0.001) ,method = "BFGS")

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
outKFS$model$Q

d <- q <- 12 #Number of diffuse initial values in the state 
w <- 3#Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15#First k autocorrelations to be used in Q-statisticlogLik <- logLik( ) dlmLL(dataUKdriversKSI, mod)
n <- 192 #Number of observations

#Maximum likelihood 
(maxLik <- logLik(fit$model, method = "BFGS")/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Maximum likelihood estimate of the initial values of the level and the seasonal components at time point t=1
(initVal <- coef(outKFS$model)[1,])

#Extracting residuals
predResid <- rstandard(outKFS) #One-step-ahead prediction residuals (standardised)
irregResid <- rstandard(outKFS, "pearson") #Auxiliary irregular  residuals (standardised)
levelResid <- rstandard(outKFS, "state") #Auxiliary level  residuals (standardised)


#Diagnostic for one-step-ahead prediction residuals (standardised)
qStat <- qStatistic(predResid, k, w)
rStat <- rStatistic(predResid, d, l)
hStat <- hStatistic(predResid, d)
nStat <- nStatistic(predResid, d)
dTable(qStat, rStat, hStat, nStat)

#Ver 3 Comparing residuals####
rm(list = setdiff(ls(), lsf.str())) 
#Model with cumulator variable
#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt"))

#We have to add one time point with NA at the begining of the time series
dataUKdriversKSI[193,] <- NA 
dataUKdriversKSI[2:nrow(dataUKdriversKSI),] <- dataUKdriversKSI[1:nrow(dataUKdriversKSI)-1,] 
dataUKdriversKSI[1,] <- NA 
#If we want to align residuals with time points we start one month earlier 12/1968 see plot
dataUKdriversKSI <- ts(dataUKdriversKSI, start = c(1968, 12), frequency=12)
#If we want to aling y hats with time points we start 1969 as the time series on 1/1969

#Fitting model

#model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree=1, Q=list(matrix(NA)), a1 = 0, P1inf = 0, P1 = 100) + SSMseasonal(period = 12, sea.type='dummy', Q=matrix(NA), a1 = rep(0, 11), P1inf = diag(rep(0, 11)), P1 = diag(rep(100, 11))) + SSMcustom(Z = 1, T = 1, R = 1, a1 = 0, P1inf = 0, P1 = 0, Q = matrix(NA), state_names = "cumulator"), H = matrix(0))

#Important: in SSMcustom P1inf = 0 and P1 = 0
model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree=1, Q=list(matrix(NA)), a1 = 0, P1inf = 1, P1 = 0) + SSMseasonal(period = 12, sea.type='dummy', Q=matrix(NA), a1 = rep(0, 11), P1inf = diag(rep(1, 11)), P1 = diag(rep(0, 11))) + SSMcustom(Z = 1, T = 1, R = 1, a1 = 0, P1inf = 0, P1 = 0, Q = matrix(NA), state_names = "cumulator"), H = matrix(0))
#Important: in SSMcustom P1inf = 0 and P1 = 0
ownupdatefn <- function(pars,model){
  
  diag(model$Q[, , 1]) <- exp(pars[1:3])#Third element of the diagonal matrix Q is the irregular disturbance for the cumulator variable
  
  model$T["cumulator", , 1] <-   model$Z[1, , 1] %*% model$T[, , 1]
  model$T <- replicate(length(dataUKdriversKSI), model$T[,,1])#Make sure that this name of data is the same as in the definition of the model above
  model$T["cumulator", "cumulator", ]  <- 0 
  #model$T["cumulator", "cumulator", 1:dim(model$T)[3] %% 3 == 1]  <- 0 #Changing to 0 the value in the row and column of the cumulator variable of matriz T for time t4, t8, t12 ...
  #model$T["cumulator", "cumulator", 1]  <- 1
  
  model$R["cumulator", , 1] <- model$Z[1, , 1] %*% model$R[, , 1]
  
  model$Z[1, colnames(model$Z)!="cumulator", 1] <- 0 # All the elements of matriz Z should be 0  except the last element that is 1 for the cumulator variable
  
  attr(model, "tv") <- as.integer(c(0, 0, 1, 0, 0)) #attr(model, "tv") <- as.integer(rep(1,5))
  #Integer vector stating whether Z,H,T,R or Q is time-varying (indicated by 1 intvand 0 otherwise).  
  #If you manually change the dimensions of the matrices youmust change this attribute also.
  #Make sure that 1 is integer!
  model
}

fit <- fitSSM(inits = log(c(0.001, 0.001, 0.001)), model = model, updatefn = ownupdatefn, method = "BFGS") #"Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

#model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree=2, Q=list(matrix(NA), matrix(NA))), H = matrix(NA))
#fit <- fitSSM(model, inits = c(0.001,0.001, 0.001) ,method = "BFGS")

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
outKFS$model$Q
predResidC <- rstandard(outKFS)

#Model without cumulator variable

#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt"))
dataUKdriversKSI <- ts(dataUKdriversKSI, start = 1969, frequency=12)

#Fitting model
model <- SSModel(dataUKdriversKSI ~ SSMtrend(1, Q=NA) + SSMseasonal(12, sea.type='dummy', Q = NA),  H=NA)

ownupdatefn <- function(pars,model,...){
  model$H[,,1] <- exp(pars[1])
  diag(model$Q[,,1]) <- exp(c(pars[2], pars[3]))
  model
}

fit <- fitSSM(inits = log(c(0.001, 0.001, 0.001)), model = model, updatefn = ownupdatefn, method = "BFGS")

#model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree=2, Q=list(matrix(NA), matrix(NA))), H = matrix(NA))
#fit <- fitSSM(model, inits = c(0.001,0.001, 0.001) ,method = "BFGS")

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
outKFS$model$Q

predResidNC <- rstandard(outKFS)


plot(predResidNC, type="l", lty="solid", col="red") # Residuals for the results without cumulator variable
lines(predResidC,  type="l", lty="dashed", col="blue") # Residuals for the case with cumulator variable
legend("topright", legend=c("Residuals without cumulator variable", "Residuals with cumulator variable"), lty=c("solid", "dashed"), col=c("red", "blue"))
# If we do not add one time index at the begining of  the time series with cumulator variable we loose the first residual



####################3

rm(list = setdiff(ls(), lsf.str())) 

UKfrontrearseatKSI <- read.table("UKfrontrearseatKSI.txt", head=TRUE)  
data <- ts(UKfrontrearseatKSI[,2:3], start = 1969,frequency=12) #Two time series: front seat passengers killed or seriously injured and rear seat passengers killed or seriously injured

seatbeltLaw <- rep(c(0.001, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
X <- cbind(seatbeltLaw, UKfrontrearseatKSI[, c(5,4)])
colnames(X)[2:3] <- c("petrolPrice", "kmDriven" )
X <- ts(X, start = 1969,frequency=12)


model <- SSModel(data ~ -1 + SSMregression(~ kmDriven, P1=100, P1inf=0) + SSMregression(~petrolPrice, P1=100, P1inf=0)  + SSMregression(~ seatbeltLaw, P1=100, P1inf=0, index=1)
                   + SSMseasonal(period=12, sea.type='dummy')
                   + SSMcustom(Z = diag(2),T = diag(2), R = matrix(1,2,1), Q = matrix(1), P1inf = diag(2)), 
                   H = matrix(NA, 2, 2), data=X)
                       

funUpdate <- function(pars, model){
  model$H[,,1] <- exp(0.5*pars[1:2])    
  model$H[1,2,1] <- model$H[2,1,1] <- tanh(pars[3])*prod(sqrt(exp(0.5*pars[1:2])))    
  model$R[28:29] <- exp(pars[4:5])    ## print(model$H[,,1])    ## print(model$R[28:29])    
  model  
  } 

fit <- fitSSM(inits = rep(0.0001, 5), model = model,updatefn = funUpdate,  method = "L-BFGS-B") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))



#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

##########WORKING###############
data("Seatbelts") # Seatbelts is in time series format
Seatbelts[, 1:7] <- log(Seatbelts[, 1:7])
data <- Seatbelts[,3:4]

X <- ts(Seatbelts[,c(5,6,8)], start = 1969,frequency=12)

#model <- SSModel(data ~ SSMregression(~ law, Q=diag(0,2), P1=1e6, P1inf=0, index=c(1, 2)) + 
                  # SSMregression(~ PetrolPrice, Q=diag(0,2), P1=1e6, P1inf=0, index=c(1, 2)) + 
                  # SSMregression(~ kms, Q=diag(0,2), P1=1e6, P1inf=0, index=c(1, 2)) +
                  # SSMtrend(degree=1, Q=matrix(NA, 2, 2)) +
                  # SSMseasonal(period=12, sea.type='dummy', Q=diag(0,2), index=c(1, 2)), H = matrix(NA, 2, 2), data=X)

model <- SSModel(data ~ SSMregression(~ law) + 
                   SSMregression(~ PetrolPrice) + 
                   SSMregression(~ kms) +
                   SSMtrend(degree=1, Q=matrix(NA, 2, 2)) +
                   SSMseasonal(period=12, sea.type='dummy'), H = matrix(NA, 2, 2), data=X)


#######################################

rm(list = setdiff(ls(), lsf.str())) 

UKfrontrearseatKSI <- read.table("UKfrontrearseatKSI.txt", head=TRUE)  
data <- ts(UKfrontrearseatKSI[,2:3], start = 1969,frequency=12) #Two time series: front seat passengers killed or seriously injured and rear seat passengers killed or seriously injured

seatbeltLaw <- rep(c(0, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
X <- cbind(seatbeltLaw, UKfrontrearseatKSI[, c(5,4)])
colnames(X)[2:3] <- c("petrolPrice", "kmDriven" )
X <- ts(X, start = 1969,frequency=12)

model <- SSModel(data ~ kmDriven + petrolPrice + seatbeltLaw + 
                   SSMtrend(degree=1, Q=list(matrix(NA, 2, 2))) +
                   SSMseasonal(period=12, sea.type='dummy'), H = matrix(NA, 2, 2), data=X)

fit <- fitSSM(inits = rep(0, 6), model = model,   method = "L-BFGS-B") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"
 
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

################################################ Does not WORK!! sometimes works

rm(list=ls())

UKfrontrearseatKSI <- read.table("UKfrontrearseatKSI.txt", head=TRUE)  
data <- ts(UKfrontrearseatKSI[,2:3], start = 1969,frequency=12) #Two time series: front seat passengers killed or seriously injured and rear seat passengers killed or seriously injured

seatbeltLaw <- rep(c(0, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
X <- cbind(seatbeltLaw, UKfrontrearseatKSI[, c(5,4)])
colnames(X)[2:3] <- c("petrolPrice", "kmDriven" )
X <- ts(X, start = 1969,frequency=12)

model <- SSModel(data ~ kmDriven + petrolPrice + seatbeltLaw + 
                   SSMtrend(degree=1, Q=list(matrix(NA, 2, 2))) +
                   SSMseasonal(period=12, sea.type='dummy'), H = matrix(NA, 2, 2), data=X)

funUpdate <- function(pars, model){
  model$H <- array(c(pars[1], pars[2], pars[2],  pars[3]) , c(2,2,1))    
  model$Q <- array(c(pars[4], pars[5], pars[5],  pars[6]) , c(2,2,1)) 
  model  
} 


fit <- fitSSM(inits = rep(0, 6), model = model,   ownupdatefn=funUpdate, method = "L-BFGS-B") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)#

######################################## WORKS
rm(list=ls())

UKfrontrearseatKSI <- read.table("UKfrontrearseatKSI.txt", head=TRUE)  
data <- ts(UKfrontrearseatKSI[,2:3], start = 1969,frequency=12) #Two time series: front seat passengers killed or seriously injured and rear seat passengers killed or seriously injured

seatbeltLaw <- rep(c(0, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
X <- cbind(seatbeltLaw, UKfrontrearseatKSI[, c(5,4)])
colnames(X)[2:3] <- c("petrolPrice", "kmDriven" )
X <- ts(X, start = 1969,frequency=12)

model <- SSModel(data ~ kmDriven + petrolPrice + seatbeltLaw + 
                   SSMtrend(degree=1, Q=list(matrix(NA, 2, 2))) +
                   SSMseasonal(period=12, sea.type='dummy'), H = matrix(NA, 2, 2), data=X)

funUpdate <- function(pars, model){
  model$H <- array(c(exp(pars[1]), pars[2], pars[2],  exp(pars[3])), c(2,2,1))    
  model$Q <- array(c(exp(pars[4]), pars[5], pars[5],  exp(pars[6])) , c(2,2,1)) 
  model  
} 


fit <- fitSSM(inits = rep(0, 6), model = model,   ownupdatefn=funUpdate, method = "L-BFGS-B") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)



#################################################################################Does not WORK!!!
rm(list=ls())

UKfrontrearseatKSI <- read.table("UKfrontrearseatKSI.txt", head=TRUE)  
data <- ts(UKfrontrearseatKSI[,2:3], start = 1969,frequency=12) #Two time series: front seat passengers killed or seriously injured and rear seat passengers killed or seriously injured

seatbeltLaw <- rep(c(0, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
X <- cbind(seatbeltLaw, UKfrontrearseatKSI[, c(5,4)])
colnames(X)[2:3] <- c("petrolPrice", "kmDriven" )
X <- ts(X, start = 1969,frequency=12)

model <- SSModel(data ~ kmDriven + petrolPrice + seatbeltLaw + 
                   SSMtrend(degree=1, Q=list(matrix(NA, 2, 2))) +
                   SSMseasonal(period=12, sea.type='dummy'), H = matrix(NA, 2, 2), data=X)

funUpdate <- function(pars, model){
  model$H[,,1] <- matrix(c(exp(pars[1]), pars[2], pars[2],  exp(pars[3])) , 2,2)    
  model$Q[,,1] <- matrix(c(exp(pars[4]), pars[5], pars[5],  exp(pars[6])) , 2,2) 
  model  
} 


fit <- fitSSM(inits = rep(0, 6), model = model,  updatefn = funUpdate, method = "L-BFGS-B") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#################################################Does not WORK!!!

rm(list=ls())

UKfrontrearseatKSI <- read.table("UKfrontrearseatKSI.txt", head=TRUE)  
data <- ts(UKfrontrearseatKSI[,2:3], start = 1969,frequency=12) #Two time series: front seat passengers killed or seriously injured and rear seat passengers killed or seriously injured

seatbeltLaw <- rep(c(0, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
X <- cbind(seatbeltLaw, UKfrontrearseatKSI[, c(5,4)])
colnames(X)[2:3] <- c("petrolPrice", "kmDriven" )
X <- ts(X, start = 1969,frequency=12)

model <- SSModel(data ~ kmDriven + petrolPrice + seatbeltLaw + 
                   SSMtrend(degree=1, Q=list(matrix(NA, 2, 2))) +
                   SSMseasonal(period=12, sea.type='dummy'), H = matrix(NA, 2, 2), data=X)

ownupdatefn <- function(pars,model){
  diag(model$H[,,1]) <- pars[1:2]
  model$H[1,2,1] <- pars[3]
  model$H[2,1,1] <- pars[3]
  diag(model$Q[,,1]) <- pars[4:5]
  model$Q[1,2,1] <- pars[6]
  model$Q[2,1,1] <- pars[6]
  model
}

fit <- fitSSM(inits = rep(0, 6), model = model,  updatefn = ownupdatefn,  method = "L-BFGS-B") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#######################################################does not WORK!!!

rm(list=ls())

UKfrontrearseatKSI <- read.table("UKfrontrearseatKSI.txt", head=TRUE)  
data <- ts(UKfrontrearseatKSI[,2:3], start = 1969,frequency=12) #Two time series: front seat passengers killed or seriously injured and rear seat passengers killed or seriously injured

seatbeltLaw <- rep(c(0.001, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
X <- cbind(seatbeltLaw, UKfrontrearseatKSI[, c(5,4)])
colnames(X)[2:3] <- c("petrolPrice", "kmDriven" )
X <- ts(X, start = 1969,frequency=12)

model <- SSModel(data ~ kmDriven + petrolPrice + seatbeltLaw + 
                   SSMtrend(degree=1, Q=list(matrix(NA, 2, 2))) +
                   SSMseasonal(period=12, sea.type='dummy'), H = matrix(NA, 2, 2), data=X)

ownupdatefn <- function(pars,model){
  diag(model$H[,,1]) <- exp(pars[1:2])
  model$H[1,2,1] <- pars[3]
  model$H[2,1,1] <- pars[3]
  diag(model$Q[,,1]) <- exp(pars[4:5])
  model$Q[1,2,1] <- pars[6]
  model$Q[2,1,1] <- pars[6]
  model
}

fit <- fitSSM(inits = rep(0, 6), model = model,  updatefn = ownupdatefn, method = "L-BFGS-B") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)




######################Does not WORK!!!!!!!!!!!!!!!!!!!!


rm(list=ls())

UKfrontrearseatKSI <- read.table("UKfrontrearseatKSI.txt", head=TRUE)  
data <- ts(UKfrontrearseatKSI[,2:3], start = 1969,frequency=12) #Two time series: front seat passengers killed or seriously injured and rear seat passengers killed or seriously injured

seatbeltLaw <- rep(c(0, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
X <- cbind(seatbeltLaw, UKfrontrearseatKSI[, c(5,4)])
colnames(X)[2:3] <- c("petrolPrice", "kmDriven" )
X <- ts(X, start = 1969,frequency=12)

model <- SSModel(data ~ kmDriven + petrolPrice + seatbeltLaw + 
                   SSMtrend(degree=1, Q=list(matrix(NA, 2, 2))) +
                   SSMseasonal(period=12, sea.type='dummy'), H = matrix(NA, 2, 2), data=X)

ownupdatefn <- function(pars, model){
  model$H <- array(c(exp(pars[1]), pars[2], pars[2],  exp(pars[3])), c(2,2,1))
  
  model$Q <- array(c(pars[4]^2, pars[4]*pars[5], pars[4]*pars[5], pars[5]^2), c(2,2,1))
    model
}



fit <- fitSSM(inits = rep(1, 5), model = model, updatefn = ownupdatefn, method = "Nelder-Mead") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

######################################################################################################
rm(list = setdiff(ls(), lsf.str())) 

UKfrontrearseatKSI <- read.table("UKfrontrearseatKSI.txt", head=TRUE)  
data <- ts(UKfrontrearseatKSI[,2:3], start = 1969,frequency=12) #Two time series: front seat passengers killed or seriously injured and rear seat passengers killed or seriously injured

seatbeltLaw <- rep(c(0, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
X <- cbind(seatbeltLaw, UKfrontrearseatKSI[, c(5,4)])
colnames(X)[2:3] <- c("petrolPrice", "kmDriven" )
X <- ts(X, start = 1969,frequency=12)

model <- SSModel(data ~ kmDriven + petrolPrice +SSMregression(~seatbeltLaw, index=1)  + 
                   SSMtrend(degree=1, Q=list(matrix(2, 2, 2))) +
                   SSMseasonal(period=12, sea.type='dummy'), H = matrix(2, 2, 2), data=X)


#model <- SSModel(data ~ SSMregression(~ kmDriven, P1=1e6, P1inf=0) + SSMregression(~petrolPrice, P1=1e6, P1inf=0)  + SSMregression(~ seatbeltLaw, P1=1e6, P1inf=0)  + 
                  # SSMtrend(degree=1, Q=list(matrix(NA, 2, 2))) +
                  # SSMseasonal(period=12, sea.type='dummy'), H = matrix(NA, 2, 2), data=X)



ownupdatefn <- function(pars, model){
  model$H <- array(c(pars[1]^2, pars[1]*pars[2], pars[1]*pars[2], pars[2]^2+pars[3]^2), c(2,2,1))
  
  
  model$Q <- array(c(pars[4]^2, pars[4]*pars[5], pars[4]*pars[5], pars[5]^2), c(2,2,1))
  model
}




fit <- fitSSM(inits = rep(0, 6), model = model,   method = "L-BFGS-B") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

eigen(Q[,,1])



#########################################################

rm(list=ls()) 

UKfrontrearseatKSI <- read.table("UKfrontrearseatKSI.txt", head=TRUE)  
data <- ts(UKfrontrearseatKSI[,2:3], start = 1969,frequency=12) #Two time series: front seat passengers killed or seriously injured and rear seat passengers killed or seriously injured

seatbeltLaw <- rep(c(0, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
X <- cbind(seatbeltLaw, UKfrontrearseatKSI[, c(5,4)])
colnames(X)[2:3] <- c("petrolPrice", "kmDriven" )
X <- ts(X, start = 1969,frequency=12)

model <- SSModel(data ~ kmDriven + petrolPrice + seatbeltLaw + 
                   SSMtrend(degree=1, Q=list(matrix(NA, 2, 2))) +
                   SSMseasonal(period=12, sea.type='dummy'), H = matrix(NA, 2, 2), data=X)

ownupdatefn <- function(pars, model){
  model$H[,,1] <- matrix(c(pars[1]^2, pars[1]*pars[2], pars[1]*pars[2], pars[2]^2+pars[3]^2), 2, 2, byrow=T)
  
  model$Q[,,1] <- matrix(c(pars[4]^2, pars[4]*pars[5], pars[4]*pars[5], pars[5]^2), 2, 2, byrow=T)
  model
}


fit <- fitSSM(inits = rep(0, 6), model = model,   method = "L-BFGS-B") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)









#########################################################

model <- SSModel(data ~ SSMregression(~ seatbeltLaw, Q=diag(0,2), P1=1e6, P1inf=0, index=c(1, 2)) + 
                   SSMregression(~ petrolPrice, Q=diag(0,2), P1=1e6, P1inf=0, index=c(1, 2)) + 
                   SSMregression(~ kmDriven, Q=diag(0,2), P1=1e6, P1inf=0, index=c(1, 2)) +
                   SSMtrend(degree=1, Q=matrix(NA, 2, 2)) +
                   SSMseasonal(period=12, sea.type='dummy', Q=diag(0,2), index=c(1, 2)), H = matrix(NA, 2, 2), data=X)



ownupdatefn <- function(pars,model){
diag(model$H[,,1]) <- exp(pars[1:2])
model$H[1,2,1] <- pars[3]
model$H[2,1,1] <- pars[3]
diag(model$Q[,,1])[7:8] <- exp(pars[4:5])
model$Q[7,8,1] <- pars[6]
model$Q[8,7,1] <- pars[6]
}

ownupdatefn <- function(pars,model){
  diag(model[["H"]][,,1]) <- exp(pars[1:2])
  model[["H"]][1,2,1] <- pars[3]
  model[["H"]][2,1,1] <- pars[3]
  diag(model[["Q"]][,,1])[7:8] <- exp(pars[4:5])
  model[["Q"]][7,8,1] <- pars[6]
  model[["Q"]][8,7,1] <- pars[6]
}


pars <- c(10, 20, 30, 40, 50, 60)

diag(model$H[,,1]) <- pars[1:2]
model$H[1,2,1] <- pars[3]
model$H[2,1,1] <- pars[3]
diag(model$Q[,,1])[7:8] <- pars[4:5]
model$Q[7,8,1] <- pars[6]
model$Q[8,7,1] <- pars[6]


fit <- fitSSM(inits = rep(0.001, 6), model = model,   method = "CG") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)
########################################################################################################


#dataNOandFinFatalities <- log(read.table("NorwayFinland.txt")[,2:3])
#colnames(dataNOandFinFatalities) <- c("Norway", "Finland")
#dataNOandFinFatalities <- ts(dataNOandFinFatalities, start = 1970, frequency=1)
#plot(dataNOandFinFatalities)







#Subchapter 2.1 Deterministic level ####################################

fn <- function(params){
  dlmModPoly(order=1, dV=exp(params[1]), dW=0)
}
fit <- dlmMLE(dataUKdriversKSI, 0, fn) #fitting a model to data
mod <- fn(fit$par)
(V <- mod$V) #Variance of the observational noise
(W <- mod$W) #Variance of the system noise
(logLik <-fit$value) #Loglikelihood

logLik2 <- dlmLL(dataUKdriversKSI, mod)
numParams <- length(fit$par)
(AIC <- 2*numParams-2*logLik) #Akaike information criterion (AIC)

filtered <- dlmFilter(dataUKdriversKSI, mod) 
smoothed <- dlmSmooth(dataUKdriversKSI, mod) # or smoothed <- dlmSmooth(filtered)
residuals <-residuals(filtered, type="standardized", sd=FALSE)  #Residuals:one-step forecast errors
plot(residuals)

plot(dataUKdriversKSI)
lines(smoothed$s)

#Independence
acf(residuals) #Autocorrelation function
maxLag <-15
indTest <- Box.test(residuals, lag = maxLag, type = "Ljung", fitdf=numParams) #Ljung-Box test examining the null hypothesis of independence in residuals; we reject the null hypotesis of independence when p<0.05
indTest$p.value                                                               # fitdf number of degrees of freedom to be subtracted if x is a series of residuals
indTestCritValue <- qchisq(0.95, numLag-numParams+1)



#Normality
qqnorm(residuals)
qqline(residuals)
normTest <- jb.norm.test(residuals) # Jarque and Bera or Shenton and Bowman test; we reject the null hypothesis of normality when p<0.05 
normTest$p.value
 
#Homoscedasticity
residualsVector <- as.matrix(residuals)
lengthOneThird<-round(length(residualsVector)/3, digits = 0)
residualsVectorB <- residualsVector[1:lengthOneThird,]
residualsVectorE <- residualsVector[(length(residualsVector)-lengthOneThird+1):length(residualsVector),]
var.test(residualsVectorB, residualsVectorE, ratio = 1, alternative = "two.sided", conf.level = 0.95) #F test examining the null hypothesis of the same variances; we reject the null hypotesis when p<0.05; (Reliability: tested series should be nomrally distributed but they are not)

#KFAS
rm(list=ls())
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt"))
dataUKdriversKSI <- ts(dataUKdriversKSI, start = 1969,frequency=12)

model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree=1, Q=list(matrix(0))), H = matrix(NA))
fit <- fitSSM(model, inits = 0.001,method = "BFGS")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))


q <- 0 #Number of diffuse initial values in the state
d <- 0 #Number of diffuse initial elements
w <- 0 #total number of disturbance variances estimated (in univariate model) or number of parameters (in multivariate model)
predResid <- rstandard(outKFS) #One-step-ahead prediction residuals (standardised)
irregResid <- rstandard(outKFS, "pearson") #Auxiliary irregular  residuals (standardised)
levelResid <- rstandard(outKFS, "state") #Auxiliary level  residuals (standardised)

q <- 1
d <- 1
#Independence
acf(predResid[-(1:d)]) #Autocorrelation function

lagtest <- c(1, 12)
acf <- list(lag = acf(predResid[-(1:d)], plot=FALSE)$acf[-1],
  critical = 2 / sqrt(length(predResid)))



maxLag <-15
indTest  <- Box.test(predResid, lag = maxLag, type = "Ljung", fitdf=q-1) #Ljung-Box test examining the null hypothesis of independence in residuals; we reject the null hypotesis of independence when p<0.05
indTest$p.value                                                               # fitdf number of degrees of freedom to be subtracted if x is a series of residuals
indTestCritValue <- qchisq(0.95, maxLag-q+1)                           # degree of freedom max lags - num parameters + 1 what means that fitdf=numParam-1


######################
QStatistic <- function(predResid, maxLag) {
  value <- Box.test(predResid, lag = maxLag, type = "Ljung", fitdf=q-1)$statistic
  criticalValue <- qchisq(0.95, maxLag-q+1) 
  list(
    value = unname(value),
    criticalValue = criticalValue
  )
  }
QStatistic <- QStatistic(predResid, maxLag)
#################

#Normality
qqnorm(predResid)
qqline(predResid)

#normTest <- jb.norm.test(predResid[-(1:d)])$statistic # Jarque and Bera or Shenton and Bowman test; we reject the null hypothesis of normality when p<0.05 
#normTest$p.value
d <- 1 #Number of diffuse initial elements (or could it be a trace of matrix Pinf)
normTestStat <- jb.norm.test(predResid[-(1:d)])$statistic
names(normTestStat) <- NULL
normTest <- list(
  stat = normTestStat,
  critical = qchisq(0.95,2)
)


#####################

NStatistic <- function(predResid, d) {
  value <- jb.norm.test(predResid[-(1:d)])$statistic
  criticalValue <- qchisq(0.95,2)
  list(
    value = unname(value),
    criticalValue = criticalValue
  )
}
NStatistic <- NStatistic(predResid, d)
######################################

#Homoscedasticity

d <- 1 #Number of diffuse initial elements (or could it be a trace of matrix Pinf)
n <- length(predResid) # Number of observations
h <- round((n-d)/3, digits = 0) #Nearest integer to (n-d)/3
valueTest <- sum(predResid[(n-h+1):n]^2) / sum(predResid[(d+1):(d+h)]^2)
 
Ftest<- list(
  size = h,
  stat = valueTest,
  upper = qf(0.975, h, h),
  lower = qf(0.025, h, h)
)
FtestStat <- ifelse(Ftest$stat > 1, Ftest$stat, 1/Ftest$stat)

####################



######################


#z <- c(2, 3) x<- sum((z)^2)
#blockSize
#vectorPredResid <- as.matrix(predResid)
#lengthOneThird<-round(length(vectorPredResid)/3, digits = 0)
#residualsVectorB <- residualsVector[1:lengthOneThird,]
#residualsVectorE <- residualsVector[(length(residualsVector)-lengthOneThird+1):length(residualsVector),]
#var.test(residualsVectorB, residualsVectorE, ratio = 1, alternative = "two.sided", conf.level = 0.95) #F test examining the null hypothesis of the same variances; we reject the null hypotesis when p<0.05; (Reliability: tested series should be nomrally distributed but they are not)


tableTemplate <- c(

  "------------------------------------------------------",    
  "                   stat     value  critical satisfied",    
  "------------------------------------------------------",    
  "independence      Q  (%2d)  %7.3f   %5.2f     %1s",  # BoxLJung, 4 args    
  "                  r  (%2d)  %7.3f  +-%4.2f     %1s", # ACF,      4 args    
  "                  r  (%2d)  %7.3f  +-%4.2f     %1s", # ACF,      4 args    
  "homoscedasticity  %-3s(%2d)  %7.3f   %5.2f     %1s",  # Homo,     5 args    
  "normality         N        %7.3f   %5.2f     %1s",    # N,        3 args    
  "------------------------------------------------------"  
  ) 

cat(    sprintf(      paste(tableTemplate, collapse="\n"),       
                      # BoxLjung, 4 args      
                      maxLag,          
                      indTest$stat,      
                      indTestCritValue,      
                      ifelse(indTest$stat < indTestCritValue, "+", "-"),     
                      # ACF, 4 args      
                      lagtest[1],       
                      acf$lag[lagtest[1]],      
                      acf$critical,      
                      ifelse(abs(acf$lag[lagtest[1]]) < acf$critical, "+", "-"),      
                      # ACF, 4 args      
                      lagtest[2],      
                      acf$lag[lagtest[2]],      
                      acf$critical,      
                      ifelse(abs(acf$lag[lagtest[2]]) < acf$critical, "+", "-"),      
                      # Homo, 5 args      
                      ifelse(Ftest$stat > 1, "H", "1/H"),      
                      Ftest$size,       
                      FtestStat,       
                      Ftest$upper,       
                      ifelse(FtestStat < Ftest$upper, "+", "-"),      
                      # N, 3 args      
                      normTest$stat,        
                      normTest$critical,      
                      ifelse(normTest$stat < normTest$critical, "+", "-")    )  )  
cat("", "\n") 

Ftest<- list(
  size = h,
  stat = valueTest,
  upper = qf(0.975, h, h),
  lower = qf(0.025, h, h)
)



#Outliers

observationDisturbances <- dataUKdriversKSI-smoothed$s
observationDisturbances <- observationDisturbances/sd(observationDisturbances)
plot(observationDisturbances)
abline(h=1.96, lt="dashed")
abline(h=-1.96, lt="dashed")


rm( fit, mod)
rm(list=ls()[ls()!= "dataUKdriversKSI"])

fn <- function(params){
  mod <- dlmModPoly(order=1) + dlmModSeas(frequency=12)
  V(mod) <- exp(params[1])
  diag(W(mod))[1:2] <- exp(params[2:3])
  return(mod)
}

rm( fit, mod)

#Figure 2.1 Deterministic level
plot(dataUKdriversKSI, xlab="",ylab = "")
lines(out.fit.ssDL.UKdriversKSI$alphahat, lt=2)
legend("topright",leg = c("log UK drivers KSI", "deterministic level"), 
       lty = c(1, 2), cex = 0.55, horiz=T)

#Subchapter 2.2 Stochastic level #########################################
ssSL.UKdriversKSI <- SSModel(logUKdriversKSI ~ SSMtrend(degree = 1, Q = list(matrix(NA)), a1=0, P1=0, P1inf=1), #
                           H = matrix(NA), data = data.UKdriversKSI) 

ssSL.UKdriversKSI <- SSModel(data.UKdriversKSI.m[,"logUKdriversKSI"] ~ SSMtrend(degree = 1, Q = list(matrix(NA)), a1=0, P1=0, P1inf=1), #
                             H = matrix(NA)) 

fit.ssSL.UKdriversKSI <- fitSSM(ssSL.UKdriversKSI, inits=c(0.000001, 0.000001), method ="BFGS")
fit.ssSL.UKdriversKSI <- fitSSM(ssSL.UKdriversKSI, inits=c(log(var(data.UKdriversKSI)), log(var(data.UKdriversKSI))), method ="BFGS")
fit.ssSL.UKdriversKSI$model$H
fit.ssSL.UKdriversKSI$model$Q
#Improtant: To get infomration for state variables and disturbances (filtered or smooth)
#we use filtering = c('state', 'disturbance'), smoothing = c('state', 'disturbance')
out.fit.ssSL.UKdriversKSI <- KFS(fit.ssSL.UKdriversKSI$model, filtering = c('state', 'disturbance'), smoothing = c('state', 'disturbance'), simplify = F)
out.fit.ssSL.UKdriversKSI$logLik

#Figure 2.3 Stochastic level
plot(data.UKdriversKSI, xlab="",ylab = "")
#To get smoothed estimates of state we use sout.fit.ssDL.UKdriversKSI$alphahat 
lines(out.fit.ssSL.UKdriversKSI$alphahat, lt=2)
legend("topright",leg = c("log UK drivers KSI", "stochastic level"), 
       lty = c(1, 2), cex = 0.55, horiz=T)

#Figure 2.4 Irregular component for local level model
#To get smoothed observation disturbances (irregular component) 
#we use sout.fit.ssDL.UKdriversKSI$epshat
irregular.component <- out.fit.ssSL.UKdriversKSI$epshat
library(reshape2)
irregular.component.m <- melt(irregular.component)$value

res2s <- rstandard(KFS(fit.ssSL.UKdriversKSI$model), type="pearson")
res2 <- residuals(KFS(fit.ssSL.UKdriversKSI$model), type="pearson") 

res1s <- rstandard(KFS(fit.ssSL.UKdriversKSI$model), type="recursive") 
res1 <- residuals(KFS(fit.ssSL.UKdriversKSI$model), type ="recursive")

res1s.m <- as.matrix(res1s) 
rownames(res1s.m) <- dates

res.s <- rstandard(KFS(fit.ssSL.UKdriversKSI$model)) 
res <- residuals(KFS(fit.ssSL.UKdriversKSI$model)) 
out.fit.ssSL.UKdriversKSI$a
par(mfrow=c(1,1))

rm(list=c("res1", "res1s", "res2", "res2s"))
plot(res1s)
plot(res1)
graphics.off()
shapiro.test(res)

out.fit.ssSL.UKdriversKSI$v

irregular.component.s <- rstudent(irregular.component)
plot(irregular.component, xlab="",ylab = "", lt=1)
abline(h=0, lt=2)
legend("topright",leg = "irregular", lty = 1, cex = 0.55, horiz=T)


plot(residuals(out.fit.ssSL.UKdriversKSI, "recursive"), xlab="",ylab = "")

#Akaike information criterion for this analysis 


#Diagnostic tests for local level model and log UK drivers KSI

#Normality test
# Doornik-Hansen introduced a multivariate version of the univariate omnibus test for normality of Shenton and Bowman (1977), based on the transformed skewness and kurtosis. 
#The Doornik-Hansen test for multivariate normality is based on the skewness and kurtosis of multivariate data that is transformed to insure independence.
#The function DH.test runs the Doornik-Hansen test for both multivariate and univariate normality. The later test follows directly from the work of Bowman and Shenton.

shapiro.test(window(res1, start=1969.02))

#Jarque and Bera or Shenton and Bowman test
jarque.bera.test(window(res1, start=1969.02)) # we reject the null hypothesis of normality when p<0.5 
shapiro.test(res1s)


shapiro.test(res2)
shapiro.test(res2s)
bptest(res1s)

#homoscedasticity
plot(function(x)(df(x, df1=2, df2=150)), xlim=c(-2,10))
plot(function(x)(df(x, df1=4, df2=10)), xlim=c(-2,10), lty="dashed", add=T)
length(res1s)
length<-round(length(window(res1s, start=1969.02))/3, digits = 0)

part1<-window(res1s, start=c(1969, 2), end=c(1969,  length))
res1sB<-window(res1s, start=1969.2)
temp<-window(res1s, start=1969.02)


length<-round(length(res1s.m)/3, digits = 0)

part1<-res1s.m[2:(length+1),1, drop=F]
length(part1)
part2<-res1.s.m[len]
#Independence
Box.test(res1s, lag = 15, type = "Ljung") 



par(mfrow=c(2,1))

acf(window(res1s, start=1969.02))
acf(window(res1, start=1969.02))
acf(res1)
if(!(require(nortest))){install.packages('nortest')}
library(nortest)
ad.test(irregular.component)
#library(reshape2)
#irregular.component.m <- melt(irregular.component)$value
#qqplot(irregular.component)
# estimate the model using MLE and diffuse initialization
library(ggplot2)
gglagplot(res1s)
ggsubseriesplot(res1s)
ggseasonplot(res1s)
ggAcf(res1s) + theme(panel.background = element_rect(fill = 'lightgrey', colour = 'white'))

library(fpp2)
fitSSKSI = fitSSM(ssKSI, inits = c(1), method = "BFGS")

ssKSI$model$H

data(Nile)

model <- SSModel( Nile ~ SSMtrend(1, Q = list(matrix(NA))) , H = matrix(NA))
model <- fitSSM( inits = c(log(var(Nile)) , log(var(Nile) ) ) ,
                 model = model, method ='BFGS' )$model
out <- KFS(model, filtering='state', smoothing=c('state','disturbance'),simplify=F)
plot(residuals(out, "recursive"), xlab="",ylab = "")



modelNile <-SSModel(Nile~SSMtrend(1,Q=list(matrix(NA))),H=matrix(NA))
modelNile <-fitSSM(inits=c(log(var(Nile)),log(var(Nile))),
                   model=modelNile,
                   method='BFGS',
                   control=list(REPORT=1,trace=0))$model
out <- KFS(modelNile,filtering='state',smoothing='state')
state.sample <- simulateSSM(modelNile, type = c("states"))

temp <- cbind(state.sample,smoothed.state)
plot.ts(temp, plot.type="single",col =c(rep("grey",1),"blue"),
        ,lwd=c(rep(1,1),2), ylab = "",xlab="")
leg <-c("realized observation error","simulated observation error")
legend("topright",leg,col=c("blue","grey"),lwd=c(3,1),cex=0.7,bty="n")



##############################################



#Pacakages used
if(!(require(KFAS))){install.packages('KFAS')}
if(!(require(normtest))){install.packages('normtest')}
library(KFAS)
library(normtest)

#Deterministic level
rm(list=ls())
setwd("~/myfiles/Koopman/Chapter2")

data <- log(read.table("UKdriversKSI.txt"))
data <- as.matrix(data, nrow=nrow(data), ncol=ncol(data))
dates <- as.character(seq(as.Date("1969-02-01"), length=nrow(data), by="month")-1)
rownames(data) <- dates
colnames(data) <- "UKdriversKSI"
plot(x=as.Date(rownames(data)), y=data[,"UKdriversKSI"],xlab="", ylab="",  type = "l")

ssmodel <- SSModel(data[,"UKdriversKSI", drop=F] ~ SSMtrend(degree = 1, Q = list(matrix(0))), 
                             H = matrix(NA)) 
fitssm <- fitSSM(ssmodel, inits=0.001, method ="BFGS")
fitssm$model$logLik
fitssm$model$H
fitssm$model$Q

#Stochastic level
rm(list=ls())
setwd("~/myfiles/Koopman/Chapter2")
data <- log(read.table("UKdriversKSI.txt"))
data <- as.matrix(data, nrow=nrow(data), ncol=ncol(data))
dates <- seq(as.Date("1969-02-01"), length=nrow(data), by="month")-1
rownames(data) <- as.character(dates)
colnames(data) <- "UKdriversKSI"
plot(x=as.Date(rownames(data)), y=data[,"UKdriversKSI"],xlab="", ylab="",  type = "l")

#ssmodel <- SSModel(data[,"UKdriversKSI"] ~ SSMtrend(degree = 1, Q = list(matrix(NA))), 
#                  H = matrix(NA)) 

ssmodel <- SSModel(data[,"UKdriversKSI"] ~ SSMtrend(degree=1, Q=list(NA)) + 
                     SSMseasonal(period=12,sea.type='dummy', Q=0), 
                                 H=NA)
ownupdatefn <- function(pars,model){
  model$H[,,1] <- exp(pars[1])
  diag(model$Q[,,1]) <- (c(exp(pars[2]), 0))
  model
}

fitssm <- fitSSM(ssmodel, inits=c(0.01, 0.01),  method ="BFGS")
logLik(fitssm$model)
fitssm$model$H
fitssm$model$Q

model <- SSModel(data[,"UKdriversKSI"]  ~  SSMtrend(degree=1, Q=list(NA)) +
                   SSMseasonal(12, sea.type = 'dummy', Q = NA),  H = NA)


ownupdatefn <- function(pars,model){
  model$H[,,1] <- exp(pars[1])
  diag(model$Q[,,1]) <- exp(c(pars[2], pars[3]))
  model
}
fitssm <- fitSSM(inits = c(0.0001, 0.0001, 0.0001), model = model, updatefn=ownupdatefn)

logLik(fitssm$model)
fitssm$model$H
fitssm$model$Q

kfs <- KFS(fitssm$model, smoothing=c("disturbance", "signal"))

#Residuals  - standardised prediction errors 
residuals <- rstandard(kfs, "recursive")

# Auxiliary residuals - standardised smoothed observation disturbances
# Inspection to detect possible outliers 
auxresiduals1 <- rstandard(kfs, "pearson")

# Auxiliary residuals - standardised smoothed level (state) disturbances
#Inspection to detec possible 
auxresiduals2 <- rstandard(kfs, "state")

par(mfrow=c(2,1))

plot(x=dates, y=auxresiduals2, xlab="", ylab="",  type = "l")
abline(h=0, lt=2)
abline(h=1.96, lt=2)
abline(h=-1.96, lt=2)
legend("topright",leg = "Structural level break t-test", lty = 1, cex = 0.55, horiz=T)

plot(x=dates, y=auxresiduals1, xlab="", ylab="",  type = "l")
abline(h=0, lt=2)
abline(h=1.96, lt=2)
abline(h=-1.96, lt=2)
legend("topright",leg = "Outlier t-test", lty = 1, cex = 0.55, horiz=T)


#############DLM

rm(list=ls())
setwd("~/myfiles/Koopman/Chapter2")

if(!(require(dlm))){install.packages('dlm')}
if(!(require(normtest))){install.packages('normtest')}
library(dlm)
library(normtest)

data <- log(read.table("UKdriversKSI.txt"))

colnames(data) <- "logUKdriversKSI"
data <- ts(data, start = 1969,frequency=12)

#dV variance of the observation noise
#dW diagonal elements of the variance matrix of the system noise
#m0 and C0 mean and variance of the initial (presample) state
#F matrix in the observation equation
#G matrix in the state equation

#Stochastic level
fn <- function(params){
  dlmModPoly(order=1, dV=exp(params[1]), dW=exp(params[2]))
}

#ML
fit <- dlmMLE(data.UKdriversKSI, rep(0, 2), fn)
temp <- dlmModPoly(1)+ dlmModSeas(12)
mod <- fn(fit$par)
(obs.error.var <- V(mod))
(state.error.var <- W(mod))

#Gibbs samler

source("dlmGibbsDIGt.R")
MCMC <- 12000
burn <- 1:6000
gibbsOut <- dlmGibbsDIG(data.UKdriversKSI, mod=dlmModPoly(1), a.y=1, b.y=1000, a.theta=10, b.theta=1000, 
                        n.sample = MCMC, thin=1, save.states = FALSE)

gibbsOut <- dlmGibbsDIGt(data.UKdriversKSI, mod=dlmModPoly(1), A_y=1, B_y=1000, A_theta=10, B_theta=1000, 
                        n.sample = MCMC, thin=1, save.states = FALSE)

model <- dlmModPoly(1)+ dlmModSeas(12)
gibbsOut <- dlmGibbsDIGt(data.UKdriversKSI, mod=model, A_y=1, B_y=1000, 
                         n.sample = MCMC, thin=1, save.states = FALSE)


mcmcMean(cbind(gibbsOut$dV[-burn], gibbsOut$dW[-burn,]))



burn <- 1 : 500
#nuRange <- c(1 : 10, seq(20, 100, by = 10))
omega_y <- ts(colMeans(gibbsOut$omega_y[-burn, ]),
              start = start(data.UKdriversKSI), freq=12)
omega_theta <- ts(apply(gibbsOut$omega_theta[,, -burn], 1 : 2,
                        mean), start = start(data.UKdriversKSI), freq = 12)

par(mfrow=c(3,1))
#layout(matrix(c(1,2,3,4), 12,1, TRUE))

#par(mar = c(3, 5, 1, 1) + 0.1)
plot(omega_y, type = "p", ylim = c(0, 1.2), pch = 16,
     xlab = "", ylab = expression(omega[list(data.UKdriversKSI, t)]), cex.lab = 1.6)
abline(h = 1, lty = "dashed")
for (i in 1 : 2)
{ 
  plot(omega_theta[,i], ylim=c(0,1.2), pch = 16,
       type = "p", xlab = "",
       ylab = bquote(omega[list(theta, t * .(i))]), cex.lab = 1.6)
  abline(h = 1, lty = "dashed")
}



#Stochastic level and slope
rm( fit, mod)

fn <- function(params){
  dlmModPoly(dV=exp(params[1]), dW=exp(params[2:3]), m0=rep(0,2), C0=1e+07*diag(2))
}

fit <- dlmMLE(data, rep(0, 3), fn)

mod <- fn(fit$par)
(obs.error.var <- V(mod))
(state.error.var <- W(mod))

#Stochastic level and deterministic slope
rm( fit, mod)

fn <- function(params){
  dlmModPoly(dV=exp(params[1]), dW=c(exp(params[2]), 0), m0=rep(0,2), C0=1e+07*diag(2))
}

fit <- dlmMLE(data, rep(0, 3), fn)

mod <- fn(fit$par)
(obs.error.var <- V(mod))
(state.error.var <- W(mod))

#Stochastic level and seasonal
rm( fit, mod)

fn <- function(params){
  mod <- dlmModPoly(order=1) + dlmModSeas(frequency=12)
  V(mod) <- exp(params[1])
  diag(W(mod))[1:2] <- exp(params[2:3])
  return(mod)
}

fit <- dlmMLE(data, rep(0,3), fn)

mod <- fn(fit$par)
(obs.error.var <- V(mod))
(state.error.var <- W(mod))

#Stochastic level and seasonal
rm( fit, mod)

fn <- function(params){
  mod <- dlmModPoly(order=1) + dlmModSeas(frequency=12)
  V(mod) <- exp(params[1])
  diag(W(mod))[1:2] <- exp(params[2:3])
  return(mod)
}

fit <- dlmMLE(data.UKdriversKSI, rep(0,3), fn)

mod <- fn(fit$par)
(obs.error.var <- V(mod))
(state.error.var <- W(mod))

filtered <- dlmFilter(data.UKdriversKSI,mod) 
res <- residuals(filtered, sd=FALSE) 
plot(res)
length(res)

#Stochastic level and deterministic seasonal
rm( fit, mod)

fn <- function(params){
  mod <- dlmModPoly(order=1, dV=exp(params[1]), dW=exp(params[2])) + 
         dlmModSeas(frequency=12, dV=exp(params[1]), dW=rep(0,11))
  return(mod)
}

fit <- dlmMLE(data, rep(0,2), fn)

mod <- fn(fit$par)
(obs.error.var <- V(mod))
(state.error.var <- W(mod))

rm( fit, mod)
fn <- function(params){
  mod <- dlmModPoly(order=1) + dlmModSeas(frequency=12)
  V(mod) <- exp(params[1])
  diag(W(mod))[1] <- exp(params[2])
  diag(W(mod))[2:12] <- rep(0,11)
  return(mod)
}

fit <- dlmMLE(data, rep(0,2), fn)

mod <- fn(fit$par)
(obs.error.var <- V(mod))
(state.error.var <- W(mod))



rm( fit, mod)
fn <- function(params){
  mod <- dlmModPoly(order=1) + dlmModSeas(frequency=12)
  V(mod) <- exp(params[1])
  diag(W(mod))[1] <- exp(params[2])
  diag(W(mod))[2] <- 0
  return(mod)
}

fit <- dlmMLE(data, rep(0,2), fn)

mod <- fn(fit$par)
(obs.error.var <- V(mod))
(state.error.var <- W(mod))

u <- rnorm(25)
myMod <- dlmModReg(u, dV=14.5)
myMod$JFF
JFF(myMod)
myMod$X

buildFun <- function(x) {
  m <- dlmModPoly(1, dV = exp(x[1]))
  m$JW<-matrix(1)
  m$X <- matrix(exp(x[2]), nc = 1, nr = length(Nile))
  j <- which(time(Nile) == 1899)
  m$X[j,1] <- m$X[j,1] * (1 + exp(x[3]))
  return(m)
}
fit <- dlmMLE(Nile, parm = c(0,0,0), build = buildFun)
fit$conv

dlmNileJump <- buildFun(fit$par)
V(dlmNileJump)
W(dlmNileJump)
X(dlmNileJump)
dlmNileJump$X[c(1, which(time(Nile) == 1889)), 1]
#9.4 An illustration of multivariate state space analysis####

#A) Model without rank restriction####
#Version 2####

rm(list = setdiff(ls(), lsf.str())) 


UKfrontrearseatKSI <- read.table("UKfrontrearseatKSI.txt", head=TRUE)  
data <- ts(UKfrontrearseatKSI[,2:3], start = 1969,frequency=12) #Two time series: front seat passengers killed or seriously injured and rear seat passengers killed or seriously injured

seatbeltLaw <- rep(c(0, 1), times=c(169, 23))  # Inroduction of the seat belt law  (intervention variable)
X <- cbind(seatbeltLaw, UKfrontrearseatKSI[, c(5,4)])
colnames(X)[2:3] <- c("petrolPrice", "kmDriven" )
X <- ts(X, start = 1969,frequency=12)

model <- SSModel(data ~ kmDriven + petrolPrice + seatbeltLaw + 
                   SSMtrend(degree=1, Q=list(matrix(NA, 2, 2))) +
                   SSMseasonal(period=12, sea.type='dummy'), H = matrix(NA, 2, 2), data=X)


fit <- fitSSM(inits = rep(0, 6), model = model, method = "L-BFGS-B") # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"

outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

(H <- fit$model$H) 
(Q <- fit$model$Q)

