 
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
    `colnames<-`(c("Log.likelihood", "Initial.value"))
  #set.seed(123)
  cat("Loop: ")
  for (j in 1:maxLoop){
    cat(paste(j, " "))
    x <- runif(1, min = 0.00001, max = 2) %>% round(3)
    fit <- fitSSM(inits = log(rep(x, w_)), model = model_, updatefn = updatefn_, method = method)
    maxLik <- (logLik(fit$model, method = method)/n) %>% round(7)
    #results[j, ] <- c(round(maxLik, 7), x)
    results[j, ] <- c(maxLik, x)
  }     
  cat("\n")
  results %>% arrange(desc(Log.likelihood), Initial.value) %>% print()
  return(results[1,2])
}
#desc(Log.likelihood)

#Function to find best initial values for optim ver. 2
initValOpt2 <- function(formula = "log(rep(x, 3))", model_ = model, updatefn_ = ownupdatefn, method = "Nelder-Mead", maxLoop = 100){
  results  <- matrix(NA, maxLoop, 2) %>% 
    data.frame() %>%
    `colnames<-`(c("Log.likelihood", "Initial.value"))
  #set.seed(123)
  cat("Loop: ")
  for (j in 1:maxLoop){
    cat(paste(j, ""))
    x <- runif(1, min = 0.00001, max = 2) %>% round(3)
    fit <- fitSSM(inits = eval(parse(text = formula)), model = model_, updatefn = updatefn_, method = method)
    maxLik <- (logLik(fit$model, method = method)/n) %>% round(7)
    results[j, ] <- c(maxLik, x)
  }     
  cat("\n")
  results %>% arrange(desc(Log.likelihood), Initial.value) %>% print()
  return(results[1,2])
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
abline(coef, col = "blue", lwd  = 2, lty = 2)
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
(maxLik <- logLik(fit$model, method = "BFGS")/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states (level)
level <- coef(outKFS)

#Initial value of smoothed level
(coef(outKFS)[1])

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
(maxLik <- logLik(fit$model, method = "BFGS")/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states (level)
level <- coef(outKFS)

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
(maxLik <- logLik(fit$model, method = "BFGS")/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states (level)
level <- coef(outKFS)

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
smoothEstStat <- coef(outKFS)

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
smoothEstStat <- coef(outKFS)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

# Trend (stochastic level + slope)
trend <- coef(outKFS, states = "level") + coef(outKFS, states = "slope")
signal(outKFS, states = "all")$signal

#Figure 3.1. Trend of stochastic linear trend model
plot(dataUKdriversKSI , xlab = "", ylab = "", lty = 1)
lines(trend, lty = 3)
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
#x <- initValOpt()
fit <- fitSSM(model, inits = log(c(0.021, 0.021, 0.021)), updatefn = ownupdatefn, method = "Nelder-Mead")
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
#x <- initValOpt()
fit <- fitSSM(model, inits = log(c(0.059, 0.059)), updatefn = ownupdatefn, method = "L-BFGS-B")
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
ownupdatefn <- function(pars, model){
  model$H[,,1] <- exp(pars[1])
  model
}
d <- q <- 12 #Number of diffuse initial values in the state 
w <- 1 #Number of estimated hyperparameters (i.e. disturbance variances)
l <- 12 #Autocorrelation at lag l to be provided by r-statistic / ACF function
k <- 15#First k autocorrelations to be used in Q-statistic
n <- 192 #Number of observations

#Fitting model and getting output
#x <- initValOpt(method = "BFGS")
fit <- fitSSM(inits = log(0.006), model = model, updatefn = ownupdatefn, method = "BFGS")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
outKFS$model$Q

#Maximum likelihood 
(maxLik <- logLik(fit$model, method = "Nelder-Mead")/n)

#Akaike information criterion (AIC)
(AIC <- (-2*logLik(fit$model)+2*(w+q))/n)

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H) 

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states (level and seasonal components)
smoothEstStat <- coef(outKFS)

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
smoothEstStat <- coef(outKFS)

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
#x <- initValOpt() #Finding best initial values for optim
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
smoothEstStat <- coef(outKFS)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])


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
fit <- fitSSM(inits = log(rep(0.004, w)), model = model, updatefn = ownupdatefn, method = "Nelder-Mead")
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
smoothEstStat <- coef(outKFS)

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
#x <- initValOpt(method = "BFGS") #Finding best initial values for optim
fit <- fitSSM(model, inits = log(rep(0.004, w)), updatefn = ownupdatefn, method = "BFGS")
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
smoothEstStat <- coef(outKFS)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])


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
#x <- initValOpt(method = "BFGS") #Finding best initial values for optim
fit <- fitSSM(model, inits = log(rep(0.008, w)), updatefn = ownupdatefn, method = "BFGS")
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
smoothEstStat <- coef(outKFS)

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
x <- initValOpt() #Finding best initial values for optim
fit <- fitSSM(model, inits = log(rep(0.015, w)), updatefn = ownupdatefn, method = "Nelder-Mead")
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
smoothEstStat <- coef(outKFS)

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
x <- initValOpt(method = "BFGS") #Finding best initial values for optim
fit <- fitSSM(model, inits = log(rep(0.009, w)), updatefn = ownupdatefn, method = "BFGS")
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
smoothEstStat <- coef(outKFS)

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
x <- initValOpt() #Finding best initial values for optim
fit <- fitSSM(model, inits = log(rep(0.022, w)), updatefn = ownupdatefn, method = "Nelder-Mead")
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
smoothEstStat <- coef(outKFS)

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
x <- initValOpt(method = "BFGS") #Finding best initial values for optim
fit <- fitSSM(model, inits = log(rep(0.003, w)), updatefn = ownupdatefn, method = "BFGS")
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
smoothEstStat <- coef(outKFS)

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
fit <- fitSSM(model, inits = log(rep(0.025, w)), updatefn = ownupdatefn, method = "Nelder-Mead")
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
smoothEstStat <- coef(outKFS)

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
fit <- fitSSM(model, inits = log(rep(0.036, w)), updatefn = ownupdatefn, method = "Nelder-Mead")
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
smoothEstStat <- coef(outKFS)

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
smoothEstStat <- coef(outKFS)

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

w <- 2 #Number of estimated hyperparameters (i.e. disturbance variances)
#n <- 192 #Number of observations

#Fitting model
fit <- fitSSM(inits = log(rep(0.009, w)), model = model, updatefn = ownupdatefn, method = "Nelder-Mead")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Level estimation error variance
levEstErVar <- outKFS$V[1, 1, ] %>% ts(start = 1969, frequency=12)

#Figure 8.1. Level estimation error variance for stochastic level and deterministic seasonal model applied to the log of UK drivers KSI
plot(levEstErVar, xlab = "", ylab = "", lty = 1)
title(main = "Figure 8.1. Level estimation error variance for stochastic level and deterministic seasonal model \n applied to the log of UK drivers KSI", 
      cex.main = 0.8)
legend("topright",leg = "level estimation error variance", 
       cex = 0.5, lty = 1, horiz = T)

#Stochastic level and sesoanl with their confidence intervals 
outPredictLev <- predict(fit$model, states = "level", interval = "confidence", se.fit = TRUE, level = 0.90)
outPredictSeas <- predict(fit$model, states = "seasonal", interval = "confidence", level = 0.90)

#Figure 8.2. Stochastic level and its 90% confidence interval for stochastic level and deterministic seasonal model applied to the log of UK drivers KSI
plot(dataUKdriversKSI, xlab = "", ylab = "", lty = 1)
lines(outPredictLev[, "fit"], lty = 3)
lines(outPredictLev[, "lwr"], lty = 3)
lines(outPredictLev[, "upr"], lty = 3)
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
 
#8.4 Filtering and prediction####

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

#Fitting model and getting output
fit <- fitSSM(model, inits = log(c(0.01, 0.01)), updatefn = ownupdatefn, method = "Nelder-Mead")
outKFS <- KFS(fit$model, filtering = "state", smoothing = "state")

#Filtered and Smoothed estimates
sigSmooth <- signal(outKFS, states = "all", filtered = FALSE)$signal
sigFilt <- signal(outKFS, states = "all", filtered = TRUE)$signal %>% window(start = c(1970, 2))

#Figure 8.5. Smoothed and filtered state of the local level model applied to Norwegian road traffic fatalities
plot(sigSmooth, xlab = "", ylab = "", lty = 1, ylim = c(5.6, 6.4))
lines(sigFilt, lty = 3)
title(main = "Figure 8.5. Smoothed and filtered state of the local level model applied to Norwegian \nroad traffic fatalities", 
      cex.main = 0.8)
legend("topright",leg = c("smoothed level", "filtered level"),
       cex = 0.6, lty = c(1, 3), horiz = T)

#One ahead prediction errors (non-standardised) with their variances 
predResid <- residuals(outKFS) %>% window(start = c(1970, 2))
predErVar <- outKFS$F[1,-1] %>% ts(start = c(1970, 2), frequency = 1)

#Figure 8.7 One-step-ahead prediction errors (top) and their variances (bottom) for the 
#local level model applied to Norwegian road traffic fatalities
# Check y scale in the first figure not the same as in the book!!!
par(mfrow = c(2, 1), mar = c(2, 2, 2, 2))
plot(predResid, xlab = "", ylab = "", lty = 3)
abline(h = 0, lty = 1)
title(main = "Figure 8.7 One-step-ahead prediction errors (top) and their variances (bottom) for \nthe local level model applied to Norwegian road traffic fatalities", 
      cex.main = 0.8)
legend("topright",leg = "predictions errors", cex = 0.6, lty = 3, horiz = T)
plot(predErVar, xlab = "", ylab = "", lty = 3)
legend("topright",leg = "prediction error variance", 
       cex = 0.6, lty = 3, horiz = T)
par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1))


#8.5 Diagnostic tests####

#A) Tests for Section 7.3####
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

#Fitting model
w <- 2#Number of estimated hyperparameters (i.e. disturbance variances)
fit <- fitSSM(model, inits = log(rep(0.036, w)), updatefn = ownupdatefn, method = "Nelder-Mead")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Residuals
predResid <- rstandard(outKFS) #One-step prediction residuals

#Figure 8.8 Standardised one-step prediction errors of model in Section 7.3
plot(predResid, xlab = "", ylab = "", lty = 3, ylim = c(-2, 2.5))
abline(h = 0, lty = 1)
title(main = "Figure 8.8 Standardised one-step prediction errors of model in Section 7.3", 
      cex.main = 0.8)
legend("topright",leg = "standardised one-step predictions errors", cex = 0.6, lty = 3, horiz = T)

#Figure 8.9. Correlogram of standardised one-step prediction errors in Figure 8.8, first 10 lags
Acf(predResid, 10, main = "", ylab = "")
title(main="Figure 8.9. Correlogram of standardised one-step prediction errors in Figure 8.8, first 10 lags", 
      cex.main=0.8)
legend("topright",leg = "ACF - standardised one-step predicition errors",cex = 0.5,lty = 1, col = "black",horiz = T)

#Figure 8.10. Histogram of standardised one-step prediction errors in Figure 8.8
hist(predResid, 
     breaks = seq(-3.5, 3.5, 0.5), 
     main = "Figure 8.10. Histogram of standardised one-step prediction errors in Figure 8.8", 
     cex.main = 0.8,
     xlab = "", ylab = "",
     prob = TRUE)
legend("topleft",leg = "N(s = 1)",cex = 0.6, lty = 3, horiz = T)
predResidNum <- predResid %>% as.numeric() %>% na.omit()
curve(dnorm(x, mean=mean(predResidNum), sd=sd(predResidNum)), add=TRUE, lty = 3) 

#B) Tests for Section 4.3####
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

#Fitting model
w <- 2 #Number of estimated hyperparameters (i.e. disturbance variances)
fit <- fitSSM(inits = log(rep(0.936, w)), model = model, updatefn = ownupdatefn, method = "Nelder-Mead")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Residuals
irregResid <- rstandard(outKFS, "pearson") #Auxiliary irregular  residuals (standardised)
levelResid <- rstandard(outKFS, "state")[, 1] #Auxiliary level  residuals (standardised)

#Figure 8.11. Standardised smoothed level disturbances (top) and standardised
#smoothed observation disturbances (bottom) for analysis of UK drivers KSI
#in Section 4.3
par(mfrow = c(2, 1), mar = c(2, 2, 2, 2))
plot(levelResid, xlab = "", ylab = "", lty = 3, ylim = c(-4, 3))
abline(h = 0, lty = 1)
abline(h = c(-2, 2), lty = 2)
title(main = "Figure 8.11. Standardised smoothed level disturbances (top) and standardised smoothed observation disturbances \n(bottom) for analysis of UK drivers KSI in Section 4.3", 
      cex.main = 0.8)
legend("topleft",leg = "Structural level break t-test", cex = 0.6, lty = 3, horiz = T)
plot(irregResid, xlab = "", ylab = "", lty = 3)
abline(h = 0, lty = 1)
abline(h = c(-2, 2), lty = 2)
legend("topleft",leg = "Outlier t-test", 
       cex = 0.6, lty = 3, horiz = T)
par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1))

#C) Tests for Section 7.3####

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

#Fitting model
w <- 2#Number of estimated hyperparameters (i.e. disturbance variances)
fit <- fitSSM(model, inits = log(rep(0.036, w)), updatefn = ownupdatefn, method = "Nelder-Mead")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Residuals
irregResid <- rstandard(outKFS, "pearson") #Auxiliary irregular  residuals (standardised)
levelResid <- rstandard(outKFS, "state")[, "level"] #Auxiliary level  residuals (standardised)

#Figure 8.11. Standardised smoothed level disturbances (top) and standardised
#smoothed observation disturbances (bottom) for analysis of UK drivers KSI
#in Section 7.3
par(mfrow = c(2, 1), mar = c(2, 2, 2, 2))
plot(levelResid, xlab = "", ylab = "", lty = 3, ylim = c(-4, 3))
abline(h = 0, lty = 1)
abline(h = c(-2, 2), lty = 2)
title(main = "Figure 8.11. Standardised smoothed level disturbances (top) and standardised smoothed observation disturbances \n(bottom) for analysis of UK drivers KSI in Section 7.3", 
      cex.main = 0.8)
legend("topright",leg = "Structural level break t-test", cex = 0.6, lty = 3, horiz = T)
plot(irregResid, xlab = "", ylab = "", lty = 3)
abline(h = 0, lty = 1)
abline(h = c(-2, 2), lty = 2)
legend("topright",leg = "Outlier t-test", 
       cex = 0.6, lty = 3, horiz = T)
par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1))


#8.6. Forecasting####

#A) Norwegian fatalities####

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

#Fitting model and getting output
fit <- fitSSM(model, inits = log(c(0.01, 0.01)), updatefn = ownupdatefn, method = "Nelder-Mead")
outKFS <- KFS(fit$model, filtering = "state", smoothing = "state")

#Forecast
outPredict1 <- predict(fit$model, states = "all", interval = "confidence", 
                      level = 0.90, filtered = TRUE) %>% window(start = 1971)
outPredict2 <- predict(fit$model, states = "all", interval = "confidence", 
                      level = 0.90, n.ahead = 5)
outPredict <- rbind(outPredict1, outPredict2) %>% 
  ts(start = start(outPredict1), frequency = frequency(outPredict1))

#Figure 8.13 Filtered level and five year forecast for Norwegian fatalities, 
#including theri 90% confidence interval
ts.plot(outPredict, lty = c(3, 2, 2), xlab = "")
lines(dataNOfatalities)
title(main = "Figure 8.13 Filtered level and five year forecast for Norwegian fatalities, including their 90% confidence interval", 
      cex.main = 0.8)
legend("topright",leg = c("log fatalities in Norway", "filtered level and forecasts"), cex = 0.6, lty = c(1, 3), horiz = T)

#B) Finnish fatalities####

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

#Fitting model and getting output
fit <- fitSSM(model, inits = log(c(0.059, 0.059)), updatefn = ownupdatefn, method = "L-BFGS-B")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Forecast
outPredict1 <- predict(fit$model, states = "all", interval = "confidence", 
                       level = 0.90, filtered = TRUE) %>% window(start = 1972)
outPredict2 <- predict(fit$model, states = "all", interval = "confidence", 
                       level = 0.90, n.ahead = 5)
outPredict <- rbind(outPredict1, outPredict2) %>% 
  ts(start = start(outPredict1), frequency = frequency(outPredict1))

#Figure 8.14 Filtered level and five year forecast for Finnish fatalities, 
#including theri 90% confidence interval
ts.plot(outPredict, lty = c(3, 2, 2), xlab = "")
lines(dataFIfatalities)
title(main = "Figure 8.14 Filtered level and five year forecast for Finnish fatalities, including their 90% confidence interval", 
      cex.main = 0.8)
legend("topright",leg = c("log fatalities in Finland", "filtered level and forecasts"), cex = 0.6, lty = c(1, 3), horiz = T)

#C) UK fatalities####
# C1) UK fatalities up to February 1983 - stochastic level + stochastic seasonal + explanatory variable####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI169 <- log(read.table("UKdriversKSI.txt")) %>% 
  ts(start = 1969, frequency=12) %>% 
  window(end = c(1983, 1))
petrolPrices169 <- read.table("logUKpetrolprice.txt")[1:169, 1]  #Explanatory variable

#Defining model
model169 <- SSModel(dataUKdriversKSI169 ~ petrolPrices169 + SSMtrend(1, Q = list(matrix(NA))) + SSMseasonal(12, sea.type ='dummy', Q = matrix(NA)),  H = matrix(NA))

ownupdatefn <- function(pars, model){
  model$H[,, 1] <- exp(pars[1])
  diag(model$Q[,, 1]) <- c(exp(pars[2:3]))
  model
}

w <- 3 #Number of estimated hyperparameters (i.e. disturbance variances)
q <- 13 #Number of the elements of the initial state vector with  difuse initialization (exact or non-exact)
n <- 169 #Number of observations

#Fitting model
#x <- initValOpt() #Finding best initial values for optim
fit169 <- fitSSM(model169, inits = log(rep(0.029, w)), updatefn = ownupdatefn, method = "Nelder-Mead")

#Maximum likelihood 
(maxLik169 <- logLik(fit169$model)/n)

#Akaike information criterion (AIC)
(AIC169 <- (-2*logLik(fit169$model)+2*(w+q))/n)


# C2) UK fatalities up to February 1983 - stochastic level + deterministic seasonal + explanatory variable####

#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data
dataUKdriversKSI169 <- log(read.table("UKdriversKSI.txt")) %>% 
  ts(start = 1969, frequency=12) %>% 
  window(end = c(1983, 1))
petrolPrices169 <- read.table("logUKpetrolprice.txt")[1:169, 1]  #Explanatory variable

#Defining model
model169 <- SSModel(dataUKdriversKSI169 ~ petrolPrices169 + SSMtrend(1, Q = list(matrix(NA))) + SSMseasonal(12, sea.type ='dummy', Q = matrix(0)),  H = matrix(NA))

ownupdatefn <- function(pars, model){
  model$H[,, 1] <- exp(pars[1])
  diag(model$Q[,, 1]) <- c(exp(pars[2]), 0)
  model
}

w <- 2 #Number of estimated hyperparameters (i.e. disturbance variances)
q <- 13 #Number of the elements of the initial state vector with  difuse initialization (exact or non-exact)
n <- 169 #Number of observations

#Fitting model
#x <- initValOpt() #Finding best initial values for optim
fit169 <- fitSSM(model169, inits = log(rep(0.004, w)), updatefn = ownupdatefn, method = "Nelder-Mead")

#Maximum likelihood 
(maxLik169 <- logLik(fit169$model)/n)

#Akaike information criterion (AIC)
(AIC169 <- (-2*logLik(fit169$model)+2*(w+q))/n)

dataUKdriversKSI23 <- rep(NA, 23) %>% 
  ts(start = c(1983, 2), frequency=12)
petrolPrices23 <- read.table("logUKpetrolprice.txt")[170:192, 1]  #Explanatory variable

newData23 <- SSModel(dataUKdriversKSI23 ~ petrolPrices23 + SSMtrend(1, Q = list(matrix(fit169$model$Q[1, 1, 1]))) + SSMseasonal(12, sea.type ='dummy', Q = matrix(0)),  H = matrix(fit169$model$H[1, 1, 1]))


outPredict1B <- predict(fit169$model, newdata = newData23, states = "all", interval = "prediction", 
                         level = 0.90) #Authors seem to use prediction interval in the book

#Figure 8.15 Forecasts for t=170,..., 192 including their 90% confidence interval 
plot(outPredict1B[, "fit"], lty = 1, xlab = "", ylab = "", ylim = c(7.1, 7.8), xaxt = "n", xlim = c(1983, 1985))
lines(outPredict1B[, 2], lty = 3)
lines(outPredict1B[, 3], lt  = 3)
title(main = "Figure 8.14 Filtered level and five year forecast for Finnish fatalities, including their 90% confidence interval", 
      cex.main = 0.8)
legend("topright",leg = "forecasts +/- 1.64SE", cex = 0.6, lty = 1, horiz = T)
axis(1, c("1983", "1984", "1985"))

outPredict1A <- predict(fit169$model, states = "all", interval = "prediction", 
                        level = 0.90, filtered = FALSE)
outPredict1 <- rbind(outPredict1A, outPredict1B) %>% 
  ts(start = start(outPredict1A), frequency = frequency(outPredict1A)) %>%
  window(start = c(1981, 12))


#Loading data
dataUKdriversKSI <- log(read.table("UKdriversKSI.txt")) %>% 
  ts(start = 1969, frequency=12)
seatbeltLaw <- as.numeric(rep(c(0, 1), times=c(169, 23)))  #Intervention variable
ptrolPrices <- read.table("logUKpetrolprice.txt")[, 1]  #Explanatory variable

#Defining model
model <- SSModel(dataUKdriversKSI ~ petrolPrices + SSMregression(~ seatbeltLaw, P1 = 100, P1inf = 0) + SSMtrend(1, Q = list(matrix(NA))) + SSMseasonal(12, sea.type ='dummy', Q = matrix(NA)),  H = matrix(NA))

ownupdatefn <- function(pars, model){
  model$H[,, 1] <- exp(pars[1])
  diag(model$Q[,, 1]) <- c(exp(pars[2]), 0)
  model
}

w <- 2 #Number of estimated hyperparameters (i.e. disturbance variances)
n <- 169 #Number of observations
#Fitting model
x <- initValOpt() #Finding best initial values for optim
fit <- fitSSM(model, inits = log(rep(0.016, w)), updatefn = ownupdatefn, method = "Nelder-Mead")
outPredict2 <- predict(fit$model, states = "all", interval = "prediction", 
                        level = 0.90) %>% window(start = c(1981, 12), filtered = FALSE)
#Figure 8.16 Last four years (1981-1984 in the time series of the log of numbers
#of drivers KSI: observed series, forecasts obtained from the analysis up to February 1983, 
#and modelled development for the complete series including an intervention variable for February 1983",
plot(window(dataUKdriversKSI, start = c(1981, 12)), xaxt = "n", xlab = "", ylab = "")
lines(outPredict1[, 1], lty = 3)
lines(outPredict2[, 1], lty = 2)
title(main = "Figure 8.16 Last four years (1981-1984 in the time series of the log of numbers\n of drivers KSI: observed series, forecasts obtained from the analysis up to February 1983, \n and modelled development for the complete series including an intervention variable for February 1983",
      cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI", "signal + forecast", "signal complete model"), 
       cex = 0.6, lty = c(1, 3, 2), horiz = T)
axis(1, c("1982", "1983", "1984", "1985"))


#8.7 Missing observations####
#Removing all objects except functions
rm(list = setdiff(ls(), lsf.str())) 

#Loading data and treating some observations as missing
dataUKdriversKSI <- read.table("UKdriversKSI.txt")[, 1] %>%
  log() %>%
  replace(c(48:62, 120:140), NA) %>%
  ts(start = 1969, frequency=12)

#Defining model
model <- SSModel(dataUKdriversKSI ~ SSMtrend(degree=1, Q=list(matrix(NA))) + SSMseasonal(12, sea.type='dummy', Q=matrix(0)),  H=matrix(NA))

ownupdatefn <- function(pars,model,...){
  model$H[,,1] <- exp(pars[1])
  diag(model$Q[,,1]) <- c(exp(pars[2]), 0)
  model
}


#Fitting model
w <- 2 #Number of estimated hyperparameters (i.e. disturbance variances)
#n <- 169 #Number of observations
#x <- initValOpt() #Finding best initial values for optim
fit <- fitSSM(inits = log(rep(0.005, w)), model = model, updatefn = ownupdatefn, method = "Nelder-Mead")
outKFS <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))

#Level and seasonal estimation error variances
levEstErVar <- outKFS$V[1, 1, ] %>% ts(start = 1969, frequency=12)
seasEstErVar <- outKFS$V[2, 2, ] %>% ts(start = 1969, frequency=12)

#Figure 8.17. Stochastic level estimation error variance for log drivers KSI
#with observations at t=48,...,62 and t=120,...,140 treated as missing
plot(levEstErVar, xlab = "", ylab = "", lty = 1)
title(main = "Figure 8.17. Stochastic level estimation error variance for log drivers KSI with observations at t=48,...,62  \n and t=120,...,140 treated as missing", 
      cex.main = 0.8)
legend("topleft",leg = "level estimation error variance", 
       cex = 0.6, lty = 1, horiz = T)

#Stochastic level and deterministic seasonal with their confidence intervals 
outPredictLev <- predict(fit$model, states = "level", interval = "confidence", se.fit = TRUE, level = 0.90)
outPredictSeas <- predict(fit$model, states = "seasonal", interval = "confidence", level = 0.90)

#Figure 8.18. Stochastic level and its 90% confidence interval for stochastic level for log drivers KSI
#with observations at t=48,...,62 and t=120,...,140 treated as missing
plot(dataUKdriversKSI, xlab = "", ylab = "", lty = 1)
lines(outPredictLev[, "fit"], lty = 3)
lines(outPredictLev[, "lwr"], lty = 3)
lines(outPredictLev[, "upr"], lty = 3)
title(main = "Figure 8.18. Stochastic level and its 90% confidence interval for stochastic level for log drivers KSI with observations at t=48,...,62 /n and t=120,...,140 treated as missing", 
      cex.main = 0.8)
legend("topright",leg = c("log UK drivers KSI", "stochastic level +/- 1.64SE"), 
       cex = 0.5, lty = c(1, 3), horiz = T)


#Figure 8.19. Seasonal estimation error variance for log drivers KSI
#with observations at t=48,...,62 and t=120,...,140 treated as missing
plot(seasEstErVar, xlab = "", ylab = "", lty = 1)
title(main = "Figure 8.17. Seasonal estimation error variance for log drivers KSI with observations at t=48,...,62  \n and t=120,...,140 treated as missing", 
      cex.main = 0.8)
legend("topleft",leg = "seasonal estimation error variance", 
       cex = 0.6, lty = 1, horiz = T)
#Results do not corrrespond to plot patern in the book. It was consulted with Mr Helske, author of KFAS
#who confirmed that resutls of KFAS correspond to resutls in other R package bssms

#Figure 8.20. Deterministic seasonal and its 90% confidence interval for t=25,..,72
plot(window(outPredictSeas[, "fit"], start = c(1970, 12), end = c(1974, 12)), xlab = "", ylab = "", lty = 3, ylim = c(-0.2, 0.3))
lines(window(outPredictSeas[, "lwr"], start = c(1970, 12), end = c(1974, 12)), lty = 3)
lines(window(outPredictSeas[, "upr"], start = c(1970, 12), end = c(1974, 12)), lty = 3)
abline(h = 0)
title(main = "Figure 8.20. Deterministic seasonal and its 90% confidence interval for t=25,..,72", 
      cex.main = 0.8)
legend("topleft",leg = "deterministic seasonal +/- 1.64SE", 
       cex = 0.6, lty = 3, horiz = T)

#Auxiliary irregular residuals (non-standardised)
irregResid <- residuals(outKFS, "pearson") 

#Figure 8.21. Irregular component
par(mfrow = c(1, 1), mar = c(4, 4, 4, 4))
plot(irregResid  , xlab = "", ylab = "", lty = 2)
abline(h = 0, lty = 1)
title(main = "Figure 8.21. Irregular component", cex.main = 0.8)
legend("topleft",leg = "irregular",cex = 0.6, lty = 2, horiz = T)


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

w <- 6 #Number of estimated hyperparameters (i.e. disturbance variances)
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

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)
 
#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states 
smoothEstStat <- coef(outKFS)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

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

ownupdatefn <- function(pars, model){
  L1T <- diag(c(exp(pars[1]), 0)) # see explanationhttps://mc-stan.org/docs/2_18/reference-manual/covariance-matrices-1.html#fn16
  L1T[upper.tri(L1T)] <- pars[2]
  model["Q"] <- crossprod(L1T) #crossprod (X,Y) = t(X) %*% Y or crossprod (X) = t(X) %*% X
  L2T <- diag(exp(pars[3:4]))
  L2T[upper.tri(L2T)] <- pars[5]
  model["H"] <- crossprod(L2T)
  model
}

w <- 5 #Number of estimated hyperparameters (i.e. disturbance variances)
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

#Maximum likelihood estimate of the irregular variance
(H <- fit$model$H)

#Maximum likelihood estimate of the state disturbance variance 
(Q <- fit$model$Q)

#Smoothed estimates of states 
smoothEstStat <- coef(outKFS)

#Initial values of the smoothed estimates of states
(initSmoothEstStat <- smoothEstStat[1,])

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



