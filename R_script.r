## Clear environment
rm(list=ls())

## Install  packages if not installed

## Load relevant packages
library(tidyr)
library(PortfolioAnalytics)
library(MASS)
library(fGarch)
library(fBasics)
library(KScorrect)
library(ADGofTest)
library(VineCopula)

## Set working directory in original data location.

## Import SSE and CAC data and create log returns
sse_data <- read.csv("SSE_weekly_2000-2018.csv")
cac_data <- read.csv("^FCHI.csv")
headings <- which(sse_data$Adj.Close=='null') # Find and store rows where the Adj.Close value is "null"
sse_data <- sse_data[-headings,] # remove null data in SSE data
cac_data <- cac_data[-headings,] # remove null data in CAC data
  
## Create log returns for SSE and CAC data
op1 <- sse_data$Adj.Close
prices1 <- as.numeric(op1)
n1 <- length(prices1)
lret1 <- log(prices1[-1]/prices1[-n1])

op2 <- cac_data$Adj.Close
prices2 <- as.numeric(op2)
n2 <- length(prices2)
lret2 <- log(prices2[-1]/prices2[-n2])

## Calculate mean and standard deviation
sse_mean <- mean(lret1)
sse_sd <- sd(lret1)
cac_mean <- mean(lret2)
cac_sd <- sd(lret2)

# Calculate the correlation between lret1 and lret2
corr <- cor(lret1,lret2,use='complete.obs') 

# Create an equally-weighted portfolio
w <- c(1/2,1/2)

# Portfolio parameters: calculate the portfolio's variance and mean:
port_vari <- ((w[1])^2)*((sse_sd)^2)+
            ((w[2])^2)*((cac_sd)^2)+
            2*w[1]*w[2]*corr*sse_sd*cac_sd
port_mean <- w[1]*sse_mean+w[2]*cac_mean


## Calculate the 99% and 95% one-week Value-at-Risk (VaR) for the portfolio:
# The VaR values give an estimate of the potential loss that the portfolio could 
# experience in a given time period (one week in this case) 
# at a specific confidence level (99% or 95%).
Port_VaR_99 <- qnorm(0.01)*sqrt(port_vari)-port_mean
Port_VaR_99

Port_VaR_95 <- qnorm(0.05)*sqrt(port_vari)-port_mean
Port_VaR_95


#A negative VaR would imply the portfolio has a high probability of making a profit,
#for example, if a one-week 5% VaR of negative $1 million implies the portfolio has
#a 95% chance of making more than $1 million over the next week. In other words, 
#there is only a 5% chance that the portfolio will lose more than $1 million 
#during that time frame. This interpretation of VaR is useful for 
#understanding the potential upside and downside risks associated with a given portfolio.

####


#### Partb

########
#Step 1: Identification
########

### Check normality
##note: SSE using data without missing values
par(mfrow=c(1,1))
plot(lret1,lret2)
# Perform the Jarque-Bera test for normality on each return series
jarqueberaTest(lret1)
jarqueberaTest(lret2)

# Create a dataframe with portfolio and log returns
##deleting the first date because log return=Rt/Rt-1
port <- data.frame(sse_data$Date[-1],lret1,lret2)
colnames(port) = c("Date","Ret1","Ret2")

### Check variance

# Check variance by finding marginal distributions
# log returns 1 
par(mfrow=c(2,2))
acf(lret1, col="green", lwd=2)
pacf(lret1, col="green", lwd=2)
acf(lret1^2, col="red", lwd=2)
par(mfrow=c(1,1))
# Returns 2
par(mfrow=c(2,2))
acf(lret2, col="green", lwd=2) 
pacf(lret2, col="green", lwd=2)
acf(lret2^2, col="red", lwd=2)
par(mfrow=c(1,1))

# The ACF and PACF plots are created for each return series (lret1 and lret2) and 
# their squared values to investigate potential autocorrelation patterns. 
# We see a spike at lag 0 from the acf graph.

# The spikes in the ACF and PACF plots at different lags provide insights into 
# the time-series properties of the log returns.
# We see 4 spikes exceeding the limits from the pacf graph.

## Plot the time series of log returns
year = as.Date(port$Date)
par(mfrow=c(1,2))
# Plot log return of idx1
plot(year,lret1,xlab="Year",ylab="Log return of idx1", main= "Time series of idx1"
     ,type="l",xaxt = "n")
axis(1,year,format(year,"%Y"),cex.axis=0.5,tick = FALSE)
# Plot log return of idx2
plot(year,lret2,xlab="Year",ylab="Log return of idx2", main= "Time series of idx2"
     ,type="l",xaxt="n")
axis(1,year,format(year,"%Y"),cex.axis=0.5,tick = FALSE)

###votality clusters

########
#Step 2: Model Estimation
########
# Fit GARCH models with different conditional distributions and AR orders
# Find the best fitting model for lret0 and lret2
# Get the final  models

# Set up GARCH models with different conditional distribution
model1_1=garchFit(~arma(0,0)+garch(1,1),data=lret1,trace=F,cond.dist="norm")
model1_2=garchFit(~arma(0,0)+garch(1,1),data=lret1,trace=F,cond.dist="std")
model1_3=garchFit(~arma(0,0)+garch(1,1),data=lret1,trace=F,cond.dist="sstd")
model1_4=garchFit(~arma(0,0)+garch(1,1),data=lret1,trace=F,cond.dist="snorm")
model1_5=garchFit(~arma(0,0)+garch(1,1),data=lret1,trace=F,cond.dist="sged")
model1_6=garchFit(~arma(0,0)+garch(1,1),data=lret1,trace=F,cond.dist="ged")
ic1=rbind(model1_1@fit$ics,model1_2@fit$ics,model1_3@fit$ics,model1_4@fit$ics,model1_5
          @fit$ics,model1_6@fit$ics)
rownames(ic1)=c("norm", "std","sstd", "snorm", "sged", "ged" )
ic1

# Find the model with the lowest BIC for SSE
min_bic1_model <- rownames(ic1)[which.min(ic1[,"BIC"])]
cat("Model with lowest BIC:", min_bic1_model, "\n")
## Student's t-distribution has the lowest BIC, thus we use this model further.

# Find the model with the lowest BIC for CAC
model2_1=garchFit(~arma(4,0)+garch(1,1),data=lret2,trace=F,cond.dist="norm")
model2_2=garchFit(~arma(4,0)+garch(1,1),data=lret2,trace=F,cond.dist="std")
model2_3=garchFit(~arma(4,0)+garch(1,1),data=lret2,trace=F,cond.dist="sstd")
model2_4=garchFit(~arma(4,0)+garch(1,1),data=lret2,trace=F,cond.dist="snorm")
model2_5=garchFit(~arma(4,0)+garch(1,1),data=lret2,trace=F,cond.dist="sged")
model2_6=garchFit(~arma(4,0)+garch(1,1),data=lret2,trace=F,cond.dist="ged")
ic2=rbind(model2_1@fit$ics,model2_2@fit$ics,model2_3@fit$ics,model2_4@fit$ics,model2_5
         @fit$ics,model2_6@fit$ics)
rownames(ic2)=c("norm", "std","sstd", "snorm", "sged", "ged" )
ic2

min_bic2_model <- rownames(ic2)[which.min(ic2[,"BIC"])]
cat("Model with lowest BIC:", min_bic2_model, "\n")
##Skew Student's t-distribution has the lowest BIC, thus we use this model further.

##calculate the AR values
model_a1=garchFit(~arma(0,0)+garch(1,1),data=lret1,trace=F,cond.dist="std")
model_b1=garchFit(~arma(1,0)+garch(1,1),data=lret1,trace=F,cond.dist="std")
model_c1=garchFit(~arma(2,0)+garch(1,1),data=lret1,trace=F,cond.dist="std")
model_d1=garchFit(~arma(3,0)+garch(1,1),data=lret1,trace=F,cond.dist="std")
model_e1=garchFit(~arma(4,0)+garch(1,1),data=lret1,trace=F,cond.dist="std")
model_f1=garchFit(~arma(5,0)+garch(1,1),data=lret1,trace=F,cond.dist="std")
model_g1=garchFit(~arma(10,0)+garch(1,1),data=lret1,trace=F,cond.dist="std")
ic.1=rbind(model_a1@fit$ics,model_b1@fit$ics,model_c1@fit$ics,model_d1@fit$ics,model_e1@fit$ics,
            model_f1@fit$ics,model_g1@fit$ics)
rownames(ic.1)=c("AR(0)","AR(1)","AR(2)","AR(3)","AR(4)","AR(5)","AR(10)")
ic.1

# Find the model with the lowest AIC and BIC
min_aic_model.1 <- rownames(ic.1)[which.min(ic.1[,"AIC"])]
min_bic_model.1 <- rownames(ic.1)[which.min(ic.1[,"BIC"])]
cat("Model with lowest AIC:", min_aic_model.1, "\n")
cat("Model with lowest BIC:", min_bic_model.1, "\n")

## In this case, AIC and BIC do not agree on the best model, we consider other factors.

### Further checking by Ljung-Box test and MSE 
# Load necessary packages
library(forecast)
library(lmtest)
library(rugarch)

# Split the data into training and testing sets
n <- length(lret1)
train_size <- floor(0.8 * n)
train_data <- lret1[1:train_size]
test_data <- lret1[(train_size + 1):n]

# Fit the models on the training set
models <- list(model_a1, model_b1, model_c1, model_d1, model_e1, model_f1, model_g1)

# Initialize variables to store the results
ljung_box_pvalues <- numeric(length(models))
mse <- numeric(length(models))

# Evaluate the models
for (i in 1:length(models)) {
  # Refit the model on the training data
  model_spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                           mean.model = list(armaOrder = c(i - 1, 0), include.mean = TRUE),
                           distribution.model = "std")
  refitted_model <- ugarchfit(spec = model_spec, data = train_data)
  
  # Ljung-Box test on residuals
  ljung_box_pvalues[i] <- Box.test(residuals(refitted_model), type = "Ljung-Box", lag = 10)$p.value
  
  # Forecast and compute the Mean Squared Error (MSE)
  forecasted_values <- ugarchforecast(refitted_model, n.ahead = length(test_data))@forecast$seriesFor
  mse[i] <- mean((test_data - forecasted_values)^2)
}

# Combine the results
model_names <- c("AR(0)", "AR(1)", "AR(2)", "AR(3)", "AR(4)", "AR(5)", "AR(10)")
results <- data.frame(Model = model_names, AIC = ic.1[,"AIC"], BIC = ic.1[,"BIC"],
                      Ljung_Box_PValue = ljung_box_pvalues, MSE = mse)
results
## We want a model with a high p-value from the Ljung-Box test (indicating no significant 
## autocorrelation in the residuals) and a low MSE value (indicating better forecasting accuracy)
## Therefore, we consider the best model based on a combination of AIC, BIC, Ljung-Box p-value, and MSE

# Define a weight for each metric
aic_weight <- 0.25
bic_weight <- 0.25
ljung_box_weight <- 0.25
mse_weight <- 0.25

# Normalize the AIC, BIC, and MSE values (the lower, the better)
norm_aic <- 1 - (results$AIC - min(results$AIC)) / (max(results$AIC) - min(results$AIC))
norm_bic <- 1 - (results$BIC - min(results$BIC)) / (max(results$BIC) - min(results$BIC))
norm_mse <- 1 - (results$MSE - min(results$MSE)) / (max(results$MSE) - min(results$MSE))

# Calculate the combined score for each model
combined_score <- aic_weight * norm_aic +
  bic_weight * norm_bic +
  ljung_box_weight * results$Ljung_Box_PValue +
  mse_weight * norm_mse

# Add the combined score to the results data frame
results$Combined_Score <- combined_score

# Find the best model based on the combined score
best_model_index <- which.max(results$Combined_Score)
best_model1 <- results[best_model_index, ]
best_model1
# Therefore, we choose the AR(4) with std for SSE data.

# Similarly, choosing the best model for CAC data
model_a2=garchFit(~arma(0,0)+garch(1,1),data=lret2,trace=F,cond.dist="sstd")
model_b2=garchFit(~arma(1,0)+garch(1,1),data=lret2,trace=F,cond.dist="sstd")
model_c2=garchFit(~arma(2,0)+garch(1,1),data=lret2,trace=F,cond.dist="sstd")
model_d2=garchFit(~arma(3,0)+garch(1,1),data=lret2,trace=F,cond.dist="sstd")
model_e2=garchFit(~arma(4,0)+garch(1,1),data=lret2,trace=F,cond.dist="sstd")
model_f2=garchFit(~arma(5,0)+garch(1,1),data=lret2,trace=F,cond.dist="sstd")
model_g2=garchFit(~arma(10,0)+garch(1,1),data=lret2,trace=F,cond.dist="sstd")
ic.2=rbind(model_a2@fit$ics,model_b2@fit$ics,model_c2@fit$ics,model_d2@fit$ics,model_e2@fit$ics,
           model_f2@fit$ics,model_g2@fit$ics)
rownames(ic.2)=c("AR(0)","AR(1)","AR(2)","AR(3)","AR(4)","AR(5)","AR(10)")
ic.2

# Find the model with the lowest AIC and BIC
min_aic_model.2 <- rownames(ic.2)[which.min(ic.2[,"AIC"])]
min_bic_model.2 <- rownames(ic.2)[which.min(ic.2[,"BIC"])]
cat("Model with lowest AIC:", min_aic_model.2, "\n")
cat("Model with lowest BIC:", min_bic_model.2, "\n")
##AIC and BIC disagree on the best model, we consider other factors.

### Further checking by Ljung-Box test and MSE 
m <- length(lret1)
train_size2 <- floor(0.8 * m)
train_data2 <- lret2[1:train_size2]
test_data2 <- lret2[(train_size2 + 1):m]

# Fit the models on the training set
models2 <- list(model_a2, model_b2, model_c2, model_d2, model_e2, model_f2, model_g2)

# Initialize variables to store the results
ljung_box_pvalues2 <- numeric(length(models2))
mse2 <- numeric(length(models2))

# Evaluate the models
for (i in 1:length(models2)) {
  # Refit the model on the training data
  model_spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                           mean.model = list(armaOrder = c(i - 1, 0), include.mean = TRUE),
                           distribution.model = "sstd")
  refitted_model <- ugarchfit(spec = model_spec, data = train_data2)
  
  # Ljung-Box test on residuals
  ljung_box_pvalues2[i] <- Box.test(residuals(refitted_model), type = "Ljung-Box", lag = 10)$p.value
  
  # Forecast and compute the Mean Squared Error (MSE)
  forecasted_values2 <- ugarchforecast(refitted_model, n.ahead = length(test_data2))@forecast$seriesFor
  mse2[i] <- mean((test_data2 - forecasted_values2)^2)
}

# Combine the results
model_names2 <- c("AR(0)", "AR(1)", "AR(2)", "AR(3)", "AR(4)", "AR(5)", "AR(10)")
results2 <- data.frame(Model = model_names2, AIC = ic.2[,"AIC"], BIC = ic.2[,"BIC"],
                      Ljung_Box_PValue = ljung_box_pvalues2, MSE = mse2)
results2
## We want a model with a high p-value from the Ljung-Box test (indicating no significant 
## autocorrelation in the residuals) and a low MSE value (indicating better forecasting accuracy)
## Therefore, we consider the best model based on a combination of AIC, BIC, Ljung-Box p-value, and MSE

# Normalize the AIC, BIC, and MSE values (the lower, the better)
norm_aic2 <- 1 - (results2$AIC - min(results2$AIC)) / (max(results2$AIC) - min(results2$AIC))
norm_bic2 <- 1 - (results2$BIC - min(results2$BIC)) / (max(results2$BIC) - min(results2$BIC))
norm_mse2 <- 1 - (results2$MSE - min(results2$MSE)) / (max(results2$MSE) - min(results2$MSE))

# Calculate the combined score for each model
combined_score2 <- aic_weight * norm_aic2 +
  bic_weight * norm_bic2 +
  ljung_box_weight * results2$Ljung_Box_PValue+
  mse_weight * norm_mse2

# Add the combined score to the results data frame
results2$Combined_Score <- combined_score2

# Find the best model based on the combined score
best_model2_index <- which.max(results2$Combined_Score)
best_model2 <- results2[best_model2_index, ]
best_model2
# Therefore, we choose the AR(2) with sstd for CAC data.


# Final selected models
model_sse = garchFit(formula=~arma(4,0)+garch(1,1),data=lret1,trace=F,cond.dist="std")
model_cac = garchFit(formula=~arma(2,0)+garch(1,1),data=lret2,trace=F,cond.dist="sstd")


## Then, to standardise the residuals for Ljung-Box tests

# Standardize the residuals of the GARCH models
resid_1 <- residuals(model_sse, standardize=TRUE)
resid_2 <- residuals(model_cac, standardize=TRUE)

# Plot ACF for the standardized residuals and their squared values
par(mfrow=c(1,2))
acf(resid_1, col="green", lwd=2)
acf(resid_1^2, col="red", lwd=2)
par(mfrow=c(1,2))
acf(resid_2, col="green", lwd=2)
acf(resid_2^2, col="red", lwd=2)
par(mfrow=c(1,1))

# Perform Ljung-Box test for autocorrelations
Box.test(resid_1, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(resid_1^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(resid_2, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(resid_2^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)
# These tests and plots help assess the goodness of fit of the GARCH models by 
# checking if there is any remaining autocorrelation in the standardized residuals 
# and their squared values. If there is no significant autocorrelation,i.e. p value > 0.05 
# it suggests that the GARCH models have adequately captured the dynamics of the data.

## PIT and plot histograms
## Should look similar to a standard uniform distribution by definition
u1 = pstd(resid_1, mean=0, sd=1,nu=coef(model_sse)['shape'])
u2 = psstd(resid_2, mean=0, sd=1,nu=coef(model_cac)['shape'],xi=coef(model_cac)[7])

par(mfrow=c(1,2))
hist(u1)
hist(u2)

# Perform Kolmogorov-Smirnov and Anderson-Darling tests for uniformity
KStest1 <- LcKS(u1, cdf = "punif")
KStest1$p.value
ADtest1 <- ad.test(u1, null="punif")
ADtest1$p.value

KStest2 <- LcKS(u2, cdf = "punif")
KStest2$p.value
ADtest2 <- ad.test(u2, null="punif")
ADtest2$p.value


## As the p-values indicate are larger than 0.05, both datasets are perfect 
## fits for the uniform distribution.

# Copula modelling
model <- BiCopSelect(u1, u2, familyset=NA, selectioncrit="AIC", indeptest=TRUE, level=0.05,se = TRUE)
model

# Fit AR-GARCH model for SSE and CAC: re-introducing autocorrelation and GARCH effects observed in data
spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(1, 0), include.mean = TRUE))
fit1 <- ugarchfit(spec1, lret1)
spec2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(1, 0), include.mean = TRUE))
fit2 <- ugarchfit(spec2, lret2)

# Perform Monte Carlo simulation using the copula model to generate correlated uniform samples
N <- 10000
set.seed(111)
u_sim <- BiCopSim(N, family = model$family, model$par, model$par2)

# Transform the uniform samples into the respective distributions
res1_sim <- qdist("std", u_sim[, 1], shape = coef(model_sse)['shape']) # Student's t distribution
res2_sim <- qdist("sstd", u_sim[, 2], shape = coef(model_cac)['shape'], skew = coef(model_cac)[7]) # Skewed Student's t distribution

# Simulate AR-GARCH residuals 
y1simulated <- fitted(ugarchsim(fit1, n.sim = N, m.sim = 1, startMethod = "sample"))
y2simulated <- fitted(ugarchsim(fit2, n.sim = N, m.sim = 1, startMethod = "sample"))

# Combine the simulated residuals with the standard normal samples
res1_correlated <- res1_sim * sd(y1simulated) + mean(y1simulated)
res2_correlated <- res2_sim * sd(y2simulated) + mean(y2simulated)

# Compute the portfolio returns
portsim <- matrix(0, nrow = N, ncol = 1)
portsim=log(1+((exp(y1simulated)-1)+(exp(y2simulated)-1))*(1/2))
varsim <- matrix(0, nrow = 1, ncol = 2)
varsim <- quantile(portsim, c(0.01, 0.05))
varsim
                  