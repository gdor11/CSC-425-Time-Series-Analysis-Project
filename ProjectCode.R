library(tseries)
library(fBasics)
library(forecast)
library(lmtest)
library(ggplot2)
library(astsa)
library(fUnitRoots)
library(lubridate)
library(scales) 
library(gridExtra)
### WEEKLY DATA###
setwd("C:\\Users\\guy.dor\\Documents\\CSC 425\\Project")
scv=read.table("WeeklyViolations.csv",header=T, sep='\t') 
head(scv)


x=scv$Violations
head(x)
xts=ts(x,frequency=52,start=c(2014,7), end = c(2017,8,20))
# create time plot
autoplot(xts)+ylab("Speeding Violations")+ggtitle("Weekly Speeding Violations")



#Check Normality
hist(x,main="Weekly Speeding Violations" , probability = TRUE, xlab = "Speeding Violations")
lines(density(x), col="blue", lwd=2)
lines(density(x, adjust=2), lty="dotted", col="darkgreen", lwd=2) 



normalTest(x,method=c("jb"))  

qqnorm(x)
qqline(x, col = 2) 


acf(x)
pacf(x)

# Dickey Fuller Tests

# tests for AR model with time trend
adfTest(xts, lags=1, type=c("ct"))
adfTest(xts, lags=2, type=c("ct"))
adfTest(xts, lags=4, type=c("ct"))
adfTest(xts, lags=17, type=c("ct"))
adfTest(xts, lags=52, type=c("ct"))
# tests for AR model with no time trend
adfTest(x, lags=1, type=c("c"))
adfTest(x, lags=2, type=c("c"))
adfTest(x, lags=4, type=c("c"))
adfTest(x, lags=17, type=c("c"))
adfTest(x, lags=52, type=c("c"))




Box.test(x,lag = 1, type="Ljung-Box", fitdf=1)

Box.test(x,lag = 3, type="Ljung-Box", fitdf=2)

Box.test(x,lag = 15, type="Ljung-Box", fitdf=2)



# APPLYING DIFFERENCING TO DATA

# compute seasonal difference for weeklydata (s=17)
sd17x=diff(x,17)
# create acf plot 
acf(as.vector(sd17x),lag.max=30, main="ACF of 17 Lag Seasonal Differenced Violations")
pacf(as.vector(sd17x),lag.max=30, main="PACF of 17 Lag Seasonal Differenced Violations")
Box.test(sd17x,lag = 1, type="Ljung-Box", fitdf=1)

Box.test(sd17x,lag = 3, type="Ljung-Box", fitdf=2)

Box.test(sd17x,lag = 15, type="Ljung-Box", fitdf=2)



hist(sd17x,main="SDX Lag 17 Speeding Violations" , probability = TRUE, xlab = "SDX Lag 17 Speeding Violations")
lines(density(sd17x), col="blue", lwd=2)
lines(density(sd17x, adjust=2), lty="dotted", col="darkgreen", lwd=2)

qqnorm(sd17x, main="Q-Q Plot on 17 Lag seasonal differenced data")
qqline(sd17x, col = 2)

normalTest(sd17x,method=c("jb"))  

#Unit-root tests on Seasonal difference
sd17xts=diff(xts,17)
plot(sd17xts)
adfTest(coredata(sd17xts), lags=2, type=c("c"))
adfTest(sd17xts, lags=2, type=c("c"))



# compute seasonal difference for weeklydata (s=52)
sdx=diff(x,52)
# create acf plot 
acf(as.vector(sdx),lag.max=30, main="ACF of Seasonal Differenced Violations")
pacf(as.vector(sdx),lag.max=30, main="PACF of Seasonal Differenced Violations")
Box.test(sdx,lag = 1, type="Ljung-Box", fitdf=1)

Box.test(sdx,lag = 3, type="Ljung-Box", fitdf=2)

Box.test(sdx,lag = 15, type="Ljung-Box", fitdf=2)



hist(sdx,main="SDX Speeding Violations" , probability = TRUE, xlab = "SDX Speeding Violations")
lines(density(sdx), col="blue", lwd=2)
lines(density(sdx, adjust=2), lty="dotted", col="darkgreen", lwd=2)

qqnorm(sdx, main="Q-Q Plot on seasonal differenced data")
qqline(sdx, col = 2)

normalTest(sdx,method=c("jb"))  

#Unit-root tests on Seasonal difference
sdxts=diff(xts,52)
plot(sdxts)
adfTest(coredata(sdxts), lags=2, type=c("c"))
adfTest(sdxts, lags=2, type=c("c"))

basicStats(sdx)


# try automated order selection
# since TS is stationary, stationary = FALSE;
auto.arima(xts)
auto.arima(xts,stepwise=FALSE,approx=FALSE)

auto.arima(sdxts)
auto.arima(sdxts,stepwise=FALSE,approx=FALSE)


auto.arima(sd17xts)
auto.arima(sd17xts,stepwise=FALSE,approx=FALSE)




m0=Arima(xts, order=c(3,1,2),seasonal=list(order=c(1,0,0),period=52), method="ML")
m1=Arima(xts, order=c(2,1,2),seasonal=list(order=c(1,0,0),period=52), method="ML")

m00=Arima(sdxts, order=c(2,0,0), method="ML")
m11=Arima(sdxts, order=c(1,0,2), method="ML")
#m11 can be translated to m12
m12=Arima(xts, order=c(1,0,2),seasonal=list(order=c(0,1,0),period=52), include.constant=TRUE, method="ML")


m000=Arima(sd17xts, order=c(3,0,1), method="ML")
m111=Arima(sd17xts, order=c(2,0,1), method="ML")

m11


m0
#Series: xts 
#ARIMA(3,1,2)(1,0,0)[52]                    

#Coefficients:
#          ar1     ar2      ar3      ma1      ma2    sar1
#       -0.0969  0.5034  -0.0383  -0.3780  -0.5682  0.2296
# s.e.   0.1681  0.0976   0.0953   0.1544   0.1506  0.1106

# sigma^2 estimated as 6614425:  log likelihood=-1454.66
# AIC=2923.32   AICc=2924.08   BIC=2944.72


m1
#Series: xts 
#ARIMA(2,1,2)(1,0,0)[52]                    

# Coefficients:
#           ar1     ar2      ma1      ma2    sar1
#       -0.1451  0.5229  -0.3429  -0.6045  0.2287
# s.e.   0.1418  0.0882   0.1423   0.1347  0.1123

# sigma^2 estimated as 6579504:  log likelihood=-1454.74
# AIC=2921.47   AICc=2922.03   BIC=2939.81


m00

#Series: sdxts 
#ARIMA(2,0,0) with non-zero mean 

#Coefficients:
#     ar1     ar2     intercept
#     0.2973  0.0603  -1879.433
#s.e. 0.0988  0.0989    447.624

#sigma^2 estimated as 9125576:  log likelihood=-998.35
#AIC=2004.7   AICc=2005.1   BIC=2015.36

m11

#Series: sdxts 
#ARIMA(1,0,2) with non-zero mean 

#Coefficients:
#       ar1      ma1     ma2      intercept
#       -0.6496  0.9994  0.4017  -1877.6381
#s.e.   0.1933   0.1885  0.0991    410.8637

#sigma^2 estimated as 8832707:  log likelihood=-996.2
#AIC=2002.4   AICc=2003   BIC=2015.72


m12

#Series: xts 
#ARIMA(1,0,2)(0,1,0)[52] with drift         

#Coefficients:
#       ar1       ma1     ma2       drift
#       -0.6496   0.9994  0.4017  -36.1084
#s.e.   0.1933    0.1885  0.0991    7.9012

#sigma^2 estimated as 8833034:  log likelihood=-996.2
#AIC=2002.4   AICc=2003   BIC=2015.72

m000

#Series: sd17xts 
#ARIMA(3,0,1) with non-zero mean 

#Coefficients:
#       ar1     ar2      ar3     ma1      intercept
#       0.1466  0.5249  -0.0675  0.3851   -462.450
#s.e.   0.2687  0.1303   0.1210  0.2593   1049.511

#sigma^2 estimated as 13671345:  log likelihood=-1356.3
#AIC=2724.6   AICc=2725.22   BIC=2742.29


m111


#Series: sd17xts 
#ARIMA(2,0,1) with non-zero mean 

#Coefficients:
#       ar1     ar2     ma1     intercept
#       0.0298  0.5635  0.4837  -460.0259
#s.e.   0.1419  0.0959  0.1573  1092.7025

#sigma^2 estimated as 13600230:  log likelihood=-1356.44
#AIC=2722.88   AICc=2723.32   BIC=2737.62




coeftest(m0)
#z test of coefficients:
  
#        Estimate Std. Error z value  Pr(>|z|)    
#  ar1  -0.096882   0.168069 -0.5764 0.5643166    
#  ar2   0.503396   0.097576  5.1590 2.483e-07 ***
#  ar3  -0.038338   0.095266 -0.4024 0.6873710    
#  ma1  -0.377974   0.154409 -2.4479 0.0143697 *  
#  ma2  -0.568167   0.150600 -3.7727 0.0001615 ***
#  sar1  0.229607   0.110587  2.0762 0.0378711 *  
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


coeftest(m1)

#z test of coefficients:
  
#          Estimate Std. Error z value  Pr(>|z|)    
#   ar1  -0.145111   0.141791 -1.0234   0.30611    
#   ar2   0.522904   0.088216  5.9276 3.074e-09 ***
#   ma1  -0.342894   0.142275 -2.4101   0.01595 *  
#   ma2  -0.604491   0.134743 -4.4863 7.248e-06 ***
#   sar1  0.228663   0.112328  2.0357   0.04178 *  
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



coeftest(m00)



#z test of coefficients:
  
#              Estimate     Std. Error z value  Pr(>|z|)    
#   ar1        2.9728e-01  9.8825e-02  3.0081  0.002628 ** 
#   ar2        6.0310e-02  9.8937e-02  0.6096  0.542138    
#   intercept  -1.8794e+03  4.4762e+02 -4.1987 2.685e-05 ***
#   ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


coeftest(m11)

#z test of coefficients:
  
#               Estimate      Std. Error    z value   Pr(>|z|)    
#   ar1        -6.4960e-01    1.9331e-01   -3.3604    0.0007782 ***
#   ma1         9.9943e-01    1.8851e-01    5.3017    1.147e-07 ***
#   ma2         4.0170e-01    9.9143e-02    4.0517    5.085e-05 ***
#  intercept   -1.8779e+03    4.1087e+02   -4.5705    4.866e-06 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

coeftest(m12)

#z test of coefficients:
  
#               Estimate    Std. Error    z value   Pr(>|z|)    
#   ar1       -0.649604     0.193325    -3.3602     0.0007789 ***
#   ma1       0.999422      0.188526    5.3012      1.150e-07 ***
#   ma2       0.401690      0.099145    4.0516      5.088e-05 ***
#   drift     -36.108429    7.901227    -4.5700     4.878e-06 ***
#  ---
 # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


coeftest(m000)

#z test of coefficients:
  
#                 Estimate    Std. Error    z value   Pr(>|z|)    
#   ar1          0.146600     0.268701      0.5456    0.5854    
#   ar2          0.524904     0.130312      4.0280    5.624e-05 ***
#   ar3         -0.067466     0.121033     -0.5574    0.5772    
#   ma1          0.385086     0.259327      1.4849    0.1376    
#   intercept   -462.450031   1049.510800  -0.4406    0.6595    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


coeftest(m111)

#z test of coefficients:
  
#               Estimate      Std. Error    z value   Pr(>|z|)    
#   ar1          0.029761     0.141869      0.2098    0.833842    
#   ar2          0.563459     0.095853      5.8784    4.143e-09 ***
#   ma1          0.483734     0.157266      3.0759    0.002099 ** 
#   intercept    -460.025917   1092.702550  -0.4210    0.673756    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#M1 and (M11) M12 Pass Coef Tests

acf(m1$resid, main="ARIMA(2,1,2)(1,0,0)[52] residual ACF")
pacf(m1$resid, main="ARIMA(2,1,2)(1,0,0)[52] residual PACF")
#  ljung box test on residuals
Box.test(m1$residuals, 26, "Ljung-Box", fitdf=5)


#Box-Ljung test

#data:  m1$residuals
#X-squared = 21.114, df = 21, p-value = 0.452

Box.test(m1$residuals, 52, "Ljung-Box", fitdf=5)
#Box-Ljung test

#data:  m1$residuals
#X-squared = 46.362, df = 47, p-value = 0.4989




acf(m12$resid, main="ARIMA(1,0,2) with drift residual ACF")
pacf(m12$resid, main="ARIMA(1,0,2) with drift residual PACF")

#  ljung box test on residuals
Box.test(m12$residuals, 26, "Ljung-Box", fitdf=4)

#Box-Ljung test

#data:  m11$residuals
#X-squared = 16.586, df = 22, p-value = 0.7858

Box.test(m12$residuals, 52, "Ljung-Box", fitdf=4)

#Box-Ljung test

#data:  m11$residuals
#X-squared = 48.039, df = 48, p-value = 0.4712







f1=forecast(m1, h=52)

autoplot(f1)+ xlab("Week")+ylab("Weekly Violations")+ggtitle("Forecasts from ARIMA(2,1,2)(1,0,0)[52]")
autoplot(m1)


f11=forecast(m11, h=52)
f12=forecast(m12, h=52)
plot(f11)
plot(f12)
autoplot(f12)+ xlab("Week")+ylab("Weekly Violations")+ggtitle("Forecasts from ARIMA(1,0,2) with drift")
autoplot(m12)




source("backtest.R")


backtest(m1, xts,h=1, orig=length(xts)*0.8)


backtest(m12, xts,h=1, orig=length(xts)*0.8)

#[1] "RMSE of out-of-sample forecasts"
#[1] 2064.699
#[1] "Mean absolute error of out-of-sample forecasts"
#[1] 1611.693
#[1] "Mean Absolute Percentage error"
#[1] 0.09048906
#[1] "Symmetric Mean Absolute Percentage error"
#[1] 0.08853523




#LogNormalize Weekly

lnx=log(x)
lnxts=ts(lnx,frequency=52,start=c(2014,7))


# create time plot
autoplot(lnxts)+ylab("Log Speeding Violations")+ggtitle("Weekly Log Speeding Violations")



#Check Normality
hist(lnx,main="Weekly Log Speeding Violations" , probability = TRUE, xlab = "Log Speeding Violations")
lines(density(lnx), col="blue", lwd=2)
lines(density(lnx, adjust=2), lty="dotted", col="darkgreen", lwd=2) 



normalTest(lnx,method=c("jb"))  

qqnorm(x)
qqline(x, col = 2) 


acf(lnx)
pacf(lnx)

# Dickey Fuller Tests

# tests for AR model with time trend
adfTest(lnxts, lags=26, type=c("ct"))
adfTest(lnxts, lags=52, type=c("ct"))
# tests for AR model with no time trend
adfTest(lnx, lags=26, type=c("c"))
adfTest(lnx, lags=52, type=c("c"))




Box.test(lnx,lag = 1, type="Ljung-Box", fitdf=1)

Box.test(lnx,lag = 3, type="Ljung-Box", fitdf=2)

Box.test(lnx,lag = 15, type="Ljung-Box", fitdf=2)

# APPLYING DIFFERENCING TO DATA
# compute regular differences
dlnx=diff(lnx) 
# create acf plot
acf(as.vector(dlnx),lag.max=32, main="ACF of DLNX Violations")
pacf(as.vector(dlnx),lag.max=32, main="PACF of DLNX Violations")

hist(dlnx,main="DLNX Speeding Violations" , probability = TRUE, xlab = "DLNX Speeding Violations")
lines(density(dlnx), col="blue", lwd=2)
lines(density(dlnx, adjust=2), lty="dotted", col="darkgreen", lwd=2) 
qqnorm(dlnx, main = "Q-Q Plot on Differenced data")
qqline(dlnx, col=2)

normalTest(dlnx,method=c("jb"))  

#Unit-root tests on first difference
dlnxts=diff(lnxts)
plot(dlnxts)
adfTest(coredata(dlnxts), lags=2, type=c("c"))
adfTest(dlnxts, lags=2, type=c("c"))




# APPLYING Seasonal DIFFERENCING TO DATA
# compute regular differences
sdlnx=diff(lnx,52) 
# create acf plot
acf(as.vector(sdlnx),lag.max=32, main="ACF of SDLNX Violations")
pacf(as.vector(sdlnx),lag.max=32, main="PACF of SDLNX Violations")

hist(sdlnx,main="SDLNX Speeding Violations" , probability = TRUE, xlab = "SDLNX Speeding Violations")
lines(density(sdlnx), col="blue", lwd=2)
lines(density(sdlnx, adjust=2), lty="dotted", col="darkgreen", lwd=2) 
qqnorm(sdlnx, main = "Q-Q Plot on Seasonal Differenced data")
qqline(sdlnx, col=2)

normalTest(sdlnx,method=c("jb"))  

#Unit-root tests on first difference
sdlnxts=diff(lnxts,52)
plot(sdlnxts)
adfTest(coredata(sdlnxts), lags=2, type=c("c"))
adfTest(sdlnxts, lags=2, type=c("c"))



#auto.arima(lnxts,stepwise=FALSE,approx=FALSE)


m2=Arima(lnxts, order=c(0,1,5), method="ML")


coeftest(m2)

acf(m2$resid, main="ARIMA(0,1,5) residual ACF")
pacf(m2$resid, main="ARIMA(0,1,5) residual PACF")
#  ljung box test on residuals
Box.test(m2$residuals, 26, "Ljung-Box", fitdf=5)

#Box-Ljung test

#data:  m2$residuals
#X-squared = 20.262, df = 21, p-value = 0.5047

Box.test(m2$residuals, 52, "Ljung-Box", fitdf=5)
#Box-Ljung test

#data:  m2$residuals
#X-squared = 46.36, df = 47, p-value = 0.499

f2=forecast(m2, h=6)

f2$mean=exp(f2$mean)
f2$upper=exp(f2$upper)
f2$lower=exp(f2$lower)
f2$x=exp(f2$x)
plot(f2, main="Forecasts from ARIMA(0,1,5)")
autoplot(f2)+ xlab("Week")+ylab("Speeding Violations")+ggtitle("Forecasts from ARIMA(0,1,5)")

backtest(m2, lnxts,h=1, orig=length(xts)*0.8)

#[1] "RMSE of out-of-sample forecasts"
#[1] 0.1126231
#[1] "Mean absolute error of out-of-sample forecasts"
#[1] 0.07869458
#[1] "Mean Absolute Percentage error"
#[1] 0.008039112
#[1] "Symmetric Mean Absolute Percentage error"
#[1] 0.008005755





















### MONTHLY DATA###
mcv=read.table("MonthlyViolations.csv",header=T, sep=',') 
head(mcv)


x=mcv$Violations
head(x)
xts=ts(x,frequency=12,start=c(2014,7))

# create time plot
plot(xts, main="Monthly Speeding Violations", xlim=c(2014.5,2017.5), ylim=c(5000,150000),  ylab="Monthly Violations")
autoplot(xts)+ylab("Speeding Violations")+ggtitle("Monthly Speeding Violations")


#Check Normality
hist(x,main="Monthly Speeding Violations" , probability = TRUE, xlab = "Speeding Violations")
lines(density(x), col="blue", lwd=2)
lines(density(x, adjust=2), lty="dotted", col="darkgreen", lwd=2) 

normalTest(x,method=c("jb"))  

qqnorm(x)
qqline(x, col = 2) 


acf(x,lag.max=13, main="ACF Montlhy Violations")
pacf(x, main="PACF Montlhy Violations")


# tests for AR model with time trend
adfTest(xts, lags=6, type=c("ct"))
adfTest(xts, lags=12, type=c("ct"))
# tests for AR model with no time trend
adfTest(x, lags=6, type=c("c"))
adfTest(x, lags=12, type=c("c"))


#DX
dx=diff(x)

# create acf plot 
acf(as.vector(dx),lag.max=13, main="ACF Differenced Violations")
pacf(as.vector(dx),lag.max=13, main="PACF Differenced Violations")

hist(dx,main="DX Speeding Violations" , probability = TRUE, xlab = "DX Speeding Violations")
lines(density(dx), col="blue", lwd=2)
lines(density(dx, adjust=2), lty="dotted", col="darkgreen", lwd=2) 

normalTest(dx,method=c("jb"))  

qqnorm(dx, main="Q-Q Plot on differenced data")
qqline(dx, col = 2)

basicStats(dx)

#Unit-root tests on first difference
dxts=diff(xts)
plot(dxts)
adfTest(coredata(dxts), lags=2, type=c("c"))
adfTest(dxts, lags=2, type=c("c"))


# SDX
sdx=diff(x,12)
# create acf plot 
acf(as.vector(sdx),lag.max=13, main="ACF of Seasonal Differenced Violations")
pacf(as.vector(sdx),lag.max=13, main="PACF of Seasonal Differenced Violations")
Box.test(sdx,lag = 1, type="Ljung-Box", fitdf=1)

Box.test(sdx,lag = 3, type="Ljung-Box", fitdf=2)

Box.test(sdx,lag = 12, type="Ljung-Box", fitdf=2)



hist(sdx,main="SDX Speeding Violations" , probability = TRUE, xlab = "SDX Speeding Violations")
lines(density(sdx), col="blue", lwd=2)
lines(density(sdx, adjust=2), lty="dotted", col="darkgreen", lwd=2) 

normalTest(sdx,method=c("jb"))  

qqnorm(sdx, main="Q-Q Plot on seasonal differenced data")
qqline(sdx, col = 2)

basicStats(sdx)

#Unit-root tests on seasonal difference
sdxts=diff(xts,12)
plot(sdxts)
adfTest(coredata(sdxts), lags=2, type=c("c"))
adfTest(sdxts, lags=2, type=c("c"))


# try automated order selection

auto.arima(xts)
auto.arima(xts,stepwise = FALSE, approx=FALSE)
#Force Seasonal Differencing to compare
auto.arima(xts,stepwise = FALSE,D=1, approx=FALSE)

m3=Arima(xts, order=c(0,1,1),seasonal=list(order=c(1,0,0),period=12), method="ML")
#Force Seasonal Differencing to compare
m4=Arima(xts, order=c(0,0,0),seasonal=list(order=c(0,1,0),period=12), method="ML")


coeftest(m3)
coeftest(m4)


acf(m3$resid)
pacf(m3$resid)

acf(m4$resid)
pacf(m4$resid)


f3=forecast(m3, h=12)
f3$mean
plot(f3)
autoplot(f3)+ xlab("Month")+ylab("Speeding Violations")+ggtitle("Forecasts from ARIMA(0,1,1)(1,0,0)[12]")
autoplot(f3, xlim=c(2017.65,2018.6))+ xlab("Month")+ylab("Speeding Violations")+ggtitle("Forecasts from ARIMA(0,1,1)(1,0,0)[12]")

f4=forecast(m4, h=12)
f4$mean
plot(f4)
autoplot(f4)+ xlab("Month")+ylab("Speeding Violations")+ggtitle("Forecasts from ARIMA(0,0,0)(0,1,0)[12]")
autoplot(f4, xlim=c(2017.65,2018.6))+ xlab("Month")+ylab("Speeding Violations")+ggtitle("Forecasts from ARIMA(0,0,0)(0,1,0)[12]")





backtest(m3, xts,h=1, orig=length(xts)*0.8)
#[1] "RMSE of out-of-sample forecasts"
#[1] 6076.338
#[1] "Mean absolute error of out-of-sample forecasts"
#[1] 3900.595
#[1] "Mean Absolute Percentage error"
#[1] 0.04977433
#[1] "Symmetric Mean Absolute Percentage error"
#[1] 0.04742143



backtest(m4, xts,h=1, orig=length(xts)*0.8)

#[1] "RMSE of out-of-sample forecasts"
#[1] 10531.4
#[1] "Mean absolute error of out-of-sample forecasts"
#[1] 9940.429
#[1] "Mean Absolute Percentage error"
#[1] 0.1178775
#[1] "Symmetric Mean Absolute Percentage error"
#[1] 0.1103653





#Violations rate of change
VrOc=diff(x)/x[-length(x)]
VrOc
r2=VrOc^2
absr=abs(VrOc)



# create acf plot 
acf(as.vector(VrOc),lag.max=13, main="ACF of Violations Rate of Change")
pacf(as.vector(VrOc),lag.max=13, main="PACF of Violations Rate of Change")

acf(as.vector(r2),lag.max=13, main="ACF of Violations Rate of Change Squared")
pacf(as.vector(r2),lag.max=13, main="PACF of Violations Rate of Change Squared")

acf(as.vector(absr),lag.max=13, main="ACF of ABS(Violations Rate of Change)")
pacf(as.vector(absr),lag.max=13, main="PACF of ABS(Violations Rate of Change)")






Box.test(VrOc,lag = 1, type="Ljung-Box", fitdf=1)

Box.test(VrOc,lag = 3, type="Ljung-Box", fitdf=2)

Box.test(VrOc,lag = 12, type="Ljung-Box", fitdf=2)



hist(VrOc,main="Speeding Violations Rate of Change" , probability = TRUE, xlab = "VrOc")
lines(density(VrOc), col="blue", lwd=2)
lines(density(VrOc, adjust=2), lty="dotted", col="darkgreen", lwd=2) 

normalTest(VrOc,method=c("jb"))  

qqnorm(VrOc, main="Q-Q Plot on Violations Rate of Change")
qqline(VrOc, col = 2)

basicStats(VrOc)
auto.arima(VrOc)
m5=Arima(VrOc, order=c(0,0,2), method="ML")
m5

coeftest(m5)

acf(m5$resid)
pacf(m5$resid)


f5=forecast(m5, h=12)
f5$mean
plot(f5)
















#### Weather Regression

library(tseries)
library(fBasics)
library(forecast)
library(lmtest)
library(ggplot2)
library(astsa)

scv=read.table("WeeklyViolations+Weather.csv",header=T, sep=",") 
head(scv)


x=scv$Violations
head(x)
xts=ts(x,frequency=52,start=c(2014,7))
xts
pairs(scv[,2:21])


dp=scv$Avg..Dew.Pt.Ave
at=scv$Avg..Temp.Ave
vl=scv$Avg..Visibility.Low
vlts=ts(vl,frequency = 52)
ccf (x,dp)
ccf (x,at)
ccf (x,vl)


ccfvalues1=ccf (x,dp)
ccfvalues2=ccf (x,at)
ccfvalues3=ccf (x,vl)
ccfvalues1
ccfvalues2
ccfvalues3

dx=diff(x) 
sdx=diff(dx,52)
head(sdx)

ddp=diff(dp) 
sddp=diff(ddp,52)
head(sddp)

dat=diff(at) 
sdat=diff(dat,52)
head(sdat)

dvl=diff(vl) 
sdvl=diff(dvl,52)
head(sdvl)



lag2.plot(sdx,sddp,10)
lag2.plot(sdx,sdat,10)
lag2.plot(sdx,sdvl,10)



#Lagged predictors. Test 0, 1, 2 or 3 lags.
weath <- cbind(dp,
                c(NA,dp[1:length(dp)-1]),
                c(NA,NA,dp[1:(length(dp)-2)]),
                c(NA,NA,NA,dp[1:(length(dp)-3)]),
                c(NA,NA,NA,NA,dp[1:(length(dp)-4)]),
                c(NA,NA,NA,NA,NA,dp[1:(length(dp)-5)]),
                c(NA,NA,NA,NA,NA,NA,dp[1:(length(dp)-6)]))
colnames(weath) <- paste("dpLag",0:6,sep="")
head(weath)

dim(weath)
# Choose optimal lag length for advertising based on AIC
# Restrict data so models use same fitting period
fit1 <- auto.arima(xts, xreg=weath[1:158,1], d=0)
fit2 <- auto.arima(xts, xreg=weath[1:158,2], d=0)
fit3 <- auto.arima(xts, xreg=weath[1:158,3], d=0)
fit4 <- auto.arima(xts, xreg=weath[1:158,4], d=0)
fit5 <- auto.arima(xts, xreg=weath[1:158,5], d=0)
fit6 <- auto.arima(xts, xreg=weath[1:158,6], d=0)
fit7 <- auto.arima(xts, xreg=weath[1:158,7], d=0)

fit1
fit2
fit3
fit4
fit5
fit6
fit7

m2=Arima(xts,order=c(1,0,0),xreg=weath[1:158,6],seasonal=list(order=c(1,0,0),period=52), method="ML")
m2
coeftest(m2)
acf(m2$residuals,na.action = na.pass)
pacf(m2$resid,na.action = na.pass)

#  ljung box test on residuals
Box.test(m2$residuals, 26, "Ljung-Box", fitdf=5) 
Box.test(m2$residuals, 52, "Ljung-Box", fitdf=5)



source("backtest.R")
backtest(m2, xts,h=1, orig=length(xts)*0.8)
length(xts)

fit8 <- auto.arima(xts, xreg=vlts, d=0,stepwise=FALSE,approx=FALSE)
fit8

m3=Arima(xts,order=c(4,0,0),xreg=vlts,seasonal=list(order=c(1,0,0),period=52), method="ML")
m3
coeftest(m3)
acf(m3$residuals,na.action = na.pass)
pacf(m3$resid,na.action = na.pass)

#  ljung box test on residuals
Box.test(m3$residuals, 26, "Ljung-Box", fitdf=5) 
Box.test(m3$residuals, 52, "Ljung-Box", fitdf=5)


backtest(m3, xts,h=1, orig=length(xts)*0.8)














