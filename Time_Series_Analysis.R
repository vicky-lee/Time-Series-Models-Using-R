# data Retrieval and libraries
setwd("C:/Users/Vicky Lee/Desktop/CSC 425/Project")

myd=read.table("ex_rate.csv", header = T, sep = ',')
head(myd)

library(rugarch)
library(tseries)
library(fBasics)
library(zoo)
library(forecast)
library(lmtest)

# missing Values and time Plot
sapply(myd, function(x) sum(is.na(x)))
ratets = zoo(myd$rate, as.Date(as.character(myd$date),format=c("%Y-%m-%d")))
rate = na.approx(ratets)
plot(rate,xlab="Time",ylab="Exchange Rate",main="U.S./U.K. Foregin Exchange Rate(# of dollars per 1 british pound)")

# ACF and PACF of daily exchange rate
acf(coredata(rate),main="ACF of Exchange Rate")
pacf(coredata(rate),main="PACF of Exchange Rate")

# log return time series, its plot, and basic statistics
rets = log(rate/lag(rate, -1))
plot(rets,xlab="Time",ylab="Rate",main="log returns of exchage rate")
basicStats(rets)

# log return exchange rates histogram
par(mfcol=c(1,1))
hist(rets, xlab="log returns of exchange rate", prob=TRUE, main="Histogram of log returns")
xfit<-seq(min(rets),max(rets), length=40)
yfit<-dnorm(xfit,mean=mean(rets),sd=sd(rets))
lines(xfit, yfit, col="blue", lwd=2)

# Plots ACF function of log returns
ret = coredata(rets);
acf(ret,main="ACF of log returns of exchange rate")
# plot ACF of squared returns to check for ARCH effect 
acf(ret^2,main="ACF of squared log returns")
# plot ACF of absolute returns to check for ARCH effect 
acf(abs(ret),main="ACF of absolute log returns")

# computes Ljung-Box test on log returns to test independence 
Box.test(coredata(rets),lag=4,type='Ljung')
Box.test(coredata(rets),lag=5,type='Ljung')
Box.test(coredata(rets),lag=7,type='Ljung')
# computes Ljung-Box test on squared log returns to test non-linear independence 
Box.test(coredata(rets^2),lag=3,type='Ljung')
Box.test(coredata(rets^2),lag=5,type='Ljung')
Box.test(coredata(rets^2),lag=7,type='Ljung')
# computes Ljung-Box test on absolute log returns to test non-linear independence 
Box.test(abs(coredata(rets)),lag=3,type='Ljung')
Box.test(abs(coredata(rets)),lag=5,type='Ljung')
Box.test(abs(coredata(rets)),lag=7,type='Ljung')

# fit ARMA(0,0)-GARCH(1,1) with normal distributed errors
garch11.spec=ugarchspec(variance.model=list(garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0)))
# estimate model 
garch11.fit=ugarchfit(spec=garch11.spec, data=rets)
garch11.fit

# fit ARMA(0,0)-GARCH(1,1) model with t-distribution
garch11.t.spec=ugarchspec(variance.model=list(garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0)), distribution.model = "std")
# estimate model 
garch11.t.fit=ugarchfit(spec=garch11.t.spec, data=rets)
garch11.t.fit

# fit ARMA(0,0)-GARCH(1,1) model with skewed t-distribution
garch11.st.spec=ugarchspec(variance.model=list(garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0)), distribution.model = "sstd")
# estimate model 
garch11.st.fit=ugarchfit(spec=garch11.st.spec, data=rets)
garch11.st.fit

# fit ARMA(0,0)-eGARCH(1,1) model with skewed t-distribution
egarch11.st.spec=ugarchspec(variance.model=list(model = "eGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0)), distribution.model = "sstd")
# estimate model 
egarch11.st.fit=ugarchfit(spec=egarch11.st.spec, data=rets)
egarch11.st.fit
plot(egarch11.st.fit)

# compute h-step ahead forecasts for h=1,2,...,10
garch11.st.fcst=ugarchforecast(garch11.st.fit, n.ahead=10)
garch11.st.fcst
plot(garch11.st.fcst)

# compute h-step ahead forecasts for h=1,2,...,10
egarch11.st.fcst=ugarchforecast(egarch11.st.fit, n.ahead=10)
egarch11.st.fcst
plot(egarch11.st.fcst)

# rolling forecasts
egarch11.st.fit=ugarchfit(spec=egarch11.st.spec, data=rets, out.sample=500)
egarch11.st.fcst=ugarchforecast(egarch11.st.fit, n.ahead=12, n.roll=450)
plot(egarch11.st.fcst)

# rolling forecasts
garch11.st.fit=ugarchfit(spec=garch11.st.spec, data=rets, out.sample=500)
garch11.st.fcst=ugarchforecast(garch11.st.fit, n.ahead=12, n.roll=450)
plot(garch11.st.fcst)

# backtesting method to compare EGARCH and GARCH models:
mod_egarch = ugarchroll(egarch11.st.spec, data = rets, n.ahead = 1,
                        n.start = 2000, refit.every = 130, refit.window = "recursive")

mod_garch = ugarchroll(garch11.st.spec, data = rets, n.ahead = 1,
                       n.start = 2000, refit.every = 130, refit.window = "recursive")

# type=VaR shows VaR at 1% level: this is the tail probability.
report(mod_egarch, type="VaR", VaR.alpha = 0.01, conf.level = 0.95)
# type="fpm" shows forecast performance measures 
# (Mean Squared Error (MSE), mean absolute error(MAE) and directional accuracy 
# of the forecasts vs realized returns(DAC)).
report(mod_egarch, type="fpm")
report(mod_garch, type="fpm")

# compute 10-step ahead forecasts 
egarch11.st.fcst=ugarchforecast(egarch11.st.fit, n.ahead=10)
stats=egarch11.st.fcst
plot(egarch11.st.fcst)

# aggregate daily into monthly exchange rate, time plot, and basic statistics
myd$date <- as.Date(myd$date)
myd$Month <- format(myd$date, format="%m")
myd$Year <- format(myd$date,format="%Y")
monthly<-aggregate( rate ~ Month + Year , myd , mean )
monthly$monthly <- as.yearmon(paste(monthly$Year,monthly$Month, sep = "-"), format=c("%Y-%m"))
rate =ts(monthly[,3], start = c(2006), frequency = 12)
ratets = zoo(monthly$rate, as.Date(monthly$monthly,format=c("%m-%Y")))

plot(rate,xlab="Time",ylab="Exchange Rate",main="Monthly U.S./U.K. Foregin Exchange Rate")
basicStats(rate)

# ACF, PACF of monthly exchange rate
acf(coredata(rate),main="ACF of monthly exchange rate")
pacf(coredata(rate),main="PACF of monthly exchange rate")

# ACF, PACF of differenced series 
dx = diff(rate)
head(dx)
acf(as.vector(dx),main="ACF of first difference")
pacf(as.vector(dx),main="PACF of first difference")
auto.arima(dx, ic =c("bic"), trace=TRUE, allowdrift=TRUE)
# Dickey Fuller test for monthly exchange rate
library(fUnitRoots)
# tests for AR model with time trend
adfTest(coredata(rate), lags=1, type=c("ct"))
adfTest(coredata(rate), lags=3, type=c("ct"))
adfTest(coredata(rate), lags=5, type=c("ct"))

# best arima based on bic
auto.arima(rate, ic =c("bic"), trace=TRUE, allowdrift=TRUE)

# fitting ARIMA(1,1,0)(1,0,0)[12] with drift  
m1=Arima(rate,order=c(1,1,0),seasonal=list(order=c(1,0,0),period=12), method="ML",include.drift=T) 
coeftest(m1)

# residual analysis of ARIMA(1,1,0)(1,0,0)[12] with drift
acf(coredata(m1$resid), main="ACF of Residuals")
pacf(coredata(m1$resid), main="PACF of Residuals")

# ljung box test on residuals
Box.test(m1$resid, 6, "Ljung-Box", fitdf=2) 
Box.test(m1$resid, 12, "Ljung-Box", fitdf=2) 

# fitting ARIMA(1,1,0)  
m2=Arima(rate,order=c(1,1,0), method="ML") 
coeftest(m2)

# residual analysis of ARIMA(1,1,0) 
acf(coredata(m2$resid), main="ACF of Residuals")
pacf(coredata(m2$resid), main="PACF of Residuals")

# ljung box test on residuals
Box.test(m2$resid, 6, "Ljung-Box", fitdf=1) 
Box.test(m2$resid, 12, "Ljung-Box", fitdf=1) 

# model validation using backtesting
source("backtest.R")
backtest(m1, rate, h=1, orig=length(rate)*0.8)
backtest(m2, rate, h=1, orig=length(rate)*0.8)

# plot forecasts of ARIMA(1,1,0)(1,0,0)[12] with drift & ARIMA(1,1,0) 
plot(forecast.Arima(m1, h=6), include=120, main="Actual vs Predicted & 6-step forecasts: SARIMA(1,1,0)(1,0,0)[12]", xlab="Date", ylab="Exchange Rate")
lines(fitted(m1), col=2)
legend('topright', legend=c("Actual","Predicted"),lty=1, col=c('black','red' ), cex=.75)
forecast(m1, h=10)
plot(forecast.Arima(m2, h=6), include=120, main="Actual vs Predicted & 6-step forecasts: ARIMA(1,1,0)", xlab="Date", ylab="Exchange Rate")
lines(fitted(m2), col=2)
legend('topright', legend=c("Actual","Predicted"),lty=1, col=c('black','red' ), cex=.75)
forecast(m2, h=10)
