library(car)
library(ggplot2)
library(GGally)
library(lubridate)
library(timeSeries)
library(xts)
source("fun.R")
####################################################################################################################
####################################################################################################################
####################################################################################################################
#########--------load and prepare data--------###########
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
load("PrecDate.RData")
load("Rhine.sts.prec_red.Rdata")
PrecDate.data <- data.frame(PrecDate)

Rhine.sts.prec_red.data <- data.frame(matrix(Rhine.sts.prec_red[,18], ncol = 1))

rain <- cbind(PrecDate.data, Rhine.sts.prec_red.data)
rain.data <- rain$Rain
colnames(rain) <- c("Year", "Rain")

#days 1 - 366
rain$Year <- strptime(rain$Year, "%Y-%m-%d")
rain.ts <- ts(rain$Rain, start = c(1971, 1, 1), frequency = 365.25)

D <- seq(as.Date("1971/01/01"), by="1 day", length.out=13880)
x <- format(D, format="%m-%d")
rain$newcolumn <- x

rain$Year <- with(rain, year(Year))
days.366 <- rep(0, 13880)

#calculate it for 1971 onwards
count <- 1
leap <- 3
index <- 1
counts_year  <- rep(0, 38)

for (i in 1:38){
  #non leap years
  if (leap %% 4 != 0){
    
    for (j in 1:365){
      
      if (j == 60){
        count <- count + 1
        days.366[index] <- count
        count <- count + 1
        index <- index + 1
      }
      else{
        days.366[index] <- count
        index <- index + 1
        count <- count + 1
      }
      
    }
    
  }
  else{
    for (k in 1:366){
      
      days.366[index] <- count
      index <- index + 1
      count <- count + 1
      
      
    }
    
  }
  leap <- leap + 1
  counts_year[i] <- index
  count <- 1
}

####################################################################################################################
####################################################################################################################
####################################################################################################################
#########--------probabilities--------###########
####################################################################################################################
####################################################################################################################
####################################################################################################################

#prepare data for probabilities
source("prob1.R")

#simple visualisation of Markvo chain proportions
par(mfrow = c(2,2))
plot(seq(1, by = 5, length = 73), r.p11, xlab = "Day of the year", ylab = "Probability of Rain", main = expression(r["11"](t)), ylim = c(0,1), cex = rain.p11/15)
plot(seq(1, by = 5, length = 73), r.p10, xlab = "Day of the year", ylab = "Probability of Rain",main = expression(r["10"](t)), ylim = c(0,1), cex = rain.p10/6)
plot(seq(1, by = 5, length = 73), r.p01, xlab = "Day of the year", ylab = "Probability of Rain",main = expression(r["01"](t)), ylim = c(0,1), cex = rain.p01/7)
plot(seq(1, by = 5, length = 73), r.p00, xlab = "Day of the year", ylab = "Probability of Rain",main = expression(r["00"](t)), ylim = c(0,1), cex = rain.p00/15)
# Binary model for rainfall occurrence ======================================

#model 1 different rainfall probabilities each day

model.a <- glm(rain.yes~factor(days.366),family=binomial)
anova(model.a)
summary(model.a)
fit.model.a <- fitted(model.a)

model.b <- glm(cbind(rain.model.yes, rain.model.no) ~ factor(c(1:366)),family=binomial)
anova(model.b)
summary(model.b)
fit.model.b <- fitted(model.b)

layout(matrix(c(2,3,1,1), 2, 2, byrow = TRUE))
plot(seq(as.Date("1971/01/01"), by="5 day", length.out=2776), fit.model.a[seq(1, by = 5, length = 2776)], ylim = c(0,1),xlab= "Year", ylab="Probability of Rain", main = "Model 1", type = "l", lwd=0.5)
plot(c(1:366), fit.model.b, type = "l", ylim = c(0,1),xlab="Day of year",ylab="Probability of rain", main = "Model 2")
plot(c(1:366), rain.model.yes/38, type = "l", ylim = c(0,1),xlab="Day of year",ylab="Probability of rain", main = "Observed proportions of rainy days")

#model 2 logistic regression
#################################################(1,1)
model.simple <- glm(cbind(rain.hist, rain.hist.no )~ 
                      
1 + sin(1 * 2*pi*c(1:366)/366) + cos(1 * 2*pi*c(1:366)/366) 
+ sin(2 * 2*pi*c(1:366)/366) + cos(2 * 2*pi*c(1:366)/366) 
+ sin(3 * 3*pi*c(1:366)/366) + cos(3 * 3*pi*c(1:366)/366) 
+ sin(4 * 3*pi*c(1:366)/366) + cos(4 * 3*pi*c(1:366)/366) 
+ sin(5 * 3*pi*c(1:366)/366) + cos(5 * 3*pi*c(1:366)/366) 
                    
                    , family = binomial(link = "logit"))

summary(model.simple)

#pvalue and AIC
1 - pchisq(sum(residuals(model.simple, type = "pearson")^2), df = 355)
model.simple$aic

anova.simple <- anova(model.simple)

fitted.simple <- fitted(model.simple)


layout(matrix(c(3,1,3,2), 2, 2, byrow = TRUE))
props <- rain.hist/(rain.hist + rain.hist.no)
plot(seq(1, by = 5, length = 73), props[seq(1, by = 5, length = 73)], ylim=c(0,1), cex = rain.hist/30, xlab="Day of year", ylab="Probability of rain", main = "Fitted probabilities against observed proportions")
lines(1:366, fitted.simple[1:366],xlab="Day of year",ylab="Probability of rain", col = "red")
plot(cooks.distance(model.simple), main = "Cook's Distance Model 3", ylab = "Cook's Distance", cex = 0.5)
abline(0.02253521, 0, col = "red")
qqPlot(normalresid(model.simple), cex = 0.5, main = "Q-Q Plot against modified standard residuals", ylab = "Sample Quantiles")

#model 3 Markov chains

#logistic regression
#################################################(1,1)
model.1 <- glm(cbind(rain.11.1,rain.11 - rain.11.1)~ 
                 
                 1 + sin(1 * 2*pi*c(1:366)/366) + cos(1 * 2*pi*c(1:366)/366) 
               + sin(2 * 2*pi*c(1:366)/366) + cos(2 * 2*pi*c(1:366)/366) +
                 sin(3 * 2*pi*c(1:366)/366) + cos(3 * 2*pi*c(1:366)/366) 
               
               , family = binomial(link = "logit"))

fit.11 <- summary(model.1)

rain.fit1 <- fitted(model.1)

anova.1 <- anova(model.1)
devres.1 <- anova.1$`Resid. Dev`[c(1,3,5,7,9,11)]

ANODev(model.1, 8)

model.1$aic


#################################################(0,1)
model.2 <- glm(cbind(rain.01.1,rain.01 - rain.01.1)~ 
                 
1 
               
               , family = binomial(link = "logit"))

model.2$aic

fit.01 <- summary(model.2)

rain.fit2 <- fitted(model.2)
anova.2 <- anova(model.2)
devres.2 <- anova.2$`Resid. Dev`[c(1,3,5,7,9,11)]

ANODev(model.2, 8)

#################################################(1,0)
model.3 <- glm(cbind(rain.10.1,rain.10 - rain.10.1)~ 
                 
1 + sin(1 * 2*pi*c(1:366)/366) + cos(1 * 2*pi*c(1:366)/366) 
+ sin(2 * 2*pi*c(1:366)/366) + cos(2 * 2*pi*c(1:366)/366) 
               
               , family = binomial(link = "logit"))

model.3$aic

fit.10 <- summary(model.3)

rain.fit3 <- fitted(model.3)

anova.3 <- anova(model.3)
devres.3 <- anova.3$`Resid. Dev`[c(1,3,5,7,9,11)]

ANODev(model.3, 8)

#################################################(0,0)
model.4 <- glm(cbind(rain.00.1,rain.00 - rain.00.1)~ 
                 
1 + sin(1 * 2*pi*c(1:366)/366) + cos(1 * 2*pi*c(1:366)/366) 
+ sin(2 * 2*pi*c(1:366)/366) + cos(2 * 2*pi*c(1:366)/366) 
+ sin(3 * 2*pi*c(1:366)/366) + cos(3 * 2*pi*c(1:366)/366) 
               
               , family = binomial(link = "logit"))

model.4$aic

anova.4 <- anova(model.4)
devres.4 <- anova.4$`Resid. Dev`[c(1,3,5,7,9,11)]


fit.00 <- summary(model.4)

rain.fit4 <- fitted(model.4)
ANODev(model.4, 8)


# normal scores plot of residuals.  
qqPlot(normalresid(model.1))
qqPlot(normalresid(model.2))
qqPlot(normalresid(model.3))
qqPlot(normalresid(model.4))


#plots of results
par(mfrow = c(2,2))
plot(rain.fit1, ylim = c(0,1), type = "l",xlab="Day of year",ylab="Probability of rain", main = expression(p["11"](t)), col = "red")
points(seq(1, by = 5, length = 73), r.p11, cex = rain.11/20)

plot(rain.fit2, ylim = c(0,1), type = "l",xlab="Day of year",ylab="Probability of rain", main = expression(p["01"](t)), col = "red")
points(seq(1, by = 5, length = 73), r.p01, cex = rain.01/10)

plot(rain.fit3, ylim = c(0,1), type = "l",xlab="Day of year",ylab="Probability of rain", main = expression(p["10"](t)), col = "red")
points(seq(1, by = 5, length = 73), r.p10, cex = rain.10/10)

plot(rain.fit4, ylim = c(0,1), type = "l",xlab="Day of year",ylab="Probability of rain", main = expression(p["00"](t)), col = "red")
points(seq(1, by = 5, length = 73), r.p00, cex = rain.p00/20)



#diagnostics
#chisq test
1 - pchisq(sum(residuals(model.1, type = "pearson")^2), df = 361)
1 - pchisq(sum(residuals(model.2, type = "pearson")^2), df = 365)
1 - pchisq(sum(residuals(model.3, type = "pearson")^2), df = 361)
1 - pchisq(sum(residuals(model.4, type = "pearson")^2), df = 359)

par(mfrow = c(2,2))
plot(cooks.distance(model.1), main = "Cook's Distance Model 1", ylab = "Cook's Distance", cex = 0.5)
abline(8 / (366 - 5), 0, col = "red")

plot(cooks.distance(model.1), main = "Cook's Distance Model 2", ylab = "Cook's Distance", cex = 0.5)
abline(8 / (366 - 1), 0, col = "red")

plot(cooks.distance(model.1), main = "Cook's Distance Model 3", ylab = "Cook's Distance", cex = 0.5)
abline(8 / (366 - 5), 0, col = "red")

plot(cooks.distance(model.1), main = "Cook's Distance Model 4", ylab = "Cook's Distance", cex = 0.5)
abline(8 / (366 - 7), 0, col = "red")

#with pooling
#logistic regression
#################################################(1,1)
model.1.pool <- glm(cbind(rain.p11.1, rain.p11 - rain.p11.1)~ 
                      
1 
                    
                    , family = binomial(link = "logit"))

fit.1.pool <- summary(model.1.pool)

rain.fit1.pool <- fitted(model.1.pool)

anova.1.pool <- anova(model.1.pool)
devres.1.pool <- anova.1.pool$`Resid. Dev`[c(1,3,5,7,9,11)]

ANODev(model.1.pool, 8)

# normal scores plot of residuals.  
rain.res1.pool <- residuals(model.1.pool)
qqPlot(normalresid(model.1.pool))


#################################################(0,1)
model.2.pool <- glm(cbind(rain.p01.1, rain.p01 - rain.p01.1)~ 
                      
1 
                    
                    , family = binomial(link = "logit"))

fit.2.pool <- summary(model.2.pool)

rain.fit2.pool <- fitted(model.2.pool)

ANODev(model.2.pool, 8)

anova.2.pool <- anova(model.2.pool)
devres.2.pool <- anova.2.pool$`Resid. Dev`[c(1,3,5,7,9,11)]

# normal scores plot of residuals.  
qqPlot(normalresid(model.2.pool))

#################################################(1,0)
model.3.pool <- glm(cbind(rain.p10.1, rain.p10 - rain.p10.1)~ 
                      
1 + sin(1 * 2*pi*c(1:73)/73) + cos(1 * 2*pi*c(1:73)/73) 
                    
                    , family = binomial(link = "logit"))

fit.3.pool <- summary(model.3.pool)

rain.fit3.pool <- fitted(model.3.pool)

ANODev(model.3.pool, 8)


anova.3.pool <- anova(model.3.pool)
devres.3.pool <- anova.3.pool$`Resid. Dev`[c(1,3,5,7,9,11)]

# normal scores plot of residuals.  
qqPlot(normalresid(model.3.pool))

#################################################(0,0)
model.4.pool <- glm(cbind(rain.p00.1, rain.p00 - rain.p00.1)~ 
                      
1 + sin(1 * 2*pi*c(1:73)/73) + cos(1 * 2*pi*c(1:73)/73) 
                    
                    , family = binomial(link = "logit"))

fit.4.pool <- summary(model.4.pool)

rain.fit4.pool <- fitted(model.4.pool)

anova.4.pool <- anova(model.4.pool)
devres.4.pool <- anova.4.pool$`Resid. Dev`[c(1,3,5,7,9,11)]

ANODev(model.4.pool, 8)

# normal scores plot of residuals.  
qqPlot(normalresid(model.3.pool))

#chisquared
#chisq test
1 - pchisq(sum(residuals(model.1.pool, type = "pearson")^2), df = 72)
1 - pchisq(sum(residuals(model.2.pool, type = "pearson")^2), df = 72)
1 - pchisq(sum(residuals(model.3.pool, type = "pearson")^2), df = 70)
1 - pchisq(sum(residuals(model.4.pool, type = "pearson")^2), df = 70)


#plotting
par(mfrow = c(2,2))
plot(seq(1, by = 5, length = 73), rain.fit1.pool[1:73], ylim = c(0,1), type = "l",xlab="Day of year",ylab="Probability of rain", main = expression(p["11"](t)), col = "red")
points(seq(1, by = 5, length = 73), r.p11, cex = rain.p11/20)

plot(seq(1, by = 5, length = 73), rain.fit2.pool[1:73], ylim = c(0,1), type = "l",xlab="Day of year",ylab="Probability of rain", main = expression(p["01"](t)), col = "red")
points(seq(1, by = 5, length = 73), r.p01, cex = rain.p01/7)

plot(seq(1, by = 5, length = 73), rain.fit3.pool[1:73], ylim = c(0,1), type = "l",xlab="Day of year",ylab="Probability of rain", main = expression(p["10"](t)), col = "red")
points(seq(1, by = 5, length = 73), r.p10, cex = rain.p10/7)


plot(seq(1, by = 5, length = 73), rain.fit4.pool[1:73], ylim = c(0,1), type = "l",xlab="Day of year",ylab="Probability of rain", main = expression(p["00"](t)), col = "red")
points(seq(1, by = 5, length = 73), r.p00, cex = rain.p00/20)


####################################################################################################################
####################################################################################################################
#########--------amount--------###########
####################################################################################################################
####################################################################################################################
####################################################################################################################
#data preparation
source("amount1.R")

#Harmonic series with Gamma Regression without markov chains
df2 <- data.frame(xbar,t.prime)

model.2 <- glm(xbar ~ 
                 
1 + sin(1 * 2*pi*c(1:366)/366) + cos(1 * 2*pi*c(1:366)/366) +  
sin(2 * 2*pi*c(1:366)/366) + cos(2 * 2*pi*c(1:366)/366) + 
sin(3 * 2*pi*c(1:366)/366) + cos(3 * 2*pi*c(1:366)/366) +
sin(4 * 2*pi*c(1:366)/366) + cos(4 * 2*pi*c(1:366)/366)
               
               ,family=Gamma(link="log"), data = df2)

model.2$aic

anova.2 <- anova(model.2,test="F")
summary.2 <- summary(model.2)
fit.2 <- fitted(model.2)
FDev(model.2, 5)

Fstats.2 <- FDev(model.2, 5)[,2]
dev.2 <- anova.2$`Resid. Dev`[c(1,3,5,7,9)]

#pooling 
df3 <- data.frame(xpool,t.pool)

model.3 <- glm(xpool ~ 
                 
1 + sin(1 * 2*pi*c(1:73)/73) + cos(1 * 2*pi*c(1:73)/73) 
+ sin(2 * 2*pi*c(1:73)/73) + cos(2 * 2*pi*c(1:73)/73) 
+ sin(3 * 2*pi*c(1:73)/73) + cos(3 * 2*pi*c(1:73)/73) 
+ sin(4 * 2*pi*c(1:73)/73) + cos(4 * 2*pi*c(1:73)/73) 
+ sin(5 * 2*pi*c(1:73)/73) + cos(5 * 2*pi*c(1:73)/73) 
               
               ,family=Gamma(link="log"), data = df3)

model.3$aic

anova.3 <- anova(model.3,test="F")
summary(model.3)
Fstats.3 <- FDev(model.3, 5)[,2]
dev.3 <- anova.3$`Resid. Dev`[c(1,3,5,7,9)]
fit.3 <- fitted(model.3)

#plot them together
par(mfrow = c(2,1))
plot(fit.2,xlab="Day of year",ylab="Mean rainfall (mm)", type = "l", col = "red", ylim = c(0 ,20), main = "Daily")
points(seq(1, by = 5, length = 73), xpool[1:73], cex = npool/125)
plot(seq(1, by = 5, length = 73), fit.3[1:73],xlab="Day of year",ylab="Mean rainfall (mm)", type = "l", col = "red", ylim = c(0 ,15), main = "5-day pooling")
points(seq(1, by = 5, length = 73), xpool[1:73], cex = npool/125)

#qqplot
par(mfrow = c(1,2))
qqPlot(normalresid(model.2), cex =0.5, main = "Daily", ylab = "Sample Quantiles")
qqPlot(normalresid(model.3), cex = 0.5, main  = "5-day Pooling", ylab = "Sample Quantiles")

#cooks

par(mfrow = c(1,2))
plot(cooks.distance(model.2), main = "Cook's Distance Model 2", ylab = "Cook's Distance", cex = 0.5)
abline(0.02272727, 0, col = "red")

plot(cooks.distance(model.3), main = "Cook's Distance Model 3", ylab = "Cook's Distance", cex = 0.5)
abline(0.02272727, 0, col = "red")

#Markov chains
############################################################
########Without pooling########################################
############################################################

df11 <- data.frame(xbar11,t.prime)

model.11 <- glm(xbar11 ~ 
                  
1 + sin(1 * 2*pi*c(1:366)/366) + cos(1 * 2*pi*c(1:366)/366) +  
sin(2 * 2*pi*c(1:366)/366) + cos(2 * 2*pi*c(1:366)/366) + 
sin(3 * 2*pi*c(1:366)/366) + cos(3 * 2*pi*c(1:366)/366) + 
sin(4 * 2*pi*c(1:366)/366) + cos(4 * 2*pi*c(1:366)/366) 
                
                ,family=Gamma(link="log"), data = df11)

anova.11 <- anova(model.11,test="F")
summary(model.11)
Fstats.11 <- FDev(model.11, 5)[,2]
dev.11 <- anova.11$`Resid. Dev`[c(1,3,5,7,9)]

fit.11 <- fitted(model.11)

############################################################(0,1)
############################################################
df01 <- data.frame(xbar01,t.prime)

model.01 <- glm(xbar01 ~ 
                  
1 + sin(1 * 2*pi*c(1:366)/366) + cos(1 * 2*pi*c(1:366)/366)
                
                ,family=Gamma(link="log"), data = df01)

anova.01 <- anova(model.01,test="F")
summary(model.01)
Fstats.01 <- FDev(model.01, 6)[,2]
dev.01 <- anova.01$`Resid. Dev`[c(1,3,5,7,9)]


fit.01 <- fitted(model.01)

############################################################(1,0)
############################################################
df10 <- data.frame(xbar10,t.prime)

model.10 <- glm(xbar10 ~ 
                  
1 + sin(1 * 2*pi*c(1:366)/366) + cos(1 * 2*pi*c(1:366)/366) 
                
                ,family=Gamma(link="log"), data = df10)

anova.10 <- anova(model.10,test="F")
summary(model.10)
Fstats.10 <- FDev(model.10, 5)[,2]
dev.10 <- anova.10$`Resid. Dev`[c(1,3,5,7,9)]

fit.10 <- fitted(model.10)

############################################################(0,0)
############################################################
df00 <- data.frame(xbar00,t.prime)

model.00 <- glm(xbar00 ~ 
                  
1 + sin(1 * 2*pi*c(1:366)/366) + cos(1 * 2*pi*c(1:366)/366) +  
sin(2 * 2*pi*c(1:366)/366) + cos(2 * 2*pi*c(1:366)/366) + 
sin(3 * 2*pi*c(1:366)/366) + cos(3 * 2*pi*c(1:366)/366)
                
                ,family=Gamma(link="log"), data = df00)

anova.00 <- anova(model.00,test="F")
summary(model.00)
Fstats.00 <- FDev(model.00, 5)[,2]
dev.00 <- anova.00$`Resid. Dev`[c(1,3,5,7,9)]

fit.00 <- fitted(model.00)

############################################################
############################################################
#qqplots

par(mfrow = c(2,2))
qqPlot(normalresid(model.11), ylab = "Sample Quantiles", main = expression(p["11"](t)), cex = 0.5)
qqPlot(normalresid(model.10), ylab = "Sample Quantiles", main = expression(p["10"](t)), cex = 0.5)
qqPlot(normalresid(model.01), ylab = "Sample Quantiles", main = expression(p["01"](t)), cex = 0.5)
qqPlot(normalresid(model.00), ylab = "Sample Quantiles", main = expression(p["00"](t)), cex = 0.5)


#plotting 
par(mfrow = c(2,2))
plot(fit.11,xlab="Day of year",ylab="Mean rainfall (mm)", type = "l", col = "red", ylim = c(0 ,17), main = expression(p["00"](t)))
points(seq(1, by = 5, length = 73), xpool11[1:73], cex = npool11/60)

plot(fit.10,xlab="Day of year",ylab="Mean rainfall (mm)", type = "l", col = "red", ylim = c(0 ,17), main = expression(p["00"](t)))
points(seq(1, by = 5, length = 73), xpool10[1:73], cex = npool10/15)

plot(fit.01,xlab="Day of year",ylab="Mean rainfall (mm)", type = "l", col = "red", ylim = c(0 ,17), main = expression(p["00"](t)))
points(seq(1, by = 5, length = 73), xpool01[1:73], cex = npool01/20)

plot(fit.00,xlab="Day of year",ylab="Mean rainfall (mm)", type = "l", col = "red", ylim = c(0 ,17), main = expression(p["00"](t)))
points(seq(1, by = 5, length = 73), xpool00[1:73], cex = npool00/20)

#cook's
par(mfrow = c(2,2))
plot(cooks.distance(model.11), main = "Cook's Distance Model 2", ylab = "Cook's Distance", cex = 0.5)
abline(0.02247191, 0, col = "red")

plot(cooks.distance(model.01), main = "Cook's Distance Model 3", ylab = "Cook's Distance", cex = 0.5)
abline(0.02247191, 0, col = "red")

plot(cooks.distance(model.10), main = "Cook's Distance Model 3", ylab = "Cook's Distance", cex = 0.5)
abline(0.02247191, 0, col = "red")

plot(cooks.distance(model.00), main = "Cook's Distance Model 3", ylab = "Cook's Distance", cex = 0.5)
abline(0.02247191, 0, col = "red")

############################################################
########With pooling########################################
############################################################
############################################################(1,1)
############################################################
df11p <- data.frame(xpool11,t.pool)

model.11p <- glm(xpool11 ~ 
                   
1 + sin(1 * 2*pi*c(1:73)/73) + cos(1 * 2*pi*c(1:73)/73) +  
sin(2 * 2*pi*c(1:73)/73) + cos(2 * 2*pi*c(1:73)/73) + 
sin(3 * 2*pi*c(1:73)/73) + cos(3 * 2*pi*c(1:73)/73)+
sin(4 * 2*pi*c(1:73)/73) + cos(4 * 2*pi*c(1:73)/73)
                 
                 ,family=Gamma(link="log"), data = df11p)

anova.11p <- anova(model.11p,test="F")
summary(model.11p)
Fstats.11p <- FDev(model.11p, 5)[,2]
dev.11p <- anova.11p$`Resid. Dev`[c(1,3,5,7,9)]

fit.11p <- fitted(model.11p)

############################################################(0,1)
############################################################
df01p <- data.frame(xpool01,t.pool)

model.01p <- glm(xpool01 ~ 
                   
1 + sin(1 * 2*pi*c(1:73)/73) + cos(1 * 2*pi*c(1:73)/73) 
                 
                 ,family=Gamma(link="log"), data = df01p)

anova.01p <- anova(model.01p,test="F")
summary(model.01p)
Fstats.01p <- FDev(model.01p, 5)[,2]
dev.01p <- anova.01p$`Resid. Dev`[c(1,3,5,7,9)]

fit.01p <- fitted(model.01p)

############################################################(1,0)
############################################################
df10p <- data.frame(xpool10,t.pool)

model.10p <- glm(xpool10 ~ 
                   
1 + sin(1 * 2*pi*c(1:73)/73) + cos(1 * 2*pi*c(1:73)/73)
                 
                 ,family=Gamma(link="log"), data = df10p)

anova.10p <- anova(model.10p,test="F")
summary(model.10p)
Fstats.10p <- FDev(model.10p, 5)[,2]
dev.10p <- anova.10p$`Resid. Dev`[c(1,3,5,7,9)]

fit.10p <- fitted(model.10p)

############################################################(0,0)
############################################################
df00p <- data.frame(xpool00,t.pool)

model.00p <- glm(xpool00 ~ 
                   
1 + sin(1 * 2*pi*c(1:73)/73) + cos(1 * 2*pi*c(1:73)/73) +  
sin(2 * 2*pi*c(1:73)/73) + cos(2 * 2*pi*c(1:73)/73) 
                 
                 ,family=Gamma(link="log"), data = df00p)

anova.00p <- anova(model.00p,test="F")
summary(model.00p)
Fstats.00p <- FDev(model.00p, 5)[,2]
dev.00p <- anova.00p$`Resid. Dev`[c(1,3,5,7,9)]

fit.00p <- fitted(model.00p)

############################################################
############################################################
#qqplot
par(mfrow = c(2,2))
qqPlot(normalresid(model.11p), ylab = "Sample Quantiles", main = expression(p["11"](t)), cex = 0.5)
qqPlot(normalresid(model.10p), ylab = "Sample Quantiles", main = expression(p["10"](t)), cex = 0.5)
qqPlot(normalresid(model.01p), ylab = "Sample Quantiles", main = expression(p["01"](t)), cex = 0.5)
qqPlot(normalresid(model.00p), ylab = "Sample Quantiles", main = expression(p["00"](t)), cex = 0.5)

#plotting 
par(mfrow = c(2,2))
plot(seq(1, by = 5, length = 73), fit.11p[1:73],xlab="Day of year",ylab="Mean rainfall (mm)", type = "l", col = "red", ylim = c(0 ,17),  main = expression(p["11"](t)))
points(seq(1, by = 5, length = 73), xpool11[1:73], cex = npool11/70)

plot(seq(1, by = 5, length = 73), fit.10p[1:73],xlab="Day of year",ylab="Mean rainfall (mm)", type = "l", col = "red", ylim = c(0 ,17),  main = expression(p["10"](t)))
points(seq(1, by = 5, length = 73), xpool10[1:73], cex = npool10/15)

plot(seq(1, by = 5, length = 73), fit.01p[1:73],xlab="Day of year",ylab="Mean rainfall (mm)", type = "l", col = "red", ylim = c(0 ,17),  main = expression(p["01"](t)))
points(seq(1, by = 5, length = 73), xpool01[1:73], cex = npool01/25)

plot(seq(1, by = 5, length = 73), fit.00p[1:73],xlab="Day of year",ylab="Mean rainfall (mm)", type = "l", col = "red", ylim = c(0 ,17), main = expression(p["00"](t)))
points(seq(1, by = 5, length = 73), xpool00[1:73], cex = npool00/20)

#cook's
par(mfrow = c(2,2))
plot(cooks.distance(model.11p), main = "Cook's Distance Model 2", ylab = "Cook's Distance", cex = 0.5)
abline(0.02247191, 0, col = "red")

plot(cooks.distance(model.01p), main = "Cook's Distance Model 3", ylab = "Cook's Distance", cex = 0.5)
abline(0.02247191, 0, col = "red")

plot(cooks.distance(model.10p), main = "Cook's Distance Model 3", ylab = "Cook's Distance", cex = 0.5)
abline(0.02247191, 0, col = "red")

plot(cooks.distance(model.00p), main = "Cook's Distance Model 3", ylab = "Cook's Distance", cex = 0.5)
abline(0.02247191, 0, col = "red")

####################################################################################################################
####################################################################################################################
#########--------basic plot of frequency and amount--------###########
####################################################################################################################
####################################################################################################################
####################################################################################################################
par(mfrow = c(1,2))
plot(rain.hist, cex = 0.5, main = "Frequency of rainy days", xlab = "Day of the year", ylab = "Frequency")
plot(xbar, cex = 0.5, main = "Amount of rainfall (mm)", xlab = "Day of the year", ylab = "Amount of rainfall")