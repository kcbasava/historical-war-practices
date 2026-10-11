setwd("~/trophy-taking")
options(scipen=999) 
library(dplyr)
library(mgcv.helper)
library(mgcv)


#read in dataset----
trophies_new <- read.csv('trophies2026.csv')
trophies_new <- trophies_new |> mutate(Intercept = 1)

# GAMs----
## PC1 linear----
tt.PC1 <- gam(
  trophies ~ 0 +
    Intercept + 
    PC1 +  te(Long, Lat, mean, bs=c('tp', 'cr'), d=c(2,1)),
  method="REML",
  family=binomial(),
  data = trophies_new
)
summary(tt.PC1) 
gam.check(tt.PC1)
confint(tt.PC1, level=0.95)


## PC1 quadratic----
tt.PC1sq <- gam(
  trophies ~ 0 +
    Intercept + 
    PC1 + I(PC1^2) + te(Long, Lat, mean, bs=c('tp', 'cr'), d=c(2,1)),
  data=trophies_new, 
  method="REML",
  family=binomial,
  fit=TRUE
)
summary(tt.PC1sq)
gam.check(tt.PC1sq)
confint(tt.PC1sq, level=0.95)


## expansion ----
tt.exp <- gam(
  trophies ~ 0 +
    Intercept + 
    expansion + te(Long, Lat, mean, bs=c('tp', 'cr'), d=c(2,1)),
  data=trophies_new, 
  method="REML",
  family=binomial
)
summary(tt.exp)
gam.check(tt.exp)
confint(tt.exp, level=0.95)


## expansion + PC1 linear----
tt.PC1exp <- gam(
  trophies ~ 0 +
    Intercept + 
    PC1 + 
    expansion + te(Long, Lat, mean, bs=c('tp', 'cr'), d=c(2,1)),
  data=trophies_new, 
  method="REML",
  family=binomial
)
summary(tt.PC1exp)
gam.check(tt.PC1exp)
confint(tt.PC1exp, level=0.95)


## expansion + PC1 quadratic----
tt.PCsqexp <- gam(
  trophies ~ 0 +
    Intercept + 
    PC1 + I(PC1^2) + 
    expansion + te(Long, Lat, mean, bs=c('tp', 'cr'), d=c(2,1)),
  data=trophies_new, 
  method="REML",
  family=binomial
)
summary(tt.PCsqexp)
gam.check(tt.PCsqexp)
confint(tt.PCsqexp, level=0.95)



