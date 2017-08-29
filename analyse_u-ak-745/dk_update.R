# Experiments with building a "stable" Dicekrigiging emulator, with 
# advice from David Sexton

setwd("/Users/dougmcneall/Documents/work/R/brazil_cssp/analyse_u-ak-745")


library(devtools)
install_github(repo = "dougmcneall/hde")

library(hde)
library(RColorBrewer)
library(fields)
library(MASS)
library(DiceKriging)
library(stepwise)

source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/imptools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/vistools.R")

yg = brewer.pal(9, "YlGn")
ryb = brewer.pal(11, "RdYlBu")
byr = rev(ryb)
rb = brewer.pal(11, "RdBu")
br = rev(rb)

#read data
frac = read.table('forest_frac.txt', header = FALSE)
lhs = read.table('lhs_u-ak745.txt', header = TRUE)

y = frac[, 2]

X = normalize(lhs)

# It turns out that a true "leave one out" function
# gives very similar answers to "trend.reestim = TRUE"
# (wth the latter being very much quicker). 
# trend.reestim = FALSE gives a much more forgiving 
# view of the fit than is justified.

true.loo = function(X,y){
  out.mean = rep(NA, length(y))
  out.sd = rep(NA, length(y))
  
  for(i in 1:nrow(X)){
    X.trunc = X[-i, ]
    y.trunc = y[-i]
    
    X.target = matrix(X[i, ], nrow = 1)
    
    fit = km(~., design = X.trunc, response = y.trunc)
    pred = predict(fit,newdata = X.target, type = 'UK')
    out.mean[i] = pred$mean
    out.sd[i = pred$sd]
  }
  return(list(mean = out.mean, sd = out.sd))
}

loo.rmse = function(loo,y){
  out = sqrt(mean((loo$mean - y)^2))
  out
}

fit0 = km(~., design = X, response = y)
loo0 = leaveOneOut.km(fit0, type = 'UK', trend.reestim = TRUE)
true.loo0 = true.loo(X = X, y = y)
loo.rmse(loo0, y)
loo.rmse(true.loo0, y)

fit1 = km(~1, design = X, response = y)
loo1 = leaveOneOut.km(fit1, type = 'UK', trend.reestim = TRUE)

plot(y, loo0$mean)
points(y, true.loo0$mean, col = 'red')
abline(0,1)


loo.rmse(loo1, y)


plot(y, loo0$mean)
abline(0,1)

fit2 = km(~., design = X, response = y, optim.method = 'gen')
loo2 = leaveOneOut.km(fit2, type = 'UK', trend.reestim = TRUE)
loo.rmse(loo2, y)
plot(fit)



fit = lm(y~X)

test = step(fit, direction = 'both',trace = 1, steps = 1000)

test = regsubsets(X,y)




