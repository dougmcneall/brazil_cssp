# test_ds_emulator
# Testing David Sexton's emulator approach.
# Conclusion: it seems to make a more robust emulator.
# An upshot of this is: I've discovered that the GP emulator
# is just returning a linear model

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

true.loo = function(X,y, lm = FALSE){
  out.mean = rep(NA, length(y))
  out.sd = rep(NA, length(y))
  
  for(i in 1:nrow(X)){
    X.trunc = X[-i, ]
    y.trunc = y[-i]
    
    X.target = matrix(X[i, ], nrow = 1)
    colnames(X.target) <- colnames(X)
    X.target.df = data.frame(X.target)
    
    
    if(lm){
      data.trunc = data.frame(y.trunc, X.trunc)
      colnames(data.trunc) = c('y', colnames(X.trunc))
      fit = lm(y~., data = data.trunc)
      pred = predict(fit, newdata = X.target.df)
      out.mean[i] = pred
    }
    else{
    
    fit = km(~., design = X.trunc, response = y.trunc)
    pred = predict(fit,newdata = X.target, type = 'UK')
    out.mean[i] = pred$mean
    out.sd[i = pred$sd]
    }
    

  }
  return(list(mean = out.mean, sd = out.sd))
}



loo.rmse = function(loo,y){
  out = sqrt(mean((loo$mean - y)^2))
  out
}


#read data
# Tropical Broadleaf Evergreen data
frac = read.table('forest_frac.txt', header = FALSE)
lhs = read.table('lhs_u-ak745.txt', header = TRUE)

y = frac[, 2]
X = normalize(lhs)

# first, build a "straight out of the box" km model, and have a look at the error statistics
m.box = km(~., design=X, response=y)

loo.box = leaveOneOut.km(m.box, type = 'UK', trend.reestim = TRUE)
loo.box.noreestim = leaveOneOut.km(m.box, type = 'UK', trend.reestim = FALSE)
plot(y, loo.box$mean)
loo.rmse(loo.box, y)
loo.rmse(loo.box.noreestim, y)

m.box.nugget = km(~., design=X, response=y, nugget = 0.01)

# source('makeccEm.R')

#xNames = as.formula(paste("~ ", paste(colnames(lhs), collapse= "+")))
xvars = colnames(lhs)

data = data.frame(response=y, x=X)
colnames(data) <- c("response", xvars)
nval = length(y)

#test = makeGPEmStep(x = X,xNames = xNames , y = y, e = 0.1)
nugget=NULL
nuggetEstim=FALSE
noiseVar=NULL
seed=NULL
trace=FALSE
maxit=500
REPORT=10
factr=1e7
pgtol=0.0
parinit=NULL
popsize=100


control_list = list(trace=trace, maxit=maxit, REPORT=REPORT, factr=factr, pgtol=pgtol, pop.size=popsize)

# Build the emulator. 
#data = data.frame(response=response, x=subinputs)
#colnames(data) <- c("response", xvars)

# Build the first emulator with a flat prior
m0 = km(y ~ 1, design=X, response=y, nugget=nugget,
        nugget.estim=nuggetEstim, noise.var=noiseVar, control=control_list)

coefs0 = m0@covariance@range.val 
print('coefs0')
print(coefs0)

start.form = as.formula(paste("y ~ ", paste(xvars, collapse= "+")))

# use BIC so that model is parsimonious and allow GP to pick up any other behaviour not
# explained by the key linear terms      
startlm = lm(start.form, data=data)
#      print('before step')
#      print(startlm)
steplm = step(startlm, direction="both", k=log(nval), trace=TRUE)
print('after step')
print(steplm)
form = as.formula(steplm)
print('Formula')
print(form)
data$response = NULL
labels = labels(terms(steplm))
labels = labels[!(labels %in% c('response'))]
if (length(labels) > 0) {
  start.form = as.formula(paste("~ ", paste(labels, collapse= "+")))
} else {
  start.form = as.formula("~ 1")
}    
print("Step has found formula:")
print(start.form)
if (!is.null(seed)) {set.seed(seed)}
nugget = 0.01
m.step = km(start.form, design=X, response=y, nugget=nugget, parinit=coefs0,
       nugget.estim=nuggetEstim, noise.var=noiseVar, control=control_list)


loo.step = leaveOneOut.km(m.step, type = 'UK', trend.reestim = TRUE)
loo.step.noreestim = leaveOneOut.km(m.step, type = 'UK', trend.reestim = FALSE)

plot(y, loo.box$mean)
points(y, loo.step$mean, col = 'red')
abline(0,1)

loo.rmse(loo.box, y)
loo.rmse(loo.box.noreestim, y)

loo.rmse(loo.step, y)
loo.rmse(loo.step.noreestim, y)


# Make David's emulator code into a (simplified) function.

twoStep = function(X, y, nugget=NULL, nuggetEstim=FALSE, noiseVar=NULL, seed=NULL, trace=FALSE, maxit=100,
                   REPORT=10, factr=1e7, pgtol=0.0, parinit=NULL, popsize=100){
  
  control_list = list(trace=trace, maxit=maxit, REPORT=REPORT, factr=factr, pgtol=pgtol, pop.size=popsize)
  
  xvars = colnames(X)
  data = data.frame(response=y, x=X)
  colnames(data) <- c("response", xvars)
  nval = length(y)
  
  # Build the first emulator with a flat prior
  m0 = km(y ~ 1, design=X, response=y, nugget=nugget,
          nugget.estim=nuggetEstim, noise.var=noiseVar, control=control_list)
  
  coefs0 = m0@covariance@range.val 
  print('coefs0')
  print(coefs0)
  
  start.form = as.formula(paste("y ~ ", paste(xvars, collapse= "+")))
  
  # use BIC so that model is parsimonious and allow GP to pick up any other behaviour not
  # explained by the key linear terms      
  startlm = lm(start.form, data=data)
  #      print('before step')
  #      print(startlm)
  steplm = step(startlm, direction="both", k=log(nval), trace=TRUE)
  print('after step')
  print(steplm)
  form = as.formula(steplm)
  print('Formula')
  print(form)
  data$response = NULL
  labels = labels(terms(steplm))
  labels = labels[!(labels %in% c('response'))]
  if (length(labels) > 0) {
    start.form = as.formula(paste("~ ", paste(labels, collapse= "+")))
  } else {
    start.form = as.formula("~ 1")
  }    
  print("Step has found formula:")
  print(start.form)
  if (!is.null(seed)) {set.seed(seed)}
  m = km(start.form, design=X, response=y, nugget=nugget, parinit=coefs0,
         nugget.estim=nuggetEstim, noise.var=noiseVar, control=control_list)
  
  return(list(x=X, y=y, nugget=nugget, nugget.estim=nuggetEstim,
              noise.var=noiseVar, emulator=m, seed=seed, coefs=m@covariance@range.val,
              trends=m@trend.coef, meanTerms=all.vars(start.form), steplm = steplm))
  
}

twostep = twoStep(X, y, nugget = 0.001, nuggetEstim = FALSE, maxit = 1000)

twostep.fit = twostep$emulator
twostep.lm = twostep$steplm

# Now do the sensitivity analysis etc. has anything changed?

n = 21
X.oaat= oaat.design(X, n, med = TRUE)
colnames(X.oaat) <- colnames(X)
X.oaat.df = data.frame(X.oaat)
y.oaat = predict(twostep.fit, newdata = X.oaat, type = 'UK')
y.oaat.ts.lm = predict(twostep.lm, newdata = X.oaat.df)
d = ncol(X) 

pdf(file = 'oaat_twostep.pdf', width = 8, height = 10)
par(mfrow = c(7,11), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = c(-0.12, 0.4)

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], y.oaat$mean[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE, lwd = 2)
  
  points(X.oaat[ix,i], y.oaat$mean[ix] + y.oaat$sd[ix] ,
       type = 'l', lwd = 2, col = 'grey')
  
  points(X.oaat[ix,i], y.oaat$mean[ix] - y.oaat$sd[ix] ,
         type = 'l', lwd = 2, col = 'grey')
  
  
  points(X.oaat[ix,i], y.oaat.ts.lm[ix], type = 'l', col = 'red')
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
dev.off()



pdf(file = 'oaat_twostep_lm.pdf', width = 8, height = 10)
par(mfrow = c(7,11), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = c(-0.12, 0.4)

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], y.oaat.ts.lm[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE)
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
dev.off()


# Construct a straight-up lm

data.all = data.frame(y=y, X=X)
colnames(data.all) <- c('y', colnames(X))
fit.lm = lm(y~., data = data.all)
y.oaat.lm = predict(fit.lm, newdata = X.oaat.df)
y.oaat.box = predict(m.box, newdata = X.oaat, type = 'UK')

# Construct a stepwise lm

steplm.fit = step(fit.lm, direction="both", k=log(length(y)), trace=TRUE)

pdf(file = 'oaat_lm.pdf', width = 10, height = 6)
par(mfrow = c(6,13), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = c(-0.12, 0.4)

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], y.oaat.lm[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE, lwd = 2)
  points(X.oaat[ix,i], y.oaat.box$mean[ix], col = 'blue', type = 'l', lty = 'dashed',
         lwd = 2)
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
plot(1:10, type = 'n', axes = FALSE)
legend('center', legend = c('km', 'lm'), col = c('black', 'blue'),
       lwd = 2, lty = c('solid', 'dashed'), bty = 'n')
                                                            
dev.off()

# -------------------------------------------------------------
# Illustrate the problem
# -------------------------------------------------------------

# 1. Construct a stepwise linear model
data = data.frame(response=y, x=X)
colnames(data) <- c("response", xvars)
nval = length(y)
d = ncol(X)


start.form = as.formula(paste("y ~ ", paste(xvars, collapse= "+")))

# use BIC so that model is parsimonious and allow GP to pick up any other behaviour not
# explained by the key linear terms      
startlm = lm(start.form, data=data)
#      print('before step')
#      print(startlm)
steplm = step(startlm, direction="both", k=log(nval), trace=TRUE)

n = 21
X.oaat= oaat.design(X, n, med = TRUE)
colnames(X.oaat) <- colnames(X)
X.oaat.df = data.frame(X.oaat)

y.startlm.oaat = predict(startlm, newdata = X.oaat.df)
y.steplm.oaat = predict(steplm, newdata = X.oaat.df)

# The stepwise linear model alters some of the sensitivities
# and sets others to zero.
pdf(file = 'oaat_lm_step.pdf', width = 10, height = 6)
par(mfrow = c(6, 13), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = c(-0.12, 0.4)

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], y.startlm.oaat[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE, lwd = 2)
  
  points(X.oaat[ix,i], y.steplm.oaat[ix],
       type = 'l', col='red', lwd = 2)
  
  mtext(1, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
plot(1:10, type = 'n', axes = FALSE)
legend('center',legend = c('lm', 'step lm'),
       lty = c('solid', 'solid'),
       col = c( 'black', 'red'),
       lwd = c(2,2), bty = 'n')
dev.off()


# We can find the steplm coefficients that it keeps, and extract the relevant
# columns from the data matrix
subterms = attr(terms(steplm),"term.labels")
subs.ix = which(colnames(X) %in% subterms)

X.subs = X[, subs.ix]
data.subs = data.frame(y, X.subs)
colnames(data.subs) = c('y', colnames(X.subs))

# Linear model with the chosen subset of inputs
sublm = lm(y~., data = data.subs)

# If we run step a second time, nothing should happen, right?
sublm.step = step(sublm, direction="both", k=log(nval), trace=TRUE)
# I think ordering changes

# Create a one-at-a-time subset input design
X.subs.oaat= oaat.design(X.subs, n, med = TRUE)
colnames(X.subs.oaat) <- colnames(X.subs)
X.subs.oaat.df = data.frame(X.subs.oaat)

y.subs.oaat = predict(sublm, newdata = X.subs.oaat.df)

subskm0 = km(y~1, design = X.subs, response = y)
subskm = km(y~., design = X.subs, response = y, parinit=subskm0@covariance@range.val)

y.subskm0.oaat = predict(subskm0, newdata = X.subs.oaat, type = 'UK')
y.subskm.oaat = predict(subskm, newdata = X.subs.oaat, type = 'UK')


d.sub = ncol(X.subs)

# This confirms that the step lm and the lm given the subset give the same
# results, and that they are both the same as the km given the subsets (it just returns a linear model)
pdf(file = 'oaat_subs.pdf', width = 10, height = 6)
par(mfrow = c(3,8), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = c(-0.12, 0.4)

for(i in 1:d.sub){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  
  # linear model built on subselected columns
  ix.subs = seq(from = ((subs.ix[i]*n) - (n-1)), to =  (subs.ix[i]*n), by = 1)
  print(ix.subs)
  
  plot(X.subs.oaat[ix,i], y.subs.oaat[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE, lwd = 3, col = 'grey')
  
  # Flat prior Gaussian process
  points(X.subs.oaat[ix,i], y.subskm0.oaat$mean[ix],
         type = 'l', col='darkred', lwd = 2)
  
  points(X.subs.oaat[ix,i], y.subskm0.oaat$mean[ix] + y.subskm0.oaat$sd[ix],
         type = 'l', col='tomato')
  
  points(X.subs.oaat[ix,i], y.subskm0.oaat$mean[ix] - y.subskm0.oaat$sd[ix],
         type = 'l', col='tomato')
  
  # Linear prior Gaussian process
  points(X.subs.oaat[ix,i], y.subskm.oaat$mean[ix],
         type = 'l', col='blue', lwd = 2, lty = 'dashed')
  
  # Non zero components of the step linear model
  points(X.oaat[ix.subs, subs.ix[i]], y.steplm.oaat[ix.subs],
         type = 'l', col='black', lwd =2, lty = 'dotted')
  
  mtext(1, text = colnames(X.subs)[i], line = 0.2, cex = 0.7)
}
plot(1:10, type = 'n', axes = FALSE)
legend('center',legend = c('subselected lm', 'step lm', 'km linear prior','km flat prior'),
       lty = c('solid', 'dotted', 'dashed', 'solid'),
       col = c('grey', 'black', 'blue', 'darkred'),
       lwd = c(3,2,2,2), bty = 'n')
dev.off()


# Leave-one-out performance of various models

lm.loo = true.loo(X = X.subs, y = y, lm = TRUE)
y.subskm0.loo = leaveOneOut.km(subskm0,trend.reestim = TRUE, type = 'UK') 
y.subskm.loo = leaveOneOut.km(subskm,trend.reestim = TRUE, type = 'UK') 

plot(y, lm.loo$mean)

plot(y, y.subskm0.loo$mean, col = 'red', xlim = c(-0.1, 0.5), ylim = c(-0.1, 0.5),
     pch = 19)
points(y, y.subskm.loo$mean, col = 'blue', pch = 19)
abline(0,1)

# It LOOKS like the flat prior is better
loo.rmse(lm.loo, y)
loo.rmse(y.subskm0.loo, y)


# Do we do better if we model log (y)??
# Or sqrt(y)?

# Not using the full model.

sqrty = sqrt(y)

# looks like this craps out
subskm0sqrty = km(sqrty~1, design = X, response = sqrty)

subskmsqrty = km(sqrty~., design = X, response = sqrty)
test = leaveOneOut.km(subskmsqrty, type = 'UK', trend.reestim = TRUE)

plot(y, (test$mean)^2)
loo.rmse(test, y)

data = data.frame(sqrt(y), X)

colnames(data) <- c('sqrty', colnames(X))
sqrty.lm = lm(sqrty~., data = data)

# If we run step a second time, nothing should happen, right?
sqrty.step = step(sqrty.lm, direction="both", k=log(nval), trace=TRUE)

test = twoStep(X, sqrty)
twostep.fit = test$emulator
twostep.lm = test$steplm

loo.ts = leaveOneOut.km(twostep.fit, type = 'UK', trend.reestim = FALSE)
sqrt(mean(loo.ts$mean - y)^2)



