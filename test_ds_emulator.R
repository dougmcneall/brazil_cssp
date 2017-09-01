# test_ds_emulator
# Testing David Sexton's emulator approach.
# Conclusion: it seems to make a more robust emulator.

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

m.box.nugget = km(~., design=X, response=y, nugget.estim = TRUE)

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
m.step = km(start.form, design=X, response=y, nugget=nugget, parinit=coefs0,
       nugget.estim=nuggetEstim, noise.var=noiseVar, control=control_list)


loo.step = leaveOneOut.km(m, type = 'UK', trend.reestim = TRUE)
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

twostep = twoStep(X, y, nuggetEstim = TRUE, maxit = 1000)

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


pdf(file = 'oaat_lm.pdf', width = 8, height = 10)
par(mfrow = c(7,11), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = c(-0.12, 0.4)

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], y.oaat.lm[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE, lwd = 2)
  points(X.oaat[ix,i], y.oaat.box$mean[ix], col = 'red', type = 'l')
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
dev.off()

# Can't get this to work
#resid = steplm.fit$residuals
#fit.steplm.resid = km(~1, design = X, response = resid )


# This shows promise - make it two-step?
# Which inputs do we select?
subs = labels(terms(steplm.fit))
subs.ix = which(colnames(X) %in% subs)

X.subs = X[, subs.ix]

test = km(~1, design = X.subs, response = y)

X.subs.oaat = oaat.design(X.subs, n, med = TRUE)
colnames(X.subs.oaat) = colnames(X.subs)
y.subs.oaat = predict(test, newdata = X.subs.oaat, type = 'UK')

d.subs = ncol(X.subs.oaat)

pdf(file = 'oaat_subs.pdf', width = 8, height = 10)
par(mfrow = c(7,11), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = c(-0.12, 0.4)

for(i in 1:d.subs){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.subs.oaat[ix,i], y.subs.oaat$mean[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE, lwd = 2)
  #points(X.oaat[ix,i], y.oaat.box$mean[ix], col = 'red', type = 'l')
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
dev.off()







