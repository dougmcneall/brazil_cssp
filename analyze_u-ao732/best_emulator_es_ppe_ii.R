# best_emulator_es_ppe_ii.R
# Establish the best emulator to use for 
# the Earth system ensemble.


# -------------------------------------------------------------------
# 0. Load functions
# -------------------------------------------------------------------
source('../per_pft.R')

load_ts_ensemble = function(fn, na.strings='-9.990000000000000000e+02', skip=1){
  dat = read.table(fn, header = FALSE, skip = skip, na.strings=na.strings)
  dat
}

anomalizeTSmatrix = function(x, ix){
  # Anomalise a matrix of timeseries
  subx = x[ ,ix]
  sweepstats = apply(subx, 1, FUN=mean)
  anom = sweep(x, 1, sweepstats, FUN = '-')
  anom
}

ts.ensemble.change = function(x, startix, endix){
  # Calculate the change in an ensemble of timeseries
  start.subx = x[ ,startix]
  start.stats = apply(start.subx, 1, FUN=mean)
  end.subx = x[ ,endix]
  end.stats = apply(end.subx, 1, FUN=mean)
  out = end.stats - start.stats
  out
}

years = 1861:2014

# Load up the data
lhs_i = read.table('data/ES_PPE_i/lhs_u-ao732.txt', header = TRUE)
lhs_ii = read.table('data/ES_PPE_ii/lhs_u-ao732a.txt', header = TRUE)

# Exclude the last 100 members of the ensemble for
# Validation of the analysis later on.
toplevel.ix = 1:400

# load up broadleaf forest
# Ensemble members P0000 to P0498, with the standard run in the 
# final row.
Amazon.area = 7.013467e+12 # m^2
frac_bl = read.table("data/ES_PPE_ii/Amazon_forest_total.txt")[toplevel.ix, -1] / Amazon.area
matplot(t(frac_bl), type = 'l', col = 'black', lty = 'solid')

frac_bl_change = anomalizeTSmatrix(frac_bl, ix = 1:10)
matplot(t(frac_bl_change), type = 'l', col = 'black', lty = 'solid')
abline(h = 0, col = 'white')

nc.precip <- nc_open("data/ES_PPE_ii/JULES-ES.0p92.vn5.0.CRUNCEPv7.P0199.Annual.Amazon.precip.global_sum.nc")
precip <- ncvar_get(nc.precip)

lhs = rbind(lhs_i, lhs_ii)[toplevel.ix, ]
X = normalize(lhs)
colnames(X) = colnames(lhs)
d = ncol(X)

fnvec = dir('data/ES_PPE_ii', pattern = 'Annual.Amazon')
fnlocvec = paste0('data/ES_PPE_ii/', fnvec)

# Extract the bit of the filename that describes the data
fs = lapply(fnvec,regexpr, pattern = 'Annual.Amazon')
fe = lapply(fnvec,regexpr, pattern = 'global_sum.txt')

fnams = rep(NA, length(fnvec))
for(i in 1:length(fnvec)){
  fnams[i] = substr(fnvec[i], attr(fs[[1]], 'match.length')+2, fe[[i]][1]-2)
}

runoff.raw = (load_ts_ensemble("data/ES_PPE_ii/Annual.Amazon.runoff.global_sum.txt"))[toplevel.ix, ]
runoff.norm = sweep(runoff.raw, 2, STATS = precip, FUN = '/')
runoff.norm.anom = anomalizeTSmatrix(runoff.norm, ix = 1:10)

pdf('graphics/ppe_ii/runoff_normalised.pdf', width =7, height = 7 )
par(mfrow = c(2,1), mar = c(4,4,2,1))
matplot(t(runoff.norm), type = 'l', col = linecols, main = 'normalised runoff')
matplot(t(runoff.norm.anom), type = 'l', col = linecols, main = 'normalised runoff anomaly')
dev.off()


runoff = (load_ts_ensemble("data/ES_PPE_ii/Annual.Amazon.runoff.global_sum.txt")/1e8)[toplevel.ix, ]
runoff.anom = anomalizeTSmatrix(runoff, ix = 1:10)

pdf('graphics/ppe_ii/runoff.pdf', width =7, height = 7 )
par(mfrow = c(2,1), mar = c(4,4,2,1))
matplot(t(runoff), type = 'l', col = linecols, main = 'runoff')
matplot(t(runoff.anom), type = 'l', col = linecols, main = 'runoff anomaly')
dev.off()

runoff.ix = which(runoff[,1] > 0.8)

X.runoff = X[runoff.ix, ]

# There's a relationship between runoff starting value and runoff change, but
# not sure about the causality.
runoff.start = runoff[runoff.ix, 1]
runoff.change = ts.ensemble.change(runoff[runoff.ix, ], 1:10, 145:154)


plot(runoff.start, runoff.change)

km.fit.runoff.start = km(~., design = X.runoff, 
                         response = runoff.start, nugget.estim = TRUE)

km.loo.runoff.start = leaveOneOut.km(km.fit.runoff.start,
                               type = 'UK', trend.reestim=TRUE)

km.fit.runoff.change = km(~., design = X.runoff, response = runoff.change)
km.loo.runoff.change = leaveOneOut.km(km.fit.runoff.change,
                               type = 'UK', trend.reestim=TRUE)

dev.new(width = 10, height = 5)
par(mfrow = c(1,2))
plot(runoff.start, km.loo.runoff.start$mean, col = 'grey', pch = 19)
abline(0,1)

plot(runoff.change, km.loo.runoff.change$mean, col = 'grey', pch = 19)
abline(0,1)


# ---------------------------------------------------------
# Lasso section
# ---------------------------------------------------------
library(glmnet)

fit.glmnet.cv = cv.glmnet(x=X.runoff,y = runoff.start)
coef(fit.glmnet.cv, s = "lambda.1se")

# The labels of the retained coefficients are here
# (retains intercept at index zero)
coef.i = (coef(fit.glmnet.cv, s = "lambda.1se"))@i

labs = labels(coef(fit.glmnet.cv, s = "lambda.1se"))[[1]]
labs = labs[-1] # remove intercept
glmnet.retained = labs[coef.i]

start.form = as.formula(paste("~ ", paste(glmnet.retained , collapse= "+")))

library(foreach)
# below an example for a computer with 2 cores, but also work with 1 core
nCores = 2
library(doParallel)

cl = makeCluster(nCores) 
registerDoParallel(cl)
# kriging model 1, with 4 starting points
km.glmnet.fit.runoff.start = km(start.form, design=X.runoff, response=runoff.start,
                                control = list(maxit = 300), multistart = 2)
stopCluster(cl)

km.glmnet.loo.runoff.start = leaveOneOut.km(km.glmnet.fit.runoff.start,
                                     type = 'UK', trend.reestim=TRUE)

plot(runoff.start, km.loo.runoff.start$mean, col = 'grey', pch = 19)
points(runoff.start, km.glmnet.loo.runoff.start$mean, col = 'tomato', pch = 19)
abline(0,1)


rmse(runoff.start,km.loo.runoff.start$mean )
rmse(runoff.start,km.glmnet.loo.runoff.start$mean )

m0 = km(y ~ 1, design=X.runoff, response=runoff.start, nugget.estim = TRUE)
coefs0 = m0@covariance@range.val 

coefs1 =
m1 = km(start.form, design=X.runoff, response=runoff.start,
        parinit=coefs0, nugget.estim = TRUE)
m1.loo = leaveOneOut.km(m1, type = 'UK', trend.reestim=TRUE)
rmse(runoff.start,m1.loo$mean)


plot(runoff.start, km.loo.runoff.start$mean, col = 'grey', pch = 19)
points(runoff.start, km.glmnet.loo.runoff.start$mean, col = 'tomato', pch = 19)
points(runoff.start, m1.loo$mean, col = 'skyblue', pch = 19)
abline(0,1)

plot(km.glmnet.loo.runoff.start$mean - m1.loo$mean)

m0 = km(y ~., design=X.norm, response=y, nugget.estim = TRUE)

# Build the first emulator with a flat prior
m0 = km(y ~ 1, design=X.norm, response=y, nugget.estim = TRUE)
m0 = km(y ~., design=X.norm, response=y, nugget.estim = TRUE)
coefs0 = m0@covariance@range.val 
m1 = km(start.form, design=X.norm, response=y,
        parinit=coefs0, nugget.estim = TRUE)
m2 = km(~., design=X.norm, response=y,
        parinit=coefs0, nugget.estim = TRUE)
nd = matrix(rep(0.5, 74), nrow = 1)
colnames(nd) = colnames(X.norm)

predict(fit.glmnet.cv, s = "lambda.1se", newx = nd)

twoStep.glmnet = function(X, y, nugget=NULL, nuggetEstim=FALSE, noiseVar=NULL, seed=NULL, trace=FALSE, maxit=100,
                   REPORT=10, factr=1e7, pgtol=0.0, parinit=NULL, popsize=100){
  # Use lasso to reduce input dimension of emulator before
  # building.
  control_list = list(trace=trace, maxit=maxit, REPORT=REPORT, factr=factr, pgtol=pgtol, pop.size=popsize)
  xvars = colnames(X)
  data = data.frame(response=y, x=X)
  colnames(data) <- c("response", xvars)
  nval = length(y)
  
  # fit a lasso by cross validation
  library(glmnet)
  fit.glmnet.cv = cv.glmnet(x=X,y=y)

  # The labels of the retained coefficients are here
  # (retains intercept at index zero)
  coef.i = (coef(fit.glmnet.cv, s = "lambda.1se"))@i
  labs = labels(coef(fit.glmnet.cv, s = "lambda.1se"))[[1]]
  labs = labs[-1] # remove intercept
  glmnet.retained = labs[coef.i]
  
  start.form = as.formula(paste("~ ", paste(glmnet.retained , collapse= "+")))
  m = km(start.form, design=X, response=y, nugget=nugget, parinit=parinit,
         nugget.estim=nuggetEstim, noise.var=noiseVar, control=control_list)
  
  return(list(x=X, y=y, nugget=nugget, nugget.estim=nuggetEstim,
              noise.var=noiseVar, emulator=m, seed=seed, coefs=m@covariance@range.val,
              trends=m@trend.coef, meanTerms=all.vars(start.form), fit.glmnet.cv=fit.glmnet.cv))
}

test = twoStep.glmnet(X = X.runoff, y = runoff.change)

ts.loo.runoff.change = leaveOneOut.km(test$emulator,
                                      type = 'UK', trend.reestim=TRUE)

plot(runoff.change, km.loo.runoff.change$mean, col = 'grey', pch = 19)
points(runoff.change, ts.loo.runoff.change$mean, col = 'tomato', pch = 19)
abline(0,1)

rmse(runoff.change, km.loo.runoff.change$mean)
rmse(runoff.change, ts.loo.runoff.change$mean)


test.sw = twoStep(X = X.runoff, y = runoff.start)




