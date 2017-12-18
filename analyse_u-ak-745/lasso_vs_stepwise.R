# lasso_vs_stepwise.R
# Testing the lasso vs stepwise regression in selecting inputs.
#

setwd("/Users/dougmcneall/Documents/work/R/brazil_cssp/analyse_u-ak-745")
## Try and build an emulator
library(devtools)

#install_github(repo = "dougmcneall/hde")

library(hde)
library(RColorBrewer)
library(fields)
library(MASS)
library(DiceKriging)
library(coefplot)

source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/imptools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/vistools.R")

yg = brewer.pal(9, "YlGn")
ryb = brewer.pal(9, "RdYlBu")
byr = rev(ryb)


# Here is the leaf area index data
setwd("/Users/dougmcneall/Documents/work/R/brazil_cssp/analyse_u-ak-745")
source("/Users/dougmcneall/Documents/work/R/brazil_cssp/per_pft.R")


filelist.global.vegfrac <- paste0('frac_area_mean/global_area_mean_PFT',0:16, '.txt')
filelist.wus.vegfrac <- paste0('frac_area_mean/WUS_area_mean_PFT',0:16, '.txt')
filelist.sam.vegfrac <- paste0('frac_area_mean/SAM_area_mean_PFT',0:16, '.txt')

filelist.global.lai <- paste0('lai_area_mean/global_area_mean_lai_PFT',0:12, '.txt')
filelist.wus.lai <- paste0('lai_area_mean/WUS_area_mean_lai_PFT',0:12, '.txt')
filelist.sam.lai <- paste0('lai_area_mean/SAM_area_mean_lai_PFT',0:12, '.txt')

c.df <- function(fn){
  out = c(read.table(fn, header = FALSE, skip = 1), recursive = TRUE)
  out
}

# creates a list of (PFTs) long, each element of which has (ensemble members) data points
global_area_means_vegfrac.list <- lapply(filelist.global.vegfrac, c.df)
global_area_means_vegfrac_standard <- as.numeric(
  readLines('frac_area_mean/global_area_mean_PFTs_standard.txt'))

wus_area_means_vegfrac.list <- lapply(filelist.wus.vegfrac, c.df)
wus_area_means_vegfrac_standard <- as.numeric(
  readLines('frac_area_mean/WUS_area_mean_PFTs_standard.txt'))

sam_area_means_vegfrac.list <- lapply(filelist.sam.vegfrac, c.df)
sam_area_means_vegfrac_standard <- as.numeric(
  readLines('frac_area_mean/SAM_area_mean_PFTs_standard.txt'))


# creates a list of (PFTs) long, each element of which has (ensemble members) data points
global_area_means_lai.list <- lapply(filelist.global.lai, c.df)
#global_area_means_standard <- as.numeric(
#  readLines('lai_area_mean/global_area_mean_lai_PFTs_standard.txt'))

wus_area_means_lai.list <- lapply(filelist.wus.lai, c.df)
#wus_area_means_standard <- as.numeric(
#  readLines('lai_area_mean/WUS_area_mean_lai_PFTs_standard.txt'))

sam_area_means_lai.list <- lapply(filelist.sam.lai, c.df)
#sam_area_means_standard <- as.numeric(
#  readLines('lai_area_mean/SAM_area_mean_lai_PFTs_standard.txt'))

# Example mean absolute error data
blet_sam_mae = c(read.table('mean_abs_anom_BLE_Trop_sam.txt', header = FALSE, skip = 1), recursive = TRUE)

lhs <- read.table('lhs_u-ak745.txt', header = TRUE)
d <- ncol(lhs)
cn <- colnames(lhs)
X.norm <- normalize(lhs)

x <- X.norm
y <- sqrt(c(sam_area_means_vegfrac.list[[2]], recursive = TRUE))
#y <- blet_sam_mae

dat = data.frame(y=y, x=X.norm)
colnames(dat) = c('y', colnames(X.norm))


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

#twostep = twoStep(X, y, nugget = 0.001, nuggetEstim = FALSE, maxit = 1000)


true.loo = function(X,y,type){
  
  out.mean = rep(NA, length(y))
  #out.sd = rep(NA, length(y))
  
  for(i in 1:nrow(X)){
    X.trunc = X[-i, ]
    y.trunc = y[-i]
    
    X.target = matrix(X[i, ], nrow = 1)
    colnames(X.target) <- colnames(X)
    X.target.df = data.frame(X.target)
    
    if(type == 'lm'){
      data.trunc = data.frame(y.trunc, X.trunc)
      colnames(data.trunc) = c('y', colnames(X.trunc))
      fit = lm(y~., data = data.trunc)
      pred = predict(fit, newdata = X.target.df)
      out.mean[i] = pred
    }
    
    else if(type == 'sw'){
      data.trunc = data.frame(y.trunc, X.trunc)
      colnames(data.trunc) = c('y', colnames(X.trunc))
      fit.lm = lm(y~., data = data.trunc)
      fit.sw = step(fit.lm, direction="both", k=log(length(y)), trace=TRUE)
      pred = predict(fit.sw, newdata = X.target.df)
      out.mean[i] = pred
    }
    
    else if(type == 'lasso'){
      fit.glmnet.cv = cv.glmnet(X.trunc,y.trunc)
      pred = predict(fit.glmnet.cv, s = "lambda.1se", newx = X.target)
      out.mean[i] = pred
    }
    
    else if(type == 'gp'){
      fit = km(~., design = X.trunc, response = y.trunc)
      pred = predict(fit,newdata = X.target, type = 'UK')
      out.mean[i] = pred$mean
      #out.sd[i = pred$sd]
    }
    
    else if(type == 'gpflat'){
      fit = km(~1, design = X.trunc, response = y.trunc)
      pred = predict(fit,newdata = X.target, type = 'UK')
      out.mean[i] = pred$mean
      #out.sd[i = pred$sd]
    }
    
    else if(type == 'twoStep'){
      fit = twoStep(X = X.trunc, y = y.trunc, nuggetEstim = TRUE, maxit = 1000)
      pred = predict(fit$emulator,newdata = X.target, type = 'UK')
      out.mean[i] = pred$mean
      #out.sd[i = pred$sd]
    }
    
  }
  return(list(mean = out.mean))
}  


# There looks to be a pretty strong bias/variance trade off going on.
# The lasso has much less variance, but is pretty heavily biased at the upper
# and lower ranges.
test.lm = true.loo(X = X.norm, y = y, type = 'lm')
test.sw = true.loo(X = X.norm, y = y, type = 'sw')
test.lasso = true.loo(X = X.norm, y = y, type = 'lasso')
test.gp = true.loo(X = X.norm, y = y, type = 'gp')
test.gpflat = true.loo(X = X.norm, y = y, type = 'gpflat')
test.twoStep = true.loo(X = X.norm, y = y, type = 'twoStep')
# How does cross validation perform compared to stepwise?
# How are the coefficients rlative to setepwise?



rmse = function(y, ydash){
  
  out = sqrt(   mean( (y-ydash)^2, na.rm = TRUE)   )
  out
}

rmse(y^2, test.lm$mean^2 )
rmse(y^2, test.sw$mean^2 )
rmse(y^2, test.lasso$mean^2 )
rmse(y^2, test.gp$mean^2 )
rmse(y^2, test.gpflat$mean^2 )
rmse(y^2, test.twoStep$mean^2 )


plot(y, test.lm$mean, col = 'black', pch = 19,
     xlim = c(0, 1), ylim = c(0,1))

points(y, test.sw$mean, col = 'red', pch = 19)
points(y, test.gp$mean, col = 'orange', pch = 19)
points(y, test.lasso$mean, col = 'blue', pch = 19)
points(y, test.gpflat$mean, col = 'purple', pch = 19)
points(y, test.twoStep$mean, col = 'grey', pch = 19)
abline(0,1)


# Stepwise section
fit.lm  = lm(y~., data = dat)
fit.steplm = step(fit.lm, direction="both", k=log(length(y)), trace=TRUE)

form = as.formula(fit.steplm)

start.form = as.formula(paste("~ ", paste(labels, collapse= "+")))
# Feed both the lasso and the Stepwise fits into a Gaussian
# process to see how they perform.

# Is "twoStep" (stepwise then GP) better than everything else?
# Are there any data transformations which might help?

# Lasso section
fit.glmnet.cv = cv.glmnet(x,y)
coef(fit.glmnet.cv, s = "lambda.1se")

# The labels of the retained coefficients are here
# (retains intercept at index zero)
coef.i = (coef(fit.glmnet.cv, s = "lambda.1se"))@i

labs = labels(coef(fit.glmnet.cv, s = "lambda.1se"))[[1]]
labs = labs[-1] # remove intercept
glmnet.retained = labs[coef.i]

start.form = as.formula(paste("~ ", paste(glmnet.retained , collapse= "+")))

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






