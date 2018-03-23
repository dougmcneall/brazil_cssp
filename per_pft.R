# per_pft.R
# A set of functions for analysing data from a 'per-pft' ensemble of JULES

library(devtools)
#install_github(repo = "dougmcneall/hde")

library(hde)
library(RColorBrewer)
library(fields)
library(MASS)
library(DiceKriging)

source("/Users/dougmcneall/Documents/work/R/packages-git/emtools.R")
source("/Users/dougmcneall/Documents/work/R/packages-git/imptools.R")
source("/Users/dougmcneall/Documents/work/R/packages-git/vistools.R")

#source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R")
#source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/imptools.R")
#source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/vistools.R")

yg = brewer.pal(9, "YlGn")
ryb = brewer.pal(11, "RdYlBu")
byr = rev(ryb)
rb = brewer.pal(11, "RdBu")
br = rev(rb)

paired = brewer.pal(11,'Paired')
linecols = c('black', paired[1], paired[2], paired[5])

# Example of adding  density plots to a pairs plot
dfunc.up <- function(x,y,...){
  require(MASS)
  require(RColorBrewer)
  
  rb = brewer.pal(9, "RdBu")
  br  = rev(rb)
  # function for plotting 2d kernel density estimates in pairs() plot.
  kde = kde2d(x,y)
  image(kde, col = rb, add = TRUE)
}


sensvar = function(oaat.pred, n, d){
  # Calculate variance as a global sensitivity meansure
  out = rep(NA,d)
  for(i in 1:d){
    ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
    out[i] = var(oaat.pred$mean[ix])
  }
  out
}

globalresponse = function(oaat.pred, n, d){
  # Calculate the effect of turning the parameterfrom lowest to highest
  out = rep(NA,d)
  for(i in 1:d){
    ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
    out[i] = (oaat.pred$mean[ix[n]] - oaat.pred$mean[ix[1]])
  }
  out
}

globalresponse.lm = function(oaat.pred, n, d){
  # Calculate the effect of turning the parameterfrom lowest to highest
  out = rep(NA,d)
  for(i in 1:d){
    ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
    out[i] = (oaat.pred[ix[n]] - oaat.pred[ix[1]])
  }
  out
}

cutoff <- function(dat, zlim){
  # Function to cut off data at thresholds for plotting
  out <- dat
  out[dat < zlim[1]] <- zlim[1]
  out[dat > zlim[2]] <- zlim[2]
  out
}

globalsens <- function(lhs, response.list){
  # run globalresponse across a list of parameters
  d <- ncol(lhs)
  lhs.norm <- normalize(lhs)
  
  n = 21
  X.oaat = oaat.design(lhs.norm, n, med = TRUE)
  colnames(X.oaat) = colnames(lhs)
  
  oaat.globalresponse = matrix(NA, nrow = d, ncol = length(response.list))
  for(i in 1: length(global_area_means.list)){
    print(i)
    try({
      y = c(response.list[i], recursive = TRUE)
      
      dat = data.frame(y=y, x=lhs.norm)
      colnames(dat) = c('y', colnames(lhs))
      
      initfit = lm(y ~ ., data = dat)
      stepfit = step(initfit, direction="both", k=log(length(y)), trace=TRUE)
      
      n = 21
      X.oaat = oaat.design(lhs.norm, n, med = TRUE)
      colnames(X.oaat) = colnames(lhs)
      oaat.pred = predict(stepfit, newdata = data.frame(X.oaat))
      oaat.globalresponse[,i] = globalresponse.lm(oaat.pred, n = n, d = d)
    },
    silent = TRUE)
    
    # try({
    #   fit = km(~., design = lhs.norm, response = c(response.list[i], recursive = TRUE) )
    #   oaat.pred = predict(fit, newdata = X.oaat, type = 'UK')
    #   oaat.globalresponse[,i] = globalresponse(oaat.pred, n = n, d = d)
    #   },
    #   silent = TRUE)
    # }
  }
  out <- oaat.globalresponse
  out
}

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

rmse = function(y, ydash){
  
  out = sqrt(   mean( (y-ydash)^2, na.rm = TRUE)   )
  out
}


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




# load data
# All the data is in frac_area_mean

surftypes = c('BLD','BLE_Trop','BLE_Temp','NLD',
              'NLE','C3G','C3C','C3P','C4G','C4C',
              'C4P','SHD','SHE','Urban','Lake','Bare Soil',
              'Ice')

pfts = c('BLD','BLE_Trop','BLE_Temp','NLD',
         'NLE','C3G','C3C','C3P','C4G','C4C',
         'C4P','SHD','SHE') # like surftypes, but without the other surfaces

