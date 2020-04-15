# Functions for es_ppe_ii_viz.rmd
  
# Example of adding  density plots to a pairs plot
dfunc.up <- function(x, y, dfunc.col = greys, ...){
  require(MASS)
  require(RColorBrewer)
  
  #rb = brewer.pal(9, "RdBu")
  #br  = rev(rb)
  
  # function for plotting 2d kernel density estimates in pairs() plot.
  kde = kde2d(x,y)
  image(kde, col = dfunc.col, add = TRUE)
}


load_ts_ensemble = function(fn, na.strings='-9.990000000000000000e+02', skip=1){
  # Load an ensemble of time series data.
  dat = read.table(fn, header = FALSE, skip = skip, na.strings=na.strings)
  dat
}


parcoord.notsilly = function (x, rx, col = 1, lty = 1, var.label = FALSE, ...) 
{
  # Function that doesn't normalise parallel coordinates plot ranges.
  # (although there are side effects)
  #   rx <- apply(x, 2L, range, na.rm = TRUE)
  #   x <- apply(x, 2L, function(x) (x - min(x, na.rm = TRUE))/(max(x, 
  #       na.rm = TRUE) - min(x, na.rm = TRUE)))
  matplot(1L:ncol(x), t(x), type = "l", col = col, lty = lty, 
          xlab = "", ylab = "", axes = FALSE, ...)
  axis(1, at = 1L:ncol(x), labels = colnames(x))
  for (i in 1L:ncol(x)) {
    lines(c(i, i), c(0, 1), col = "grey70")
    if (var.label) 
      text(c(i, i), c(0, 1), labels = format(rx[, i], digits = 3), 
           xpd = NA, offset = 0.3, pos = c(1, 3), cex = 0.7)
  }
  invisible()
}


makeTransparent<-function(someColor, alpha=100)
  # Transparent colours for plotting
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

allin = function(x, mins, maxes){
  # are all the elements of a vector in range?
  all(x > mins & x < maxes)
}

normalize.na = function(X, wrt = NULL){
  # Normalize with respect to a matrix.
  # This version handles NAs
  
  f <- function(X){
    (X-min(X, na.rm = TRUE))/(max(X, na.rm = TRUE)-min(X, na.rm = TRUE))
  }
  
  # test to see if we have a matrix, array or data frame
  if(length(dim(X))==2){
    out <- apply(X,2,f)
  }
  
  else{	
    out <- f(X)
  }
  
  if(is.null(wrt) == FALSE){
    # if argument wrt is given
    
    n <- nrow(X)
    mmins <- t(kronecker(apply(wrt,2,min, na.rm = TRUE),t(rep(1,n))))
    mmaxs <- t(kronecker(apply(wrt,2,max, na.rm = TRUE),t(rep(1,n))))
    
    out <- (X-mmins)/(mmaxs-mmins)
    
  }
  
  out
}

reset <- function() {
  # Allows annotation of graphs, resets axes
  par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
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


twoStep.sens = function(X, y, n=21, predtype = 'UK', nugget=NULL, nuggetEstim=FALSE, noiseVar=NULL, seed=NULL, trace=FALSE, maxit=100,
                        REPORT=10, factr=1e7, pgtol=0.0, parinit=NULL, popsize=100){
  # Sensitivity analysis with twoStep emulator. 
  # Calculates the variance of the output varied one at a time across each input.
  d = ncol(X)
  X.norm = normalize(X)
  X.oaat = oaat.design(X.norm, n, med = TRUE)
  colnames(X.oaat) = colnames(X)
  
  twoStep.em = twoStep.glmnet(X=X, y=y, nugget=nugget, nuggetEstim=nuggetEstim, noiseVar=noiseVar,
                              seed=seed, trace=trace, maxit=maxit,
                              REPORT=REPORT, factr=factr, pgtol=pgtol,
                              parinit=parinit, popsize=popsize)
  
  oaat.pred = predict(twoStep.em$emulator, newdata = X.oaat, type = predtype)
  
  sens = sensvar(oaat.pred = oaat.pred, n=n, d=d)
  out = sens
  out
}

anomalizeTSmatrix = function(x, ix){
  # Anomalise a timeseries matrix
  subx = x[ ,ix]
  sweepstats = apply(subx, 1, FUN=mean)
  anom = sweep(x, 1, sweepstats, FUN = '-')
  anom
}

ts.ensemble.change = function(x, startix, endix){
  # Calculate how much a timeseries changes
  start.subx = x[ ,startix]
  start.stats = apply(start.subx, 1, FUN=mean)
  
  end.subx = x[ ,endix]
  end.stats = apply(end.subx, 1, FUN=mean)
  
  out = end.stats - start.stats
  out
}

ensTShist <- function(x, dat,grid = TRUE,colvec,histcol, mainvec,...){
  # Plot comparison ensemble time series
  # add a histogram on the end
  #source('/home/h01/hadda/code/R/useful/dougplotpars.R')  
  #par(dougpar_web)
  par(mar = c(5,5,4,0), mgp = c(3.5,1,0))
  nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),widths = c(10,2), TRUE)
  
  matplot(x, t(dat), type = 'l',
          lty = 1,
          col = colvec,
          axes = FALSE,
          ...
  )
  axis(1)
  axis(2)
  # matlines(colnames(dat2), t(dat2), type = 'l',
  #           lty = 1,
  #           col = colvec[2]
  #  )
  if(grid) {grid(lty ='dashed',col = 'grey')}
  mtext(side = 3, line = 1, adj = 0, mainvec, cex = 1.5, col = 'black')
  #  legend('topleft', legvec,
  #         fill = colvec,
  #         text.col = colvec,
  #         bg = 'white',
  #         border = par()$fg,
  #         cex = 1.5
  #  )
  # Add the histograms
  datran <- range(dat, na.rm = TRUE)
  breaks <- seq(from = datran[1], to = datran[2], length = 20)
  datHist <- hist( dat[,ncol(dat)],breaks = breaks, plot = FALSE)
  #dat2Hist <- hist( dat2[,ncol(dat2)],breaks = breaks,  plot = FALSE)
  
  xlim = c(0, max(datHist$counts, na.rm = TRUE))
  par(mar = c(5,0,4,1), fg = 'white')
  barplot(datHist$counts, horiz = TRUE, col = histcol, space = 0, axes = FALSE, xlim = xlim)
  #barplot(dat2Hist$counts, horiz = TRUE, col = colvec[2], space = 0, axes = FALSE, xlim = xlim)
}
