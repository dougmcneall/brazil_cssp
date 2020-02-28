# es_ppe_ii_history_match_common.R
# History matching for the ESPPEii ensemble.
# This script loads data and functions used in common in the other analyses.
# Uses data in  data/ES_PPE_ii_test/, which removed a bug where
# rows of the data were repeated rather than NAs being installed
# 

#setwd('analyze_u-ao732')
source('../per_pft.R')

#------------------------------------------------------------------------------
# Local functions
#------------------------------------------------------------------------------

load_ts_ensemble = function(fn, na.strings='-9.990000000000000000e+02', skip=1){
  # Load an ensemble of time series data.
  dat = read.table(fn, header = FALSE, skip = skip, na.strings=na.strings)
  dat
}

rx = rbind(rep(0,32), rep(1,32))

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

normalize.wrt <- function(a,b){
  ## n.wrt(x)
  ##
  ## A function that takes a matrix a and normalizes each column,
  ## from zero to 1, depending upon the minimum and maximum of
  ## matrix b.
  ##
  ## D.McNeall 27th November 2006
  ## I think this came from RKSH!
  ##
  n <- nrow(a)
  mmins <- t(kronecker(apply(b,2,min),t(rep(1,n))))
  mmaxs <- t(kronecker(apply(b,2,max),t(rep(1,n))))
  
  (a-mmins)/(mmaxs-mmins)
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

# -----------------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------------

years = 1861:2014
ysec = 60*60*24*365
norm.vec = c(1e12, 1e12, 1e12/ysec , 1e12, 1e12, 1e9)

# Load up the data
lhs_i = read.table('data/ES_PPE_i/lhs_u-ao732.txt', header = TRUE)
lhs_ii = read.table('data/ES_PPE_ii/lhs_u-ao732a.txt', header = TRUE)

toplevel.ix = 1:499

lhs = rbind(lhs_i, lhs_ii)[toplevel.ix, ]

X = normalize(lhs)
colnames(X) = colnames(lhs)
d = ncol(X)

# --------------------------------------------------------------------------------
# Apply constraints to the input space by history matching with global data
#
# npp 35-80 GtC
# nbp > 0
# cVeg 300 - 800 GtC
# cSoil 750 - 3000 GtC
# --------------------------------------------------------------------------------
#fnvec = c("Annual.cs_gb.global_sum.txt",
#          "Annual.cv.global_sum.txt",
#          "Annual.gpp_gb.global_sum.txt",
#          "Annual.nbp.global_sum.txt",
#          "Annual.npp_n_gb.global_sum.txt",
#          "Annual.runoff.global_sum.txt")

fnallvec = dir('data/ES_PPE_ii_test/', pattern = 'Annual')
# WARNING - hard coded hack to sort
fidx = grep("Annual.(?!Amazon).*", fnallvec, perl=TRUE)
fnvec_interim = fnallvec[fidx]
fidx2 = grep("sum.(?!standard).*", fnvec_interim, perl=TRUE)
fnvec = fnvec_interim[fidx2]
fnlocvec = paste0('data/ES_PPE_ii_test/', fnvec)

# Extract the bit of the filename that describes the data
fs = lapply(fnvec,regexpr, pattern = 'Annual')
fe = lapply(fnvec,regexpr, pattern = 'global_sum.txt')

fnams = rep(NA, length(fnvec))
for(i in 1:length(fnvec)){
  fnams[i] = substr(fnvec[i], attr(fs[[1]], 'match.length')+2, fe[[i]][1]-2)
}

# Constrain on runoff first
datmat.raw = matrix(nrow = nrow(X), ncol = length(fnlocvec))
for(i in 1:length(fnlocvec)){
  dat = load_ts_ensemble(fnlocvec[i])[toplevel.ix, ]
  dat.modern = dat[ ,135:154]
  mean.modern = apply(dat.modern, 1, mean)
  datmat.raw[ , i] = mean.modern
}
colnames(datmat.raw) = fnams

dat.norm = sweep(datmat.raw, 2, norm.vec, FUN = '/')
p = ncol(dat.norm)
