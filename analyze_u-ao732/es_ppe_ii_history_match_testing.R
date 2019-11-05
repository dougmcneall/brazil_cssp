# es_ppe_ii_history_match_testing.R
# A more formal history match for the ESPPEii ensemble.
# Uses data in  data/ES_PPE_ii_test/, which removed a bug where
# rows of the data were repeated rather than NAs being installed
# 

setwd('analyze_u-ao732')
source('../per_pft.R')

#------------------------------------------------------------------------------
# Local functions
#------------------------------------------------------------------------------

load_ts_ensemble = function(fn, na.strings='-9.990000000000000000e+02', skip=1){
  dat = read.table(fn, header = FALSE, skip = skip, na.strings=na.strings)
  dat
}

rx = rbind(rep(0,32), rep(1,32))

parcoord.notsilly = function (x, rx, col = 1, lty = 1, var.label = FALSE, ...) 
{
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
  par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}

load_ts_ensemble = function(fn, na.strings='-9.990000000000000000e+02', skip=1){
  dat = read.table(fn, header = FALSE, skip = skip, na.strings=na.strings)
  dat
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
  subx = x[ ,ix]
  sweepstats = apply(subx, 1, FUN=mean)
  anom = sweep(x, 1, sweepstats, FUN = '-')
  anom
}

ts.ensemble.change = function(x, startix, endix){
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

# My suspicion is that the duplicated rows should be NA, and need to
# be removed from analysis as they will affect the emulator.

# Start by removing only the duplicated (e.g. later) rows, and keeping
# the initial rows



dev.new(width = 9, height = 10)
par(mfrow = c(3,2))
for(i in 1:6){
  hist(dat.norm[,i], main = fnams[i])
}

# This gives a logical vector of all the repated rows.
test = duplicated(datmat.raw) | duplicated(datmat.raw[nrow(datmat.raw):1, ])[nrow(datmat.raw):1]
#test = duplicated(dat) | duplicated(dat[nrow(datmat):1, ])[nrow(dat):1]


# ----------------------------------------------------------
# Level zero constraint
# ----------------------------------------------------------

allix = 1:(nrow(dat.norm))

# Constrain with global runoff and nbp, so that the emulator
# is not negatively affected by really bad points
level0.ix = which(dat.norm[,'runoff'] >0.5 & dat.norm[,'nbp'] > -10)
dat.level0  = dat.norm[level0.ix, ]
X.level0 = X[level0.ix, ]
lhs.level0 = lhs[level0.ix, ]

nlevel0.ix = setdiff(allix, level0.ix)
X.nlevel0 = X[nlevel0.ix, ]
lhs.nlevel0 = lhs[nlevel0.ix, ]

# These are the ensemble members which do not run (produce NA)
na.ix = which(is.na(dat.norm[,'runoff']))
X.na = X[na.ix, ]
lhs.na = lhs[na.ix, ]



lhs.rx = round(apply(lhs,2,range),1)


# Visualise the constrained space ...
mins  = apply(X.level0, 2, min)
maxes = apply(X.level0, 2, max)


# expand this to visualise the level 1 constraints
dev.new(width = 9, height = 10)
par(mfrow = c(3,2))
for(i in 1:6){
  hist(dat.norm[level0.ix,i], main = fnams[i])
}


# Expand this to include all outputs (it only presently does one)
dev.new(width = 10, height = 7)
par(mfrow = c(4, 8), mar = c(1,1,1,1))
for(i in 1:d){
  plot(lhs[ ,i],dat.norm[ ,6], axes = FALSE, xlab = '', ylab = '')
}




# Parallel Coordinates plot of NROY and ruled out members, level 0
#dev.new(width = 20, height = 9)
pdf(file = 'graphics/pcp_level0.pdf', width = 20, height = 9)
par(mfrow = c(2,1), las = 2, mar = c(7,4,4,1), cex.axis = 0.8)

parcoord.notsilly(X.level0, rx, col = makeTransparent('black', 100), ylim = c(0,1),
         main = 'level 0 NROY ensemble members', var.label = TRUE)


                  
parcoord.notsilly(X.nlevel0, rx, col = makeTransparent('black', 100), ylim = c(0,1),
         main = 'level 0 ruled out ensemble members', var.label = TRUE)
                  
parcoord.notsilly(X.na,rx, col = makeTransparent('red', 255), add = TRUE)

reset()
legend('left',
       legend = 'Did Not Run',
       inset = 0.1,
       col = 'red',
       cex = 1.2,
       lwd = 2,
       horiz = TRUE)

dev.off()



# Pairs plot of level 0 constraint
# (This is hard to interpret)
dev.new(width = 10, height = 10)
pairs(X.level0, gap = 0, lower.panel = NULL,
      labels = 1:d,
      cex.lab = 0.8,
      xlim = c(0,1), ylim = c(0,1),
      col = makeTransparent('black', 50),
      pch = 21,
      xaxt = 'n', yaxt = 'n'
      )
par(xpd = NA)
#text(0.2, 0.6, labels = paste0(colnames(lhs)[1:5],'\n'))

legend('left', legend = paste(1:d, colnames(lhs)), cex = 0.9, bty = 'n')


# Pairs plot of level 0 constraint excluded members is much easier to interpret
#dev.new(width = 10, height = 10)
pdf(file = 'graphics/pairs_nlevel0.pdf', width = 10, height = 10)
pairs(X.nlevel0, gap = 0, lower.panel = NULL,
      labels = 1:d,
      cex.lab = 0.8,
      xlim = c(0,1), ylim = c(0,1),
      col = makeTransparent('black', 50),
      pch = 21,
      xaxt = 'n', yaxt = 'n'
      )
par(xpd = NA)
#text(0.2, 0.6, labels = paste0(colnames(lhs)[1:5],'\n'))

legend('left', legend = paste(1:d, colnames(lhs)), cex = 0.9, bty = 'n')
mtext(side = 1, text = 'Level 0 ruled out ensemble members', cex = 2)
dev.off()



#dev.new(width = 20, height = 9)
#par(mfrow = c(2,1), las = 2, mar = c(7,4,4,1), cex.axis = 0.8)

#parcoord(lhs.level0, col = makeTransparent('black', 100), ylim = lhs.rx,
#         main = 'level 0 NROY ensemble members', var.label = TRUE)
                  
#parcoord(lhs.nlevel0, col = makeTransparent('black', 100),
#         main = 'level 0 ruled out ensemble members')
                  
#parcoord(lhs.na, col = makeTransparent('red', 255), add = TRUE)

#reset()#legend('left',
#       legend = 'Did Not Run',
#       inset = 0.1,
#       col = 'red',
#       cex = 1.2,
#       lwd = 2,
#       horiz = TRUE)

# Both 
dev.new(width = 20, height = 5)
par(las = 2, mar = c(7,4,4,1), cex.axis = 0.8)
parcoord(X.level0, col = makeTransparent('black', 100), ylim = c(0,1),
         main = 'level 0 NROY ensemble members', var.label = TRUE)

parcoord(X.nlevel0, col = makeTransparent('tomato2', 100), add = TRUE)

# The minimum and maximum multiplication factors that produce runs that
# comply with the level 0 constraints are therefore as follows

level0.min = apply(X.level0, 2, min)
level0.max = apply(X.level0, 2, max)

level0.min.f = apply(lhs.level0, 2, min)
level0.max.f = apply(lhs.level0, 2, max)



# Histogram of level 1 constraints
hcol = 'darkgrey'
lcol = 'black'
pdf(file = 'graphics/l1_constraint_hists.pdf', width = 8, height = 8)
#dev.new(width = 8, height = 8)
par(mfrow = c(3,2), fg = 'darkgrey', las = 1)

hist(dat.norm[,'runoff'], col = hcol, main = 'Runoff', xlab = 'Sv')
polygon(x = c(0.5, 100, 100, 0.5), y = c(0, 0, 1000, 1000),
        col = makeTransparent('tomato2'),
        border = makeTransparent('tomato2'))
hist(dat.norm[,'nbp'], col = hcol, main = 'NBP', xlab = 'GtC/year')
polygon(x = c(-10, 100, 100, -10), y = c(0, 0, 1000, 1000),
        col = makeTransparent('tomato2'),
        border = makeTransparent('tomato2'))
hist(dat.norm[,'cs_gb'], col = hcol, main = 'Soil Carbon', xlab = 'GtC')
polygon(x = c(750, 3000, 3000, 750), y = c(0, 0, 1000, 1000),
        col = makeTransparent('tomato2'),
        border = makeTransparent('tomato2'))

hist(dat.norm[,'cv'], col = hcol, main = 'Vegetation Carbon', xlab = 'GtC')
polygon(x = c(300, 800, 800, 300), y = c(0, 0, 1000, 1000),
        col = makeTransparent('tomato2'),
       border =  makeTransparent('tomato2'))
hist(dat.norm[,'npp_n_gb'], col = hcol , main = 'NPP', xlab = 'GtC/year')
polygon(x = c(35, 80, 80, 35), y = c(0, 0, 1000, 1000),
        col = makeTransparent('tomato2'),
        border = makeTransparent('tomato2'))
dev.off()



# Apply level 1 constraints
level1.ix = which(dat.norm[,'cs_gb'] > 750 & dat.norm[,'cs_gb'] < 3000 &
  dat.norm[,'cv'] > 300 & dat.norm[,'cv'] < 800 & 
  dat.norm[,'npp_n_gb'] > 35 &
  dat.norm[,'npp_n_gb'] < 80 &
  dat.norm[,'runoff'] >0.5 &
  dat.norm[,'nbp'] > -10
  )

X.level1 = X[level1.ix, ]

# Ranges of the inputs after a level 1 constraint
rn.l1 = apply(X.level1,2, range)

nlevel1.ix = setdiff(allix, level1.ix)
X.nlevel1 = X[nlevel1.ix, ]

dat.level1 = dat.norm[level1.ix, ]

# Parallel Coordinates plot of NROY and ruled out members, level 1
pdf(file = 'graphics/pcp_level1.pdf', width = 20, height = 9)
par(mfrow = c(2,1), las = 2, mar = c(7,4,4,1), cex.axis = 0.8)
parcoord.notsilly(X.level1, rx = rx, col = makeTransparent('black', 100), ylim = c(0,1), var.label = TRUE,
         main = 'level 1 NROY ensemble members')

parcoord.notsilly(X.nlevel1, rx = rx, col = makeTransparent('black', 100), ylim = c(0,1),
         main = 'level 1 ruled out ensemble members', var.label = TRUE)
dev.off()


# Pairs plot of level 1 constraint
#dev.new(width = 10, height = 10)
pdf(file = 'graphics/pairs_level1.pdf', width = 10, height = 10)
pairs(X.level1, gap = 0, lower.panel = NULL,
      labels = 1:d,
      cex.lab = 0.8,
      xlim = c(0,1), ylim = c(0,1),
      col = makeTransparent('black', 100),
      pch = 21,
      xaxt = 'n', yaxt = 'n'
      )
par(xpd = NA)

legend('left', legend = paste(1:d, colnames(lhs)), cex = 0.9, bty = 'n')
mtext(side = 1, text = 'Level 1 NROY ensemble members', cex = 2)
dev.off()


# Pairs plot of ensemble members removed by level 1 constraint
# Almost no information in this!
dev.new(width = 10, height = 10)
pairs(X.nlevel1, gap = 0, lower.panel = NULL,
      labels = 1:d,
      cex.lab = 0.8,
      xlim = c(0,1), ylim = c(0,1),
      col = makeTransparent('black', 10),
      pch = 21,
      xaxt = 'n', yaxt = 'n'
      )
par(xpd = NA)
#text(0.2, 0.6, labels = paste0(colnames(lhs)[1:5],'\n'))

legend('left', legend = paste(1:d, colnames(lhs)), cex = 0.9, bty = 'n')




# level 1 constrained marginal hgistograms
dev.new(width = 12, height = 11)
par(mfrow = c(8,4), mar = c(4,3,2,1), oma = c(0,0,2,0))
for(i in 1:ncol(X)){

hist(X.level1[,i], xlim = c(0,1), main = colnames(X)[i], xlab = '', ylab = '',
     breaks = seq(from = 0, to = 1, by = 0.05), col = 'grey', border = 'grey')
abline(v = rn.l1[1,i])
abline(v = rn.l1[2,i])
}


dev.new(width = 12, height = 11)
par(mfrow = c(8,4), mar = c(4,3,2,1), oma = c(0,0,2,0))
for(i in 1:ncol(X)){

hist(X.nlevel1[,i], xlim = c(0,1), main = colnames(X)[i], xlab = '', ylab = '',
     breaks = seq(from = 0, to = 1, by = 0.1), col = 'grey', border = 'grey')
}

# -----------------------------------------------------------------
# How good is the emulator when using the level0 constraint?
# Use a leave-one-out metric
# -----------------------------------------------------------------


level0.emlist = vector('list',length(fnams))
level0.loolist = vector('list',length(fnams))

for(i in 1:ncol(dat.level0)){
  
  em = twoStep.glmnet(X = X.level0, y = dat.level0[,i])
  level0.emlist[[i]] = em
  loo = leaveOneOut.km(model = em$emulator, type = 'UK', trend.reestim=FALSE)
  level0.loolist[[i]] = loo 
}

# leave one out error of a raw GP (rather than glmnet dimension reduction)
level0.rawemlist = vector('list',length(fnams))
level0.rawloolist = vector('list',length(fnams))

for(i in 1:ncol(dat.level0)){
  
  em = km(~.,design = X.level0, response = dat.level0[,i], nugget.estim = TRUE)
  
  level0.rawemlist[[i]] = em
  loo = leaveOneOut.km(model = em, type = 'UK', trend.reestim=FALSE)
  level0.rawloolist[[i]] = loo
}


# Plot the leave-one-out predictions against the observed model output
# for all the outputs.
dev.new(width = 15, height = 10)
par(mfrow = c(2,3))

for(i in 1:p){

lwr = level0.loolist[[i]]$mean-(2*level0.loolist[[i]]$sd)
upr = level0.loolist[[i]]$mean+(2*level0.loolist[[i]]$sd)
ylim = range(c(lwr, upr))

plot(dat.level0[,i], level0.loolist[[i]]$mean, pty = 'n', ylim = ylim,
     main = colnames(dat.level0)[i], xlab = 'observed', ylab = 'predicted')

segments(dat.level0[,i], lwr,
         dat.level0[,i],upr,
         col = makeTransparent('black', 100)
         )
points(dat.level0[,i], level0.loolist[[i]]$mean, pch = 20, col = 'black')
abline(0,1)

}


# Plot the leave-one-out predictions against the observed model output
# for all the outputs.
dev.new(width = 15, height = 10)
par(mfrow = c(2,3))

for(i in 1:p){

lwr = level0.rawloolist[[i]]$mean-(2*level0.rawloolist[[i]]$sd)
upr = level0.rawloolist[[i]]$mean+(2*level0.rawloolist[[i]]$sd)
ylim = range(c(lwr, upr))

plot(dat.level0[,i], level0.rawloolist[[i]]$mean, pty = 'n', ylim = ylim,
     main = colnames(dat.level0)[i], xlab = 'observed', ylab = 'predicted')

segments(dat.level0[,i], lwr,
         dat.level0[,i],upr,
         col = makeTransparent('black', 100)
         )
points(dat.level0[,i], level0.rawloolist[[i]]$mean, pch = 20, col = 'black')
abline(0,1)

}

# It looks as though vegetation carbon is a little underestimated lower down.
# It might be good to have a look lower down.
# There is some structure in the errors, so maybe worth looking at (e.g.)
# vegetation carbon emulator.

dev.new(width = 15, height = 10)
par(mfrow = c(2,3))

for(i in 1:ncol(dat.level0)){

  ix.bymean  = sort(level0.loolist[[i]]$mean,index.return = TRUE)
  err = level0.loolist[[i]]$mean - dat.level0[,i]
  datlength = length(err)
  
  lwr = err - (2*level0.loolist[[i]]$sd)
  upr = err + (2*level0.loolist[[i]]$sd)
  
  ylim = range(c(lwr, upr))
  
  plot(1:datlength, err[ix.bymean$ix], pty = 'n',
       ylim = ylim,
       xlab = 'index',
       ylab = 'error',
       main = colnames(dat.level0)[i]
       )
  
  segments(1:datlength, lwr[ix.bymean$ix],
           1:datlength, upr[ix.bymean$ix],
           col = makeTransparent('black', 100)
           )
  
  points(1:datlength, err[ix.bymean$ix], pch = 19)
  abline(h = 0, col = 'lightgrey')  
}

# How accurate are the emulators using level zero data?

em.mae.level0 = rep(NA, ncol(dat.level0))
raw.em.mae.level0 = rep(NA, ncol(dat.level0))
for(i in 1:ncol(dat.level0)){

em.mae.level0[i] = mean(abs(level0.loolist[[i]]$mean - dat.level0[,i]))
raw.em.mae.level0[i] = mean(abs(level0.rawloolist[[i]]$mean - dat.level0[,i]))

}

mean.level0 = apply(dat.level0, 2, mean)
prop.mae.level0 = (em.mae.level0 / c(mean.level0)) * 100
prop.raw.mae.level0 = (raw.em.mae.level0 / c(mean.level0)) * 100


# ---------------------------------------------------------------
# Now, how accurate are the level1 constrained emulators?
# ---------------------------------------------------------------

level1.emlist = vector('list',length(fnams))
level1.loolist = vector('list',length(fnams))

for(i in 1:p){
  print(i) 
  em = twoStep.glmnet(X = X.level1, y = dat.level1[,i])
  level1.emlist[[i]] = em
  loo = leaveOneOut.km(model = em$emulator, type = 'UK', trend.reestim=FALSE)
  level1.loolist[[i]] = loo 
}


# having a problem here - glmnet isn't keeping any of the inputs
# Not enough runs?
# Revert to a standard km model.
# (could we try stepwise?)

level1.rawemlist = vector('list',length(fnams))
level1.rawloolist = vector('list',length(fnams))

for(i in 1:ncol(dat.level1)){
  print(i)
  em = km(~., design = X.level1, response = dat.level1[,i])
  level1.rawemlist[[i]] = em
  loo = leaveOneOut.km(model = em, type = 'UK', trend.reestim=FALSE)
  level1.rawloolist[[i]] = loo  
}


em.mae.level1 = rep(NA, ncol(dat.level1))
for(i in 1:ncol(dat.level1)){
em.mae.level1[i] = mean(abs(level1.rawloolist[[i]]$mean - dat.level1[,i]))
}

# The 
mean.level1 = apply(dat.level1, 2, mean)
prop.mae.level1 = (em.mae.level1 / c(mean.level0)) * 100



# Plot the leave-one-out predictions against the observed model output
# for all the outputs.
dev.new(width = 15, height = 10)
par(mfrow = c(2,3))

for(i in 1:ncol(dat.level1)){

lwr = level1.rawloolist[[i]]$mean-(2*level1.rawloolist[[i]]$sd)
upr = level1.rawloolist[[i]]$mean+(2*level1.rawloolist[[i]]$sd)
ylim = range(c(lwr, upr))

plot(dat.level1[,i], level1.rawloolist[[i]]$mean, pty = 'n', ylim = ylim,
     main = colnames(dat.level1)[i], xlab = 'observed', ylab = 'predicted')

segments(dat.level1[,i], lwr,
         dat.level1[,i],upr,
         col = makeTransparent('black', 100)
         )
points(dat.level1[,i], level1.rawloolist[[i]]$mean, pch = 20, col = 'black')
abline(0,1)

}

#

dev.new(width = 15, height = 10)
par(mfrow = c(2,3))

for(i in 1:ncol(dat.level1)){

  ix.bymean  = sort(level1.rawloolist[[i]]$mean,index.return = TRUE)
  err = level1.rawloolist[[i]]$mean - dat.level1[,i]
  datlength = length(err)
  
  lwr = err - (2*level1.rawloolist[[i]]$sd)
  upr = err + (2*level1.rawloolist[[i]]$sd)
  
  ylim = range(c(lwr, upr))
  
  plot(1:datlength, err[ix.bymean$ix], pty = 'n',
       ylim = ylim,
       xlab = 'index',
       ylab = 'error',
       main = colnames(dat.level1)[i]
       
       )
  
  segments(1:datlength, lwr[ix.bymean$ix],
           1:datlength, upr[ix.bymean$ix],
           col = makeTransparent('black', 100)
           )
  
  points(1:datlength, err[ix.bymean$ix], pch = 19)
  abline(h = 0, col = 'lightgrey')  
}


# How well does the level0 emulator predict the level 1 data?



# what is the subset of level0.ix that is level1.ix?
level1.in.level0.ix = which(level0.ix %in% level1.ix)


l0atl1.meanlist = vector('list',length(fnams))
l0atl1.sdlist = vector('list',length(fnams))
l0atl1.errlist = vector('list',length(fnams))


for(i in 1:ncol(dat.level1)){

  # level0 loo at level 1 locations.
l0atl1.meanlist[[i]] = level0.loolist[[i]]$mean[level1.in.level0.ix]
l0atl1.sdlist[[i]]   = level0.loolist[[i]]$sd[level1.in.level0.ix]
  
}
# get these from the level0 loo data and compare.

em.mae.l10atl1 = rep(NA, ncol(dat.level1))

for(i in 1:ncol(dat.level1)){
em.mae.l10atl1[i] = mean(abs(l0atl1.meanlist[[i]] - dat.level0[level1.in.level0.ix, i]))
}

cbind(dat.level0[level1.in.level0.ix, i], dat.level1[,i])

# Normalize all of the leave-one-out errors to the range of the ensemble output.
# Chosen level 0 range to normalise to, as feels fairer.
raw.l0.range = apply(dat.level0, 2, range, na.rm = TRUE)
abs.l0.range = abs(raw.l0.range[2,] - raw.l0.range[1,])


# It appears that the level 1 emulator is better than the level 0
# emulator at the same inputs

em.mae.l0.norm = (em.mae.level0 / abs.l0.range) * 100
em.mae.l10atl1.norm = (em.mae.l10atl1 / abs.l0.range) * 100
em.mae.l1.norm = (em.mae.level1 / abs.l0.range) * 100


dev.new()
par(las = 0, mar = c(7,5,2,1))
plot(1:p, em.mae.l0.norm, pch = 19, ylim = c(0,10), type = 'o',
     axes = FALSE, xlab = '', ylab = 'mean absolute error (%)',
     main = 'Leave-one-out emulator error at different constraints'
     )
axis(1, at = 1:p, labels = colnames(dat.norm), las = 2)
axis(2)

points(1:p, em.mae.l10atl1.norm, col = 'red', pch = 19 , type = 'o')
points(1:p, em.mae.l1.norm, col = 'blue', pch = 19, type = 'o')

text(x = c(4, 4, 4), y = c(6,9,3), col = c('black', 'red', 'blue'),
     labels = c('level 0','level 0 at level1 points', 'level 1'))


# --------------------------------------------------------------------
# twoStep sensitivity analysis
# --------------------------------------------------------------------


# Ensure the outputs are all on the same scale
dat.level0.norm = normalize(dat.level0) # there is some space removed, but we'll ignore it.
sensmat.level0 = matrix(NA, ncol=d, nrow=p)

# Run over all outputs. Constrain to the "runoff approved" input space
for(i in 1:p){
  dat.const = dat.level0.norm[, i]
  ts.sens = twoStep.sens(X=X.level0, y = dat.const)
  sensmat.level0[i, ] = ts.sens
}


abssum.sensmat.level0 = apply(abs(sensmat.level0), 2, sum)

# Sensitivity matrix of runoff-constrained inputs
#pdf(file = 'graphics/ppe_ii/sensmat_level0.pdf', width = 9, height = 9)
dev.new(width = 9, height = 9)
par(mfrow = c(2,1), mar = c(8,7,3,2))
image(t(sensmat.level0), col = blues, axes = FALSE)
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/(ndat-1) ),
     labels = datnames, las = 1)

par(mar = c(8,7,0,2))
plot(1:d, abssum.sensmat.level0, axes = FALSE, xlab = '', ylab = 'Sum abs. sensitivity', pch = 19)
segments(x0 = 1:d, y0 = rep(0,d), x1 = 1:d, y1 = abssum.sensmat.level0)
axis(1,labels = colnames(lhs), at = 1:d, las=3, cex.axis = 0.8 )
axis(2, las = 1) 
#dev.off()

abssum.sensmat.level0.sort = sort(abssum.sensmat.level0, decreasing = TRUE, index.return = TRUE)

# Sorted summary sensitivity of runoff to inputs
#pdf(file = 'graphics/ppe_ii/ordered_oaat_SA_level0.pdf', width = 7, height = 5)
dev.new(width = 7, height = 5)
par(mar = c(8,4,3,1))
plot(1:d, abssum.sensmat.level0.sort$x, axes = FALSE, pch = 19,
     xlab = '', ylab = 'OAAT Sensitivity Index')
segments(x0 = 1:d, y0 = rep(0,d), x1 = 1:d, y1 = abssum.sensmat.level0.sort$x)
axis(1, at = 1:d,  labels = colnames(X.level0)[abssum.sensmat.level0.sort$ix], las = 3, cex.axis = 0.8)
axis(2,las =1)
#dev.off()



# A linear model could provide a simple sensitivity analysis for situations where
# it's difficult to build a Lasso

lm.sensmat.level0 = matrix(NA, ncol=d, nrow=p)
# Run over all outputs. Constrain to the "runoff approved" input space
for(i in 1:p){
  dat.const = dat.level0.norm[, i]
  lm.sens = lm(dat.const~X.level0)
  lm.sensmat.level0[i, ] = tail(coef(lm.sens), -1) #remove intercept
}

abssum.lm.sensmat.level0 = apply(abs(lm.sensmat.level0), 2, sum)

abssum.lm.sensmat.level0.sort = sort(abssum.lm.sensmat.level0, decreasing = TRUE, index.return = TRUE)

dev.new(width = 7, height = 5)
par(mar = c(8,4,3,1))
plot(1:d, abssum.lm.sensmat.level0.sort$x, axes = FALSE, pch = 19,
     xlab = '', ylab = 'LM Sensitivity Index')
segments(x0 = 1:d, y0 = rep(0,d), x1 = 1:d, y1 = abssum.lm.sensmat.level0.sort$x)
axis(1, at = 1:d,  labels = colnames(X.level0)[abssum.lm.sensmat.level0.sort$ix], las = 3, cex.axis = 0.8)
axis(2,las =1)


# Although the sizes of the sensitivities show some similarities, the ranks
# look quite different

dev.new()
plot(abssum.sensmat.level0, abssum.lm.sensmat.level0,
     xlab = 'LOO twoStep', ylab = 'Linear model coefficients')
abline(0,1)

dev.new()
plot(abssum.sensmat.level0.sort$ix, abssum.lm.sensmat.level0.sort$ix)
abline(0,1)


# OK, we know the
# Making decisions on choosing parameter sets. Are there any margins we
# can rule out?




#leaveOneOut.km(model, type, trend.reestim=FALSE)

#test = twoStep.glmnet(X = X.level0, y = dat.level0[,1])
  
#leaveOneOut.km(model, type, trend.reestim=FALSE)


# What are the error rates for the various emulators?
#level0.loo.cs = true.loo(X.level0, dat.level0[,1], type = 'twoStep')


#level1.loo.cs = true.loo(X.level1, dat.level1[,1], type = 'twoStep')


#level1.loo.cv = true.loo(X.level1, dat.level1[,2], type = 'twoStep')
#level1.loo.gpp = true.loo(X.level1, dat.level1[,3], type = 'twoStep')
#level1.loo.nbp = true.loo(X.level1, dat.level1[,4], type = 'twoStep')
#level1.loo.npp = true.loo(X.level1, dat.level1[,5], type = 'twoStep')
#level1.loo.runoff = true.loo(X.level1, dat.level1[,6], type = 'twoStep')

dev.new()
plot(dat.level1[,1], level1.loo.cs$mean,
     xlim = c(0,max(dat.level1[,1]))
       )
abline(0,1)

# -------------------------------------------------------------------------
# Level1 constraint using a level0 emulator.
#
# -------------------------------------------------------------------------

nsamp.unif = 100000
X.unif = samp.unif(nsamp.unif, mins = mins, maxes = maxes)

y.unif = matrix(nrow = nsamp.unif, ncol = ncol(dat.level0))
colnames(y.unif) = colnames(dat.norm)

global.emlist = vector('list',length(fnams))

for(i in 1:ncol(y.unif)){
  em = twoStep.glmnet(X = X.level0, y = dat.level0[,i])
  global.emlist[[i]] = em
  pred = predict(em$emulator, newdata = X.unif, type = 'UK')
  y.unif[,i] = pred$mean
}

# Try using a level1 emulator
# It turns out that you don't exclude as much input space if you use a
# 
#for(i in 1:ncol(y.unif)){ 
#  em = km(~., design = X.level1, response = dat.level1[,i])
#  global.emlist[[i]] = em
#  pred = predict(em, newdata = X.unif, type = 'UK')
#  y.unif[,i] = pred$mean
#}


ix.kept = which(
  y.unif[,'cs_gb'] > 750 & y.unif[,'cs_gb'] < 3000 &
  y.unif[,'cv'] > 300 & y.unif[,'cv'] < 800 & 
  y.unif[,'npp_n_gb'] > 35 & y.unif[,'npp_n_gb'] < 80 &
  y.unif[,'runoff'] >0.5 &
  y.unif[,'nbp'] > -10)
X.kept = X.unif[ix.kept, ]

# we've removed 80% of our prior input space
(nrow(X.kept) / nsamp.unif) * 100


ix.rejected = setdiff(1:nsamp.unif, ix.kept)
X.rejected = X.unif[ix.rejected, ]


#pdf(file = 'pcpTest.pdf', width = 20, height = 4)
dev.new(width = 20, height = 9)
par(las = 2, cex.axis = 0.8, mfrow = c(2,1))

parcoord.notsilly(X.kept[1:5000,], col = makeTransparent('black', 5), ylim = c(0,1))
parcoord.notsilly(X.rejected[1:5000,], col = makeTransparent('black', 5), ylim = c(0,1))

#dev.off()


dev.new(width = 12, height = 12)
pairs(X.kept[1:5000,] , pch = '.', gap = 0, upper.panel = NULL,
      xlim = c(0,1), ylim = c(0,1),
      col = makeTransparent('black', 20)
      )
 

blues = brewer.pal(9, 'Blues')
pdf(file = 'graphics/pairs_dens_level0_km.pdf', width = 10, height = 10)
#dev.new(width = 10, height = 10)
par(oma = c(0,0,0,3))
test = pairs(X.kept,
             labels = 1:d,
             gap = 0, lower.panel = NULL, xlim = c(0,1), ylim = c(0,1),
             panel = dfunc.up,
             cex.labels = 1,
             col.axis = 'white',
  dfunc.col = blues)

image.plot(legend.only = TRUE,
           zlim = c(0,1),
           col = blues,
           legend.args = list(text = 'Density of model runs matching the criteria', side = 3, line = 1),
           horizontal = TRUE
)
#par(xpd = NA)
#text(0.2, 0.6, labels = paste0(colnames(lhs)[1:5],'\n'))
legend('left', legend = paste(1:d, colnames(lhs)), cex = 0.9, bty = 'n')

dev.off()


# Next, run the one-at-a-time sensitivity measures again

# The idea is that we do a sensitivity analysis, looking only within the NROY
# space. So, we build an emulator using all or non-bonkers points, then 
# sample one-at-a-time, but without including NROY points.
#***

# Need a sensitivity analysis function that
# 1) Build multiple emulators
# 2) Allows constraints




#test = matrix(1:6, ncol = 2)
#apply(test, 1, FUN = allin, mins = c(0,0), maxes = c(2,5))

#test = matrix(1:9, ncol = 3)
#apply(test, 1, FUN = allin, mins = c(0,0,0), maxes = c(3,5,8))

constrained.oaat = function(X, Y, n.oaat = 21, mins, maxes, hold = NULL,
                            predtype = 'UK',
                            nugget=NULL, nuggetEstim=FALSE, noiseVar=NULL,
                            seed=NULL, trace=FALSE, maxit=100,
                            REPORT=10, factr=1e7, pgtol=0.0, parinit=NULL, 
                            popsize=100)
  {
  # X        ...   Experiment design
  # Y        ...   Matrix of model output, with variables in columns
  # maxes    ...   Vector maximum tolerable value of variables corresponding
  #                to columns of Y
  # mins     ...   Vector minimum tolerable value of variables corresponding
  #                to columns of Y
  
  # generate oaat design
  d = ncol(X)
  X.norm = normalize(X)
  X.oaat = oaat.design(X.norm, n = n.oaat, med = TRUE, hold = hold)
  colnames(X.oaat) = colnames(X)
  
  # generate ncol(Y) emulators
  p = ncol(Y)
  pred.mean = matrix(NA, ncol = p, nrow = nrow(X.oaat))
  pred.sd = matrix(NA, ncol = p, nrow = nrow(X.oaat))
  
  for(i in 1:p){
    
    y = Y[, i]
    
    em = twoStep.glmnet(X=X.norm, y=y, nugget=nugget, nuggetEstim=nuggetEstim,
                                  noiseVar=noiseVar,
                                  seed=seed, trace=trace, maxit=maxit,
                                 REPORT=REPORT, factr=factr, pgtol=pgtol,
                                 parinit=parinit, popsize=popsize)
    
    oaat.pred = predict(em$emulator, newdata = X.oaat, type = predtype)
    
    # produce the whole oaat emulator output
    pred.mean[, i] = oaat.pred$mean
    pred.sd[, i] = oaat.pred$sd
    
  }
  
  ix.kept = apply(pred.mean, 1, FUN = allin, mins = mins, maxes = maxes)
  
  # Replace out-of-bound rows in X with NA, to make plotting
  # the final results easier
  X.oaat.constr = X.oaat
  X.oaat.constr[ix.kept==FALSE, ] <- NA
  
  pred.constr = pred.mean
  pred.constr[ix.kept==FALSE, ] <- NA
  
  pred.sd.constr = pred.sd
  pred.sd.constr[ix.kept==FALSE, ] <- NA
  
  # keep only the oaat emulator output within constraints
  return(list(X.oaat.constr = X.oaat.constr,
              ix.kept = ix.kept,
              pred.constr = pred.constr,
              pred.sd.constr = pred.sd.constr,
              X.oaat = X.oaat,
              pred.mean = pred.mean,
              pred.sd = pred.sd)
         )
}

#                 cs   cv   gpp  nbp  npp runoff
mins.constr   = c(750, 300, min(dat.norm[,'gpp_gb'], na.rm = TRUE), -10, 35, 0.5)
maxes.constr  = c(3000, 800, max(dat.norm[,'gpp_gb'], na.rm = TRUE),
                  max(dat.norm[,'nbp'], na.rm = TRUE),80,
                  max(dat.norm[,'runoff'], na.rm = TRUE))


#level1.ix = which(dat.norm[,'cs_gb'] > 750 & dat.norm[,'cs_gb'] < 3000 &
#  dat.norm[,'cv'] > 300 & dat.norm[,'cv'] < 800 & 
#  dat.norm[,'npp_n_gb'] > 35 &
#  dat.norm[,'npp_n_gb'] < 80 &
#  dat.norm[,'runoff'] >0.5 &
#  dat.norm[,'nbp'] > -10
#  )


# Y

normalize.na = function(X, wrt = NULL){ 
  
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

# Express the standard in terms of lhs, then normalised matrices
X.stan = c(rep(1, d-1),0)
lhs.range = round( apply(lhs,2,range),1)

# Standard parameters when compared to the initial latin hypercube
X.stan.wrt.lhs = normalize(matrix(X.stan, nrow = 1), wrt = lhs.range)

# Normalize BEFORE putting it in to the SA
# Keep everything in relation to original design


# Need to normalize the constraints too
# Normalize everything compared to the initial data (dat.norm)
mins.norm = normalize.na(matrix(mins.constr, nrow = 1), wrt = dat.norm)
maxes.norm = normalize.na(matrix(maxes.constr, nrow = 1), wrt = dat.norm)


Y.norm = normalize.na(dat.level0, wrt = dat.norm)

glob.const.oaat = constrained.oaat(X = X.level0,
  Y = Y.norm,
  n.oaat = 21,
  mins = mins.norm,
  maxes = maxes.norm,
  hold = X.stan.wrt.lhs
  )


# It's not constraining at the moment: why not?
# Its the corners that are ruled out!!


n = 21
# Here is the full sensitivity analysis

linecols.ext = c('black', paired)
ylim = c(0,1)
#pdf(file = 'graphics/ppe_ii/global_constrained_oaat.pdf', width = 9, height = 9)
dev.new(width = 9, height = 9)
par(mfrow = c(4,8), mar = c(2,3,2,0.3), oma = c(0.5,0.5, 3, 0.5))
ndat = ncol(glob.const.oaat$pred.mean)
for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  
  y.oaat = glob.const.oaat$pred.mean

  plot(glob.const.oaat$X.oaat[ix,i], y.oaat[ix],
       type = 'n',
       ylab= '', ylim = ylim, axes = FALSE,
       main = '',
       xlab = '')
  
  for(j in 1:ndat){
    y.oaat = glob.const.oaat$pred.mean[ix,j]
    lines(glob.const.oaat$X.oaat[ix,i],y.oaat, col = linecols.ext[j], lwd = 2) 
  }
  axis(1, col = 'grey', col.axis = 'grey', las = 1)
  axis(2, col = 'grey', col.axis = 'grey', las = 1)
  mtext(3, text = colnames(lhs)[i], line = 0.2, cex = 0.7)  
}
reset()
legend('top',
       legend = colnames(dat.norm), 
       col = linecols.ext,
       cex = 0.8,
       lwd = 2,
       horiz = TRUE)
#dev.off()


linecols.ext = c('black', paired)
ylim = c(0,1)
#pdf(file = 'graphics/ppe_ii/global_constrained_oaat_2.pdf', width = 9, height = 9)
dev.new(width = 9, height = 9)
par(mfrow = c(4,8), mar = c(2,3,2,0.3), oma = c(0.5,0.5, 3, 0.5))
ndat = ncol(glob.const.oaat$pred.constr)
for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  
  y.oaat = glob.const.oaat$pred.constr
  
  plot(c(0,1), c(0,1),
       type = 'n',
       ylab= '', ylim = ylim, axes = FALSE,
       main = '',
       xlab = '')
  
  for(j in 1:ndat){
    y.oaat = glob.const.oaat$pred.constr[ix,j]
    lines(glob.const.oaat$X.oaat.constr[ix,i],y.oaat, col = linecols.ext[j], lwd = 2) 
  }
  axis(1, col = 'grey', col.axis = 'grey', las = 1)
  axis(2, col = 'grey', col.axis = 'grey', las = 1)
  mtext(3, text = colnames(lhs)[i], line = 0.2, cex = 0.7)  
}
reset()
legend('top',
       legend = colnames(dat.level0), 
       col = linecols.ext,
       cex = 0.8,
       lwd = 2,
       horiz = TRUE)






