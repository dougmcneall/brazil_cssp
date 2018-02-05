# timeseries_u-ao732.R
# Analysis of timeseries of u-ao732, including
# carbon stores, npp, gpp and runoff, both globally and in the
# Amazon.

source('../per_pft.R')

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
  
  twoStep.em = twoStep(X=X, y=y, nugget=nugget, nuggetEstim=nuggetEstim, noiseVar=noiseVar,
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

years = 1861:2014

# Load up the data
lhs = read.table('data/lhs_u-ao732.txt', header = TRUE)
stanparam = read.table('data/stanparms_u-ao732.txt', header = TRUE)
stanparam = matrix(c(rep(1, 31),0.01), nrow = 1)
stanparam.norm = normalize.wrt(stanparam, lhs)

X = normalize(lhs)
colnames(X) = colnames(lhs)
d = ncol(X)

fnvec = dir('data', pattern = 'Annual.Amazon')
fnlocvec = paste0('data/', fnvec)

# Extract the bit of the filename that describes the data
fs = lapply(fnvec,regexpr, pattern = 'Annual.Amazon')
fe = lapply(fnvec,regexpr, pattern = 'global_sum.txt')

fnams = rep(NA, length(fnvec))
for(i in 1:length(fnvec)){
  fnams[i] = substr(fnvec[i], attr(fs[[1]], 'match.length')+2, fe[[i]][1]-2)
}

greys = brewer.pal(9, "Greys")
blues = brewer.pal(9, "Blues")

# -----------------------------------------------------------------------------------
# Analysis starts here
#
# -----------------------------------------------------------------------------------

# Plot the timeseries
pdf(width = 7, height = 12, file = 'graphics/tsplot.pdf')
par(mfrow = c(6, 2), mar = c(2,2,2,2))
for(i in 1:length(fnlocvec)){
  
  dat = load_ts_ensemble(fnlocvec[i])
  dat.anom = anomalizeTSmatrix(dat, 1:30)
  
  matplot(years, t(dat), type = 'l', col = linecols, lty = 'solid', 
          lwd = 0.5, axes = FALSE,
          main = fnams[i])
  axis(1, col = 'grey')
  axis(2, col = 'grey')
  
  matplot(years, t(dat.anom), type = 'l', col = linecols,
          lty = 'solid',lwd = 0.5, axes = FALSE)
  axis(1, col = 'grey')
  axis(2, col = 'grey')
    
}
dev.off()

# Sensitivity analysis to starting value
startvalue.sensmat = matrix(NA, ncol=d, nrow=length(fnlocvec))
colnames(startvalue.sensmat) = colnames(lhs)

# Run over all outputs
for(i in 1:length(fnlocvec)){
  
  dat = load_ts_ensemble(fnlocvec[i])
  dat.start = dat[,1]
  
  ts.sens = twoStep.sens(X=X, y = dat.start)
  sens.norm = ts.sens/max(ts.sens)
  startvalue.sensmat[i, ] = sens.norm
}



pdf(file = 'graphics/sensitivity_matrix_fullspace.pdf', width = 9, height = 9)
par(mfrow = c(2,1), mar = c(8,7,3,2))
image(t(startvalue.sensmat), col = blues, axes = FALSE)
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/(length(fnams)-1) ),
     labels = fnams, las = 1)

abssum.startvalue.sensmat = apply(abs(startvalue.sensmat), 2, sum)

par(mar = c(8,7,0,2))
plot(1:d, abssum.startvalue.sensmat, axes = FALSE, xlab = '', ylab = 'Sum abs. sensitivity', pch = 19)
segments(x0 = 1:d, y0 = rep(0,d), x1 = 1:d, y1 = abssum.startvalue.sensmat)
axis(1,labels = colnames(lhs), at = 1:d, las=3, cex.axis = 0.8 )
axis(2, las = 1) 
dev.off()




change.sensmat = matrix(NA, ncol=d, nrow=length(fnlocvec))
colnames(change.sensmat) = colnames(lhs)

for(i in 1:length(fnlocvec)){
  
  dat = load_ts_ensemble(fnlocvec[i])
  dat.change = ts.ensemble.change(dat, 1:30, 125:154)
  
  ts.sens = twoStep.sens(X=X, y = dat.change)
  sens.norm = ts.sens/max(ts.sens)
  change.sensmat[i, ] = sens.norm
}

abssum.change.sensmat = apply(abs(change.sensmat), 2, sum)

pdf(file = 'graphics/sensitivity_tschange_fullspace.pdf', width = 9, height = 9)
par(mfrow = c(2,1), mar = c(8,7,3,2))
image(t(change.sensmat), col = blues, axes = FALSE)
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/(length(fnams)-1) ),
     labels = fnams, las = 1)

par(mar = c(8,7,0,2))
plot(1:d, abssum.change.sensmat, axes = FALSE, xlab = '', ylab = 'Sum abs. sensitivity', pch = 19)
segments(x0 = 1:d, y0 = rep(0,d), x1 = 1:d, y1 = abssum.change.sensmat)
axis(1,labels = colnames(lhs), at = 1:d, las=3, cex.axis = 0.8 )
axis(2, las = 1) 

dev.off()


# remove parts of the parameter space where the model produces bad output,
# and then re-do the sensitivity plots.

# First, remove everything with terrible runoff.
#



runoff = load_ts_ensemble(fnlocvec[6])/1e8
runoff.ix = which(runoff[,1] > 0.8)

X.runoff = X[runoff.ix, ]

runoffconst.sensmat = matrix(NA, ncol=d, nrow=length(fnlocvec))
colnames(runoffconst.sensmat) = colnames(lhs)

# Run over all outputs
for(i in 1:length(fnlocvec)){
  
  dat = load_ts_ensemble(fnlocvec[i])
  dat.start = dat[,1]
  dat.const = dat.start[runoff.ix]
  
  ts.sens = twoStep.sens(X=X.runoff, y = dat.const)
  sens.norm = ts.sens/max(ts.sens)
  runoffconst.sensmat[i, ] = sens.norm
}

abssum.runoffconst.sensmat = apply(abs(runoffconst.sensmat), 2, sum)

pdf(file = 'graphics/sensitivity_matrix_runoff_constrained.pdf', width = 9, height = 9)
par(mfrow = c(2,1), mar = c(8,7,3,2))
image(t(runoffconst.sensmat), col = blues, axes = FALSE)
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/(length(fnams)-1) ),
     labels = fnams, las = 1)

par(mar = c(8,7,0,2))
plot(1:d, abssum.runoffconst.sensmat, axes = FALSE, xlab = '', ylab = 'Sum abs. sensitivity', pch = 19)
segments(x0 = 1:d, y0 = rep(0,d), x1 = 1:d, y1 = abssum.runoffconst.sensmat)
axis(1,labels = colnames(lhs), at = 1:d, las=3, cex.axis = 0.8 )
axis(2, las = 1) 
dev.off()


runoffconst.sort = sort(abssum.runoffconst.sensmat, decreasing = TRUE, index.return = TRUE)

pdf(file = 'graphics/ordered_oaat_SA_runoff_constrained.pdf', width = 7, height = 5)
par(mar = c(8,4,3,1))
plot(1:d, runoffconst.sort$x, axes = FALSE, pch = 19,
     xlab = '', ylab = 'OAAT Sensitivity Index')
segments(x0 = 1:d, y0 = rep(0,d), x1 = 1:d, y1 = runoffconst.sort$x)
axis(1, at = 1:d,  labels = names(runoffconst.sort$x), las = 3, cex.axis = 0.8)
axis(2,las =1)
dev.off()


# one-at-time sensitivity plots, containing all outputs
n = 21
X.oaat.runoff = oaat.design(X.runoff, n = n)
colnames(X.oaat.runoff) = colnames(lhs)

oaat.mat = matrix(NA, ncol = length(fnlocvec), nrow = nrow(X.oaat.runoff)) 

for(i in 1:length(fnlocvec)){
  
  dat = load_ts_ensemble(fnlocvec[i])
  dat.start = dat[,1]
  dat.const = dat.start[runoff.ix]
  
  em = twoStep(X = X.runoff, y = dat.const)
  pred = predict(em$emulator, newdata = X.oaat.runoff, type = 'UK')

  oaat.mat[, i] = pred$mean

}

oaat.norm = normalize(oaat.mat)

linecols.ext = c('black', paired)

ylim = c(0,1)
pdf(file = 'graphics/runoff_constrained_amazon_oaat.pdf', width = 9, height = 9)
par(mfrow = c(4,8), mar = c(2,3,2,0.3), oma = c(0.5,0.5, 3, 0.5))

for(i in 1:d){
  
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  y.oaat = oaat.norm[,1]

  plot(X.oaat.runoff[ix,i], y.oaat[ix],
       type = 'n',
       ylab= '', ylim = ylim, axes = FALSE,
       main = '',
       xlab = '')
  
  for(j in 1:length(fnlocvec)){
    
    y.oaat = oaat.norm[ix,j]
    lines(X.oaat.runoff[ix,i],y.oaat, col = linecols.ext[j], lwd = 2) 
  }
  
  axis(1, col = 'grey', col.axis = 'grey', las = 1)
  axis(2, col = 'grey', col.axis = 'grey', las = 1)
  mtext(3, text = colnames(lhs)[i], line = 0.2, cex = 0.7)  

}

reset()
legend('top',
       legend = fnams, 
       col = linecols.ext,
       lwd = 2,
       horiz = TRUE)
  
dev.off()

# --------------------------------------------------------------------------------
# Now apply constraints to the global data and apply it to the local.
#
#
#
# npp 35-80 GtC
# nbp > 0
# cVeg 300 - 800 GtC
# cSoil 750 - 3000 GtC
# --------------------------------------------------------------------------------
#dir('data', pattern = 'Annual/.(?!Amazon).*')
fnallvec = dir('data', pattern = 'Annual')
# WARNING - hard coded hack to sort
fidx = grep("Annual.(?!Amazon).*", fnallvec, perl=TRUE)
fnvec = fnallvec[fidx]
fnlocvec = paste0('data/', fnvec)

# Extract the bit of the filename that describes the data
fs = lapply(fnvec,regexpr, pattern = 'Annual')
fe = lapply(fnvec,regexpr, pattern = 'global_sum.txt')

fnams = rep(NA, length(fnvec))
for(i in 1:length(fnvec)){
  fnams[i] = substr(fnvec[i], attr(fs[[1]], 'match.length')+2, fe[[i]][1]-2)
}



datmat = matrix(nrow = nrow(lhs), ncol = length(fnlocvec))

for(i in 1:length(fnlocvec)){
  dat = load_ts_ensemble(fnlocvec[i])
  dat.modern = dat[,135:154]
  mean.modern = apply(dat.modern, 1, mean)
  datmat[ , i] = mean.modern
}
colnames(datmat) = fnams

norm.vec = c(1e12, 1e12, 1e6, 1e12, 1e12, 1e9)

dat.norm = sweep(datmat, 2, norm.vec, FUN = '/')

dev.new(width = 9, height = 10)
par(mfrow = c(3,2))
for(i in 1:6){
  hist(dat.norm[,i], main = fnams[i])
}

# Constrain with runoff - removes 62 members
globrunoff.ix = dat.norm[,'runoff'] >0.5
dat.globrunoff  = dat.norm[globrunoff.ix, ]
X.globrunoff = X[globrunoff.ix, ]

dev.new(width = 8, height = 8)
pairs(X.globrunoff[,1:10], xlim = c(0,1), ylim =c(0,1), pch = 19, gap = 0)

dev.new(width = 9, height = 10)
par(mfrow = c(3,2))
for(i in 1:6){
  hist(dat.runoffconst[,i], main = fnams[i])
}

#constrain with nbp - removes 209 members (not failed runoff)
dat.nbpconst  = dat.norm[dat.norm[,'nbp'] >= 0, ]
dev.new(width = 9, height = 10)
par(mfrow = c(3,2))
for(i in 1:6){
  hist(dat.nbpconst[,i], main = fnams[i])
}

# constrain with npp - removes 172 members (including failed runoff)
dat.nppconst  = dat.norm[dat.norm[,'npp_n_gb'] > 35 & dat.norm[,'npp_n_gb'] <80, ]
dev.new(width = 9, height = 10)
par(mfrow = c(3,2))
for(i in 1:6){
  hist(dat.nppconst[,i], main = fnams[i])
}


# constrain with cv - removes 255 members (including failed runoff)
dat.cvconst  = dat.norm[dat.norm[,'cv'] > 300 & dat.norm[,'cv'] <800, ]
dev.new(width = 9, height = 10)
par(mfrow = c(3,2))
for(i in 1:6){
  hist(dat.cvconst[,i], main = fnams[i])
}


# constrain with cs - removes 149 members (not failed runoff)
dat.csconst  = dat.norm[dat.norm[,'cs_gb'] > 750 & dat.norm[,'cs_gb'] <3000, ]
dev.new(width = 9, height = 10)
par(mfrow = c(3,2))
for(i in 1:6){
  hist(dat.csconst[,i], main = fnams[i])
}

# sensitivity using npp as a constraint
npp.ix = which(dat.norm[,'npp_n_gb'] > 35 & dat.norm[,'npp_n_gb'] <80)
X.npp = X[npp.ix, ]

nppconst.sensmat = matrix(NA, ncol=d, nrow=length(fnlocvec))
colnames(nppconst.sensmat) = colnames(lhs)

# Run over all outputs
for(i in 1:ncol(dat.norm)){
  
  dat.const = dat.norm[npp.ix, i]
  
  ts.sens = twoStep.sens(X=X.npp, y = dat.const)
  sens.norm = ts.sens/max(ts.sens)
  nppconst.sensmat[i, ] = sens.norm
}

abssum.nppconst.sensmat = apply(abs(nppconst.sensmat), 2, sum)

pdf(file = 'graphics/sensitivity_matrix_npp_constrained.pdf', width = 9, height = 9)
par(mfrow = c(2,1), mar = c(8,7,3,2))
image(t(nppconst.sensmat), col = blues, axes = FALSE)
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/(length(fnams)-1) ),
     labels = fnams, las = 1)

par(mar = c(8,7,0,2))
plot(1:d, abssum.nppconst.sensmat, axes = FALSE, xlab = '', ylab = 'Sum abs. sensitivity', pch = 19)
segments(x0 = 1:d, y0 = rep(0,d), x1 = 1:d, y1 = abssum.nppconst.sensmat)
axis(1,labels = colnames(lhs), at = 1:d, las=3, cex.axis = 0.8 )
axis(2, las = 1) 
dev.off()



n = 21
X.oaat.npp = oaat.design(X.npp, n = n)
colnames(X.oaat.npp) = colnames(lhs)

oaat.mat.npp = matrix(NA, ncol = length(fnlocvec), nrow = nrow(X.oaat.npp)) 

for(i in 1:length(fnlocvec)){
  
  dat.const = dat.norm[npp.ix,i]
  em = twoStep(X = X.npp, y = dat.const)
  pred = predict(em$emulator, newdata = X.oaat.npp, type = 'UK')

  oaat.mat.npp[, i] = pred$mean

}

oaat.npp.norm = normalize(oaat.mat.npp)
linecols.ext = c('black', paired)

ylim = c(0,1)
pdf(file = 'graphics/npp_constrained_global_oaat.pdf', width = 9, height = 9)
par(mfrow = c(4,8), mar = c(2,3,2,0.3), oma = c(0.5,0.5, 3, 0.5))

for(i in 1:d){
  
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  y.oaat = oaat.npp.norm[,1]

  plot(X.oaat.npp[ix,i], y.oaat[ix],
       type = 'n',
       ylab= '', ylim = ylim, axes = FALSE,
       main = '',
       xlab = '')
  
  for(j in 1:length(fnlocvec)){
    
    y.oaat = oaat.npp.norm[ix,j]
    lines(X.oaat.npp[ix,i],y.oaat, col = linecols.ext[j], lwd = 2) 
  }
  
  axis(1, col = 'grey', col.axis = 'grey', las = 1)
  axis(2, col = 'grey', col.axis = 'grey', las = 1)
  mtext(3, text = colnames(lhs)[i], line = 0.2, cex = 0.7)  

}

reset()
legend('top',
       legend = fnams, 
       col = linecols.ext,
       lwd = 2,
       horiz = TRUE)
  
dev.off()


# Build an emulator with the runoff-constrained ensemble


X.oaat.globrunoff = oaat.design(X.globrunoff, med = FALSE, hold = stanparam.norm, n = 21)

#produce an oaat output matrix
y.oaat.globrunoff  = matrix(nrow = nrow(X.oaat.globrunoff),
                            ncol = ncol(dat.globrunoff ))

colnames(y.oaat.globrunoff) = colnames(dat.globrunoff)

for(i in 1:ncol(y.oaat.globrunoff)){
  
  em = twoStep(X = X.globrunoff, y = dat.globrunoff[,i])
  y.oaat = predict(em$emulator, newdata = X.oaat.globrunoff, type = 'UK')
  y.oaat.globrunoff[,i] = y.oaat$mean
  
}


# The idea is to compare the sensitivity when constrained by each global thing.

# Visualise the constrained space ...
mins  = apply(X.globrunoff, 2, min)
maxes = apply(X.globrunoff, 2, max)

nsamp.unif = 100000
X.unif = samp.unif(nsamp.unif, mins = mins, maxes = maxes)

y.unif = matrix(nrow = nsamp.unif, ncol = ncol(dat.globrunoff))
colnames(y.unif) = colnames(dat.globrunoff)

for(i in 1:ncol(y.unif)){
  em = twoStep(X = X.globrunoff, y = dat.globrunoff[,i])
  pred = predict(em$emulator, newdata = X.unif, type = 'UK')
  y.unif[,i] = pred$mean
  
}

dev.new(width = 7, height = 7)
par(mfrow = c(2,3))

for(i in 1:ncol(y.unif)){
  hist(y.unif[,i], main = fnams[i])
  
}

# Now apply constraints to the output matrix
# npp 35-80 GtC
# nbp > 0
# cVeg 300 - 800 GtC
# cSoil 750 - 3000 GtC

ix.kept = which(y.unif[,'cs_gb'] > 750 & y.unif[,'cs_gb'] < 3000 & y.unif[,'cv'] > 300 & y.unif[,'cv'] < 800 & y.unif[,'npp_n_gb'] > 35 & y.unif[,'npp_n_gb'] < 80)
X.kept = X.unif[ix.kept, ]

rb = brewer.pal(9, "RdBu")
br = rev(rb)

pdf(file = 'graphics/pairs_dens_all_constraints.pdf', width = 10, height = 10)
par(oma = c(0,0,0,3))

# Emulate all input space and keep only those inputs which match a 
# criteria (such as having absolute error below a threshold)

test = pairs(X.kept,
             labels = 1:d,
             gap = 0, lower.panel = NULL, xlim = c(0,1), ylim = c(0,1),
             panel = dfunc.up,
             cex.labels = 1,
             col.axis = 'white')

image.plot(legend.only = TRUE,
           zlim = c(0,1),
           col = rb,
           legend.args = list(text = 'Density of model runs matching the criteria', side = 3, line = 1),
           horizontal = TRUE
)
#par(xpd = NA)
#text(0.2, 0.6, labels = paste0(colnames(lhs)[1:5],'\n'))
legend('left', legend = paste(1:d, colnames(lhs)), cex = 0.9, bty = 'n')

dev.off()

pdf(file = 'graphics/hist_kept_all_constraints.pdf', width = 8, height = 6)
#dev.new(width = 8, height = 6)
par(mfrow =c(4,8), fg = 'white', mar = c(3,1,3,1))
for(i in 1:d){
  
  hist(X.kept[, i], xlim = c(0,1), axes = FALSE, col = 'black',
       xlab = '', ylab = '', main = '')
  abline(v = stanparam.norm[i], col = 'red')
  #axis(1, col = 'black')
  mtext(side = 3, text = colnames(X.kept)[i], line = 1, col = 'black', cex = 0.8)
  
}
dev.off()






