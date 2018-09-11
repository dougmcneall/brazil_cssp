# es_ppe_ii.R
# Analyse the full 500 (ish member ensemble)
#

# 1. Check how much difference an increasing
# CO2 concentration makes, particularly for runoff.

# 2. How well does the emulator work, for various
# outputs? How much does the extra 100 members help
# with accuracy?

# 3. Constrain to all outputs (initially those fram Andy before
# we do any History Matching), and *then* run the sensitivity
# Analysis. Remove space where the model fails before building
# any emulators.

# 4. Is there a relationship between parameters and *trend*
# for the "broadly plausible" set of models? If so, we might be
# able to constrain the future values even more than just through the
# starting values.

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

years = 1861:2014
ysec = 60*60*24*365
norm.vec = c(1e12, 1e12, 1e12/ysec , 1e12, 1e12, 1e9)

# Load up the data
lhs_i = read.table('data/ES_PPE_i/lhs_u-ao732.txt', header = TRUE)
lhs_ii = read.table('data/ES_PPE_ii/lhs_u-ao732a.txt', header = TRUE)

# Exclude the last 100 members of the ensemble for
# Validation of the analysis later on.
toplevel.ix = 1:400

lhs = rbind(lhs_i, lhs_ii)[toplevel.ix, ]

#stanparam = read.table('data/stanparms_u-ao732.txt', header = TRUE)
#stanparam = matrix(c(rep(1, 31),0.01), nrow = 1)
#stanparam.norm = normalize.wrt(stanparam, lhs)

X = normalize(lhs)
colnames(X) = colnames(lhs)
d = ncol(X)

# load Amazon broadleaf forest fraction
# Ensemble members P0000 to P0498, with the standard run in the 
# final row.
Amazon.area = 7.013467e+12 # m^2
frac_bl = read.table("data/ES_PPE_ii/Amazon_forest_total.txt")[toplevel.ix, -1] / Amazon.area
frac_bl_anom = anomalizeTSmatrix(frac_bl, ix = 1:10)
frac_bl_change = ts.ensemble.change(frac_bl, startix = 1:30, endix = 125:154)

# Histogram of modern forest fraction
par(mfrow = c(2,1))
hist(apply(frac_bl[, 125:154], 1, FUN = 'mean'))
hist(frac_bl_change)

# load Amazon precipitation
nc.precip <- nc_open("data/ES_PPE_ii/JULES-ES.0p92.vn5.0.CRUNCEPv7.P0199.Annual.Amazon.precip.global_sum.nc")
precip <- ncvar_get(nc.precip)

# Plot forest fraction
pdf('graphics/ppe_ii/forest_fraction.pdf', width = 5, height = 7)
par(las = 1)
ensTShist(years, frac_bl, 
          colvec = linecols[c(1,3,4)], 
          histcol = linecols[3],
          ylab = 'fraction', xlab = '',
          grid = FALSE, 
          mainvec = 'Amazon BL Forest Fraction')
dev.off()


# Plot forest fraction anomaly
pdf('graphics/ppe_ii/forest_fraction_anomaly.pdf', width = 5, height = 7)
par(las = 1)
ensTShist(years, frac_bl_anom, 
          colvec = linecols[c(1,3,4)], 
          histcol = linecols[3],
          ylab = 'fraction', xlab = '',
          grid = FALSE, 
          mainvec = 'Amazon BL Forest Fraction Anomaly')
dev.off()

# Plot normalised runoff
runoff.raw = (load_ts_ensemble("data/ES_PPE_ii/Annual.Amazon.runoff.global_sum.txt"))[toplevel.ix, ]
runoff.norm = sweep(runoff.raw, 2, STATS = precip, FUN = '/')
runoff.norm.anom = anomalizeTSmatrix(runoff.norm, ix = 1:10)
runoff.norm.change = ts.ensemble.change(runoff.norm, startix = 1:30, endix = 125:154)


pdf('graphics/ppe_ii/runoff_normalised.pdf', width = 8, height = 5 )
par(las = 1)
ensTShist(years, runoff.norm, 
          colvec = linecols[c(1,3,4)], 
          histcol = linecols[3],
          ylab = 'fraction of precipitation', xlab = '',
          grid = FALSE, 
          mainvec = 'Amazon Normalized Runoff')

dev.off()

pdf('graphics/ppe_ii/runoff_normalised_anomaly.pdf', width = 8, height = 5 )
par(las = 1)
ensTShist(years, runoff.norm.anom, 
          colvec = linecols[c(1,3,4)], 
          histcol = linecols[3],
          ylab = 'fraction of precipitation', xlab = '',
          grid = FALSE, 
          mainvec = 'Amazon Normalized Runoff Anomaly')

dev.off()

#fnvec = dir('data/ES_PPE_ii', pattern = 'Annual.Amazon')
fnvec = c("Annual.Amazon.cs_gb.global_sum.txt",
          "Annual.Amazon.cv.global_sum.txt",
          "Annual.Amazon.gpp_gb.global_sum.txt",
          "Annual.Amazon.nbp.global_sum.txt",
          "Annual.Amazon.npp_n_gb.global_sum.txt",
          "Annual.Amazon.runoff.global_sum.txt"
          )

fnlocvec = paste0('data/ES_PPE_ii/', fnvec)

# Extract the bit of the filename that describes the data
fs = lapply(fnvec,regexpr, pattern = 'Annual.Amazon')
fe = lapply(fnvec,regexpr, pattern = 'global_sum.txt')

fnams = rep(NA, length(fnvec))
for(i in 1:length(fnvec)){
  fnams[i] = substr(fnvec[i], attr(fs[[1]], 'match.length')+2, fe[[i]][1]-2)
}

nams = c(
  "Soil Carbon",
  "Vegetation Carbon",
  "Gross Primary Production",
  "Net Biospheric Production",
  "Net Primary Production",
  "Runoff"
)

greys = brewer.pal(9, "Greys")
blues = brewer.pal(9, "Blues")

ylabs = c(
  'GtC',
  'GtC',
  'GtC/year',
  'GtC/year',
  'GtC/year',
  'Sv'
)
# Plot the timeseries
pdf(width = 7, height = 12, file = 'graphics/ppe_ii/tsplot_Amazon.pdf')
par(mfrow = c(6, 2), mar = c(5,4,2,2), las = 1)
for(i in 1:length(fnlocvec)){
  
  dat = load_ts_ensemble(fnlocvec[i])/norm.vec[i]
  dat.anom = anomalizeTSmatrix(dat, 1:30)
  
  matplot(years, t(dat), type = 'l', col = linecols[c(1,3,4)], lty = 'solid', 
          lwd = 0.5, axes = FALSE,
          main = '',
          ylab = ylabs[i],
          xlab = '')
  axis(1, at = seq(from = 1860, to = 2020, by = 30))
  axis(2)
  mtext(side = 3, adj = 0, line = 0.5, text = nams[i])
  
  
  matplot(years, t(dat.anom), type = 'l', col = linecols[c(1,3,4)],
          lty = 'solid',lwd = 0.5, axes = FALSE,
          ylab = ylabs[i],
          xlab = '')
  axis(1, at = seq(from = 1860, to = 2020, by = 30))
  axis(2)
  
}
dev.off()

# ----------------------------------------------------------------------
# Remove parts of parameter space where the model does really
# badly at simulating runoff at the start of the run.
# Use this to get some basic idea about how well the 
# emulator works.
# ----------------------------------------------------------------------
runoff = (load_ts_ensemble("data/ES_PPE_ii/Annual.Amazon.runoff.global_sum.txt")/1e8)[toplevel.ix, ]
runoff.ix = which(runoff[,1] > 0.8)

X.runoff = X[runoff.ix, ]

# There's a relationship between runoff starting value and runoff change, but
# not sure about the causality.
runoff.start = runoff[runoff.ix, 1]
runoff.change = ts.ensemble.change(runoff[runoff.ix, ], 1:10, 145:154)
plot(runoff.start, runoff.change)



# datmat contains the starting values of the timeseries
nruns = length(toplevel.ix) # number of runs
ndat = 8 # number of data streams 

# Sensitivity analysis of space left over once we've removed the
# non-performing models
runoffconst.sensmat = matrix(NA, ncol=d, nrow=ndat)
colnames(runoffconst.sensmat) = colnames(lhs)

datmat = matrix(NA, nrow = nruns, ncol =  ndat)
datnames = c(fnams,'n_runoff', 'frac_bl')
colnames(datmat) = datnames

for(i in 1:length(fnlocvec)){
  datmat[,i] = load_ts_ensemble(fnlocvec[i])[toplevel.ix, 1]
}
datmat[,7 ] = runoff.norm[, 1]
datmat[, 8] = frac_bl[,1]

# Run over all outputs. Constrain to the "runoff approved" input space
for(i in 1:ncol(datmat)){
  dat.const = datmat[runoff.ix, i]
  ts.sens = twoStep.sens(X=X.runoff, y = dat.const)
  sens.norm = ts.sens/max(ts.sens)
  runoffconst.sensmat[i, ] = sens.norm
}

abssum.runoffconst.sensmat = apply(abs(runoffconst.sensmat), 2, sum)

# Sensitivity matrix of runoff-constrained inputs
pdf(file = 'graphics/ppe_ii/sensitivity_matrix_runoff_constrained.pdf', width = 9, height = 9)
par(mfrow = c(2,1), mar = c(8,7,3,2))
image(t(runoffconst.sensmat), col = blues, axes = FALSE)
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/(ndat-1) ),
     labels = datnames, las = 1)

par(mar = c(8,7,0,2))
plot(1:d, abssum.runoffconst.sensmat, axes = FALSE, xlab = '', ylab = 'Sum abs. sensitivity', pch = 19)
segments(x0 = 1:d, y0 = rep(0,d), x1 = 1:d, y1 = abssum.runoffconst.sensmat)
axis(1,labels = colnames(lhs), at = 1:d, las=3, cex.axis = 0.8 )
axis(2, las = 1) 
dev.off()

runoffconst.sort = sort(abssum.runoffconst.sensmat, decreasing = TRUE, index.return = TRUE)

# Sorted summary sensitivity of runoff to inputs
pdf(file = 'graphics/ppe_ii/ordered_oaat_SA_runoff_constrained.pdf', width = 7, height = 5)
par(mar = c(8,4,3,1))
plot(1:d, runoffconst.sort$x, axes = FALSE, pch = 19,
     xlab = '', ylab = 'OAAT Sensitivity Index')
segments(x0 = 1:d, y0 = rep(0,d), x1 = 1:d, y1 = runoffconst.sort$x)
axis(1, at = 1:d,  labels = names(runoffconst.sort$x), las = 3, cex.axis = 0.8)
axis(2,las =1)
dev.off()


# One-at-a-time sensitivity plots of all outputs
n = 21
X.oaat.runoff = oaat.design(X.runoff, n = n)
colnames(X.oaat.runoff) = colnames(lhs)

oaat.mat = matrix(NA, ncol = ndat, nrow = nrow(X.oaat.runoff)) 

for(i in 1:ndat){
  dat.const = datmat[runoff.ix,i]
  em = twoStep.glmnet(X = X.runoff, y = dat.const)
  pred = predict(em$emulator, newdata = X.oaat.runoff, type = 'UK')
  oaat.mat[, i] = pred$mean
}

oaat.norm = normalize(oaat.mat)

linecols.ext = c('black', paired)

ylim = c(0,1)
pdf(file = 'graphics/ppe_ii/runoff_constrained_amazon_oaat.pdf', width = 9, height = 9)
par(mfrow = c(4,8), mar = c(2,3,2,0.3), oma = c(0.5,0.5, 3, 0.5))

for(i in 1:d){
  
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  y.oaat = oaat.norm[,1]
  
  plot(X.oaat.runoff[ix,i], y.oaat[ix],
       type = 'n',
       ylab= '', ylim = ylim, axes = FALSE,
       main = '',
       xlab = '')
  
  for(j in 1:ndat){
    
    y.oaat = oaat.norm[ix,j]
    lines(X.oaat.runoff[ix,i],y.oaat, col = linecols.ext[j], lwd = 2) 
  }
  
  axis(1, col = 'grey', col.axis = 'grey', las = 1)
  axis(2, col = 'grey', col.axis = 'grey', las = 1)
  mtext(3, text = colnames(lhs)[i], line = 0.2, cex = 0.7)  
  
}

reset()
legend('top',
       legend = datnames, 
       col = linecols.ext,
       lwd = 2,
       horiz = TRUE,
       cex = 0.8)

dev.off()


# Now do the same for the *change* in the variable 
# over the time period
datchangemat = matrix(NA, nrow = nruns, ncol =  ndat)
colnames(datchangemat) = datnames

for(i in 1:length(fnlocvec)){
  dat = load_ts_ensemble(fnlocvec[i])[toplevel.ix, ]
  dat.change = ts.ensemble.change(dat, 1:10, 145:154)
  datchangemat[,i] = dat.change
}

datchangemat[,7 ] = runoff.norm.change
datchangemat[, 8] = frac_bl_change

for(i in 1:ndat){
  dat.change = datchangemat[runoff.ix, i]
  em = twoStep.glmnet(X = X.runoff, y = dat.change)
  pred = predict(em$emulator, newdata = X.oaat.runoff, type = 'UK')
  oaat.mat[, i] = pred$mean
}

oaat.norm = normalize(oaat.mat)

linecols.ext = c('black', paired)
ylim = c(0,1)
pdf(file = 'graphics/ppe_ii/runoff_constrained_amazon_change_oaat.pdf', width = 9, height = 9)
par(mfrow = c(4,8), mar = c(2,3,2,0.3), oma = c(0.5,0.5, 3, 0.5))

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  y.oaat = oaat.norm[,1]
  
  plot(X.oaat.runoff[ix,i], y.oaat[ix],
       type = 'n',
       ylab= '', ylim = ylim, axes = FALSE,
       main = '',
       xlab = '')
  
  for(j in 1:ndat){
    y.oaat = oaat.norm[ix,j]
    lines(X.oaat.runoff[ix,i],y.oaat, col = linecols.ext[j], lwd = 2) 
  }
  axis(1, col = 'grey', col.axis = 'grey', las = 1)
  axis(2, col = 'grey', col.axis = 'grey', las = 1)
  mtext(3, text = colnames(lhs)[i], line = 0.2, cex = 0.7)  
}
reset()
legend('top',
       legend = datnames, 
       col = linecols.ext,
       cex = 0.8,
       lwd = 2,
       horiz = TRUE)
dev.off()

# Need sensitivity and ordered sensitivity of runoff change
# Sensitivity analysis of space left over once we've removed the
# non-performing models
runoffchange.sensmat = matrix(NA, ncol=d, nrow=ndat)
colnames(runoffchange.sensmat) = colnames(lhs)

# Run over all outputs
for(i in 1:ndat){
  
  dat.change = datchangemat[runoff.ix, i]
  
  ts.sens = twoStep.sens(X=X.runoff, y = dat.change)
  sens.norm = ts.sens/max(ts.sens)
  runoffchange.sensmat[i, ] = sens.norm
}

abssum.runoffchange.sensmat = apply(abs(runoffchange.sensmat), 2, sum, na.rm = TRUE)

# Sensitivity matrix of runoff-constrained inputs
pdf(file = 'graphics/ppe_ii/sensitivity_matrix_runoff_change_constrained.pdf', width = 9, height = 9)
par(mfrow = c(2,1), mar = c(8,7,3,2))
image(t(runoffchange.sensmat), col = blues, axes = FALSE)
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/(ndat-1) ),
     labels = datnames, las = 1)

par(mar = c(8,7,0,2))
plot(1:d, abssum.runoffchange.sensmat, axes = FALSE, xlab = '', ylab = 'Sum abs. sensitivity', pch = 19)
segments(x0 = 1:d, y0 = rep(0,d), x1 = 1:d, y1 = abssum.runoffchange.sensmat)
axis(1,labels = colnames(lhs), at = 1:d, las=3, cex.axis = 0.8 )
axis(2, las = 1) 
dev.off()

runoffchange.sort = sort(abssum.runoffchange.sensmat, decreasing = TRUE, index.return = TRUE)

# Sorted summary sensitivity of runoff to inputs
pdf(file = 'graphics/ppe_ii/ordered_oaat_SA_runoff_change_constrained.pdf', width = 7, height = 5)
par(mar = c(8,4,3,1))
plot(1:d, runoffchange.sort$x, axes = FALSE, pch = 19,
     xlab = '', ylab = 'OAAT Sensitivity Index')
segments(x0 = 1:d, y0 = rep(0,d), x1 = 1:d, y1 = runoffchange.sort$x)
axis(1, at = 1:d,  labels = names(runoffchange.sort$x), las = 3, cex.axis = 0.8)
axis(2,las =1)
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
#fnvec = c("Annual.cs_gb.global_sum.txt",
#          "Annual.cv.global_sum.txt",
#          "Annual.gpp_gb.global_sum.txt",
#          "Annual.nbp.global_sum.txt",
#          "Annual.npp_n_gb.global_sum.txt",
#          "Annual.runoff.global_sum.txt")

fnallvec = dir('data/ES_PPE_ii/', pattern = 'Annual')
# WARNING - hard coded hack to sort
fidx = grep("Annual.(?!Amazon).*", fnallvec, perl=TRUE)
fnvec_interim = fnallvec[fidx]
fidx2 = grep("sum.(?!standard).*", fnvec_interim, perl=TRUE)
fnvec = fnvec_interim[fidx2]
fnlocvec = paste0('data/ES_PPE_ii/', fnvec)

# Extract the bit of the filename that describes the data
fs = lapply(fnvec,regexpr, pattern = 'Annual')
fe = lapply(fnvec,regexpr, pattern = 'global_sum.txt')

fnams = rep(NA, length(fnvec))
for(i in 1:length(fnvec)){
  fnams[i] = substr(fnvec[i], attr(fs[[1]], 'match.length')+2, fe[[i]][1]-2)
}

# Constrain on runoff first
datmat = matrix(nrow = nrow(X), ncol = length(fnlocvec))
for(i in 1:length(fnlocvec)){
  dat = load_ts_ensemble(fnlocvec[i])[toplevel.ix, ]
  dat.modern = dat[ ,135:154]
  mean.modern = apply(dat.modern, 1, mean)
  datmat[ , i] = mean.modern
}
colnames(datmat) = fnams


dat.norm = sweep(datmat, 2, norm.vec, FUN = '/')

dev.new(width = 9, height = 10)
par(mfrow = c(3,2))
for(i in 1:6){
  hist(dat.norm[,i], main = fnams[i])
}

# Constrain with global runoff
globrunoff.ix = which(dat.norm[,'runoff'] >0.5 & dat.norm[,'nbp'] > -10)
dat.globrunoff  = dat.norm[globrunoff.ix, ]
X.globrunoff = X[globrunoff.ix, ]

# Visualise the constrained space ...
mins  = apply(X.globrunoff, 2, min)
maxes = apply(X.globrunoff, 2, max)

nsamp.unif = 100000
X.unif = samp.unif(nsamp.unif, mins = mins, maxes = maxes)

y.unif = matrix(nrow = nsamp.unif, ncol = ncol(dat.globrunoff))
colnames(y.unif) = colnames(dat.globrunoff)

global.emlist = vector('list',length(fnams))

for(i in 1:ncol(y.unif)){
  em = twoStep.glmnet(X = X.globrunoff, y = dat.globrunoff[,i])
  global.emlist[[i]] = em
  pred = predict(em$emulator, newdata = X.unif, type = 'UK')
  y.unif[,i] = pred$mean
}

ix.kept = which(y.unif[,'cs_gb'] > 750 & y.unif[,'cs_gb'] < 3000 & y.unif[,'cv'] > 300 & y.unif[,'cv'] < 800 & y.unif[,'npp_n_gb'] > 35 & y.unif[,'npp_n_gb'] < 80)
X.kept = X.unif[ix.kept, ]

# we've removed 80% of our prior input space
(nrow(X.kept) / nsamp.unif) * 100

# Any marginal constraint? (not really)
apply(X.kept, 2, min)
apply(X.kept, 2, max)


bl_frac_modern = frac_bl[globrunoff.ix, 154]

em = twoStep.glmnet(X = X.globrunoff, y = bl_frac_modern)
pred = predict(em$emulator, newdata = X.unif, type = 'UK')

ix.kept2 = which(y.unif[,'cs_gb'] > 750 & y.unif[,'cs_gb'] < 3000 & 
                  y.unif[,'cv'] > 300 & y.unif[,'cv'] < 800 &
                  y.unif[,'npp_n_gb'] > 35 & y.unif[,'npp_n_gb'] < 80 &
                  pred$mean > 0.5
                  )
X.kept2 = X.unif[ix.kept2, ]
(nrow(X.kept2) / nsamp.unif) * 100

# Histograms of the constraint outputs
#note: always pass alpha on the 0-255 scale
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

hcol = 'grey'
lcol = 'black'
pdf(file = 'graphics/ppe_ii/constraint_hists.pdf', width = 8, height = 8)
par(mfrow = c(3,2), fg = 'white', las = 1)

hist(dat.norm[,'runoff'], col = hcol, main = 'Runoff', xlab = 'Sv')
polygon(x = c(0.5, 100, 100, 0.5), y = c(0, 0, 1000, 1000), col = makeTransparent('tomato2'))
hist(dat.norm[,'nbp'], col = hcol, main = 'NBP', xlab = 'GtC/year')
polygon(x = c(-10, 100, 100, -10), y = c(0, 0, 1000, 1000), col = makeTransparent('tomato2'))
hist(dat.norm[,'cs_gb'], col = hcol, main = 'Soil Carbon', xlab = 'GtC')
polygon(x = c(750, 3000, 3000, 750), y = c(0, 0, 1000, 1000), col = makeTransparent('tomato2'))

hist(dat.norm[,'cv'], col = hcol, main = 'Vegetation Carbon', xlab = 'GtC')
polygon(x = c(300, 800, 800, 300), y = c(0, 0, 1000, 1000), col = makeTransparent('tomato2'))
hist(dat.norm[,'npp_n_gb'], col = hcol , main = 'NPP', xlab = 'GtC/year')
polygon(x = c(35, 80, 80, 35), y = c(0, 0, 1000, 1000), col = makeTransparent('tomato2'))
hist(bl_frac_modern, col = hcol, main = 'Amazon Forest Fraction', xlab = 'fraction')
polygon(x = c(0.5, 1, 1, 0.5), y = c(0, 0, 1000, 1000), col = makeTransparent('tomato2'))
dev.off()



dev.new(width = 7, height = 7)
par(mfrow = c(2,3))

for(i in 1:ncol(y.unif)){
  hist(y.unif[ix.kept2,i], main = fnams[i])
}

rb = brewer.pal(9, "RdBu")
br = rev(rb)

pdf(file = 'graphics/ppe_ii/pairs_dens_all_constraints2.pdf', width = 10, height = 10)
par(oma = c(0,0,0,3))
test = pairs(X.kept2,
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

pdf("graphics/ppe_ii/constrained_inputs_hists.pdf", width = 9, height = 6)
par(mfrow = c(4,8), mar = c(3,3,2,0.3), oma = c(0.5,0.5, 3, 0.5), fg = 'white')
for(i in 1:d){
  hist(X.kept2[,i], xlim = c(0,1), col = 'darkgrey', axes = FALSE, main = colnames(lhs)[i])
  axis(1, col = 'black')
  #axis(2)
}
dev.off()






