# timeseries_u-ao732.R
# Analysis of timeseries of u-ao732, including
# carbon stores, npp, gpp and runoff, both globally and in the
# Amazon.

source('../per_pft.R')


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

# Golbal numbers - what would these look like 
# npp 35-80 GtC
# nbp > 0
# cVeg 300 - 800 GtC
# cSoil 750 - 3000 GtC

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








