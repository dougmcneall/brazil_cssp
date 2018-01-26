# timeseries_u-ao732.R
# Analysis of timeseries of u-ao732, including
# carbon stores, npp, gpp and runoff, both globally and in the
# Amazon.

source('../per_pft.R')

years = 1861:2014

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

lhs = read.table('data/lhs_u-ao732.txt', header = TRUE)
X = normalize(lhs)
colnames(X) = colnames(lhs)
d = ncol(X)

fnvec = dir('data', pattern = 'Annual.Amls()azon')
fnlocvec = paste0('data/', fnvec)

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


# Extract the bit of the filename that describes the data
fs = lapply(fnvec,regexpr, pattern = 'Annual.Amazon')
fe = lapply(fnvec,regexpr, pattern = 'global_sum.txt')

fnams = rep(NA, length(fnvec))
for(i in 1:length(fnvec)){
  fnams[i] = substr(fnvec[i], attr(fs[[1]], 'match.length')+2, fe[[i]][1]-2)
}

greys = brewer.pal(9, "Greys")
blues = brewer.pal(9, "Blues")

pdf(file = 'graphics/sensitivity_matrix_fullspace.pdf', width = 9, height = 5)
par(mar = c(8,7,4,7))
image(t(startvalue.sensmat), col = blues, axes = FALSE)
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/(length(fnams)-1) ),
     labels = fnams, las = 1)

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

pdf(file = 'graphics/sensitivity_tschange_fullspace.pdf', width = 9, height = 5)
par(mar = c(8,7,4,7))
image(t(change.sensmat), col = blues, axes = FALSE)
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/(length(fnams)-1) ),
     labels = fnams, las = 1)

dev.off()


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









