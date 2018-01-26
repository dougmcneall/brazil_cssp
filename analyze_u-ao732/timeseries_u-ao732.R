# timeseries_u-ao732.R
# Analysis of timeseries of u-ao732, including
# carbon stores, npp, gpp and runoff, both globally and in the
# Amazon.

source('../per_pft.R')

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

lhs = read.table('data/lhs_u-ao732.txt', header = TRUE)
X = normalize(lhs)
colnames(X) = colnames(lhs)
d = ncol(X)

fnvec = dir('data', pattern = 'Annual.Amazon')
fnlocvec = paste0('data/', fnvec)

cs = load_ts_ensemble(fnlocvec[1])/1e13
cs.start = cs[,1]


test = twoStep.sens(X=X, y = cs.start)
# normalise output from the SA

startvalue.sensmat = matrix(NA, ncol=d, nrow=length(fnlocvec))
colnames(startvalue.sensmat) = colnames(lhs)

# Include changes
for(i in 1:length(fnlocvec)){
  
  dat = load_ts_ensemble(fnlocvec[i])
  dat.start = dat[,1]
  
  ts.sens = twoStep.sens(X=X, y = dat.start)
  sens.norm = ts.sens/max(ts.sens)
  startvalue.sensmat[i, ] = sens.norm
}









filelist.global <- paste0('frac_area_mean/global_area_mean_PFT',0:16, '.txt')
filelist.wus <- paste0('frac_area_mean/WUS_area_mean_PFT',0:16, '.txt')
filelist.sam <- paste0('frac_area_mean/SAM_area_mean_PFT',0:16, '.txt')

c.df <- function(fn){
  out = c(read.table(fn, header = FALSE, skip = 1), recursive = TRUE)
  out
}

# creates a list of (PFTs) long, each element of which has (ensemble members) data points
global_area_means.list <- lapply(filelist.global, c.df)
global_area_means_standard <- as.numeric(
  readLines('frac_area_mean/global_area_mean_PFTs_standard.txt'))