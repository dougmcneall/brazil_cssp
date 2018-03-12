# es_ppe_ii.R
# Analyse the full 500 (ish member ensemble)

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
lhs_i = read.table('data/ES_PPE_i/lhs_u-ao732.txt', header = TRUE)
lhs_ii = read.table('data/ES_PPE_ii/lhs_u-ao732a.txt', header = TRUE)

lhs = rbind(lhs_i, lhs_ii)

#stanparam = read.table('data/stanparms_u-ao732.txt', header = TRUE)
#stanparam = matrix(c(rep(1, 31),0.01), nrow = 1)
#stanparam.norm = normalize.wrt(stanparam, lhs)

X = normalize(lhs)
colnames(X) = colnames(lhs)
d = ncol(X)

fnvec = dir('data/ES_PPE_ii', pattern = 'Annual.Amazon')
fnlocvec = paste0('data/ES_PPE_ii/', fnvec)

dat = load_ts_ensemble(fnlocvec[4])


pdf(width = 10, height = 10, file = 'test.pdf')
par(mfrow =c(6,6), mar = c(2,2,2,2))
for(i in 1:d){
  
  plot(lhs[-500,i],dat[,1], pch = '.', axes = FALSE)
  
}
dev.off()


