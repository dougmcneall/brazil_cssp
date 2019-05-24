# es_ppe_ii_history_match.R
# A more formal history match for the ESPPEii ensemble.
# Uses all 500 members

setwd('analyze_u-ao732')
source('../per_pft.R')

load_ts_ensemble = function(fn, na.strings='-9.990000000000000000e+02', skip=1){
  dat = read.table(fn, header = FALSE, skip = skip, na.strings=na.strings)
  dat
}

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

# Express the "standard" runs (factor = 1) in terms of the
# latin hypercube design
X.stan = matrix(rep(1,32), nrow =1)
X.stan.norm = normalize(X.stan, wrt = lhs)

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

# ----------------------------------------------------------
# Level zero constraint
# ----------------------------------------------------------

# Constrain with global runoff and nbp, so that the emulator
# is not negatively affected by really bad points
level0.ix = which(dat.norm[,'runoff'] >0.5 & dat.norm[,'nbp'] > -10)
dat.level0  = dat.norm[level0.ix, ]
X.level0 = X[level0.ix, ]

# Visualise the constrained space ...
mins  = apply(X.level0, 2, min)
maxes = apply(X.level0, 2, max)

dev.new(width = 9, height = 10)
par(mfrow = c(3,2))
for(i in 1:6){
  hist(dat.norm[level0.ix,i], main = fnams[i])
}

dev.new(width = 10, height = 7)
par(mfrow = c(4, 8), mar = c(1,1,1,1))
for(i in 1:d){
  plot(lhs[level0.ix,i], dat.norm[level0.ix,2], axes = FALSE, xlab = '', ylab = '')
}

# Just looking at which inputs in the ensemble remain after constraint will tell us about
# which inputs are compatible with the constraints. 
ix.X.level1 = which(dat.norm[,'cs_gb'] > 750 & dat.norm[,'cs_gb'] < 3000 &
                  dat.norm[,'cv'] > 300 & dat.norm[,'cv'] < 800 & 
                  dat.norm[,'npp_n_gb'] > 35 &
                  dat.norm[,'npp_n_gb'] < 80)

X.level1 = X[ix.X.level1, ]

dev.new(width = 10, height = 10)
pairs(X.level1, gap = 0, xlim = c(0,1), ylim = c(0,1), lower.panel = NULL,
      pch = '.')


# Build emulators and do the constraint more thouroughly.
nsamp.unif = 99999
# The last row is the "standard" set of parameters
X.unif = rbind( samp.unif(nsamp.unif, mins = mins, maxes = maxes), X.stan.norm)

y.unif = matrix(nrow = nrow(X.unif), ncol = ncol(dat.level0))
colnames(y.unif) = colnames(dat.level0)

global.emlist = vector('list',length(fnams))


for(i in 1:ncol(y.unif)){
  em = twoStep.glmnet(X = X.level0, y = dat.level0[,i])
  global.emlist[[i]] = em
  pred = predict(em$emulator, newdata = X.unif, type = 'UK')
  y.unif[,i] = pred$mean
}

ix.kept = which(y.unif[,'cs_gb'] > 750 & y.unif[,'cs_gb'] < 3000 &
                  y.unif[,'cv'] > 300 & y.unif[,'cv'] < 800 & 
                  y.unif[,'npp_n_gb'] > 35 &
                  y.unif[,'npp_n_gb'] < 80)
X.kept = X.unif[ix.kept, ]

# we've removed 80% of our prior input space
(nrow(X.kept) / nsamp.unif) * 100

# for comparison, what does the emulator think the the "standard"
# parameters would produce?

dev.new(width = 9, height = 10)
par(mfrow = c(3,2))
for(i in 1:6){
  hist(dat.norm[level0.ix,i], main = fnams[i])
  rug(tail(y.unif,1)[, i], col = 'red', lwd = 2)
}

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
#pdf(file = 'graphics/ppe_ii/constraint_hists_standard.pdf', width = 8, height = 8)
dev.new()
par(mfrow = c(3,2), fg = 'white', las = 1)

hist(dat.norm[level0.ix,'runoff'], col = hcol, main = 'Runoff', xlab = 'Sv')
polygon(x = c(0.5, 100, 100, 0.5), y = c(0, 0, 1000, 1000), 
        col = makeTransparent('tomato2', alpha = 80))
rug(tail(y.unif,1)[,'runoff'], col = 'red', lwd = 3)

hist(dat.norm[level0.ix,'nbp'], col = hcol, main = 'NBP', xlab = 'GtC/year')
polygon(x = c(-10, 100, 100, -10), y = c(0, 0, 1000, 1000),
        col = makeTransparent('tomato2', alpha = 80))
rug(tail(y.unif,1)[,'nbp'], col = 'red', lwd = 3)

hist(dat.norm[level0.ix,'cs_gb'], col = hcol, main = 'Soil Carbon', xlab = 'GtC')

polygon(x = c(750, 3000, 3000, 750), y = c(0, 0, 1000, 1000),
        col = makeTransparent('tomato2', alpha = 80))
# AR5 numbers
polygon(x = c(1500, 2400, 2400, 1500), y = c(0, 0, 1000, 1000),
        col = makeTransparent('skyblue2', alpha = 80))

rug(tail(y.unif,1)[,'cs_gb'], col = 'red', lwd = 3)

hist(dat.norm[level0.ix,'cv'], col = hcol, main = 'Vegetation Carbon', xlab = 'GtC')
polygon(x = c(300, 800, 800, 300), y = c(0, 0, 1000, 1000),
        col = makeTransparent('tomato2', alpha = 80))

polygon(x = c(450, 650, 650, 450), y = c(0, 0, 1000, 1000),
        col = makeTransparent('skyblue2', alpha = 80))

rug(tail(y.unif,1)[,'cv'], col = 'red', lwd = 3)

hist(dat.norm[level0.ix,'npp_n_gb'], col = hcol , main = 'NPP', xlab = 'GtC/year')
polygon(x = c(35, 80, 80, 35), y = c(0, 0, 1000, 1000),
        col = makeTransparent('tomato2', alpha = 80))
rug(tail(y.unif,1)[,'npp_n_gb'], col = 'red', lwd = 3)


#hist(bl_frac_modern, col = hcol, main = 'Amazon Forest Fraction', xlab = 'fraction')
#polygon(x = c(0.5, 1, 1, 0.5), y = c(0, 0, 1000, 1000), col = makeTransparent('tomato2'))
#dev.off()

# Le Quere (2018) say that Ciais (2013) say that carbon stocks are:
# soil 1500-2400 GtC
# Veg 450-650 GtC
# Although it isn't clear what level of uncertainty that represents.
# what would that do to our input space?

ix.kept.AR5 = which(y.unif[,'cs_gb'] > 1500 & y.unif[,'cs_gb'] < 2400 &
                  y.unif[,'cv'] > 450 & y.unif[,'cv'] < 650 & 
                  y.unif[,'npp_n_gb'] > 35 &
                  y.unif[,'npp_n_gb'] < 80)
X.kept.AR5 = X.unif[ix.kept.AR5, ]

# we've removed 98.7% of our prior input space, including our standard set of parameters
# - chiefly by requiring a higher soil carbon that JULES is willing to simulate.
(nrow(X.kept.AR5) / nsamp.unif) * 100


# It's pretty clear that a problem with soil carbon is going to drive the 
# acceptance or rejection of input space, to a large degree. In that case, we have
# two options: 1) Accept the new space, (and also reject the standard parameters), 
# 2) Add a model discrepancy term and related uncertainty.


# First, it woudl be useful to have a simple sensitivity analysis for soil carbon.
X.oaat = oaat.design(X.level0, n=21, med = TRUE)
colnames(X.oaat) = colnames(X.level0)

twoStep.em = twoStep.glmnet(X=X.level0, y=dat.level0[,'cs_gb'])
oaat.pred = predict(twoStep.em$emulator, newdata = X.oaat, type = 'UK')

n = 21
dev.new(width = 8, height = 6)
par(mfrow = c(4,8), mar = c(2,3,2,0.3), oma = c(0.5,0.5, 3, 0.5))

y.oaat = oaat.pred$mean
y.upper = oaat.pred$mean+(2*oaat.pred$sd)
y.lower = oaat.pred$mean-(2*oaat.pred$sd)
ylim = range(y.oaat)

for(i in 1:d){
  
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  
  plot(X.oaat[ix,i], y.oaat[ix], type = 'l',
       ylab= '',ylim = ylim, axes = FALSE,
       main = '',
       xlab = '')
  lines(X.oaat[ix,i], y.upper[ix], col = 'grey')
  lines(X.oaat[ix,i],y.lower[ix], col = 'grey')
  
  axis(1, col = 'grey', col.axis = 'grey', las = 1)
  axis(2, col = 'grey', col.axis = 'grey', las = 1)
  mtext(3, text = colnames(lhs)[i], line = 0.2, cex = 0.7)  
  
}

soil.sens = sensvar(oaat.pred = oaat.pred, n=21, d=ncol(X.oaat))
soil.sens.sort = sort(soil.sens, decreasing = TRUE, index.return = TRUE)

#pdf(file = 'graphics/ppe_ii/soil_sens_level0.pdf', width = 7, height = 5)
dev.new(width = 7, height = 5)
par(mar = c(8,4,3,1))
plot(1:d, soil.sens.sort$x, axes = FALSE, pch = 19,
     xlab = '', ylab = 'OAAT Sensitivity Index')
segments(x0 = 1:d, y0 = rep(0,d), x1 = 1:d, y1 = soil.sens.sort$x)
axis(1, at = 1:d,  labels = colnames(X)[soil.sens.sort$ix], las = 3, cex.axis = 0.8)
axis(2,las =1)
#dev.off()

# We find that the soil carbon is most sensitive to kaps_roth
# Type:	real(4)
# Default:	3.22e-7, 9.65e-9, 2.12e-8, 6.43e-10
# Specific soil respiration rate for the RothC submodel for each soil carbon pool.
# Only used if using the TRIFFID vegetation model (l_triffid = TRUE),
# in which case soil carbon is modelled using four pools
# (biomass, humus, decomposable plant material, resistant plant material).
#
# we half-and-double kaps_roth in the ensemble.
#
# Soil carbon is second-most-sensitive to n_inorg_turnover
# Type:	real
# Default:	1.0

# Parameter controlling the lifetime of the inorganic N pool.
# A value of 1 implies the whole pool will turnover in 360 days.

# third most sensitive is alpha_io
# Type:	real(npft)
# Default:	None
#Quantum efficiency (mol CO2 per mol PAR photons).


# For NPP, Cramer et al (1999) found 44.4 - 66.3 PgC a year in models
# https://www.pik-potsdam.de/members/cramer/publications/edited-books/potsdam95/Cramer_1999b_GCB.pdf

# mean 53 range 40.5 - 78 in literature from Melillo (1993)
# http://www.as.wvu.edu/biology/bio463/Melillo%20et%20al%201993TEM%20NPP%20Estimations.pdf




