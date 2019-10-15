# es_ppe_ii_history_match.R
# A more formal history match for the ESPPEii ensemble.
# Uses all 500 members

setwd('analyze_u-ao732')
source('../per_pft.R')

load_ts_ensemble = function(fn, na.strings='-9.990000000000000000e+02', skip=1){
  dat = read.table(fn, header = FALSE, skip = skip, na.strings=na.strings)
  dat
}

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

years = 1861:2014
ysec = 60*60*24*365
norm.vec = c(1e12, 1e12, 1e12/ysec , 1e12, 1e12, 1e9)

# Load up the data
lhs_i = read.table('data/ES_PPE_i/lhs_u-ao732.txt', header = TRUE)
lhs_ii = read.table('data/ES_PPE_ii/lhs_u-ao732a.txt', header = TRUE)

toplevel.ix = 1:499

lhs = rbind(lhs_i, lhs_ii)[toplevel.ix, ]

X.raw = normalize(lhs)
colnames(X.raw) = colnames(lhs)
d = ncol(X.raw)

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
datmat.raw = matrix(nrow = nrow(X.raw), ncol = length(fnlocvec))
for(i in 1:length(fnlocvec)){
  dat = load_ts_ensemble(fnlocvec[i])[toplevel.ix, ]
  dat.modern = dat[ ,135:154]
  mean.modern = apply(dat.modern, 1, mean)
  datmat.raw[ , i] = mean.modern
}
<<<<<<< HEAD
colnames(datmat.raw) = fnams

#remove.ix = c(1, which(duplicated(datmat.raw[,1])))

#datmat = datmat.raw[-remove.ix, ]
#X = X.raw[-remove.ix, ]
X = X.raw

dat.norm = sweep(datmat.raw, 2, norm.vec, FUN = '/')
p = ncol(dat.norm)

# My suspicion is that the duplicated rows should be NA, and need to
# be removed from analysis as they will affect the emulator.

# Start by removing only the duplicated (e.g. later) rows, and keeping
# the initial rows


=======
colnames(datmat) = fnams
dat.norm = sweep(datmat, 2, norm.vec, FUN = '/')
>>>>>>> 6652f158f40d90f925a47092f6bb521ff6b92c07

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

nlevel0.ix = setdiff(allix, level0.ix)
X.nlevel0 = X[nlevel0.ix, ]


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

<<<<<<< HEAD
# Parallel Coordinates plot of NROY and ruled out members, level 0
dev.new(width = 20, height = 9)
par(mfrow = c(2,1), las = 2, mar = c(7,4,4,1), cex.axis = 0.8)
parcoord(X.level0, col = makeTransparent('black', 100), ylim = c(0,1),
         main = 'level 0 NROY ensemble members')

parcoord(X.nlevel0, col = makeTransparent('black', 100), ylim = c(0,1),
         main = 'level 0 ruled out ensemble members')


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

# Parallel Coordinates plot of NROY and ruled out members, level 0
dev.new(width = 20, height = 9)
par(mfrow = c(2,1), las = 2, mar = c(7,4,4,1), cex.axis = 0.8)
parcoord(X.level1, col = makeTransparent('black', 100), ylim = c(0,1),
         main = 'level 1 NROY ensemble members')

parcoord(X.nlevel1, col = makeTransparent('black', 100), ylim = c(0,1),
         main = 'level 1 ruled out ensemble members')



# Pairs plot of level 1 constraint
dev.new(width = 10, height = 10)
pairs(X.level1, gap = 0, lower.panel = NULL,
      labels = 1:d,
      xlim = c(0,1), ylim = c(0,1),
      col = makeTransparent('black', 150),
      pch = 21
      )
par(xpd = NA)
#text(0.2, 0.6, labels = paste0(colnames(lhs)[1:5],'\n'))

legend('left', legend = paste(1:d, colnames(lhs)), cex = 0.9, bty = 'n')



dev.new(width = 12, height = 11)
par(mfrow = c(8,4), mar = c(4,3,2,1), oma = c(0,0,2,0))
for(i in 1:ncol(X)){
#barplot(X.level1[,i], xlim = c(0,1), col = 'black')
#plot(table(cut(X.level1[,i], breaks = seq(from = 0, to = 1, length.out = 7))),
#     ylab = '', xlab = colnames(X)[i])

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


# OK, we know the
# Making decisions on choosing parameter sets. Are there any margins we
# can rule out?

apply(X.level1,2, range)





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



=======
# Just looking at which inputs in the ensemble remain after constraint will tell us about
# which inputs are compatible with the constraints. 
ix.X.level1 = which(dat.norm[,'cs_gb'] > 750 & dat.norm[,'cs_gb'] < 3000 &
                  dat.norm[,'cv'] > 300 & dat.norm[,'cv'] < 800 & 
                  dat.norm[,'npp_n_gb'] > 35 &
                  dat.norm[,'npp_n_gb'] < 80)

X.level1 = X[ix.X.level1, ]
>>>>>>> 6652f158f40d90f925a47092f6bb521ff6b92c07

dev.new(width = 10, height = 10)
pairs(X.level1, gap = 0, xlim = c(0,1), ylim = c(0,1), lower.panel = NULL,
      pch = '.')

<<<<<<< HEAD
y.unif = matrix(nrow = nsamp.unif, ncol = ncol(dat.level0))
colnames(y.unif) = colnames(dat.norm)
=======

# Build emulators and do the constraint more thouroughly.
nsamp.unif = 99999
# The last row is the "standard" set of parameters
X.unif = rbind( samp.unif(nsamp.unif, mins = mins, maxes = maxes), X.stan.norm)

y.unif = matrix(nrow = nrow(X.unif), ncol = ncol(dat.level0))
colnames(y.unif) = colnames(dat.level0)
>>>>>>> 6652f158f40d90f925a47092f6bb521ff6b92c07

global.emlist = vector('list',length(fnams))


for(i in 1:ncol(y.unif)){
  em = twoStep.glmnet(X = X.level0, y = dat.level0[,i])
  global.emlist[[i]] = em
  pred = predict(em$emulator, newdata = X.unif, type = 'UK')
  y.unif[,i] = pred$mean
}

ix.kept = which(
  y.unif[,'cs_gb'] > 750 & y.unif[,'cs_gb'] < 3000 &
  y.unif[,'cv'] > 300 & y.unif[,'cv'] < 800 & 
  y.unif[,'npp_n_gb'] > 35 & y.unif[,'npp_n_gb'] < 80 &
  y.unif[,'runoff'] >0.5 &
  y.unif[,'nbp'] > -10)
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


ix.rejected = setdiff(1:nsamp.unif, ix.kept)
X.rejected = X.unif[ix.rejected, ]


#pdf(file = 'pcpTest.pdf', width = 20, height = 4)
dev.new(width = 20, height = 9)
par(las = 2, cex.axis = 0.8, mfrow = c(2,1))

parcoord(X.kept, col = makeTransparent('black', 10), ylim = c(0,1))
parcoord(X.rejected[1:2000,], col = makeTransparent('black', 10), ylim = c(0,1))

#dev.off()


dev.new(width = 12, height = 12)
pairs(X.kept[1:5000,] , pch = '.', gap = 0, upper.panel = NULL,
      xlim = c(0,1), ylim = c(0,1),
      col = makeTransparent('black', 20)
      )
 

blues = brewer.pal(9, 'Blues')
#pdf(file = 'graphics/ppe_ii/pairs_dens_all_constraints0.pdf', width = 10, height = 10)
dev.new(width = 10, height = 10)
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

#dev.off()


# Next, run the one-at-a-time sensitivity measures again






