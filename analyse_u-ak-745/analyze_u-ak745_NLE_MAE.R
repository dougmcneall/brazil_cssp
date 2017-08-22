setwd("/Users/dougmcneall/Documents/work/R/brazil_cssp/analyse_u-ak-745")
## Try and build an emulator
library(devtools)

install_github(repo = "dougmcneall/hde")

library(hde)
library(RColorBrewer)
library(fields)
library(MASS)
library(DiceKriging)

source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/imptools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/vistools.R")

yg = brewer.pal(9, "YlGn")
ryb = brewer.pal(9, "RdYlBu")
byr = rev(ryb)

# To first order, if we grow more NLE forest, we minimise absolute error
NLE = c(read.table('WUS_NLE_mean.txt', header = FALSE, skip = 1), recursive = TRUE)
NLE.mae.wus = c(read.table('test_absanom.txt', header = FALSE, skip = 1), recursive = TRUE)


# How do we minimise absolute error?
# What are the most important parameters, and how do we alter them to
# minimise error?
surftypes = c('BLD','BLE_Trop','BLE_Temp','NLD',
              'NLE','C3G','C3C','C3P','C4G','C4C',
              'C4P','SHD','SHE','Urban','Lake','Bare Soil',
              'Ice')

# ------------------------------------------------------
# Rank the parameters
# ------------------------------------------------------
# Build emulators
lhs = read.table('lhs_u-ak745.txt', header = TRUE)
d = ncol(lhs)
cn = colnames(lhs)
lhs.norm = normalize(lhs)


NLE.mae.wus.fit = km(~., design = lhs.norm, response = NLE.mae.wus)

n = 21
X.oaat= oaat.design(lhs.norm, n, med = TRUE)
NLE.oaat = predict(NLE.mae.wus.fit, newdata = X.oaat, type = 'UK')

dev.new(width = 8, height = 10)
par(mfrow = c(7,11), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = c(0.11, 0.15)

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], NLE.oaat$mean[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE)
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
#dev.off

sensvar = function(oaat.pred, n, d){
  # Calculate variance as a global sensitivity meansure
  out = rep(NA,d)
  for(i in 1:d){
    ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
    out[i] = var(oaat.pred$mean[ix])
  }
  out
}

NLE.sensvar = sensvar(NLE.oaat, n = n, d = d)

ix.sorted = sort(NLE.sensvar, decreasing = TRUE, index.return = TRUE)$ix

dev.new()
dotchart((NLE.sensvar[ix.sorted])[1:10], labels = (colnames(lhs)[ix.sorted])[1:10])


# seems to be a bug with the 'two at a time'
# x.taat = taat.design(X = lhs.norm, 21)

x.unif = samp.unif(100000, mins = rep(0, ncol(lhs)), maxes = rep(1, ncol(lhs)) )
pred.unif = predict(NLE.mae.wus.fit, newdata = x.unif, type = 'UK')

quilt.plot(x = x.unif[ , ix.sorted[1]],
           y = x.unif[, ix.sorted[2]],
           z = pred.unif$mean,
           xlab = colnames(lhs)[ix.sorted[1]],
           ylab = colnames(lhs)[ix.sorted[2]]
           )

quilt.plot(x = x.unif[ , ix.sorted[3]],
           y = x.unif[, ix.sorted[4]],
           z = pred.unif$mean,
           xlab = colnames(lhs)[ix.sorted[3]],
           ylab = colnames(lhs)[ix.sorted[4]]
)

quilt.plot(x = x.unif[ , ix.sorted[5]],
           y = x.unif[, ix.sorted[6]],
           z = pred.unif$mean,
           xlab = colnames(lhs)[ix.sorted[5]],
           ylab = colnames(lhs)[ix.sorted[6]]
)







