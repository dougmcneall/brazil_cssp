# find_vegfrac_optimum_u-ak-745.R
# Which directions should we move the ensemble
# in, if we want to minimise error in vegetation fraction?
# Using Brazil and Western US as example.
setwd("/Users/dougmcneall/Documents/work/R/brazil_cssp/analyse_u-ak-745")
## Try and build an emulator
library(devtools)

#install_github(repo = "dougmcneall/hde")

library(hde)
library(RColorBrewer)
library(fields)
library(MASS)
library(DiceKriging)
library(coefplot)

source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/imptools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/vistools.R")

yg = brewer.pal(9, "YlGn")
ryb = brewer.pal(9, "RdYlBu")
byr = rev(ryb)

# Load the target data - mean absolute (pointwise) error of 
# veg fraction for Needleleaf Evergreen and Tropical Broadleaf Evergreen, globally and
# regionall (Western US and South America, respectively)
nle_wus_mae = c(read.table('mean_abs_anom_NLE_wus.txt', header = FALSE, skip = 1), recursive = TRUE)
blet_sam_mae = c(read.table('mean_abs_anom_BLE_Trop_sam.txt', header = FALSE, skip = 1), recursive = TRUE)

nle_glob_mae = c(read.table('mean_abs_anom_NLE_glob.txt', header = FALSE, skip = 1), recursive = TRUE)
blet_glob_mae = c(read.table('mean_abs_anom_BLE_Trop_glob.txt', header = FALSE, skip = 1), recursive = TRUE)

# Load and normalize the input data
lhs = read.table('lhs_u-ak745.txt', header = TRUE)
d = ncol(lhs)
cn = colnames(lhs)
X.norm = normalize(lhs)

# This is the cost function, gets fed to optim 
fn.step = function(newdata, cn, stepfit){
  newdata.df  = data.frame(matrix(newdata, nrow = 1))
  colnames(newdata.df) = cn
  out = predict(stepfit, newdata = newdata.df)
  out
}

# # initial values for optim in the middle of the design
 startin.mat <- matrix(rep(0.5, ncol(X.norm)), nrow = 1)
 startin <- data.frame(startin.mat)
 colnames(startin) <- colnames(X.norm)
# 
# # Find the values which minimise needleleaf absolute error
# best.X <- optim(par = startin,
#                 fn = fn.step,
#                 method = "L-BFGS-B",
#                 lower = rep(0,ncol(X.norm)),
#                 upper = rep(1,ncol(X.norm)),
#                 control = list(maxit = 2000),
#                 cn = colnames(X.norm)
# )
# 
# best.X$par[best.X$par!=0.5]


best.corner = function(X.norm, y, fn){
  # Find the ["best"] corner of the data that minimises output (e.g. MAE)
  
   dat = data.frame(y=y, x=X.norm)
  colnames(dat) = c('y', colnames(X.norm))
  
  initfit = lm(y ~ ., data = dat)
  stepfit = step(initfit, direction="both", k=log(length(y)), trace=TRUE)
  
  startin.mat <- matrix(rep(0.5, ncol(X.norm)), nrow = 1)
  startin <- data.frame(startin.mat)
  colnames(startin) <- colnames(X.norm)
  
  best.X <- optim(par = startin,
                  fn = fn,
                  method = "L-BFGS-B",
                  lower = rep(0,ncol(X.norm)),
                  upper = rep(1,ncol(X.norm)),
                  control = list(maxit = 2000),
                  cn = colnames(X.norm),
                  stepfit = stepfit
  )
  
  nd = data.frame(matrix(best.X$par, nrow = 1))
  colnames(nd) = colnames(X.norm)
  best.y = predict(stepfit, newdata = nd)
  
  return(list(best.X = best.X, best.y = best.y, stepfit = stepfit, y=y))
  
}

best.nle.wus = best.corner(X.norm = X.norm, y = nle_wus_mae, fn = fn.step)
best.nle.glob = best.corner(X.norm = X.norm, y = nle_glob_mae, fn = fn.step)
best.blet.sam = best.corner(X.norm = X.norm, y = blet_sam_mae, fn = fn.step)
best.blet.glob = best.corner(X.norm = X.norm, y = blet_glob_mae, fn = fn.step)

# Visualise the "corners" that you should explore, to minimise the 
# error.
pdf(file = 'corners_NLE.pdf', width = 2.5, height = 7)
par(mar = c(2,5,3,1),las = 1)
plot(best.nle.wus$best.X$par,1:length(best.nle.wus$best.X$par),
     axes = FALSE, xlab = '', ylab = '',
     pty = 'n')
abline(h = 1:74, lty = 'dashed', col = 'grey')

points(best.nle.wus$best.X$par,1:length(best.nle.wus$best.X$par),
     pch = 20,
     col = 'black')

points(best.nle.glob$best.X$par,1:length(best.nle.glob$best.X$par),
       pch  = 20,
       col = 'red')

axis(2, at = 1:74, labels = colnames(lhs), cex.axis = 0.4)
par(xpd = TRUE)
legend(x = 0, y = 80, 
       legend = c('NLE WUS', 'NLE Global'), 
       pch = 19, col = c('black','red'),
       pt.cex = c(1,0.8),
       ncol = 2,
       cex = 0.5,
       bty = 'n'
)

dev.off()

pdf(file = 'corners_BLET.pdf', width = 2.5, height = 7)
par(mar = c(2,5,3,1),las = 1)
plot(best.blet.sam$best.X$par,1:length(best.blet.sam$best.X$par),
     axes = FALSE, xlab = '', ylab = '',
     pty = 'n')
abline(h = 1:74, lty = 'dashed', col = 'grey')

points(best.blet.sam$best.X$par,1:length(best.blet.sam$best.X$par),
       pch = 19, col = 'orange'
       )

points(best.blet.glob$best.X$par,1:length(best.blet.glob$best.X$par),
       pch  = 19,
       cex = 0.8,
       col = 'dodgerblue')

axis(2, at = 1:74, labels = colnames(lhs), cex.axis = 0.4)
par(xpd = TRUE)
legend(x = 0, y = 80, 
       legend = c('BLET SAM', 'BLET Global'), 
       pch = 19, col = c('orange','dodgerblue'),
       pt.cex = c(1,0.8),
       ncol = 2,
       cex = 0.5,
       bty = 'n'
)
dev.off()


pdf(file = 'corners.pdf', width = 2.5, height = 7)
par(mar = c(2,5,3,1),las = 1)
plot(best.nle.wus$best.X$par,1:length(best.nle.wus$best.X$par),
     axes = FALSE, xlab = '', ylab = '',
     pty = 'n',
     pch = 19)

abline(h = 1:74, lty = 'dashed', col = 'grey')

points(best.nle.wus$best.X$par,1:length(best.nle.wus$best.X$par),
       pch = 19)

points(best.blet.sam$best.X$par,1:length(best.blet.sam$best.X$par),
       pch  = 19, cex = 0.8,
       col = 'orange')

axis(2, at = 1:74, labels = colnames(lhs), cex.axis = 0.4)

points(best.nle.glob$best.X$par,1:length(best.nle.glob$best.X$par),
       col = 'red', pch = 19, cex = 0.6
)

points(best.blet.glob$best.X$par,1:length(best.blet.glob$best.X$par),
       col = 'dodgerblue',
       pch = 19, cex = 0.4)

par(xpd = TRUE)
legend(x = 0, y = 80, 
       legend = c('NLE WUS', 'BLET SAM', 'NLE Global', 'BLET Global'), 
       pch = 19, col = c('black','orange', 'red', 'dodgerblue'),
       pt.cex = c(1,0.8,0.6,0.4),
       ncol = 2,
       cex = 0.5,
       bty = 'n'
       )
dev.off()


# Summarise the "optimum" output against the ensemble, and
# plot the regression coefficients of the stepwise model.

coefplot.doug = function(fit, ...){
  # Plot coefficients and uncertainty from a linear model.
  y = coef(fit)
  segs = confint(fit)
  
  plot(y, 1:length(y), pty = 'n', axes = FALSE, xlab = '', ylab = '',
      xlim = range(segs) )
  abline(v=0, col = 'grey')
  points(y, 1:length(y), ...)
  segments(x0 = segs[,1], y0 = 1:length(y), x1 = segs[,2], y1 = 1:length(y))
  axis(1)
  axis(2, labels = names(y), at = 1:length(y), las = 1, ...)
}


pdf(file = 'blet_optimum_output.pdf', width = 7, height = 7)
par(fg = 'white', mfrow = c(2,2))
hist(best.blet.sam$y, xlim = c(0, max(best.blet.sam$y)),
     main = 'Regional Broadleaf Tropical Evergreen',
     xlab = 'Mean abs. fraction error',
     col = 'grey',
     axes = FALSE)
par(fg = 'black')
axis(1)
axis(2)

rug(best.blet.sam$best.y, lwd = 3, col = 'red')
legend('topleft', pch = '|', pt.lwd = 5, col = 'red', legend = 'optimum', bty = 'n')
par(mar = c(5.1, 6.1, 4.1, 2.1))
coefplot.doug(best.blet.sam$stepfit, pch = 19, cex.axis = 0.8)

par(fg = 'white',mar = c(5.1, 4.1, 4.1, 2.1))
hist(best.blet.glob$y, xlim = c(0, max(best.blet.glob$y)),
     xlab = 'Mean abs. fraction error',
     col = 'grey',
     main = 'Global Broadleaf Tropical Evergreen',
     axes = FALSE)
par(fg = 'black')
axis(1)
axis(2)

rug(best.blet.glob$best.y, lwd = 3, col = 'red')
par(mar = c(5.1, 6.1, 4.1, 2.1))
coefplot.doug(best.blet.glob$stepfit, pch = 19, cex.axis = 0.8)
dev.off()



# NLE
pdf(file = 'nle_optimum_output.pdf', width = 7, height = 7)
par(fg = 'white', mfrow = c(2,2))
hist(best.nle.wus$y, xlim = c(0, max(best.nle.wus$y)),
     main = 'Regional Needleleaf Evergreen',
     xlab = 'Mean abs. fraction error',
     col = 'grey',
     axes = FALSE)
par(fg = 'black')
axis(1)
axis(2)

rug(best.nle.wus$best.y, lwd = 3, col = 'red')
legend('topleft', pch = '|', pt.lwd = 5, col = 'red', legend = 'optimum', bty = 'n')
par(mar = c(5.1, 6.1, 4.1, 2.1))
coefplot.doug(best.nle.wus$stepfit, pch = 19, cex.axis = 0.8)

par(fg = 'white',mar = c(5.1, 4.1, 4.1, 2.1))
hist(best.nle.glob$y, xlim = c(-0.02, max(best.nle.glob$y)),
     xlab = 'Mean abs. fraction error',
     col = 'grey',
     main = 'Global Needleleaf Evergreen',
     axes = FALSE)
abline(v = 0, col = 'grey')
par(fg = 'black')
axis(1)
axis(2)

rug(best.nle.glob$best.y, lwd = 3, col = 'red')
par(mar = c(5.1, 6.1, 4.1, 2.1))
coefplot.doug(best.nle.glob$stepfit, pch = 19, cex.axis = 0.8)
dev.off()

# Visualise the response to a few of the most important inputs


# Sample from the stepwise emulator
stepfit.sam = best.blet.sam$stepfit

X.unif = data.frame(samp.unif(100000, mins = rep(0, ncol(X.norm)), maxes = rep(1, ncol(X.norm))))
colnames(X.unif) = colnames(X.norm)

pred.unif.sam = predict(stepfit.sam, newdata = X.unif)


# Which are the largest model coefficients? (i.e. most important inputs)
# labels that we want to pull out from the input
important.X = names(coef(stepfit.sam))[order(abs(coef(stepfit.sam)), decreasing = TRUE)]
important.X = important.X[-1] # although we don't want the intercept.

ix = match(important.X , colnames(X.norm))
colnames(X.norm)[ix]

pdf(file = 'quilts_sam.pdf', width = 8, height = 8)
par(mfrow = c(2,2), oma = c(1,1,1,1), mar = c(5,5,4,2))

quilt.plot(x = X.unif[ , ix[1]],
           y = X.unif[ , ix[2]],
           z = pred.unif.sam,
           xlab = colnames(X.norm)[ix[1]],
           ylab = colnames(X.norm)[ix[2]],
           col = byr
)

quilt.plot(x = X.unif[ , ix[3]],
           y = X.unif[ , ix[4]],
           z = pred.unif.sam,
           xlab = colnames(X.norm)[ix[3]],
           ylab = colnames(X.norm)[ix[4]],
           col = byr
)

quilt.plot(x = X.unif[ , ix[4]],
           y = X.unif[ , ix[5]],
           z = pred.unif.sam,
           xlab = colnames(X.norm)[ix[4]],
           ylab = colnames(X.norm)[ix[5]],
           col = byr
)


dev.off()






