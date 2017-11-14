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

# Load the target data
nle_wus_mae = c(read.table('mean_abs_anom_NLE_wus.txt', header = FALSE, skip = 1), recursive = TRUE)
blet_sam_mae = c(read.table('mean_abs_anom_BLE_Trop_sam.txt', header = FALSE, skip = 1), recursive = TRUE)

nle_glob_mae = c(read.table('mean_abs_anom_NLE_glob.txt', header = FALSE, skip = 1), recursive = TRUE)
blet_glob_mae = c(read.table('mean_abs_anom_BLE_Trop_glob.txt', header = FALSE, skip = 1), recursive = TRUE)


# Load the input data
# Build emulators
lhs = read.table('lhs_u-ak745.txt', header = TRUE)
d = ncol(lhs)
cn = colnames(lhs)
X.norm = normalize(lhs)


# y = c(nle_wus_mae, recursive = TRUE)
# dat = data.frame(y=y, x=X.norm)
# colnames(dat) = c('y', colnames(lhs))
# 
# initfit = lm(y ~ ., data = dat)
# stepfit = step(initfit, direction="both", k=log(length(y)), trace=TRUE)


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
  # Find the ["best"] corner of the data that minimises output
  
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
       pch = 19, col = 'darkgrey'
       )

points(best.blet.glob$best.X$par,1:length(best.blet.glob$best.X$par),
       pch  = 19,
       cex = 0.8,
       col = 'dodgerblue')

axis(2, at = 1:74, labels = colnames(lhs), cex.axis = 0.4)
par(xpd = TRUE)
legend(x = 0, y = 80, 
       legend = c('BLET SAM', 'BLET Global'), 
       pch = 19, col = c('darkgrey','dodgerblue'),
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
       col = 'darkgrey')

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
       pch = 19, col = c('black','darkgrey', 'red', 'dodgerblue'),
       pt.cex = c(1,0.8,0.6,0.4),
       ncol = 2,
       cex = 0.5,
       bty = 'n'
       )
dev.off()


# with the South American Broadleaf Evergreen Tropical, we actually get a very small
# model back from the stepwise regression

# Fit stepwise models as a sensitivity analysis

dat = data.frame(y=c(blet_sam_mae, recursive = TRUE), x=X.norm)
colnames(dat) = c('y', colnames(lhs))
initfit = lm(y ~ ., data = dat)
stepfit = step(initfit, direction="both", k=log(length(y)), trace=TRUE)



# However, it looks like we can't get close to "zero" error
hist(best.blet.sam$y, xlim = c(0, max(best.blet.sam$y)), main = 'SAM BLE')
rug(best.blet.sam$best.y, lwd = 3, col = 'red')

coefplot(best.blet.sam$stepfit)

# with global Broadleaf Evergreen Tropical, we get a larger model back
y = blet_glob_mae
dat = data.frame(y=y, x=X.norm)
colnames(dat) = c('y', colnames(lhs))

initfit = lm(y ~ ., data = dat)
stepfit = step(initfit, direction="both", k=log(length(y)), trace=TRUE)
# Much larger model for the global fit
coefplot(best.blet.glob$stepfit)
# As in the regional, we can roughly half the MAE
hist(y, xlim = c(0, max(y)))
rug(best.blet.glob$best.y, lwd = 3, col = 'red')


# Regional NLE
y = c(nle_wus_mae, recursive = TRUE)
dat = data.frame(y=y, x=X.norm)
colnames(dat) = c('y', colnames(lhs))

initfit = lm(y ~ ., data = dat)
stepfit = step(initfit, direction="both", k=log(length(y)), trace=TRUE)

coefplot(stepfit)
# As in the regional, we can roughly half the MAE
hist(y, xlim = c(0, max(y)))
rug(best.nle.wus$best.y, lwd = 3, col = 'red')

# Global NLE
y = c(nle_glob_mae, recursive = TRUE)
dat = data.frame(y=y, x=X.norm)
colnames(dat) = c('y', colnames(lhs))

initfit = lm(y ~ ., data = dat)
stepfit = step(initfit, direction="both", k=log(length(y)), trace=TRUE)

coefplot(stepfit)
# As in the regional, we can roughly half the MAE
hist(y, xlim = c(0, max(y)))
rug(best.nle.glob$best.y, lwd = 3, col = 'red')




