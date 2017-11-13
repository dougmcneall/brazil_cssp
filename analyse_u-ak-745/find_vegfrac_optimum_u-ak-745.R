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


best.corner = function(X.norm, dat, fn){
  # Find the ["best"] corner of the data that minimises output
  
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
  
  return(list(best.X = best.X, best.y = best.y, stepfit))
  
}

y = c(nle_wus_mae, recursive = TRUE)
dat = data.frame(y=y, x=X.norm)
colnames(dat) = c('y', colnames(lhs))
best.nle.wus = best.corner(X.norm = X.norm, dat = dat, fn = fn.step)

y = c(blet_sam_mae, recursive = TRUE)
dat = data.frame(y=y, x=X.norm)
colnames(dat) = c('y', colnames(lhs))

best.blet.sam = best.corner(X.norm = X.norm, dat = dat, fn = fn.step)


y = c(nle_glob_mae, recursive = TRUE)
dat = data.frame(y=y, x=X.norm)
colnames(dat) = c('y', colnames(lhs))
best.nle.glob = best.corner(X.norm = X.norm, dat = dat, fn = fn.step)

y = c(blet_glob_mae, recursive = TRUE)
dat = data.frame(y=y, x=X.norm)
colnames(dat) = c('y', colnames(lhs))

best.blet.glob = best.corner(X.norm = X.norm, dat = dat, fn = fn.step)


best.mat = rbind(best.nle.wus$best.X$par, best.blet.sam$best.X$par,
                 best.nle.glob$best.X$par, best.blet.glob$best.X$par)

parcoord(best.mat)


