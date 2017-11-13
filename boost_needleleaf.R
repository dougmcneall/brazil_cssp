# boost_needleleaf.R
# What parameter values would we choose to get more
# needleleaf forest globally and in the American West?

setwd("/Users/dougmcneall/Documents/work/R/brazil_cssp/analyse_u-ak-745")
source("/Users/dougmcneall/Documents/work/R/brazil_cssp/per_pft.R")

# load data
# All the data is in frac_area_mean

#surftypes = c('BLD','BLE_Trop','BLE_Temp','NLD',
#              'NLE','C3G','C3C','C3P','C4G','C4C',
#              'C4P','SHD','SHE','Urban','Lake','Bare Soil',
#              'Ice')

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

wus_area_means.list <- lapply(filelist.wus, c.df)
wus_area_means_standard <- as.numeric(
  readLines('frac_area_mean/WUS_area_mean_PFTs_standard.txt'))

sam_area_means.list <- lapply(filelist.sam, c.df)
sam_area_means_standard <- as.numeric(
  readLines('frac_area_mean/SAM_area_mean_PFTs_standard.txt'))


# ------------------------------------------------------
# Build emulators
# ------------------------------------------------------
lhs <- read.table('lhs_u-ak745.txt', header = TRUE)
d <- ncol(lhs)
cn <- colnames(lhs)
lhs.norm <- normalize(lhs)


# At the moment, let's just look at the area mean western US evergreen needleleaf
nle.wus <- wus_area_means.list[[5]]

# How do we maximise? Build a stepwise linear model to get rid of the less important inputs.

y = c(nle.wus, recursive = TRUE)

dat = data.frame(y=y, x=lhs.norm)
colnames(dat) = c('y', colnames(lhs))

initfit = lm(y ~ ., data = dat)
stepfit = step(initfit, direction="both", k=log(length(y)), trace=TRUE)



n = 21
X.oaat = oaat.design(lhs.norm, n, med = TRUE)
colnames(X.oaat) = colnames(lhs)
oaat.pred = predict(stepfit, newdata = data.frame(X.oaat))
#oaat.globalresponse[,i] = globalresponse.lm(oaat.pred, n = n, d = d)

dev.new(width = 8, height = 10)
par(mfrow = c(7,11), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = c(0, 0.1)

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], oaat.pred[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE)
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
#dev.off

# optimise to find the parameters
fn.step = function(newdata){
  newdata.df  = data.frame(matrix(newdata, nrow = 1))
  colnames(newdata.df) = colnames(lhs)
  out = predict(stepfit, newdata = newdata.df)
  -out
}

startin.mat <- matrix(rep(0.5, 74), nrow = 1)
startin <- data.frame(startin.mat)
colnames(startin) <- colnames(lhs)

testpred <- fn.step(rep(0.5, 74))

test <- optim(par = startin,
              fn = fn.step,
              method = "L-BFGS-B",
              lower = rep(0,74),
              upper = rep(1,74),
              control = list(maxit = 2000)
              )

test$par[test$par!=0.5]

plot(test$par)

test2 <- optim(par = startin,
              fn = fn.step,
              control = list(maxit = 2000)
              )

nd = data.frame(matrix(test$par, nrow = 1))
colnames(nd) = colnames(lhs)
predict(stepfit, newdata = nd)
hist(y)





