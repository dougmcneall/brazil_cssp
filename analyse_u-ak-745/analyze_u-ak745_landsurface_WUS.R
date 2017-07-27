# setwd("/Users/dougmcneall/Documents/work/R/brazil_cssp/analyse_u-ak-745")
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

#read data
NLE = c(read.table('WUS_NLE_mean.txt', header = FALSE, skip = 1), recursive = TRUE)
NLD = c(read.table('WUS_NLD_mean.txt', header = FALSE, skip = 1), recursive = TRUE)
C3G = c(read.table('WUS_C3G_mean.txt', header = FALSE, skip = 1), recursive = TRUE)
BLE_Temp = c(read.table('WUS_BLE_TEMP_mean.txt', header = FALSE, skip = 1), recursive = TRUE)
SHE = c(read.table('WUS_SHE_mean.txt', header = FALSE, skip = 1), recursive = TRUE)

NLE.standard = 0.105954818915
NLD.standard = 0.0119113730396
C3G.standard = 0.0406237504424
BLE_Temp.standard = 0.0135605978606



breaks = seq(from = 0, to =0.2, by = 0.02)
xlim = c(0,0.2)
dev.new(width=7, height=7)
#pdf(file = 'fractypes_hists.pdf', width = 7, height = 7)
par(mfrow = c(3,2), fg = 'white', las = 1)

lwd = 2
hist(NLE, xlim = xlim, col = 'grey', axes = FALSE, breaks = breaks, main = 'NLE')
axis(1, col = 'black'); axis(2, col = 'black')
rug(NLE.standard, col = 'red', lwd = lwd)

hist(NLD, xlim = xlim, col = 'grey', axes = FALSE, breaks = breaks, main = 'NLD')
axis(1, col = 'black'); axis(2, col = 'black')
rug(NLD.standard, xlim = xlim, col = 'red', lwd = lwd)

hist(C3G, xlim = xlim, col = 'grey', axes = FALSE, breaks = breaks, main = 'C3G')
axis(1, col = 'black'); axis(2, col = 'black')
rug(C3G.standard, col = 'red', lwd = lwd)

hist(BLE_Temp, xlim = xlim, col = 'grey', axes = FALSE, breaks = breaks, main = 'BLE_Temp')
axis(1, col = 'black'); axis(2, col = 'black')
rug(BLE_Temp.standard, col = 'red', lwd = lwd)

hist(SHE, xlim = xlim, col = 'grey', axes = FALSE, breaks = breaks, main = 'SHE')
axis(1, col = 'black'); axis(2, col = 'black')


#dev.off()

dat = cbind(NLE, NLD, C3G, BLE_Temp, SHE)
names(dat) = c('NLE', 'NLD', 'C3G', 'BLE_Temp', 'SHE')
#pdf(file = 'fractypes_pairs.pdf',width = 5, height = 5)
pairs(dat, upper.panel = NULL, xlim = c(0,0.2), ylim = c(0,0.2))
#dev.off()



codet = function(x1, x2){
  fit = lm(x1~x2)
  summary(fit)$r.squared
}

regcoef = function(x1, x2){
  fit = lm(x1~x2)
  summary(fit)$coef
}

codet(NLE, NLD)
regcoef(NLE, NLD)

codet(C3G, NLD)
regcoef(C3G, NLD)

codet(BLE_Temp, NLD)
regcoef(BLE_Temp, NLD)

codet(SHE, NLD)
regcoef(SHE, NLD)


codet(SHE, BLE_Temp)
regcoef(SHE, BLE_Temp)

codet(SHE, C3G)
regcoef(SHE, C3G)




