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
SHE.standard = 0.111781059264

# 0.00837554171292
# 0.00508592070801
# 0.0135605978606
# 0.0119113730396
# 0.105954818915
# 0.0406237504424
# 8.57981544253e-05
# 3.98540348509e-05
# 0.0121790902756
# 3.78826704019e-05
# 0.000118142083638
# 0.00149984928681
# 0.111781059264

#obs = read.table('lccci_wus_fractypes_mean.txt', head = FALSE)
#obs_names = read.table('lccci_wus_fractypes_mean_names.txt', head = FALSE)


obs = as.numeric(readLines('lccci_wus_fractypes_mean.txt'))
obs_names = readLines('lccci_wus_fractypes_mean_names.txt')

obs_tab = data.frame(pft = c(obs_names), frac = c(obs))

chosen.pfts = c('NLE', 'NLD', 'C3G', 'BLE_Temp', 'SHE')

obs.ix = match(chosen.pfts, obs_names)
obs.vec = obs[obs.ix]

standard.vec = c(NLE.standard, NLD.standard, C3G.standard, BLE_Temp.standard, SHE.standard)

pdf(file = 'frac_comparison.pdf', width = 7, height = 5)
par(las = 1)
plot(1:5, standard.vec, axes = FALSE, ylim = c(0,0.3), pch = 19,
     main = 'Land surface fractions in the Western US',
     xlab = 'PFT', ylab = 'Mean fraction')
points(obs.vec, col = 'red', pch = 19)
axis(1, at =1:5, labels = chosen.pfts )
axis(2)
legend('topleft', pch = 19, col = c('black','red'),
       legend = c('JULES standard parameters', 'observations'),
       box.col = 'grey90', bg = 'grey90')
dev.off()




breaks = seq(from = 0, to =0.2, by = 0.02)
xlim = c(0,0.2)
dev.new(width=7, height=7)
pdf(file = 'fractypes_hists.pdf', width = 7, height = 7)
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
rug(SHE.standard, col = 'red', lwd = lwd)

dev.off()


dat = cbind(NLE, NLD, C3G, BLE_Temp, SHE)
names(dat) = c('NLE', 'NLD', 'C3G', 'BLE_Temp', 'SHE')
pdf(file = 'fractypes_pairs.pdf',width = 5, height = 5)
pairs(dat, upper.panel = NULL, xlim = c(0,0.2), ylim = c(0,0.2))
dev.off()


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

# ------------------------------------------------------
# Build emulators
# ------------------------------------------------------
lhs = read.table('lhs_u-ak745.txt', header = TRUE)
d = ncol(lhs)
cn = colnames(lhs)
lhs.norm = normalize(lhs)


NLE.fit = km(~., design = lhs.norm, response = NLE)
NLD.fit = km(~., design = lhs.norm, response = NLD)
C3G.fit = km(~., design = lhs.norm, response = C3G)
BLE_Temp.fit = km(~., design = lhs.norm, response = BLE_Temp)
SHE.fit = km(~., design = lhs.norm, response = SHE)

n = 21
X.oaat= oaat.design(lhs.norm, n, med = TRUE)


NLE.oaat = predict(NLE.fit, newdata = X.oaat, type = 'UK')
NLD.oaat = predict(NLD.fit, newdata = X.oaat, type = 'UK')
C3G.oaat = predict(C3G.fit, newdata = X.oaat, type = 'UK')
SHE.oaat = predict(SHE.fit, newdata = X.oaat, type = 'UK')
BLE_Temp.oaat = predict(BLE_Temp.fit, newdata = X.oaat, type = 'UK')


dev.new(width = 8, height = 10)
par(mfrow = c(7,11), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = c(0, 0.5)

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

dev.new(width = 8, height = 10)
par(mfrow = c(7,11), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = c(0, 0.5)

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], NLD.oaat$mean[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE)
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
#dev.off

dev.new(width = 8, height = 10)
par(mfrow = c(7,11), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = c(0, 0.5)

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], C3G.oaat$mean[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE)
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
#dev.off


#pdf(file = 'oaat.pdf', width = 8, height = 10)
dev.new(width = 8, height = 10)
par(mfrow = c(7,11), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = c(-0.12, 0.4)

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], SHE.oaat$mean[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE)
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
#dev.off()

#pdf(file = 'oaat.pdf', width = 8, height = 10)
dev.new(width = 8, height = 10)
par(mfrow = c(7,11), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = c(-0.12, 0.4)

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], BLE_Temp.oaat$mean[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE)
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
#dev.off()

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
NLD.sensvar = sensvar(NLD.oaat, n = n, d = d)
C3G.sensvar = sensvar(C3G.oaat, n = n, d = d)
SHE.sensvar = sensvar(SHE.oaat, n = n, d = d)
BLE_Temp.sensvar = sensvar(BLE_Temp.oaat, n = n, d = d)

plot(1:d, NLE.sensvar, ylim = c(0,0.0005))
points(1:d, NLD.sensvar, col = 'red')

sensmat = as.matrix(cbind(NLE.sensvar, NLD.sensvar, C3G.sensvar, SHE.sensvar, BLE_Temp.sensvar))



sensmat.norm = normalize(sensmat)
pdf(file = 'sensitivity_summary_WUS.pdf', width = 12, height = 4)
par(mar = c(8,7,4,2))
image(sensmat.norm, col = yg, axes = FALSE, main = 'Parameter sensitivity in the Western US (normalised)')
axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels = c('NLE (5)', 'NLD (4)', 'C3G (6)', 'SHE (13)', 'BLE_Temp (3)'), las = 1)
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
dev.off()


image(t(test[d:1, ]), col = yg)


globalresponse = function(oaat.pred, n, d){
  # Calculate the effect of turning the parameterfrom lowest to highest
  out = rep(NA,d)
  for(i in 1:d){
    ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
    out[i] = (oaat.pred$mean[ix[n]] - oaat.pred$mean[ix[1]])
  }
  out
}

NLE.response = globalresponse(NLE.oaat, n = n, d = d)
NLD.response= globalresponse(NLD.oaat, n = n, d = d)
C3G.response = globalresponse(C3G.oaat, n = n, d = d)
SHE.response = globalresponse(SHE.oaat, n = n, d = d)
BLE_Temp.response = globalresponse(BLE_Temp.oaat, n = n, d = d)

responsemat = as.matrix(cbind(NLE.response, NLD.response, C3G.response, SHE.response, BLE_Temp.response))

responsemat[responsemat > 0.08] <- 0.08
responsemat[responsemat < -0.08] <- -0.08

pdf(file = 'response_summary_WUS.pdf', width = 12, height = 3.5)
par(mar = c(8,7,4,7))
image(responsemat, col = byr, axes = FALSE, main = 'Response to parameters in the Western US ')
axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels = c('NLE (5)', 'NLD (4)', 'C3G (6)', 'SHE (13)', 'BLE_Temp (3)'), las = 1)
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)

image.plot(responsemat, col = byr,legend.only = TRUE)

dev.off()





