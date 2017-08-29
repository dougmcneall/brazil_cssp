setwd("/Users/dougmcneall/Documents/work/R/brazil_cssp/analyse_u-ak-745")

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
ryb = brewer.pal(11, "RdYlBu")
byr = rev(ryb)
rb = brewer.pal(11, "RdBu")
br = rev(rb)


sensvar = function(oaat.pred, n, d){
  # Calculate variance as a global sensitivity meansure
  out = rep(NA,d)
  for(i in 1:d){
    ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
    out[i] = var(oaat.pred$mean[ix])
  }
  out
}

globalresponse = function(oaat.pred, n, d){
  # Calculate the effect of turning the parameterfrom lowest to highest
  out = rep(NA,d)
  for(i in 1:d){
    ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
    out[i] = (oaat.pred$mean[ix[n]] - oaat.pred$mean[ix[1]])
  }
  out
}

cutoff <- function(dat, zlim){
  # Function to cut off data at thresholds for plotting
  out <- dat
  out[dat < zlim[1]] <- zlim[1]
  out[dat > zlim[2]] <- zlim[2]
  out
}

globalsens <- function(lhs, response.list){
  # run globalresponse across a list of parameters
  d <- ncol(lhs)
  lhs.norm <- normalize(lhs)
  
  n = 21
  X.oaat = oaat.design(lhs.norm, n, med = TRUE)
  
  oaat.globalresponse = matrix(NA, nrow = d, ncol = length(response.list))
  for(i in 1: length(global_area_means.list)){
    
    try({
    fit = km(~., design = lhs.norm, response = c(response.list[i], recursive = TRUE) )
    oaat.pred = predict(fit, newdata = X.oaat, type = 'UK')
    oaat.globalresponse[,i] = globalresponse(oaat.pred, n = n, d = d)
    },
    silent = TRUE)
    
  }
  
  out <- oaat.globalresponse
  out
}

# load data
# All the data is in frac_area_mean

surftypes = c('BLD','BLE_Trop','BLE_Temp','NLD',
              'NLE','C3G','C3C','C3P','C4G','C4C',
              'C4P','SHD','SHE','Urban','Lake','Bare Soil',
              'Ice')

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

global.sensmat <- globalsens(lhs, global_area_means.list)
wus.sensmat <- globalsens(lhs, wus_area_means.list)
sam.sensmat <- globalsens(lhs, sam_area_means.list)

global.sensmat.co <- cutoff(global.sensmat, c(-0.09, 0.09))

pdf(file = 'response_summary_global.pdf', width = 12, height = 6)
par(mar = c(8,7,4,7))
image(global.sensmat.co, col = br, axes = FALSE)
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/16) , labels = surftypes, las = 1)

image.plot(global.sensmat.co, col = br,legend.only = TRUE)
dev.off()

wus.sensmat.co <- cutoff(wus.sensmat, c(-0.09, 0.09))

pdf(file = 'response_summary_wus.pdf', width = 12, height = 6)
par(mar = c(8,7,4,7))
image(wus.sensmat.co, col = br, axes = FALSE)
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/16) , labels = surftypes, las = 1)

image.plot(wus.sensmat.co, col = br,legend.only = TRUE)
dev.off()


sam.sensmat.co <- cutoff(sam.sensmat, c(-0.09, 0.09))
pdf(file = 'response_summary_sam.pdf', width = 12, height = 6)
par(mar = c(8,7,4,7))
image(sam.sensmat.co, col = br, axes = FALSE)
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/16) , labels = surftypes, las = 1)

image.plot(sam.sensmat.co, col = br,legend.only = TRUE)
dev.off()


pdf(file = 'response_summary_all.pdf', width = 12, height = 10)
par(mar = c(8,7,4,10), mfrow = c(3,1))

image(wus.sensmat.co, col = br, axes = FALSE, main = 'Western USA')
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/16) , labels = surftypes, las = 1)

image(sam.sensmat.co, col = br, axes = FALSE, main = 'South America')
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/16) , labels = surftypes, las = 1)

image(global.sensmat.co, col = br, axes = FALSE, main = 'Global')
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/16) , labels = surftypes, las = 1)
image.plot(wus.sensmat.co, col = br,legend.only = TRUE)

dev.off()

# Plot summaries of the sensitivity measures - here, summing the absolute
# values of the sensitivity across all land surface types.

wus.sens.summary <- apply(abs(wus.sensmat),1, sum, na.rm = TRUE)
wus.sens.ix <- order(wus.sens.summary, decreasing = TRUE)

pdf(width = 12, height = 5, file = 'sens_summary_wus.pdf')
par(mar = c(7,5,3,2))
plot(wus.sens.summary[wus.sens.ix], axes = FALSE, 
     xlab = '', ylab = 'Sensitivity summary', ylim = c(0,0.3))
axis(1, labels=(colnames(lhs))[wus.sens.ix], at = 1:length(wus.sens.summary), las = 3, cex.axis = 0.8 )
axis(2, las = 1)
dev.off()

sam.sens.summary <- apply(abs(sam.sensmat),1, sum, na.rm = TRUE)
sam.sens.ix <- order(sam.sens.summary, decreasing = TRUE)

pdf(width = 12, height = 5, file = 'sens_summary_sam.pdf')
par(mar = c(7,5,3,2))
plot(sam.sens.summary[sam.sens.ix], axes = FALSE, 
     xlab = '', ylab = 'Sensitivity summary', ylim = c(0,1))
axis(1, labels=(colnames(lhs))[sam.sens.ix], at = 1:length(sam.sens.summary), las = 3, cex.axis = 0.8 )
axis(2, las = 1)
dev.off()

global.sens.summary <- apply(abs(global.sensmat),1, sum, na.rm = TRUE)
global.sens.ix <- order(global.sens.summary, decreasing = TRUE)

pdf(width = 12, height = 5, file = 'global_summary_sam.pdf')
par(mar = c(7,5,3,2))
plot(global.sens.summary[global.sens.ix], axes = FALSE, 
     xlab = '', ylab = 'Sensitivity summary', ylim = c(0,0.5))
axis(1, labels=(colnames(lhs))[global.sens.ix], at = 1:length(global.sens.summary), las = 3, cex.axis = 0.8 )
axis(2, las = 1)
dev.off()





