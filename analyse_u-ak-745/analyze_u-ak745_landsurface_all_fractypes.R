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
ryb = brewer.pal(9, "RdYlBu")
byr = rev(ryb)


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


# All the data is in frac_area_mean

# load data

# Try land surface types in rows of a data frame

filelist.global <- paste0('frac_area_mean/global_area_mean_PFT',0:16, '.txt')

c.df <- function(fn){
    out = c(read.table(fn, header = FALSE, skip = 1), recursive = TRUE)
  out
}

# creates a list of (PFTs) long, each element of which has (ensemble members) data points
global_area_means.list <- lapply(filelist.global, c.df)

# ------------------------------------------------------
# Build emulators
# ------------------------------------------------------
lhs <- read.table('lhs_u-ak745.txt', header = TRUE)
d <- ncol(lhs)
cn <- colnames(lhs)
lhs.norm <- normalize(lhs)

n = 21
X.oaat = oaat.design(lhs.norm, n, med = TRUE)

oaat.sensvar = matrix(NA, nrow = d, ncol = 17)
oaat.globalresponse = matrix(NA, nrow = d, ncol = 17)

for(i in 1: length(global_area_means.list)){
  
  fit = km(~., design = lhs.norm, response = c(global_area_means.list[i], recursive = TRUE) )
  oaat.pred = predict(fit, newdata = X.oaat, type = 'UK')
  oaat.sensvar[,i] = sensvar(oaat.pred, n = n, d = d)
  oaat.globalresponse[,i] = globalresponse(oaat.pred, n = n, d = d)
  
}

surftypes = c('BLD','BLE_Trop','BLE_Temp','NLD','NLE','C3G','C3C','C3P','C4G','C4C','C4P','SHD','SHE','Urban','Lake','Bare Soil','Ice')
image(oaat.globalresponse, col = byr, axes = FALSE)
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/16) , labels = surftypes, las = 1)

# build emulators
# test and summarise sensitivity
# global land surface





