# per_pft.R
# A set of functions for analysing data from a 'per-pft' ensemble of JULES

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

globalresponse.lm = function(oaat.pred, n, d){
  # Calculate the effect of turning the parameterfrom lowest to highest
  out = rep(NA,d)
  for(i in 1:d){
    ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
    out[i] = (oaat.pred[ix[n]] - oaat.pred[ix[1]])
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
  colnames(X.oaat) = colnames(lhs)
  
  oaat.globalresponse = matrix(NA, nrow = d, ncol = length(response.list))
  for(i in 1: length(global_area_means.list)){
    print(i)
    try({
      y = c(response.list[i], recursive = TRUE)
      
      dat = data.frame(y=y, x=lhs.norm)
      colnames(dat) = c('y', colnames(lhs))
      
      initfit = lm(y ~ ., data = dat)
      stepfit = step(initfit, direction="both", k=log(length(y)), trace=TRUE)
      
      n = 21
      X.oaat = oaat.design(lhs.norm, n, med = TRUE)
      colnames(X.oaat) = colnames(lhs)
      oaat.pred = predict(stepfit, newdata = data.frame(X.oaat))
      oaat.globalresponse[,i] = globalresponse.lm(oaat.pred, n = n, d = d)
    },
    silent = TRUE)
    
    # try({
    #   fit = km(~., design = lhs.norm, response = c(response.list[i], recursive = TRUE) )
    #   oaat.pred = predict(fit, newdata = X.oaat, type = 'UK')
    #   oaat.globalresponse[,i] = globalresponse(oaat.pred, n = n, d = d)
    #   },
    #   silent = TRUE)
    # }
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

pfts = c('BLD','BLE_Trop','BLE_Temp','NLD',
         'NLE','C3G','C3C','C3P','C4G','C4C',
         'C4P','SHD','SHE') # like surftypes, but without the other surfaces

