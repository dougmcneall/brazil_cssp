# screen_test.R
# testing methods for screening input parameters
# in the test ensemble of JULES

## Try and build an emulator
library(devtools)

install_github(repo = "dougmcneall/hde")

library(hde)
library(RColorBrewer)
library(fields)
library(DiceEval)
library(MASS)
library(DiceKriging)

source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/imptools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/vistools.R")

yg = brewer.pal(9, "YlGn")
ryb = brewer.pal(9, "RdYlBu")
byr = rev(ryb)

#read data
frac = read.table('forest_frac.txt', header = FALSE)
lhs = read.table('lhs_u-ak745.txt', header = TRUE)

forest_frac = frac[, 2]

d = ncol(lhs)
n = nrow(lhs)
cn = colnames(lhs)

lhs.norm = normalize(lhs)

fit = km(~., design = lhs.norm, response = forest_frac)

newdata = matrix(rep(0.5,d), nrow = 1)
colnames(newdata) = cn

pred = predict(fit, newdata = newdata, type = 'UK')


# How does the emulator error change as we enlarge the
# size of the design? Where are we on the curve? Use a 
# true-holdout method.


onestep = function(x, y, reps, startsize){

  n = nrow(x)
  d = ncol(x)
  cn = colnames(x)
  err_matrix = matrix(NA, ncol = reps, nrow = n)
  
  for(i in startsize:(n-1)){
    
    for(j in 1:reps){
  
    ix = sample(1:n, size = i, replace = FALSE)
    x_trunc = x[head(ix, -1), ]
    y_trunc = y[head(ix, -1)]
    
    x_target = matrix(x[tail(ix, 1), ], nrow = 1)
    colnames(x_target) = cn
    
    y_target = y[tail(ix, 1)]
    fit_trunc = km(~., design = x_trunc, response = y_trunc)
    
    pred = predict(fit_trunc, newdata = x_target, type = 'UK')
    err_matrix[(i+1),j] = pred$mean - y_target
    
    }
  }
  err_matrix
}

osa = onestep(x = lhs.norm, y = forest_frac, reps = 33, startsize = 77)

osa_mean = apply(abs(osa), 1, mean)
plot(osa_mean)


library(glmnet)

glmnetfit = glmnet(x = lhs.norm, y = forest_frac)

cvfit = cv.glmnet(x = lhs.norm, y = forest_frac)

