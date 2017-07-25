# add_inputs_screen.R
# What is the impact on the emulator of building it with
# a reduced subset of inputs? We'll try and find the best
# set to build the emulator with.

# Result: We compare the reduction in MAE as we add ordered (by importance) inputs
# with randomly chosen inputs. The ordered inputs have lower absoliute error, and reduces
# to an error floor more rapidly, after adding around 40-50 inputs.
# This is a larger number of inputs to hit bottom than I expected!

library(DiceKriging)
library(hde)
library(MASS)
library(sensitivity)
library(lhs)

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


# Use a variance based one-at-a-time metric
# to rank the inputs in order of importance.
n = 21
X.oaat= oaat.design(lhs.norm, n, med = TRUE)
y.oaat = predict(fit, newdata = X.oaat, type = 'UK')

sens.var = rep(NA,d)

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  sens.var[i] = var(y.oaat$mean[ix])
}

#ix.sorted = sort(sens.range, decreasing = TRUE, index.return = TRUE)$ix
#sens.range[ix.sorted]

ix.sorted = sort(sens.var, decreasing = TRUE, index.return = TRUE)$ix
plot(sens.var[ix.sorted])

# full.order.var = cn[ix.sorted]

# How does the predictive (or model fit) error change
# as we add input variables to the emulator?

mae = rep(NA, d)
for(i in 2:d){
  
  lhs.trunc = lhs.norm[, ix.sorted[1:i]]
  fit = km(~., design = lhs.trunc, response = forest_frac)
  loo =  leaveOneOut.km(fit, type = 'UK')
  err = loo$mean - forest_frac
  mae[i] = mean(abs(err))
}

mae.rand = rep(NA, d)
for(i in 2:d){
  
  lhs.trunc = lhs.norm[, sample(1:d, i, replace = FALSE)]
  fit = km(~., design = lhs.trunc, response = forest_frac)
  loo =  leaveOneOut.km(fit, type = 'UK')
  err = loo$mean - forest_frac
  mae.rand[i] = mean(abs(err))
}


dev.new(width = 4, height = 6)
par(las = 1)
plot(2:74, mae.rand[2:74], xlab = 'No. included inputs', ylab = 'Mean absolute error',
     pch = 19, col = 'black',
     ylim = c(0, max(mae.rand, na.rm = TRUE))
      )
points(mae[2:74], col = 'red', pch = 19)
legend('topright', c('Random', 'Ordered'), pch = 19, col = c('black', 'red'))

# which model has the minimum absolute error?

wmin = which.min(mae)

# repeat the whole analysis

x.trunc1 = lhs.norm[, ix.sorted[1:wmin]]

fit = km(~., design = x.trunc1, response = forest_frac)
n = 21
X.oaat= oaat.design(x.trunc1, n, med = TRUE)
y.oaat = predict(fit, newdata = X.oaat, type = 'UK')

sens.var = rep(NA,d)

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  sens.var[i] = var(y.oaat$mean[ix])
}

ix.sorted = sort(sens.var, decreasing = TRUE, index.return = TRUE)$ix
plot(sens.var[ix.sorted])

mae = rep(NA, ncol(x.trunc1))
for(i in 2:d){
  
  lhs.trunc = x.trunc1[, ix.sorted[1:i]]
  fit = km(~., design = lhs.trunc, response = forest_frac)
  loo =  leaveOneOut.km(fit, type = 'UK')
  err = loo$mean - forest_frac
  mae[i] = mean(abs(err))
}

# Now we are down to ~ 50 that make a difference to the emulator fit!

x.trunc2 = x.trunc1[ ,ix.sorted[1:50]]

fit = km(~., design = x.trunc2, response = forest_frac)
n = 21
X.oaat= oaat.design(x.trunc2, n, med = TRUE)
y.oaat = predict(fit, newdata = X.oaat, type = 'UK')


par(mfrow = c(5,10), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = c(-0.12, 0.4)

for(i in 1:ncol(x.trunc2)){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], y.oaat$mean[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE)
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(x.trunc2)[i], line = 0.2, cex = 0.7)
}

sens.var = rep(NA,ncol(x.trunc2))

for(i in 1:ncol(x.trunc2)){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  sens.var[i] = var(y.oaat$mean[ix])
}

plot(sens.var)







