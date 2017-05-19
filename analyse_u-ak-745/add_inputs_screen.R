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

par(las = 1)
plot(mae.rand, xlab = 'No. included inputs', ylab = 'Mean absolute error',
     pch = 19, col = 'black',
     ylim = c(0, max(mae.rand, na.rm = TRUE))
      )
points(mae, col = 'red', pch = 19)
legend('topright', c('Random', 'Ordered'), pch = 19, col = c('black', 'red'))





