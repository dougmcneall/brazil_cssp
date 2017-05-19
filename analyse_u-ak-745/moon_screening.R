# moon_screening.R
# Testing screening algorithms with the functions from Moon (2010, 2012)

library(DiceKriging)
library(hde)
library(MASS)
library(sensitivity)
library(lhs)

source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/imptools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/vistools.R")

moon10hdc1 <- function(xx)
{
  #########################################################################
  #
  # MOON (2010) HIGH-DIMENSIONALITY FUNCTION, C-1
  # This function is a modification of the function moon10hd.r.
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2, ..., x20)
  #
  ##########################################################################
  
  x1  <- xx[1]
  x7  <- xx[7]
  x12 <- xx[12]
  x18 <- xx[18]
  x19 <- xx[19]
  
  term1 <- -19.71*x1*x18 + 23.72*x1*x19
  term2 <- -13.34*x19^2 + 28.99*x7*x12
  
  y <- term1 + term2
  return(y)
}

# Generate a MMLHS sample from the test function.
# Include a dummy variable
d = 50
xx = maximinLHS(100, d)
xt = xx[, 1:(d-1)]
y = apply(xt, 1, FUN = moon10hdc1)
colnames(xx) <- paste0('x',1:d)

# A plot of each input in turn, plotted against the output
#dev.new(width = 10, height = 5)
#par(mfrow = c(3,7), mar = c(2,2,2,2))

#for(i in 1:d){
#  plot(xx[, i ], y, axes = FALSE, xlab = '', ylab = '')
#}


# # Build a GP emulator of the output using DiceKriging
 fit = km(~., design = xx, response = y)
# 
# # Have a look at the one-at-a-time response, comparing the
# # true function against the emulator
 n = 21
 X.oaat= oaat.design(xx, n, med = TRUE)
 colnames(X.oaat) <- colnames(xx)
 dk.oaat = predict(fit, newdata = X.oaat, type = 'UK')
 y.oaat = apply(X.oaat, 1, FUN = moon10hdc1)
# 
# 
# # plot the emulated oat output and true function
par(mfrow = c(10,5), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = range(y.oaat)
# 
for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], y.oaat[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE)
  lines(X.oaat[ix,i], dk.oaat$mean[ix], col = 'red')
  mtext(1, text = colnames(xx)[i], line = 0.2, cex = 0.7)
}

# summarise the oaat sensitivity by measuring the variance of the
# output across each input (do this for the true function too)
sens.var = rep(NA,d)
for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  sens.var[i] = var(dk.oaat$mean[ix])
}

# # do this for the true function too
 active.ix = c(1,7,12,18,19)
 dummy.ix = d
 colvec = rep(1,d)
 colvec[active.ix] = 2
 colvec[dummy.ix] = 4
 dev.new()
 plot(sens.var, col = colvec, pch = 19)


subvar = function(x, y, subsize, cutoff, reps){
  
  # choose subsets of the full input matrix,
  # and summarise the sensitivity as we go.
  
  # pick a subset of inputs
  # build a gp
  # rank the effect size
  # record the rank, keep the top ones and iterate
  
  d = ncol(x)
  # subranks contains the index of the input which has the
  # largest effect (measured by variance of output)
  subranks = matrix(data = NA, nrow = reps, ncol = cutoff)
  
  # subvars is the measure of the effect size for the [cutoff] inputs
  # with the largest effect size
  subvars  = matrix(data = NA, nrow = reps, ncol = cutoff)
  
  for(i in 1:reps){
    
    sub.ix = sample(1:d, size = subsize)
    x.trunc = x[, sub.ix]
    
    fit = km(~., design = x.trunc, response = y)
    
    n = 21
    x.oaat= oaat.design(x.trunc, n, med = TRUE)
    y.oaat = predict(fit, newdata = x.oaat, type = 'UK')
    
    sens.var = rep(NA,subsize)
    
    for(j in 1:subsize){
      ix = seq(from = ((j*n) - (n-1)), to =  (j*n), by = 1)
      sens.var[j] = var(y.oaat$mean[ix])
    }
    ix.sorted = sort(sens.var, decreasing = TRUE, index.return = TRUE)$ix
    subranks[i, ] = (sub.ix[ix.sorted])[1:cutoff]
    subvars[i, ] =  (sens.var[ix.sorted])[1:cutoff]
  }
  
  return(list(subranks = subranks, subvars = subvars))
  
}

test = subvar(x = xx, y = y, subsize = 10, cutoff = 3, reps = 100)

varsum = rep(NA, d)
varmean = rep(NA, d)
varvar = rep(NA, d)
varsamp = rep(NA, d)
for(i in 1:d){
  
  varsum[i] = sum(test$subvars[test$subrank == i], na.rm = TRUE)
  varmean[i] = mean(test$subvars[test$subrank == i], na.rm = TRUE)
  
  varvar[i] = var( test$subvars[test$subrank == i])
}

plot(varsum, col = colvec, pch = 19)

# I think we ought to check against FAST99. How well does it do?
library(sensitivity)
xfast = fast99(model = NULL, factors = colnames(xx), n = 5000,
               q = "qunif", q.arg = list(min = 0, max = 1))
fast.pred = predict(fit, newdata = xfast$X, type = 'UK')
fast <- tell(xfast, fast.pred$mean)
#pdf(file = 'fast.pdf', width = 12, height = 6)
par(las = 2, mar = c(10,4,2,1))
plot(fast)
#dev.off()

fast.summary <- print(fast)



r = 10000
varsum.samp = rep(NA, r)
# resample the sum
# This doesn't seem to adequately represent the range
# that the (non known) zero coefficients have.
for(i in 1:r){
  varsum.samp[i] = sum(sample(test$subvars[test$subrank == 21], replace = TRUE), na.rm = TRUE)
}

vsrange = range(varsum.samp)
plot(varsum, col = colvec, pch = 19)
segments(x0 = 21, y0 = vsrange[1], x1 = 21, y1 = vsrange[2], col = 4)

#plot(varvar, col = colvec, pch = 19)


varsort = sort(varsum, decreasing = TRUE, index.return = TRUE)


xx.trunc = xx[, varsort$ix[1:10]] 

test2  = subvar(x = xx.trunc, y = y, subsize = 5, cutoff = 2, reps = 500)

varsum = rep(NA, ncol(xx.trunc))
for(i in 1:ncol(xx.trunc)){
  varsum[i] = sum(test2$subvars[test2$subrank == i], na.rm = TRUE)
}

# Ideas to follow up on
# Different forms of sensitivity measure (e.g. FAST, others?)
# More than oat measures (total effect?)
# Include a step that compares against dummy variables
# How does the number of e.g. samples, replications etc impact
# the quality of the screening?
# Increase the number of non-active inputs
# Try out different functions (see below)


morretal06 <- function(xx, k1=2)
{
  ##########################################################################
  #
  # MORRIS ET AL. (2006) FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., x30)
  # k1 = number of arguments with an effect (optional), with default value
  #      2  
  #
  ##########################################################################
  
  alpha <- sqrt(12) - 6*sqrt(0.1)*sqrt(k1-1)
  beta <- 12 * sqrt(0.1) * sqrt(k1-1)
  
  xi <- xx[1:k1]
  ximat <- matrix(rep(xi,times=k1), k1, k1, byrow=TRUE)
  ximatlow <- ximat
  ximatlow[!upper.tri(ximatlow)] <- 0
  
  inner <- rowSums(xi*ximatlow)
  outer <- sum(xi + beta*inner)
  
  y <- alpha * outer
  return(y)
}

n = 100
d = 300
# Generate a MMLHS sample from the test function
xx = maximinLHS(n, d)

active = 5
active.ix = 1:active

y = apply(xx, 1, FUN = morretal06, k1 = active)
colnames(xx) <- paste0('x',1:d)
test = subvar(x = xx, y = y, subsize = 30, cutoff = 2, reps = 100)

varsum = rep(NA, d)
for(i in 1:d){
  varsum[i] = sum( test$subvars[test$subrank == i])
}

colvec = rep(1, d)
colvec[active.ix] = 3

plot(varsum, col = colvec)

# Add a different (non active) dummy variable for each repliction, and keep inputs
# whose sensitivity falls outside some threshold of the distribution
# of dummy variables


