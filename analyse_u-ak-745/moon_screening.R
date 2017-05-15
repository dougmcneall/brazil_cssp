# moon_screening.R
# Testing screening algorithms with the functions from Moon (2010, 2012)

library(DiceKriging)
library(hde)
library(DiceEval)
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

xx = maximinLHS(500, 20)
colnames(xx) <- paste0('x',1:d)

y = apply(xx, 1, FUN = moon10hdc1)
d = ncol(xx)

#dev.new(width = 10, height = 5)
#par(mfrow = c(3,7), mar = c(2,2,2,2))

#for(i in 1:d){
#  plot(xx[, i ], y, axes = FALSE, xlab = '', ylab = '')
#}


# Build a GP emulator of the output


fit = km(~., design = xx, response = y)

n = 21
X.oaat= oaat.design(xx, n, med = TRUE)
colnames(X.oaat) <- colnames(xx)
dk.oaat = predict(fit, newdata = X.oaat, type = 'UK')
y.oaat = apply(X.oaat, 1, FUN = moon10hdc1)

# Run the oaat output through the function so that we can compare the
# emulator


# plot the emulated oat output
par(mfrow = c(5,4), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = range(y.oaat)

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  
  
  plot(X.oaat[ix,i], y.oaat[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE)
  lines(X.oaat[ix,i], dk.oaat$mean[ix], col = 'red')
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(xx)[i], line = 0.2, cex = 0.7)
}



sens.var = rep(NA,d)
for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  # sens.range[i] = diff(range(y.oaat$mean[ix]))
  sens.var[i] = var(dk.oaat$mean[ix])
}


active.ix = c(1,7,12,18,19)
colvec = rep(1,d)
colvec[active.ix] = 2
plot(sens.var, col = colvec, pch = 19)


subvar = function(x, y, subsize, cutoff, reps){
  
  # choose subsets of the full input matrix,
  # and summarise the sensitivity as we go.
  
  # pick a subset of inputs
  # build a gp
  # rank the effect size
  # record the rank, keep the top ones and iterate
  
  d = ncol(x)
  subranks = matrix(data = NA, nrow = reps, ncol = cutoff)
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

test = subvar(x = xx, y = y, subsize = 5, cutoff = 2, reps = 200)

varsum = rep(NA, d)
for(i in 1:d){
  
  varsum[i] = sum( test$subvars[test$subrank == i])
}

varsort = sort(varsum, decreasing = TRUE, index.return = TRUE)



