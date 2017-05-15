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

pred.gp = leaveOneOut.km(fit, type = 'UK')



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

# test = onestep(x = lhs.norm, y = forest_frac, reps = 5, startsize = 95)


# Use glmnet
library(glmnet)

glmnetfit = glmnet(x = lhs.norm, y = forest_frac)

cvfit = cv.glmnet(x = lhs.norm, y = forest_frac)
coef(cvfit, s = "lambda.min")

pred.glmnet = predict(cvfit, newx = lhs.norm, s = "lambda.min")


# OK, what happens to the prediction error of the Gaussian Process fit,
# if we just keep the non-zero elements, as identified by the glmfit

glmnet.coef = coef(cvfit, s = "lambda.min")
coef_ix = glmnet.coef@i

km.fit_trunc = km(~., design = lhs.norm[, coef_ix[coef_ix > 0]], response = forest_frac)
pred.gp_trunc = leaveOneOut.km(km.fit_trunc, type = 'UK')


plot(forest_frac, pred.gp$mean, xlim = c(0,0.5), ylim = c(0,0.5), pch = 19)
points(forest_frac, pred.glmnet, col = 'red', pch = 19)
points(forest_frac, pred.gp_trunc$mean, col = 'blue', pch = 19)
abline(0,1)

# The GP emulator pretty comprehensively outperforms the best glmnet
gp.err = pred.gp$mean - forest_frac
gp_trunc.err = pred.gp_trunc$mean - forest_frac
glmnet.err = pred.glmnet  - forest_frac

mean(abs(gp.err))
mean(abs(gp_trunc.err))

# It looks like screening inputs with glmnet and then 
# building a GP leads to approximately the same (higher)
# LOO error as just using the glmnet


n = 21
X.oaat= oaat.design(lhs.norm, n, med = TRUE)
y.oaat = predict(fit, newdata = X.oaat, type = 'UK')

sens.range = rep(NA,d)
sens.var = rep(NA,d)

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  sens.range[i] = diff(range(y.oaat$mean[ix]))
  sens.var[i] = var(y.oaat$mean[ix])
}

ix.sorted = sort(sens.range, decreasing = TRUE, index.return = TRUE)$ix
sens.range[ix.sorted]

ix.sorted = sort(sens.var, decreasing = TRUE, index.return = TRUE)$ix
plot(sens.var[ix.sorted])

full.order.var = cn[ix.sorted]

# It might be better to use variance as a measure of sensitivity.



km.fit_trunc = km(~., design = lhs.norm[, ix.sorted[1:50]], response = forest_frac)
pred.gp_trunc = leaveOneOut.km(km.fit_trunc, type = 'UK')
gp_trunc.err = pred.gp_trunc$mean - forest_frac
mean(abs(gp_trunc.err))

# ----------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------

# Here's a plan for screening with GP emulators post-hoc (i.e. after LH designed runs)

# 1) Sample manageable set of design columns
# 2) Build GP emulator and calculate effects
# 3) Keep any where effects are > 0
# 4) Iterate until some metric
# 5) Build GP and calculate effect sizes with final 

# It's basically subset regression?
# OK, how about randomly sampling and then keeping track of the ranks?

sens.range = rep(NA,d)

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  sens.range[i] = diff(range(y.oaat$mean[ix]))
}

ix.sorted = sort(sens.range, decreasing = TRUE, index.return = TRUE)$ix


subrank = function(x, y, subsize, cutoff, reps){
  
  # pick a subset of inputs
  # build a gp
  # rank the effect size
  # record the rank, keep the top ones and iterate
  
  d = ncol(x)
  subranks = matrix(data = NA, nrow = reps, ncol = cutoff)
  
  for(i in 1:reps){
    
    sub.ix = sample(1:d, size = subsize)
    x.trunc = x[, sub.ix]
    
    fit = km(~., design = x.trunc, response = y)
    
    n = 21
    x.oaat= oaat.design(x.trunc, n, med = TRUE)
    y.oaat = predict(fit, newdata = x.oaat, type = 'UK')
    
    sens.range = rep(NA,subsize)
    
    for(j in 1:subsize){
      ix = seq(from = ((j*n) - (n-1)), to =  (j*n), by = 1)
      sens.range[j] = diff(range(y.oaat$mean[ix]))
    }
    ix.sorted = sort(sens.range, decreasing = TRUE, index.return = TRUE)$ix
   subranks[i, ] = (sub.ix[ix.sorted])[1:cutoff]
  }
  
  subranks
  
}


test = subrank(x = lhs.norm, y = forest_frac, subsize = 15, cutoff = 5, reps = 300)


cn = colnames(lhs.norm)
tab = table(test)
tabsort = sort(tab, decreasing = TRUE)

ordered.list = cn[as.numeric(names(tabsort))]

# lets choose the top 12

chosen.cols = ordered.list[1:20]

x.chosen = lhs.norm[ , cn%in%chosen.cols]
fit.chosen = km(~., design = x.chosen, response = forest_frac)

n = 21
d = ncol(x.chosen)
X.oaat= oaat.design(x.chosen, n, med = TRUE)
y.oaat = predict(fit.chosen, newdata = X.oaat, type = 'UK')

#pdf(file = 'oaat.pdf', width = 8, height = 10)
par(mfrow = c(5,4), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))
ylim = c(-0.12, 0.4)

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], y.oaat$mean[ix],
       type = 'l',
       ylab= '', ylim = c(0,0.5), axes = FALSE)
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(x.chosen)[i], line = 0.2, cex = 0.7)
}
#dev.off()



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
  

test = subvar(x = lhs.norm, y = forest_frac, subsize = 30, cutoff = 5, reps = 300)


varsum = rep(NA, 74)
for(i in 1:74){
  
  varsum[i] = sum( test$subvars[test$subrank == i])
}

varsort = sort(varsum, decreasing = TRUE, index.return = TRUE)

# Compare with a FAST sens
library(sensitivity)
xfast = fast99(model = NULL, factors = colnames(lhs.norm), n = 5000,
            q = "qunif", q.arg = list(min = 0, max = 1))
fast.pred = predict(fit, newdata = xfast$X, type = 'UK')
fast <- tell(xfast, fast.pred$mean)
#pdf(file = 'fast.pdf', width = 12, height = 6)
par(las = 2, mar = c(10,4,2,1))
plot(fast)
#dev.off()

fast.summary <- print(fast)

fast.ix = sort(fast.summary[, 1], index.return = TRUE, decreasing = TRUE)$ix

dotchart(rev(fast.summary[fast.ix,1]), labels = cn[rev(fast.ix)])




dotchart(rev(sens.range[ix.sorted]), labels = rev(colnames(lhs)[ix.sorted]), cex = 0.7, pch = 19,
         xlab = 'One-at-a-time response range')

cn(fast.ix)[1:10]



plot(varsort$x)
points(fast.summary[varsort$ix, 1], col = 'red')


plot(fast.summary[,1], varsum)

# it looks like our sorting algorithm doesn't work very well :(
plot(varsum, sens.var)

# Interesting! It looks like the FAST algorithm produces a (first order) sensitivity measure
# that is linearly proportional to the variance across the inputs, and non-linearly
# proportional to the range across the inputs.
plot(fast.summary[,1], sens.range)
plot(fast.summary[,1], sens.var, col = 'red')

plot(sens.var, sens.range)










