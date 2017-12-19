# initial_analysis_u-ao732

source('../per_pft.R')

lhs = read.table('data/lhs_u-ao732.txt', header = TRUE)
X = normalize(lhs)
colnames(X) = colnames(lhs)


npp = c(read.table('data/global_sum_present_npp.txt',
                   na.strings = '-9.990000000000000000e+02',
                   header = FALSE), recursive = TRUE)[-1] / 1e12

nbp = c(read.table('data/global_sum_present_nbp.txt',
                   na.strings = '-9.990000000000000000e+02',
                   header = FALSE), recursive = TRUE)[-1] / 1e12

cv = c(read.table('global_sum_present_cv.txt',
             na.strings = '-9.990000000000000000e+02',
             header = FALSE), recursive = TRUE)[-1] / 1e12

cs = c(read.table('global_sum_present_cs.txt',
                  na.strings = '-9.990000000000000000e+02',
                  header = FALSE), recursive = TRUE)[-1] / 1e12

# I think we can get rid of quite a lot of parameter space
# using only the model runs

# Raw data histograms

pdf('graphics/raw_histograms.pdf', width = 6, height = 6)
par(mfrow = c(2,2), fg = 'white', las = 1)

hist(npp, col = 'tomato', axes = FALSE)
axis(1);axis(2)

hist(nbp, col = 'cadetblue3', axes = FALSE)
axis(1);axis(2)

hist(cv, col = 'slategray4', axes = FALSE)
axis(1);axis(2)

hist(cs, col = 'coral', axes = FALSE)
axis(1);axis(2)
dev.off()


# Now a marginal scatter plot

pdf('graphics/marginal_scatter_npp.pdf', width = 15, height = 10)
par(mfrow =c(4,8), mar = c(3,1,1,1))
for(i in 1:ncol(lhs)){
   plot(lhs[, i], npp,
       col = 'black',
       pch = 20,
       xlab = '',
       ylab = '',
       axes = FALSE)
  mtext(1, text = colnames(lhs[i]), line = 1)
}
dev.off()

pdf('graphics/marginal_scatter_nbp.pdf', width = 15, height = 10)
par(mfrow =c(4,8), mar = c(3,1,1,1))
for(i in 1:ncol(lhs)){
  plot(lhs[, i], nbp,
       col = 'black',
       pch = 20,
       xlab = '',
       ylab = '',
       ylim = c(-1,1),
       axes = FALSE)
  mtext(1, text = colnames(lhs[i]), line = 1)
}
dev.off()


pdf('graphics/marginal_scatter_cv.pdf', width = 15, height = 10)
par(mfrow =c(4,8), mar = c(3,1,1,1))
for(i in 1:ncol(lhs)){
  plot(lhs[, i], cv,
       col = 'black',
       pch = 20,
       xlab = '',
       ylab = '',
       axes = FALSE)
  mtext(1, text = colnames(lhs[i]), line = 1)
}
dev.off()

pdf('graphics/marginal_scatter_cs.pdf', width = 15, height = 10)
par(mfrow =c(4,8), mar = c(3,1,1,1))
for(i in 1:ncol(lhs)){
  plot(lhs[, i], cs,
       col = 'black',
       pch = 20,
       xlab = '',
       ylab = '',
       axes = FALSE)
  mtext(1, text = colnames(lhs[i]), line = 1)
}
dev.off()
  

gp.loo = true.loo()


# km doesn't like NAs - set to zero or exclude?
npp.naix  = which(is.na(npp))

npp.ex = npp[-npp.naix]
X.ex = X[-npp.naix, ]

fit.km.npp = km(~., design = X.ex, response = npp.ex)
fit.km.npp.sqrt = km(~., design = X.ex, response = sqrt(npp.ex))


loo.km.npp = leaveOneOut.km(fit.km.npp, type='UK', trend.reestim=TRUE)
loo.km.npp.sqrt = leaveOneOut.km(fit.km.npp.sqrt, type='UK', trend.reestim=TRUE)

rmse(npp.ex, loo.km.npp$mean)
rmse(npp.ex, loo.km.npp.sqrt$mean^2)

plot(npp.ex, loo.km.npp$mean, ylim = c(-20,120))
points(npp.ex, loo.km.npp.sqrt$mean^2, col = 'red')
abline(0,1)

# Actual hold-out gaussian process loo
km.true.loo = true.loo(X = X.ex, y = sqrt(npp.ex), type = 'gp')

# Rather nicely, the true.loo is very similar to the
# trend.reestim=TRUE method for LeaveOneOut.km (at least in mean)
plot(npp.ex, loo.km.npp$mean, ylim = c(-20,120))
points(npp.ex, loo.km.npp.sqrt$mean^2, col = 'red')
points(npp.ex, km.true.loo$mean^2, col = 'blue')
abline(0,1)



