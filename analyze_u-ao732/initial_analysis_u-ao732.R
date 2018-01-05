# initial_analysis_u-ao732

source('../per_pft.R')

paired = brewer.pal(11,'Paired')

linecols = c('black', paired[1], paired[2], paired[5])



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

# Andy Wiltshire suggests some initial values to help rule out parts 
# of parameter space.

# npp 35-80 GtC
# nbp > 0
# cVeg 300 - 800 GtC
# cSoil 750 - 3000 GtC

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
  

#gp.loo = true.loo()

exclude.na = function(y){
  na.ix  = which(is.na(y))
  ex = y[-na.ix]
  return(list(ex = ex, na.ix = na.ix))
}

npp.naex =  exclude.na(npp)
# km doesn't like NAs - set to zero or exclude?
npp.ex = npp.naex$ex
X.ex = X[-npp.naex$na.ix, ]

fit.km.npp = km(~., design = X.ex, response = npp.ex)
fit.km.npp.sqrt = km(~., design = X.ex, response = sqrt(npp.ex))


loo.km.npp = leaveOneOut.km(fit.km.npp, type='UK', trend.reestim=TRUE)
loo.km.npp.sqrt = leaveOneOut.km(fit.km.npp.sqrt, type='UK', trend.reestim=TRUE)

# Actual hold-out gaussian process loo
# Warning - SLOW
km.true.loo = true.loo(X = X.ex, y = sqrt(npp.ex), type = 'gp')

# Rather nicely, the true.loo is very similar to the
# trend.reestim=TRUE method for LeaveOneOut.km (at least in mean)
plot(npp.ex, loo.km.npp$mean, ylim = c(-20,120))
points(npp.ex, loo.km.npp.sqrt$mean^2, col = 'red')
points(npp.ex, km.true.loo$mean^2, col = 'blue')
abline(0,1)
legend('topleft', legend = c('NPP', 'sqrt(NPP)', 'true loo sqrt(NPP)'),
       col = c('black', 'blue', 'red'),
       pch = 21
       )

npp.twoStep = twoStep(X = X.ex, y = sqrt(npp.ex))

#abline(0,1)

twostep.pred = leaveOneOut.km(npp.twoStep$emulator, type='UK', trend.reestim=TRUE)

points(npp.ex, twostep.pred$mean^2, col = 'orange')

# npp GP rmse
rmse(npp.ex, loo.km.npp$mean)
# npp sqrt(GP) rmse
rmse(npp.ex, loo.km.npp.sqrt$mean^2)

# npp sqrt(twostep) rmse
rmse(npp.ex, twostep.pred$mean^2)

plot(loo.km.npp.sqrt$mean^2, km.true.loo$mean^2)
abline(0,1)


# Which inputs does the twoStep keep?

# One-at-a-time Sensitivity plots

# Here's what the steplm leaves in:
terms(npp.twoStep$steplm)
attr(terms(npp.twoStep$steplm), 'term.labels')

n = 21
d = ncol(X)
X.oaat = oaat.design(X, n, med = TRUE)


npp.oaat = predict(npp.twoStep$emulator, newdata = X.oaat, type = 'UK')
y.oaat = npp.oaat$mean^2
ylim = c(0, 70)
#pdf(file = 'npp_oaat.pdf', width = 8, height = 10)
dev.new()
par(mfrow = c(4,8), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], y.oaat[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE, col = 'black')
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
#dev.off()


nbp.constr.ix = which(nbp < -10)
nbp.constr = nbp[-nbp.constr.ix]
X.nbp.constr = X[-nbp.constr.ix, ]

nbp.naex =  exclude.na(nbp.constr)
nbp.ex = nbp.naex$ex
X.nbp.ex = X.nbp.constr[-nbp.naex$na.ix, ]

# Looks like the extra low input in nbp isn't being managed well.
# Exclude it too.
nbp.twoStep = twoStep(X = X.nbp.ex, nbp.ex)
nbp.twostep.pred = leaveOneOut.km(nbp.twoStep$emulator,
                                  type='UK', trend.reestim=TRUE)

attr(terms(nbp.twoStep$steplm), 'term.labels')

nbp.oaat = predict(nbp.twoStep$emulator, newdata = X.oaat, type = 'UK')
y.oaat = nbp.oaat$mean
ylim = c(-1, 1)
pdf(file = 'nbp_oaat.pdf', width = 8, height = 10)
dev.new()
par(mfrow = c(4,8), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], y.oaat[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE, col = 'black')
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
dev.off()


# vegetation carbon
cv.naex =  exclude.na(cv)
cv.ex = cv.naex$ex
X.cv.ex = X[-cv.naex$na.ix, ]

cv.twoStep = twoStep(X = X.cv.ex, sqrt(cv.ex))

cv.oaat = predict(cv.twoStep$emulator, newdata = X.oaat, type = 'UK')

y.oaat = cv.oaat$mean^2
ylim = range(y.oaat)
#pdf(file = 'nbp_oaat.pdf', width = 8, height = 10)
dev.new()
par(mfrow = c(4,8), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], y.oaat[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE, col = 'black')
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
#dev.off()

attr(terms(cv.twoStep$steplm), 'term.labels')


# Soil carbon
cs.naex =  exclude.na(cs)
cs.ex = cs.naex$ex
X.cs.ex = X[-cs.naex$na.ix, ]

cs.twoStep = twoStep(X = X.cs.ex, sqrt(cs.ex))

cs.oaat = predict(cs.twoStep$emulator, newdata = X.oaat, type = 'UK')

y.oaat = cs.oaat$mean^2
ylim = range(y.oaat)
#pdf(file = 'nbp_oaat.pdf', width = 8, height = 10)
dev.new()
par(mfrow = c(4,8), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 0.5, 0.5))

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], y.oaat[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE, col = 'black')
  #abline(v = 0, col = 'grey')
  #abline(h = 0.4, col = 'grey')
  mtext(1, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
#dev.off()
attr(terms(cs.twoStep$steplm), 'term.labels')

which(names(lhs) %in% attr(terms(cs.twoStep$steplm), 'term.labels'))


factors.matrix = matrix(0, ncol = d, nrow = 4)

factors.matrix[1, ] = names(lhs) %in% attr(terms(npp.twoStep$steplm), 'term.labels')
factors.matrix[2, ] = names(lhs) %in% attr(terms(nbp.twoStep$steplm), 'term.labels')
factors.matrix[3, ] = names(lhs) %in% attr(terms(cv.twoStep$steplm), 'term.labels')
factors.matrix[4, ] = names(lhs) %in% attr(terms(cs.twoStep$steplm), 'term.labels')

factors.sum = apply(factors.matrix,2, sum)

pdf(file = 'graphics/steplm_chosen_factors.pdf',width = 8, height = 5)
par(mfrow = c(2,1), las = 0, mar = c(1,5,3,1))
image(t(factors.matrix), col = yg, axes = FALSE)
axis(2, at = seq(from = 0, to = 1, by = 1/3),
     labels = c('npp', 'nbp','cv', 'cs'),
     las = 1)
par(las = 2, mar = c(8,5,1,1))
barplot(factors.sum, xaxs = 'i', names.arg=  colnames(lhs), las = 3, ylab = 'effect sum' )

dev.off()



# Create an "output matrix" so that we can plot everything on one sensitivity
# analysis plot.

all.out = normalize(
  matrix(
  cbind(npp.oaat$mean^2, nbp.oaat$mean, cv.oaat$mean^2, cs.oaat$mean^2),
  ncol = 4
                )
)


pdf(file = 'graphics/oaat_all.pdf', width = 7, height = 6)
par(mfrow = c(4,8), mar = c(1,0.3,0.3,0.3), oma = c(0.5,0.5, 3, 0.5))

# set up the plots
for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], all.out[ix, 1],
       ylab= '', ylim = c(0,1), type = 'n', axes = FALSE)
  for(j in 1:4){
    points(X.oaat[ix,i], all.out[ix, j],
      type = 'l',
      lwd = 2,
      col = linecols[j])
      }
  mtext(1, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
#par(xpd = TRUE)

reset <- function() {
  par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}
reset()
legend('top',
       legend = c('npp','nbp', 'cv', 'cs'), 
       col = linecols[1:4],
       lwd = 2,
       horiz = TRUE)
#legend("top", legend=c("A", "B"), fill=c("red", "blue"), ncol=2, bty="n")
dev.off()







