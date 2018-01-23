# runoff_u-ao732.R

source('../per_pft.R')

#paired = brewer.pal(11,'Paired')
#linecols = c('black', paired[1], paired[2], paired[5])

lhs = read.table('data/lhs_u-ao732.txt', header = TRUE)
X = normalize(lhs)
colnames(X) = colnames(lhs)

amazon_runoff_historical = read.table('data/amazon_historical_runoff.txt',
                              na.strings = '-9.990000000000000000e+02',
                              header = FALSE)[-1, ]/1e8
                                        
anomalizeTSmatrix = function(x, ix){
  subx = x[ ,ix]
  sweepstats = apply(subx, 1, FUN=mean)
  anom = sweep(x, 1, sweepstats, FUN = '-')
  anom
}

years = 1861:2014
pdf(file = 'graphics/amazon_runoff_timeseries.pdf')
layout(matrix(c(1,1,1,1,2,2), 3, 2, byrow = TRUE))
 ## show the regions that have been allocated to each plot
#layout.show(2)
par(las = 1, mar = c(2,4,2,1))
matplot(t(amazon_runoff_historical),
        type = 'l', col = linecols, lty = 1,
        axes = FALSE,
        ylab = 'Amazon Runoff (kg/s)/1E8')
axis(1, at = seq(from = 1860, to = 2020, by = 20))
axis(2)
par(las = 1, mar = c(5,4,1,1))
amazon_runoff_anom = anomalizeTSmatrix(amazon_runoff_historical, 1:30)
matplot(years, t(amazon_runoff_anom), type = 'l', col = linecols, lty = 1,
        axes = FALSE,
        ylab = 'Amazon Runoff Anomaly (kg/s)/1E8')
axis(1, at = seq(from = 1860, to = 2020, by = 20))
axis(2)
dev.off()


amazon_runoff_init = amazon_runoff_historical[,1]

# A significant fraction of runs are close to zero runoff in the Amazon,
# which is handy
pdf(file='graphics/amazon_runoff_init_hist.pdf')
par(fg = 'white', las = 1)
hist(amazon_runoff_init, breaks = 20, xlab = 'Amazon Runoff in 1860 (kg/s)/1E8',
     col = 'grey', main = '')
axis(1, col = 'black')
axis(2, col = 'black')
dev.off()





# Where are the parameters that produce near-zero runoff
# in the amazon?
X.bt = X[amazon_runoff_init<1, ]

# It's really hard to see any patterns in a parallel coordinates plot:
parcoord(X.bt)

# It's easier to see patterns in a pairs plot
pdf(file = 'graphics/amazon_runoff_init_pairs.pdf', width = 10, height = 10)
pairs(X.bt, gap = 0, pch = 20, cex = 0.6, cex.labels = 0.5, lower.panel = NULL)
dev.off()


pdf('graphics/marginal_scatter_amazon_runoff_init.pdf', width = 15, height = 10)
par(mfrow =c(4,8), mar = c(3,1,1,1))
for(i in 1:ncol(lhs)){
  plot(lhs[, i], amazon_runoff_init,
       col = 'black',
       pch = 20,
       xlab = '',
       ylab = '',
       axes = FALSE)
  mtext(1, text = colnames(lhs[i]), line = 1)
}
dev.off()

# Build a two-step Gaussian Process emulator of the 
# inital runoff value
amazon_runoff_init.twoStep = twoStep(X = X, y = sqrt(amazon_runoff_init))

# A stepwise algorithm finds that only b_wl_io and f0_io are important
# for building the emulator

plot(amazon_runoff_init.twoStep$emulator)

twostep.pred = leaveOneOut.km(amazon_runoff_init.twoStep$emulator, 
                              type='UK', trend.reestim=TRUE)

plot(amazon_runoff_init,twostep.pred$mean^2)
abline(0,1)

# Here's what the steplm leaves in:
terms(amazon_runoff_init.twoStep$steplm)
attr(terms(amazon_runoff_init.twoStep$steplm), 'term.labels')

n = 21
d = ncol(X)
X.oaat = oaat.design(X, n, med = TRUE)


amazon_runoff_init.oaat = predict(amazon_runoff_init.twoStep$emulator, 
                                  newdata = X.oaat, type = 'UK')
y.oaat = amazon_runoff_init.oaat$mean^2

ylim = c(0,4)
pdf(file = 'graphics/amazon_runoff_init_oaat.pdf', width = 9, height = 6)
par(mfrow = c(4,8), mar = c(2,3,2,0.3), oma = c(0.5,0.5, 0.5, 0.5))

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], y.oaat[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE, col = 'black',
       main = '',
       xlab = '')
  axis(1, col = 'grey', col.axis = 'grey', las = 1)
  axis(2, col = 'grey', col.axis = 'grey', las = 1)
  mtext(3, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
dev.off()

X.unif = samp.unif(10000, mins = rep(0, d), maxes = rep(1,d))
colnames(X.unif) <- colnames(lhs)

y.unif = predict(amazon_runoff_init.twoStep$emulator, 
                 newdata = X.unif, type = 'UK')$mean^2

ix.kept = which(y.unif^2 < 1)

X.kept = X.unif[ix.kept, ]

pdf(file = 'graphics/pairs_dens_amazon_runoff_init.pdf', width = 8, height = 8)
par(oma = c(0,0,0,3))

# Emulate all input space and keep only those inputs which match a 
# criteria (such as having absolute error below a threshold)

test = pairs(X.kept,
             gap = 0, lower.panel = NULL, xlim = c(0,1), ylim = c(0,1),
             panel = dfunc.up)

image.plot(legend.only = TRUE,
           zlim = c(0,1),
           col = rb,
           legend.args = list(text = 'Density of model runs matching the criteria', side = 3, line = 1),
           horizontal = TRUE
)
dev.off()

# How about the change over time?

hist_runoff = apply(amazon_runoff_historical[,1:30], 1, mean)

modern_runoff = apply(amazon_runoff_historical[,125:154], 1, mean)

runoff_change = modern_runoff - hist_runoff

hist(runoff_change)

pdf('graphics/marginal_scatter_amazon_runoff_change.pdf', width = 15, height = 10)
par(mfrow =c(4,8), mar = c(3,1,1,1))
for(i in 1:ncol(lhs)){
  plot(lhs[, i], runoff_change,
       col = 'black',
       pch = 20,
       xlab = '',
       ylab = '',
       axes = FALSE)
  mtext(1, text = colnames(lhs[i]), line = 1)
}
dev.off()

runoff_change.twoStep = twoStep(X = X, y = (runoff_change))

# A stepwise algorithm finds that only b_wl_io and f0_io are important
# for building the emulator

runoff_changetwostep.pred = leaveOneOut.km(runoff_change.twoStep$emulator, 
                              type='UK', trend.reestim=TRUE)

amazon_runoff_change.oaat = predict(runoff_change.twoStep$emulator, 
                                  newdata = X.oaat, type = 'UK')
y.oaat = amazon_runoff_change.oaat$mean

ylim = c(-0.1, 0.2)
pdf(file = 'graphics/amazon_runoff_change_oaat.pdf', width = 9, height = 6)
par(mfrow = c(4,8), mar = c(2,3,2,0.3), oma = c(0.5,0.5, 0.5, 0.5))

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], y.oaat[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE, col = 'black',
       main = '',
       xlab = '')
  axis(1, col = 'grey', col.axis = 'grey', las = 1)
  axis(2, col = 'grey', col.axis = 'grey', las = 1)
  mtext(3, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
dev.off()


# How do the sensitivities change if you rule out the
# "Implausible" space?

## throw away nans and zeros
ix.plausible = which(amazon_runoff_init>1)
X.plausible = X[ix.plausible, ]
y.plausible =  runoff_change[ix.plausible]
runoff_change_plausible.twoStep = twoStep(X = X.plausible, y = y.plausible)
amazon_runoff_change_plausible.oaat = predict(runoff_change_plausible.twoStep$emulator, 
                                    newdata = X.oaat, type = 'UK')
y.plausible.oaat = amazon_runoff_change_plausible.oaat$mean
# Another way to do this is to reject the input space that we know is no good
# - in this case, b_wl_io below ~ 0.35 and f0_io above 0.8

ix.mankeep = which(X[, which( colnames(X) == "b_wl_io")] > 0.35 & X[, which( colnames(X) == "f0_io")] <0.8)
X.mankeep = X[ ix.mankeep, ] 
y.mankeep=  runoff_change[ix.mankeep]
runoff_change_mankeep.twoStep = twoStep(X = X.mankeep, y = y.mankeep)
amazon_runoff_change_mankeep.oaat = predict(runoff_change_mankeep.twoStep$emulator, 
                                            newdata = X.oaat, type = 'UK')

y.mankeep.oaat = amazon_runoff_change_mankeep.oaat$mean
ylim = range(y.oaat)

pdf(file = 'graphics/amazon_runoff_change_plausible_oaat.pdf', width = 9, height = 6)
par(mfrow = c(4,8), mar = c(2,3,2,0.3), oma = c(0.5,0.5, 0.5, 0.5))

for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], y.plausible.oaat[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE, col = 'black',
       main = '',
       xlab = '')
  lines(X.oaat[ix,i], y.mankeep.oaat[ix], col = 'red')
  axis(1, col = 'grey', col.axis = 'grey', las = 1)
  axis(2, col = 'grey', col.axis = 'grey', las = 1)
  mtext(3, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
dev.off()

sensvar.plausible = sensvar(amazon_runoff_change_plausible.oaat, n = 21, d = d)
sensvar.mankeep = sensvar(amazon_runoff_change_mankeep.oaat, n = 21, d = d)

plot(1:d,sensvar.plausible )
points(1:d, sensvar.mankeep, col = 'red')


pdf(file = 'graphics/amazon_runoff_change_mankeep_oaat.pdf', width = 9, height = 6)
par(mfrow = c(4,8), mar = c(2,3,2,0.3), oma = c(0.5,0.5, 0.5, 0.5))
for(i in 1:d){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], y.oaat[ix],
       type = 'l',
       ylab= '', ylim = ylim, axes = FALSE, col = 'black',
       main = '',
       xlab = '')
  axis(1, col = 'grey', col.axis = 'grey', las = 1)
  axis(2, col = 'grey', col.axis = 'grey', las = 1)
  mtext(3, text = colnames(lhs)[i], line = 0.2, cex = 0.7)
}
dev.off()




###

oaat.twoStep = function(design, n, med = TRUE, hold = NULL, emulator){
  # emulator is output from twoStep
  
  X.oaat = oaat.design(emulator$x, n, med = TRUE)
  y.oaat = predict(emulator$emulator, newdata = X.oaat)
  
  return(list(X.oaat = X.oaat, y.oaat = y.ooat))
  
  
}

plot.oaat(X.oaat, y.oaat, rows, cols, ...){
  
  ylim = range(y.oaat)
  par(mfrow = c(rows,cols), mar = c(2,3,2,0.3), oma = c(0.5,0.5, 0.5, 0.5))
  
  for(i in 1:d){
    ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
    plot(X.oaat[ix,i], y.oaat[ix],
         type = 'l',
         ylab= '', ylim = ylim, axes = FALSE, col = 'black',
         main = '',
         xlab = '')
    axis(1, col = 'grey', col.axis = 'grey', las = 1)
    axis(2, col = 'grey', col.axis = 'grey', las = 1)
    mtext(3, text = colnames(X.oaat)[i], line = 0.2, cex = 0.7)
  }
  dev.off()
  
  
}




