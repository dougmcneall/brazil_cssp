# famous_bias.R
# R script to illustrate fixing a bias in famous by using
# temperature and precipitation as inputs to the emulator.
# dougmcneall@gmail.com

# ------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------
library(DiceKriging)
library(RColorBrewer)
library(MASS)
library(fields)
library(parallel)

load('famous_forest_fraction.RData')
load('famous_agg.RData')

source('https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R')
source('https://raw.githubusercontent.com/dougmcneall/packages-git/master/imptools.R')
source('https://raw.githubusercontent.com/dougmcneall/packages-git/master/vistools.R')

# pallettes
rb <- brewer.pal(11, "RdBu")
ryg <- brewer.pal(11, "RdYlGn")
pbg <- brewer.pal(9, "PuBuGn")
bg <- brewer.pal(9, "BuGn")
yg <- brewer.pal(9, "YlGn")
byr <- rev(brewer.pal(11,'RdYlBu'))
br <- rev(rb)
blues <-  brewer.pal(9,'Blues')
rblues <-  rev(blues)

greens <-  brewer.pal(9,'Greens')
ygb <- brewer.pal(9, "YlGnBu")
brbg <- brewer.pal(11, "BrBG")
yob <- brewer.pal(9, "YlOrBr")
yor <- brewer.pal(9, "YlOrRd")
acc <- brewer.pal(8,'Paired')

col.amaz <- acc[1]
col.namerica <- acc[2]
col.seasia <- acc[3]
col.congo <- acc[4]
col.global <- acc[5]


pch.global <- 3
pch.amaz <- 1
pch.congo <- 2
pch.seasia <- 5
pch.namerica <- 4

dfunc.up <- function(x,y,...){
  require(MASS)
  require(RColorBrewer)
  
  br <- brewer.pal(9, 'Blues')
  # function for plotting 2d kernel density estimates in pairs() plot.
  kde <- kde2d(x,y)
  image(kde, col = br, add = TRUE)
  
}


X = famous_agg[, 2:8]
X.norm = normalize(X)
X.stan.norm <- normalize(matrix(X.standard, nrow = 1), wrt = X)
colnames(X.stan.norm) = colnames(X)

# Matrix of (repeated) input parameters, modelled temperatures,
# precipitation and dorest fractions.

X_AMAZON = famous_agg[ ,c(2,3,4,5,6,7,8,15,19)] # AMAZON
X_SEASIA = famous_agg[, c(2,3,4,5,6,7,8, 17,21)] # SEASIA
X_CONGO = famous_agg[, c(2,3,4,5,6,7,8,16,20)] # CONGO


colnames(X_AMAZON) = c(colnames(famous_agg)[2:8], 'MOD_TEMP', 'MOD_PRECIP')
colnames(X_SEASIA) = c(colnames(famous_agg)[2:8], 'MOD_TEMP', 'MOD_PRECIP')
colnames(X_CONGO) = c(colnames(famous_agg)[2:8], 'MOD_TEMP', 'MOD_PRECIP')

X_tropics = rbind(X_AMAZON, X_SEASIA, X_CONGO)

X_tropics_norm = normalize(X_tropics)

Y_tropics = c(famous_agg$AMAZ_MOD_FRAC, famous_agg$SEASIA_MOD_FRAC, famous_agg$CONGO_MOD_FRAC)

tropics_fit = km(~., design = X_tropics_norm, response=Y_tropics)


# normalize amazon obs
tp.amaz.norm <- normalize(
  matrix(c(temps_obs$AMAZ_OBS_TEMP+273.15, precips_obs$AMAZ_OBS_PRECIP),nrow=1),
  wrt=X_tropics[, 8:9]
)
#colnames(tp.amaz.norm) <- c('OBS_TEMP', 'OBS_PRECIP')
model_climate_names = c('MOD_TEMP', 'MOD_PRECIP')
colnames(tp.amaz.norm) <- model_climate_names

tp.seasia.norm <- normalize(
  matrix(c(temps_obs$SEASIA_OBS_TEMP+273.15, precips_obs$SEASIA_OBS_PRECIP),nrow=1),
  wrt=X_tropics[, 8:9]
)
#colnames(tp.seasia.norm) <- c('OBS_TEMP', 'OBS_PRECIP')
colnames(tp.seasia.norm) <- model_climate_names

tp.congo.norm <- normalize(
  matrix(c(temps_obs$CONGO_OBS_TEMP+273.15, precips_obs$CONGO_OBS_PRECIP),nrow=1),
  wrt=X_tropics[, 8:9]
)
#colnames(tp.congo.norm) <- c('OBS_TEMP', 'OBS_PRECIP')
colnames(tp.congo.norm) <- model_climate_names

amaz.x  <- cbind(X.stan.norm, tp.amaz.norm)
congo.x <- cbind(X.stan.norm, tp.congo.norm)
seasia.x <- cbind(X.stan.norm, tp.seasia.norm)

pred.amaz.bc <- predict(tropics_fit, newdata=amaz.x, type='UK')
pred.congo.bc <- predict(tropics_fit, newdata=congo.x, type='UK')
pred.seasia.bc <- predict(tropics_fit, newdata=seasia.x, type='UK')

# Emulator fit with no temp/precip
fit.amazon <- km(~., design=X.norm, response=famous_agg$AMAZ_MOD_FRAC)
fit.congo <- km(~., design=X.norm, response=famous_agg$CONGO_MOD_FRAC)
fit.seasia <- km(~., design=X.norm, response=famous_agg$SEASIA_MOD_FRAC)


#fit.namerica <- km(~., design=X.norm, response=NAMERICA_MOD_FRAC)
#fit.global <- km(~., design=X.norm, response=GLOB_MOD_FRAC)

standard.amazon <- predict(fit.amazon, newdata=X.stan.norm, type='UK')
standard.congo <- predict(fit.congo, newdata=X.stan.norm, type='UK')
standard.seasia <- predict(fit.seasia, newdata=X.stan.norm, type='UK')

#standard.namerica <- predict(fit.namerica, newdata=X.stan.norm, type='UK')
#standard.global <- predict(fit.global, newdata=X.stan.norm, type='UK')

obs_amazon = obs[,'AMAZON']
obs_congo = obs[,'CONGO']
obs_seasia = obs[,'SEASIA']

# It looks like bias correction via T and P improves the Amazon,
# slightly worsens the Congo, and very slightly improves SE Asia. 
obs_amazon - standard.amazon$mean
obs_amazon - pred.amaz.bc$mean

obs_congo - standard.congo$mean
obs_congo - pred.congo.bc$mean

obs_seasia - standard.seasia$mean
obs_seasia- pred.seasia.bc$mean

# --------------------------------------------------------------------
# Sensitivity analysis, including temperature and precipitation
# --------------------------------------------------------------------
library(sensitivity)

n = 21
X.oat = oaat.design(X_tropics_norm, n = n, hold = amaz.x)
colnames(X.oat) <- colnames(X_tropics_norm)

fit.sens = km(~.,design = X_tropics_norm, response = Y_tropics)

pred.sens = predict(fit.sens, newdata = X.oat, type = 'UK')

col.chosen = col.amaz
col.transp = adjustcolor(col.chosen, alpha = 0.5)

#dev.new(width = 5, height = 5)
pdf(width = 7, height = 6, file = 'graphics/sensitivity_TP_amazon.pdf') 
par(mfrow = c(2,5), las = 1, mar = c(5,0.5,2,0.5), oma = c(0,5,0,0), fg = 'grey')
for(i in 1: ncol(X_tropics_norm)){
  
  ix <- seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oat[ix, i], pred.sens$mean[ix], ylim = c(0,1), xlab = colnames(X.oat)[i], type = 'n', axes = FALSE)
  axis(1)
  if (i==1 | i==6 ) {axis(2)
    mtext(side = 2, line = 3.5, text = 'FOREST FRACTION', las = 0, col = 'black')
  }
  
  polygon(x = c(X.oat[ix, i], rev(X.oat[ix, i])),
          y = c( (pred.sens$mean[ix] - pred.sens$sd[ix]), rev(pred.sens$mean[ix] + pred.sens$sd[ix])),
          col = col.transp, border = col.transp)
  
  lines(X.oat[ix, i], pred.sens$mean[ix], ylim = c(0,1), xlab = colnames(X.oat)[i], col = col.chosen )
}
dev.off()

inputs.set <- function(X, y, thres, obs, obs.sd = 0, disc = 0, disc.sd = 0, n = 100000, abt = FALSE){ 
  # find a set of inputs that are consistent with a particular
  # set of implausibility (either below or above)
  
  X.mins <- apply(X,2,min)
  X.maxes <- apply(X,2,max)
  X.unif <- samp.unif(n, mins = X.mins, maxes = X.maxes)
  colnames(X.unif) <- colnames(X)
  
  fit <- km(~., design = X, response = y, control = list(trace = FALSE))
  pred <- predict(fit, newdata = X.unif, type = 'UK')
  pred.impl <- impl(em = pred$mean, em.sd = pred$sd,
                    disc = disc, obs = obs, disc.sd = disc.sd, obs.sd = obs.sd)
  
  if(abt){
    # choose those above the threshold 
    ix.bt <- pred.impl > thres
  }
  
  else{
    ix.bt <- pred.impl < thres
  }
  
  X.out <- X.unif[ix.bt, ]
  
  return(list(X.out = X.out, fit = fit, X.unif = X.unif, pred = pred,pred.impl = pred.impl))   
  
}


# Find the set of plausible inputs, when temperature and precip are included in the inputs
plausible.amazon.bc <- inputs.set(X = X_tropics_norm, y = Y_tropics,thres = 3,
                                  obs = obs_amazon,
                                  obs.sd = 0,
                                  disc = 0,
                                  disc.sd = 0.01,
                                  n = 100000,
                                  abt = FALSE)


# Emulated surface of forest fraction given T & P
#dev.new(width = 5, height =5)
pdf(file = 'graphics/emulated_fraction_vs_temp_precip.pdf',width = 7, height = 7)
quilt.plot(plausible.amazon.bc$X.unif[,8], plausible.amazon.bc$X.unif[, 9], plausible.amazon.bc$pred$mean, 
           col = byr, xlab = 'Normalised regional temperature', ylab = 'Normalised regional precipitation',
           legend.args = list(text = "forest\nfraction",col="black", cex=1.2, side=3, line=1))
cex = 1.5
points(X_tropics_norm[1:100,8], X_tropics_norm[1:100,9], col = 'black', bg = col.amaz, pch = 21, cex = cex)
points(X_tropics_norm[101:200,8], X_tropics_norm[101:200,9], col = 'black', bg = col.seasia, pch = 21, cex = cex)
points(X_tropics_norm[201:300,8], X_tropics_norm[201:300,9], col = 'black', bg = col.congo, pch = 21, cex = cex)

points(tp.amaz.norm, col = 'black', pch = 21, cex = 2.5, bg = col.amaz, lwd = 2.5)
points(tp.congo.norm, col = 'black', pch = 21, cex = 2.5, bg = col.congo, lwd = 2.5)
points(tp.seasia.norm, col = 'black', pch = 21, cex = 2.5, bg = col.seasia, lwd = 2.5)

text(tp.amaz.norm, 'Amazon', pos = 4, font = 2)
text(tp.congo.norm, 'Central Africa', pos = 4, font = 2)
text(tp.seasia.norm, 'SE Asia', pos = 4, font = 2)
dev.off()

# Plot colour as the third dimension, to compare model runs with
# mean emulated surface.
col3rd = function(n, pal, z){
  # produce a set of colours that match the values of z
  # Use for colour in 3rd dimension on a scatterplot.
  cRP  = colorRampPalette(pal)
  cols = cRP(n)
  out = cols[(z - min(z))/diff(range(z))*n + 1]
  out
}

zcolor = col3rd(n=3, pal=byr, z = Y_tropics) 
plot(X_tropics[, 8], X_tropics[, 9], col = 'black', bg = zcolor, pch = 21, cex = 2)

# Normalize the colours to the background
# the first 300 points are Y_tropics
allz = c(Y_tropics,plausible.amazon.bc$pred$mean)
zcolor = col3rd(n=11, pal=byr, z = allz) 
plot(X_tropics[, 8], X_tropics[, 9], col = 'black', bg = zcolor, pch = 21, cex = 2)

pdf(file = 'graphics/emulated_fraction_vs_temp_precip_pcolcor.pdf',width = 7, height = 7)
quilt.plot(plausible.amazon.bc$X.unif[,8], plausible.amazon.bc$X.unif[, 9], plausible.amazon.bc$pred$mean, col = byr, xlab = 'Normalised regional temperature', ylab = 'Normalised regional precipitation', legend.args = list(text = "forest\nfraction",col="black", cex=1.2, side=3, line=1))
cex = 1.5
lwd = 2
points(X_tropics_norm[1:100,8], X_tropics_norm[1:100,9], col = col.amaz, bg = zcolor[1:100], pch = 21, cex = cex, lwd = lwd)
points(X_tropics_norm[101:200,8], X_tropics_norm[101:200,9], col = col.seasia, bg = zcolor[101:200], pch = 21, cex = cex, lwd = lwd)
points(X_tropics_norm[201:300,8], X_tropics_norm[201:300,9], col = col.congo, bg = zcolor[201:300], pch = 21, cex = cex, lwd = lwd)

points(tp.amaz.norm, col = 'black', pch = 21, cex = 2.5, bg = col.amaz, lwd = 2.5)
points(tp.seasia.norm, col = 'black', pch = 21, cex = 2.5, bg = col.seasia, lwd = 2.5)
points(tp.congo.norm, col = 'black', pch = 21, cex = 2.5, bg = col.congo, lwd = 2.5)

text(tp.amaz.norm, 'Amazon', pos = 4, font = 2)
text(tp.congo.norm, 'Central Africa', pos = 4, font = 2)
text(tp.seasia.norm, 'SE Asia', pos = 4, font = 2)
dev.off()

# ------------------------------------------------------------------------
# We worked out the response surface for y = f(X, T, P).
# Now bias correct T and P in all of the ensemble members.
# ------------------------------------------------------------------------

# Fit the whole data set
fit.tropics = km(~.,design = X_tropics_norm, response = Y_tropics)

# Bias correct the model runs using the observed (normalised)
# temperature and precipitation.
bc.x  = rbind(cbind(X.norm, matrix(tp.amaz.norm, nrow = 100, ncol = 2, byrow = TRUE)),
              cbind(X.norm, matrix(tp.seasia.norm, nrow = 100, ncol = 2, byrow = TRUE)),
              cbind(X.norm, matrix(tp.congo.norm, nrow = 100, ncol = 2, byrow = TRUE))
)
colnames(bc.x) <- colnames(X_tropics_norm)

model_runs.bc = predict(fit.tropics, newdata = bc.x , type = 'UK')

plot(Y_tropics,model_runs.bc$mean, xlim = c(0,1), ylim = c(0,1))
abline(0,1)

dev.new()
par(mfrow = c(2,3))
hist(Y_tropics[1:100] - obs_amazon, xlim = c(-1,1))
hist(Y_tropics[101:200] - obs_seasia, xlim = c(-1,1))
hist(Y_tropics[201:300] - obs_congo, xlim = c(-1,1))

hist(model_runs.bc$mean[1:100] - obs_amazon, xlim = c(-1,1))
hist(model_runs.bc$mean[101:200] - obs_seasia, xlim = c(-1,1))
hist(model_runs.bc$mean[201:300] - obs_congo, xlim = c(-1,1))

# How much does the forest fraction improve?
mean(abs(Y_tropics[1:100] - obs_amazon))
mean(abs(model_runs.bc$mean[1:100] - obs_amazon))

mean(abs(Y_tropics[101:200] - obs_seasia))
mean(abs(model_runs.bc$mean[101:200] - obs_seasia))

mean(abs(Y_tropics[201:300] - obs_congo))
mean(abs(model_runs.bc$mean[201:300] - obs_congo))

pdf(width = 7, height = 5, file='graphics/bias_corrected_fractions.pdf')
par(las=1)
plot(c(1,2,3), c(obs_amazon, obs_congo, obs_seasia), xlim=c(0.5,3.5), ylim=c(0,1), pch=19,
     col = c(col.amaz, col.congo, col.seasia), cex=1.5, axes=FALSE, xlab='', ylab='Forest fraction',
     pty = 'n')

points(rep(1.1, 100), Y_tropics[1:100], col = 'grey', pch = '-')
points(rep(1.2, 100), model_runs.bc$mean[1:100], col = 'pink', pch = '-')

points(rep(2.1, 100), Y_tropics[101:200], col = 'grey', pch = '-')
points(rep(2.2, 100), model_runs.bc$mean[101:200], col = 'pink', pch = '-')

points(rep(3.1, 100), Y_tropics[201:300], col = 'grey', pch = '-')
points(rep(3.2, 100), model_runs.bc$mean[201:300], col = 'pink', pch = '-')

points(c(1,2,3), c(obs_amazon, obs_congo, obs_seasia), pch=19,
     col = c(col.amaz, col.congo, col.seasia), cex=1.5)

points(c(1.1,2.1,3.1), c(standard.amazon$mean, standard.congo$mean,standard.seasia$mean), pch=19)

segments(1.1, standard.amazon$mean - standard.amazon$sd,
         1.1, standard.amazon$mean + standard.amazon$sd, col='black')

segments(2.1, standard.congo$mean - standard.congo$sd,
         2.1, standard.congo$mean + standard.congo$sd, col='black')

segments(3.1, standard.seasia$mean - standard.seasia$sd,
         3.1, standard.seasia$mean + standard.seasia$sd, col='black')

points(c(1.2,2.2,3.2),c(pred.amaz.bc$mean, pred.congo.bc$mean, pred.seasia.bc$mean),  col='red3', pch=19)

segments(1.2, pred.amaz.bc$mean - pred.amaz.bc$sd,
         1.2, pred.amaz.bc$mean + standard.amazon$sd, col='red3')

segments(2.2, pred.congo.bc$mean - standard.congo$sd,
         2.2, pred.congo.bc$mean + standard.congo$sd, col='red3')

segments(3.2, pred.seasia.bc$mean - pred.seasia.bc$sd,
         3.2, pred.seasia.bc$mean + pred.seasia.bc$sd, col='red3')

axis(1, labels = c('Amazon', 'Africa', 'SE Asia'), at = c(1.1,2.1,3.1))
axis(2)

text(1, obs_amazon, 'observation', pos=2, col='grey', cex=0.9)
text(1.1, standard.amazon$mean, 'default\nparameters', pos=2, cex=0.9, col='grey')
text(1.2, pred.amaz.bc$mean, 'bias\ncorrected', pos=4, col='red3', cex=0.9)
legend('topright', legend = c('model runs', 'bias corrected model runs'),
                              col = c('grey', 'pink'), pch = '-', bty = 'n'
       )
dev.off()

# --------------------------------------------------------
# Find points which are NROY for all three systems,
# When T and P are held at observed values and not 
# --------------------------------------------------------
# Find the set of plausible inputs, when temperature and precip are included in the inputs

# create a 'core' X
n = 100000
X.mins = apply(X.norm,2,min)
X.maxes = apply(X.norm,2,max)
X.unif = samp.unif(n, mins = X.mins, maxes = X.maxes)

# append the observed T and P
X.unif.amaz = cbind(X.unif, 
                    matrix(tp.amaz.norm, nrow = n,ncol = 2, byrow = TRUE))
colnames(X.unif.amaz) <- colnames(X_tropics_norm)
X.unif.seasia = cbind(X.unif, 
                    matrix(tp.seasia.norm, nrow = n,ncol = 2, byrow = TRUE))
colnames(X.unif.seasia) <- colnames(X_tropics_norm)
X.unif.congo = cbind(X.unif, 
                      matrix(tp.congo.norm, nrow = n,ncol = 2, byrow = TRUE))
colnames(X.unif.congo) <- colnames(X_tropics_norm)

disc = 0
obs.sd = 0
disc.sd = 0.01
thres = 3

pred.unif.amaz = predict(fit.amazon, newdata = X.unif, type = 'UK')
amaz.impl = impl(em = pred.unif.amaz$mean, em.sd = pred.unif.amaz$sd,
                 disc = disc, obs = obs_amazon,
                 disc.sd = disc.sd,
                 obs.sd = obs.sd)
nroy.ix.amaz = which(amaz.impl < thres)

pdf(width = 7, height = 7, file = 'graphics/best_inputs_amazon.pdf')
pairs(X.unif[nroy.ix.amaz, ], panel = dfunc.up, gap = 0, upper.panel = NULL)
dev.off()


pred.unif.amaz.bc = predict(fit.tropics, newdata = X.unif.amaz, type = 'UK')
amaz.impl.bc = impl(em = pred.unif.amaz.bc$mean, em.sd = pred.unif.amaz.bc$sd,
                 disc = disc, obs = obs_amazon,
                 disc.sd = disc.sd,
                 obs.sd = obs.sd)
nroy.ix.amaz.bc = which(amaz.impl.bc < thres)

pdf(width = 7, height = 7, file = 'graphics/best_inputs_amazon_bc.pdf')
pairs(X.unif[nroy.ix.amaz.bc, ], panel = dfunc.up, gap = 0, upper.panel = NULL)
dev.off()


# Now SE Asia
pred.unif.seasia = predict(fit.seasia, newdata = X.unif, type = 'UK')
seasia.impl = impl(em = pred.unif.seasia$mean, em.sd = pred.unif.seasia$sd,
                 disc = disc, obs = obs_seasia,
                 disc.sd = disc.sd,
                 obs.sd = obs.sd)
nroy.ix.seasia = which(seasia.impl < thres)

pdf(width = 7, height = 7, file = 'graphics/best_inputs_seasia.pdf')
pairs(X.unif[nroy.ix.seasia, ], panel = dfunc.up, gap = 0, upper.panel = NULL)
dev.off()

pred.unif.seasia.bc = predict(fit.tropics, newdata = X.unif.seasia, type = 'UK')
seasia.impl.bc = impl(em = pred.unif.seasia.bc$mean, em.sd = pred.unif.seasia.bc$sd,
                    disc = disc, obs = obs_seasia,
                    disc.sd = disc.sd,
                    obs.sd = obs.sd)
nroy.ix.seasia.bc = which(seasia.impl.bc < thres)

pdf(width = 7, height = 7, file = 'graphics/best_inputs_seasia_bc.pdf')
pairs(X.unif[nroy.ix.seasia.bc, ], panel = dfunc.up, gap = 0, upper.panel = NULL)
dev.off()

# Now Congo
pred.unif.congo = predict(fit.congo, newdata = X.unif, type = 'UK')
congo.impl = impl(em = pred.unif.congo$mean, em.sd = pred.unif.congo$sd,
                   disc = disc, obs = obs_congo,
                   disc.sd = disc.sd,
                   obs.sd = obs.sd)
nroy.ix.congo = which(congo.impl < thres)

pdf(width = 7, height = 7, file = 'graphics/best_inputs_congo.pdf')
pairs(X.unif[nroy.ix.congo, ], panel = dfunc.up, gap = 0, upper.panel = NULL)
dev.off()

pred.unif.congo.bc = predict(fit.tropics, newdata = X.unif.congo, type = 'UK')
congo.impl.bc = impl(em = pred.unif.congo.bc$mean, em.sd = pred.unif.congo.bc$sd,
                      disc = disc, obs = obs_congo,
                      disc.sd = disc.sd,
                      obs.sd = obs.sd)
nroy.ix.congo.bc = which(congo.impl.bc < thres)

pdf(width = 7, height = 7, file = 'graphics/best_inputs_congo_bc.pdf')
pairs(X.unif[nroy.ix.congo.bc, ], panel = dfunc.up, gap = 0, upper.panel = NULL)
dev.off()


# what part of parameter space matches everything?
# Very low uncertainty here.

nroy.bc.ix = intersect(intersect(nroy.ix.amaz.bc,nroy.ix.seasia.bc ), nroy.ix.congo.bc)
nroy.nobc.ix = intersect(intersect(nroy.ix.amaz,nroy.ix.seasia ), nroy.ix.congo)

(length(nroy.nobc.ix)/n) *100 # 2% of space is NROY with no bias correction.
(length(nroy.bc.ix)/n) *100 # 28% of space is NROY with bias correction.

pdf(width = 7, height = 7, file = 'graphics/best_inputs_all_bc.pdf')
pairs(X.unif[nroy.bc.ix, ], panel = dfunc.up, gap = 0, upper.panel = NULL)
dev.off()

pdf(width = 7, height = 7, file = 'graphics/best_inputs_all_nobc.pdf')
pairs(X.unif[nroy.nobc.ix, ], panel = dfunc.up, gap = 0, upper.panel = NULL)
dev.off()


# What fraction of parameter space is NROY for all three observations? 




