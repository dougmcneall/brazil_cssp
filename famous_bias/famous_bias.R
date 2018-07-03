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
cbPal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


col.amaz <- cbPal[1]
col.seasia <- cbPal[2]
col.congo <- cbPal[3]
col.global <- cbPal[4]
col.namerica <- cbPal[5]

col.tropics = c(rep(col.amaz, 100), rep(col.seasia,100), rep(col.congo,100))


#col.amaz <- acc[1]
#col.namerica <- acc[2]
#col.seasia <- acc[3]
#col.congo <- acc[4]
#col.global <- acc[5]


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

shadowtext <- function(x, y=NULL, labels, 
                       col='white', bg='black',
                       theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {
  # put text on a plot - white with a black background
  xy <- xy.coords(x,y)
  xo <- r*strwidth('A')
  yo <- r*strheight('A')
  for (i in theta) {
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, 
          labels, col=bg, ... )
  }
  text(xy$x, xy$y, labels, col=col, ... )
}

# Plot colour as the third dimension, to compare model runs with
# mean emulated surface.
col3rd = function(n, pal, z){
  # produce a set of colours that match the values of z
  # Use for colour in 3rd dimension on a scatterplot.
  cRP  = colorRampPalette(pal)
  cols = cRP(n)
  
  #out = cols[(z - min(z))/diff(range(z))*n + 1 ]
  out = cols[cut(z, breaks = n) ]
  out
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

# ----------------------------------------------------------------------------------
# Emulator diagnostics
#
# ----------------------------------------------------------------------------------

# Plot the emulator against the true leave-one-out prediction
# This usually gives a very similar response to trend.reestim = TRUE
# in leaveOneOut.km

true.loo = function(X,y){
  out.mean = rep(NA, length(y))
  out.sd = rep(NA, length(y))
  
  for(i in 1:nrow(X)){
    X.trunc = X[-i, ]
    y.trunc = y[-i]
    
    X.target = matrix(X[i, ], nrow = 1)
    colnames(X.target) <- colnames(X)
    X.target.df = data.frame(X.target)
    
      fit = km(~., design = X.trunc, response = y.trunc)
      pred = predict(fit,newdata = X.target, type = 'UK')
      out.mean[i] = pred$mean
      out.sd[i] = pred$sd
  }
  return(list(mean = out.mean, sd = out.sd))
}  

true.loo.all = true.loo(X = X_tropics_norm, y = Y_tropics)

# Mean absolute error is about 0.03 or 3%
print(paste('With T/P mean absolute cross validation error = ', mean(abs(true.loo.all$mean - Y_tropics))))


pdf(width = 6, height = 6, file = 'graphics/true_loo_all.pdf' )
xlim = c(-0.05, 1.05)
ylim = c(-0.05, 1.05)
par(las =1)
plot(Y_tropics, true.loo.all$mean, pch = 20,
     xlab = 'observation', ylab = 'prediction',
     col = col.tropics,
     xlim = xlim, 
     ylim = ylim,
    bty = 'n',
    axes = FALSE,
    xaxs = 'i', yaxs = 'i')

segments(x0 = Y_tropics, y0 = true.loo.all$mean - (2*true.loo.all$sd),
          x1 = Y_tropics, y1 = true.loo.all$mean +(2*true.loo.all$sd),
         col = col.tropics)
axis(1, pos = 0, col = 'grey')
axis(2, pos = 0, col = 'grey')
abline(0,1, col = 'grey')
legend('top', legend = c('Amazon', 'Asia', 'Africa'),
       pch = 20, col = c(col.amaz, col.seasia, col.congo),
       bty = 'n')
dev.off()


# Rank histograms for checking the uncertainty?
# The principle behind the rank histogram is quite simple. 
# Ideally, one property that is desired from an EF is reliable probabilities;
# if ensemble relative frequency suggests P percent probability of occurrence,
# the event truly ought to have P probability of occurring. 
# For this probability to be reliable, the set of ensemble member forecast values
# at a given point and the true state (the verification) ought to be able to be 
# considered random samples from the same probability distribution.
# This reliability then implies in turn that if an n-member ensemble and the
# verification are pooled into a vector and sorted from lowest to highest,
# then the verification is equally likely to occur in each of the n + 1 possible ranks. 
# If the rank of the verification is tallied and the process repeated over many
# independent sample points, a uniform histogram over the possible ranks should result.
# From https://journals.ametsoc.org/doi/full/10.1175/1520-0493%282001%29129%3C0550%3AIORHFV%3E2.0.CO%3B2


loo.rankhist = function(obs, pred.mean, pred.sd, n = 1000){
  # a version of the rank histogram     
  obs.ranks = rep(NA, length(obs))
  
  for(i in 1:length(obs)){
    ranks = rank(c(obs[i], rnorm(n = n, mean = pred.mean[i], sd = pred.sd[i])))
    obs.ranks[i] = ranks[1]
  }
  
  out = obs.ranks/(n+1)
  out
}

true.loo.rankhist = loo.rankhist(obs = Y_tropics, pred.mean = true.loo.all$mean, pred.sd = true.loo.all$sd, n = 500)

pdf(file = 'graphics/rankhist.pdf', width = 6, height =4)
par(las = 1, fg = 'white')
hist(true.loo.rankhist, col = 'grey', 
     axes = FALSE, main = '',
     xlab = 'Rank'
     )
axis(1, col = 'black')
axis(2, col = 'black')
dev.off()

# Prediction errors normalised by their standard deviations
# approximately follow a normal distribution - pehaps the tails are a little off.
true.loo.err.norm = (true.loo.all$mean - Y_tropics) / true.loo.all$sd
pdf(file = 'graphics/normalisedQQ.pdf')
par(las = 1)
qqnorm(true.loo.err.norm)
abline(0,1)
dev.off()


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


plot(X_tropics[, 8], X_tropics[, 9], col = 'black', bg = zcolor, pch = 21, cex = 2)

# Normalize the colours to the background
# the first 300 points are Y_tropics

allz = c(Y_tropics,obs_amazon,obs_seasia, obs_congo, plausible.amazon.bc$pred$mean)
zcolor = col3rd(n=9, pal=viridis(7), z = allz) 

pdf(file = 'graphics/emulated_fraction_vs_temp_precip_pcolcor.pdf',width = 7, height = 7)
par(las = 1)
quilt.plot(plausible.amazon.bc$X.unif[,8],
           plausible.amazon.bc$X.unif[, 9], plausible.amazon.bc$pred$mean,
           col = viridis(7),
           xlab = 'Normalised regional temperature', 
           ylab = 'Normalised regional precipitation', 
           legend.args = list(text = "forest\nfraction",
                              col="black", cex=1.2, side=3, line=1))
cex = 1.4
lwd = 1.5
points(X_tropics_norm[1:100,8], X_tropics_norm[1:100,9], 
       col = 'black', bg = zcolor[1:100], pch = 21, cex = cex, lwd = lwd)
points(X_tropics_norm[101:200,8], X_tropics_norm[101:200,9], 
       col = 'black', bg = zcolor[101:200], pch = 22, cex = cex, lwd = lwd)
points(X_tropics_norm[201:300,8], X_tropics_norm[201:300,9], 
       col = 'black', bg = zcolor[201:300], pch = 24, cex = cex, lwd = lwd)

points(tp.amaz.norm, col = 'black', pch = 21, cex = 2.5, bg = zcolor[301], lwd = 2)
points(tp.seasia.norm, col = 'black', pch = 22, cex = 2.5, bg = zcolor[302], lwd = 2)
points(tp.congo.norm, col = 'black', pch = 24, cex = 2.5, bg = zcolor[303], lwd = 2)

shadowtext(tp.amaz.norm[1],tp.amaz.norm[2], 'Amazon', pos = 4, font = 2,r =0.2)
shadowtext(tp.congo.norm[1], tp.congo.norm[2], 'Central Africa', pos = 4, font = 2, r = 0.2)
shadowtext(tp.seasia.norm[1], tp.seasia.norm[2], 'SE Asia', pos = 4, font = 2, r = 0.2)
dev.off()

# No emulated surface in this version


pdf(file = 'graphics/fraction_vs_temp_precip_pcolcor.pdf',width = 8, height = 7)
par(las = 1, fg = 'black', mar = c(5,6,3,7))
Y_obs = c(Y_tropics, obs_amazon, obs_seasia, obs_congo)
Y_obs_ix = 1:300
zcolor = col3rd(n=9, pal= viridis(9), z = Y_obs) 
pr =1e5
plot(X_tropics[, 8]-273.15, X_tropics[, 9]*pr, col = 'black', pch = 21, cex = 2, type = 'n',
     xlab = expression(paste('Regional temperature (', degree,'C)')),
     ylab = expression(paste('Regional Precipitation x',10^5,' kgm'^-2,'s'^-1))
)
                       


cex = 1.4
lwd = 1.5

points(X_tropics[1:100,8]-273.15, X_tropics[1:100,9]*pr, col = 'black', bg = zcolor[1:100], pch = 21, cex = cex, lwd = lwd)
points(X_tropics[101:200,8]-273.15, X_tropics[101:200,9]*pr, col = 'black', bg = zcolor[101:200], pch = 22, cex = cex, lwd = lwd)
points(X_tropics[201:300,8]-273.15, X_tropics[201:300,9]*pr, col = 'black', bg = zcolor[201:300], pch = 24, cex = cex, lwd = lwd)

points(temps_obs$AMAZ_OBS_TEMP, precips_obs$AMAZ_OBS_PRECIP*pr, col = 'black', pch = 21, cex = 2.5, bg = zcolor[301], lwd = 2)
points(temps_obs$SEASIA_OBS_TEMP, precips_obs$SEASIA_OBS_PRECIP*pr, col = 'black', pch = 22, cex = 2.5, bg = zcolor[302], lwd = 2)
points(temps_obs$CONGO_OBS_TEMP, precips_obs$CONGO_OBS_PRECIP*pr, col = 'black', pch = 24, cex = 2.5, bg = zcolor[303], lwd = 2)

shadowtext(temps_obs$AMAZ_OBS_TEMP,precips_obs$AMAZ_OBS_PRECIP*pr, 'Amazon', pos = 4, font = 2,r =0.2)
shadowtext(temps_obs$SEASIA_OBS_TEMP, precips_obs$SEASIA_OBS_PRECIP*pr, 'Central Africa', pos = 4, font = 2, r = 0.2)
shadowtext(temps_obs$CONGO_OBS_TEMP, precips_obs$CONGO_OBS_PRECIP*pr, 'SE Asia', pos = 4, font = 2, r = 0.2)
image.plot(z = Y_obs, legend.only = TRUE, col = viridis(9), horizontal = FALSE,  legend.args = list(text = "forest\nfraction",col="black", cex=1.2, side=3, line=1))
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


# --------------------------------------------------------------
# Here's a fun thing - suggested by Michael Goldstein - 
# What value does the model add over just using T and P to
# fit the data?
# -------------------------------------------------------------

# Produce genuine LOO for all these, put them together and compare with 
# true.loo
fit.x.amaz = km(~., design = X.norm, response=famous_agg$AMAZ_MOD_FRAC)
fit.x.seasia = km(~., design = X.norm, response=famous_agg$SEASIA_MOD_FRAC)
fit.x.congo = km(~., design = X.norm, response=famous_agg$CONGO_MOD_FRAC)

# This is much quicker than adding them all together!
true.loo.amaz = true.loo(X = X.norm, y = famous_agg$AMAZ_MOD_FRAC)
true.loo.seasia = true.loo(X = X.norm, y = famous_agg$SEASIA_MOD_FRAC)
true.loo.congo = true.loo(X = X.norm, y = famous_agg$CONGO_MOD_FRAC)

true.loo.X.mean = c(true.loo.amaz$mean, true.loo.seasia$mean, true.loo.congo$mean)
true.loo.X.sd = c(true.loo.amaz$sd, true.loo.seasia$sd, true.loo.congo$sd)


# Mean absolute error is about 0.06 or 6%
print(paste('Just X mean absolute cross validation error =', mean(abs(true.loo.X.mean - Y_tropics))))


plot(Y_tropics, true.loo.X.mean)
pdf(width = 6, height = 6, file = 'graphics/true_loo_X.pdf' )
xlim = c(-0.05, 1.05)
ylim = c(-0.05, 1.05)
par(las =1)
plot(Y_tropics, true.loo.X.mean, pch = 20,
     xlab = 'observation', ylab = 'prediction',
     col = col.tropics,
     xlim = xlim, 
     ylim = ylim,
     bty = 'n',
     axes = FALSE,
     xaxs = 'i', yaxs = 'i')

segments(x0 = Y_tropics, y0 = true.loo.X.mean - (2*true.loo.X.sd),
         x1 = Y_tropics, y1 = true.loo.X.mean +(2*true.loo.X.sd),
         col = col.tropics)
axis(1, pos = 0, col = 'grey')
axis(2, pos = 0, col = 'grey')
abline(0,1, col = 'grey')
legend('top', legend = c('Amazon', 'Asia', 'Africa'),
       pch = 20, col = c(col.amaz, col.seasia, col.congo),
       bty = 'n')
dev.off()


# Mean absolute error is about 0.03 or 3%
print(paste('mean absolute cross validation error = ', mean(abs(true.loo.all$mean - Y_tropics))))




# This is much quicker than adding them all together!
true.loo.tp.amaz = true.loo(X = X_tropics_norm[1:100, 8:9], y = famous_agg$AMAZ_MOD_FRAC)
true.loo.tp.seasia = true.loo(X = X_tropics_norm[1:100, 8:9], y = famous_agg$SEASIA_MOD_FRAC)
true.loo.tp.congo = true.loo(X = X_tropics_norm[1:100, 8:9], y = famous_agg$CONGO_MOD_FRAC)

true.loo.tp.mean = c(true.loo.tp.amaz$mean, true.loo.tp.seasia$mean, true.loo.tp.congo$mean)
true.loo.tp.sd = c(true.loo.tp.amaz$sd, true.loo.tp.seasia$sd, true.loo.tp.congo$sd)
true.loo.tp.err  = Y_tropics

print(paste('mean absolute cross validation error = ', mean(abs(true.loo.tp.mean - Y_tropics))))
pdf(width = 6, height = 6, file = 'graphics/true_loo_tp.pdf' )
xlim = c(-0.05, 1.05)
ylim = c(-0.05, 1.05)
par(las =1)
plot(Y_tropics, true.loo.tp.mean, pch = 20,
     xlab = 'observation', ylab = 'prediction',
     col = col.tropics,
     xlim = xlim, 
     ylim = ylim,
     bty = 'n',
     axes = FALSE,
     xaxs = 'i', yaxs = 'i')

segments(x0 = Y_tropics, y0 = true.loo.X.mean - (2*true.loo.tp.sd),
         x1 = Y_tropics, y1 = true.loo.X.mean +(2*true.loo.tp.sd),
         col = col.tropics)
axis(1, pos = 0, col = 'grey')
axis(2, pos = 0, col = 'grey')
abline(0,1, col = 'grey')
legend('top', legend = c('Amazon', 'Asia', 'Africa'),
       pch = 20, col = c(col.amaz, col.seasia, col.congo),
       bty = 'n')
dev.off()






# fit using just temperature andprecip
fit.tp  = km(~., design = X_tropics_norm[, 8:9], response=Y_tropics)
pred.tp = leaveOneOut.km(fit.tp, type="UK", trend.reestim=TRUE)

plot(Y_tropics, pred.tp$mean )
fit.tp.resid = pred.tp$mean - Y_tropics

# Have to split these out per-forest
fit.resid.amazon = km(~., design = X.norm, response=fit.tp.resid[1:100])
fit.resid.seasia = km(~., design = X.norm, response=fit.tp.resid[101:200])
fit.resid.congo = km(~., design = X.norm, response=fit.tp.resid[201:300])

plot(fit.resid.congo)


