# extended_famous_bias.R
# An extended look at fixing a bias in famous by using
# temperature and precipitation as inputs to the emulator.
# Now with more regions.
# dougmcneall@gmail.com

# ------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------
library(DiceKriging)
library(RColorBrewer)
library(MASS)
library(fields)
library(parallel)

source('https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R')
source('https://raw.githubusercontent.com/dougmcneall/packages-git/master/imptools.R')
source('https://raw.githubusercontent.com/dougmcneall/packages-git/master/vistools.R')

load('famous_forest_fraction_v2.RData')
load('famous_forest_fraction.RData') # for standard X at the moment

X = famous_agg_v2[, 2:8]
X.norm = normalize(X)
X.stan.norm <- normalize(matrix(X.standard, nrow = 1), wrt = X)
colnames(X.stan.norm) = colnames(X)


cbPal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

col.amaz <- cbPal[1]
col.seasia <- cbPal[2]
col.congo <- cbPal[3]

col.namerica <- cbPal[4]
col.fmec <- cbPal[5]
col.eurasia <- cbPal[6]
col.global <- cbPal[7]


X_AMAZON = famous_agg_v2[ ,c('F0', 'LAI_MIN', 'NL0', 'R_GROW', 'TUPP',
                             'Q10', 'V_CRIT_ALPHA', 'AMAZ_MOD_TEMP',
                             'AMAZ_MOD_PRECIP')] # AMAZON
X_SEASIA = famous_agg_v2[ ,c('F0', 'LAI_MIN', 'NL0', 'R_GROW', 'TUPP',
                             'Q10', 'V_CRIT_ALPHA', 'SEASIA_MOD_TEMP',
                             'SEASIA_MOD_PRECIP')] # SEASIA
X_CONGO = famous_agg_v2[ ,c('F0', 'LAI_MIN', 'NL0', 'R_GROW', 'TUPP',
                           'Q10', 'V_CRIT_ALPHA', 'CONGO_MOD_TEMP',
                           'CONGO_MOD_PRECIP') ] # CONGO
X_NAMERICA = famous_agg_v2[,c('F0', 'LAI_MIN', 'NL0', 'R_GROW', 'TUPP',
                            'Q10', 'V_CRIT_ALPHA', 'NAMERICA_MOD_TEMP',
                           'NAMERICA_MOD_PRECIP') ] # NAMERICA
X_FMEC = famous_agg_v2[ ,c('F0', 'LAI_MIN', 'NL0', 'R_GROW', 'TUPP',
                              'Q10', 'V_CRIT_ALPHA', 'FMEC_MOD_TEMP',
                              'FMEC_MOD_PRECIP') ] # FMEC
X_EURASIA = famous_agg_v2[ ,c('F0', 'LAI_MIN', 'NL0', 'R_GROW', 'TUPP',
                          'Q10', 'V_CRIT_ALPHA', 'EURASIA_MOD_TEMP',
                          'EURASIA_MOD_PRECIP') ] # 
Xnames = c('F0', 'LAI_MIN', 'NL0', 'R_GROW', 'TUPP',
                    'Q10', 'V_CRIT_ALPHA', 'MOD_TEMP',
                    'MOD_PRECIP')
names(X_AMAZON) = Xnames
names(X_SEASIA) = Xnames
names(X_CONGO) = Xnames
names(X_NAMERICA) = Xnames
names(X_FMEC) = Xnames
names(X_EURASIA) = Xnames

X_forests = rbind(X_AMAZON, X_SEASIA, X_CONGO, X_NAMERICA, X_FMEC, X_EURASIA)

Y_forests = c(famous_agg_v2$AMAZ_MOD_FRAC,
              famous_agg_v2$SEASIA_MOD_FRAC,
              famous_agg_v2$CONGO_MOD_FRAC,
              famous_agg_v2$NAMERICA_MOD_FRAC,
              famous_agg_v2$FMEC_MOD_FRAC,
              famous_agg_v2$EURASIA_MOD_FRAC)

Y_forests.cols = c(rep(col.amaz, 100),
                   rep(col.seasia, 100),
                   rep(col.congo, 100),
                   rep(col.namerica, 100),
                   rep(col.fmec, 100),
                   rep(col.eurasia, 100))

X_forests_norm = normalize(X_forests)



# Forests in Temperature/precipitation phase space
plot(X_forests[,'MOD_TEMP'] - 273.15, X_forests[, 'MOD_PRECIP'],
     xlab = 'Temperature (deg C)', ylab = 'Precipitation', col = Y_forests.cols)

points(temps_obs_v2[1:6], precips_obs_v2[1:6],
       col = 'black', bg =  cbPal[1:6], pch = 21 )



# Where is everything compared to its observations?
dev.new(width = 5, height = 6)
par(mfrow = c(3,2))
breaks = seq(from = 0, to = 1, by = 0.05)
hist(famous_agg_v2$AMAZ_MOD_FRAC, col = col.amaz, xlim = c(0,1), breaks = breaks)
rug(forest_fraction_obs_v2$AMAZ_OBS_FRAC, lwd = 3)

hist(famous_agg_v2$SEASIA_MOD_FRAC, col = col.seasia, xlim = c(0,1), breaks = breaks)
rug(forest_fraction_obs_v2$SEASIA_OBS_FRAC, lwd = 3)

hist(famous_agg_v2$CONGO_MOD_FRAC, col = col.congo, xlim = c(0,1), breaks = breaks)
rug(forest_fraction_obs_v2$CONGO_OBS_FRAC, lwd = 3)

hist(famous_agg_v2$NAMERICA_MOD_FRAC, col = col.namerica, xlim = c(0,1), breaks = breaks)
rug(forest_fraction_obs_v2$NAMERICA_OBS_FRAC, lwd = 3)

hist(famous_agg_v2$FMEC_MOD_FRAC, col = col.fmec, xlim = c(0,1), breaks = breaks)
rug(forest_fraction_obs_v2$FMEC_OBS_FRAC, lwd = 3)

hist(famous_agg_v2$EURASIA_MOD_FRAC, col = col.eurasia, xlim = c(0,1), breaks = breaks)
rug(forest_fraction_obs_v2$ERUASIA_OBS_FRAC, lwd = 3)


# Looks like this works pretty well 
forests_fit = km(~., design = X_forests_norm, response=Y_forests)


# Bias correct and see if we do better at the standard parameters
model_climate_names = c('MOD_TEMP', 'MOD_PRECIP')

# normalize amazon obs
tp.amaz.norm <- normalize(
  matrix(c(temps_obs_v2$AMAZ_OBS_TEMP+273.15, precips_obs_v2$AMAZ_OBS_PRECIP),nrow=1),
  wrt=X_forests[, 8:9]
)
colnames(tp.amaz.norm) <- model_climate_names

tp.seasia.norm <- normalize(
  matrix(c(temps_obs_v2$SEASIA_OBS_TEMP+273.15, precips_obs_v2$SEASIA_OBS_PRECIP),nrow=1),
  wrt=X_forests[, 8:9]
)
colnames(tp.seasia.norm) <- model_climate_names

tp.congo.norm <- normalize(
  matrix(c(temps_obs_v2$CONGO_OBS_TEMP+273.15, precips_obs_v2$CONGO_OBS_PRECIP),nrow=1),
  wrt=X_forests[, 8:9]
)
colnames(tp.congo.norm) <- model_climate_names

tp.namerica.norm <- normalize(
  matrix(c(temps_obs_v2$NAMERICA_OBS_TEMP+273.15, precips_obs_v2$NAMERICA_OBS_PRECIP),nrow=1),
  wrt=X_forests[, 8:9]
)
colnames(tp.namerica.norm) <- model_climate_names

tp.fmec.norm <- normalize(
  matrix(c(temps_obs_v2$FMEC_OBS_TEMP+273.15, precips_obs_v2$FMEC_OBS_PRECIP),nrow=1),
  wrt=X_forests[, 8:9]
)
colnames(tp.fmec.norm) <- model_climate_names

tp.eurasia.norm <- normalize(
  matrix(c(temps_obs_v2$ERUASIA_OBS_TEMP+273.15, precips_obs_v2$ERUASIA_OBS_PRECIP),nrow=1),
  wrt=X_forests[, 8:9]
)
colnames(tp.eurasia.norm) <- model_climate_names


# Default parameters attached and observed temperature and precip
amaz.x  <- cbind(X.stan.norm, tp.amaz.norm)
seasia.x <- cbind(X.stan.norm, tp.seasia.norm)
congo.x <- cbind(X.stan.norm, tp.congo.norm)
namerica.x <- cbind(X.stan.norm, tp.namerica.norm)
fmec.x <- cbind(X.stan.norm, tp.fmec.norm)
eurasia.x <- cbind(X.stan.norm, tp.eurasia.norm)


# Emulator predicted bias corrected at default parameters
pred.amaz.bc <- predict(forests_fit, newdata=amaz.x, type='UK')
pred.seasia.bc <- predict(forests_fit, newdata=seasia.x, type='UK')
pred.congo.bc <- predict(forests_fit, newdata=congo.x, type='UK')
pred.namerica.bc <- predict(forests_fit, newdata=namerica.x, type='UK')
pred.fmec.bc <- predict(forests_fit, newdata=fmec.x, type='UK')
pred.eurasia.bc <- predict(forests_fit, newdata=eurasia.x, type='UK')


# Emulator fit with no temp/precip
fit.amazon <- km(~., design=X.norm, response=famous_agg_v2$AMAZ_MOD_FRAC)
fit.seasia <- km(~., design=X.norm, response=famous_agg_v2$SEASIA_MOD_FRAC)
fit.congo <- km(~., design=X.norm, response=famous_agg_v2$CONGO_MOD_FRAC)
fit.namerica <- km(~., design=X.norm, response=famous_agg_v2$NAMERICA_MOD_FRAC)
fit.fmec <- km(~., design=X.norm, response=famous_agg_v2$FMEC_MOD_FRAC)
fit.eurasia <- km(~., design=X.norm, response=famous_agg_v2$EURASIA_MOD_FRAC)


standard.amazon <- predict(fit.amazon, newdata=X.stan.norm, type='UK')
standard.seasia <- predict(fit.seasia, newdata=X.stan.norm, type='UK')
standard.congo <- predict(fit.congo, newdata=X.stan.norm, type='UK')
standard.namerica <- predict(fit.namerica, newdata=X.stan.norm, type='UK')
standard.fmec <- predict(fit.fmec, newdata=X.stan.norm, type='UK')
standard.eurasia <- predict(fit.eurasia, newdata=X.stan.norm, type='UK')


pdf(file = 'graphics/extended_dotchart_fractions.pdf', width = 6, height = 7)
par(mar = c(5,7,2,8))
plot(c(forest_fraction_obs_v2$AMAZ_OBS_FRAC,
       forest_fraction_obs_v2$SEASIA_OBS_FRAC,
       forest_fraction_obs_v2$CONGO_OBS_FRAC,
       forest_fraction_obs_v2$NAMERICA_OBS_FRAC,
       forest_fraction_obs_v2$FMEC_OBS_FRAC,
       forest_fraction_obs_v2$ERUASIA_OBS_FRAC),
     c(1,2,3,4,5,6), ylim=c(0.5,7), xlim=c(0,1), pch=19,
     col = 'blue', cex=1.5, ylab='', xlab='Forest fraction',
     axes = FALSE, type = 'n', bty = 'l',
     panel.first = rect(c(0,0,0),c(0.5,2.5,4.5),c(1,1,1),c(1.5,3.5,5.5), 
                        col = 'grey95', border = 'grey95'),
     xaxs = 'i', yaxs = 'i'
)

abline(h = c(0.75, 1, 1.25, 1.75, 2, 2.25, 2.75, 3, 3.25, 3.75, 4, 4.25, 4.75,5,5.25, 5.75, 6, 6.25), col = 'grey', lty = 'dashed')

# observed points
points(c(forest_fraction_obs_v2$AMAZ_OBS_FRAC,
         forest_fraction_obs_v2$SEASIA_OBS_FRAC,
         forest_fraction_obs_v2$CONGO_OBS_FRAC,
         forest_fraction_obs_v2$NAMERICA_OBS_FRAC,
         forest_fraction_obs_v2$FMEC_OBS_FRAC,
         forest_fraction_obs_v2$ERUASIA_OBS_FRAC),c(1.25,2.25,3.25, 4.25, 5.25, 6.25), pch=19,
       col = cbPal[1:6], cex=1.5)


# standard model run points
points(c(standard.amazon$mean,standard.seasia$mean,  standard.congo$mean,
         standard.namerica$mean, standard.fmec$mean, standard.eurasia$mean),
       1:6, pch=17, col = cbPal[1:6], cex = 1.5)

segments(standard.amazon$mean - standard.amazon$sd,1,
         standard.amazon$mean + standard.amazon$sd,1, col=col.amaz, lwd = 1.5)

segments(standard.seasia$mean - standard.seasia$sd,2,
         standard.seasia$mean + standard.seasia$sd,2, col=col.seasia, lwd = 1.5)

segments(standard.congo$mean - standard.congo$sd,3,
         standard.congo$mean + standard.congo$sd,3, col=col.congo, lwd =1.5)

segments(standard.namerica$mean - standard.namerica$sd,4,
         standard.namerica$mean + standard.namerica$sd,4, col=col.namerica, lwd =1.5)

segments(standard.fmec$mean - standard.fmec$sd,5,
         standard.fmec$mean + standard.fmec$sd,5, col=col.fmec, lwd =1.5)

segments(standard.eurasia$mean - standard.eurasia$sd,6,
         standard.eurasia$mean + standard.eurasia$sd,6, col=col.eurasia, lwd =1.5)


# bias corrected points
points(c(pred.amaz.bc$mean, pred.seasia.bc$mean,  pred.congo.bc$mean,
         pred.namerica.bc$mean,pred.fmec.bc$mean, pred.eurasia.bc$mean),
       c(0.75,1.75,2.75, 3.75, 4.75, 5.75), col=cbPal[1:6], pch=15, cex = 1.5)

segments( pred.amaz.bc$mean - pred.amaz.bc$sd, 0.75,
          pred.amaz.bc$mean + pred.amaz.bc$sd, 0.75, col=col.amaz, lwd = 1.5)

segments(pred.seasia.bc$mean - pred.seasia.bc$sd,1.75,
         pred.seasia.bc$mean + pred.seasia.bc$sd,1.75, col=col.seasia, lwd = 1.5)

segments(pred.congo.bc$mean - pred.congo.bc$sd,2.75,
         pred.congo.bc$mean + pred.congo.bc$sd,2.75,  col=col.congo, lwd = 1.5)

segments(pred.namerica.bc$mean - pred.namerica.bc$sd,3.75,
         pred.namerica.bc$mean + pred.namerica.bc$sd,3.75,  col=col.namerica, lwd = 1.5)

segments(pred.fmec.bc$mean - pred.fmec.bc$sd,4.75,
         pred.fmec.bc$mean + pred.fmec.bc$sd,4.75,  col=col.fmec, lwd = 1.5)

segments(pred.eurasia.bc$mean - pred.eurasia.bc$sd,5.75,
         pred.eurasia.bc$mean + pred.eurasia.bc$sd,5.75,  col=col.eurasia, lwd = 1.5)

axis(2, labels = c('Amazon', 'SE Asia', 'Africa','N America', 'FMEC', 'Eurasia'), cex.axis = 1.3, at = 1:6, col = NA, las = 1)
axis(1)
par(xpd = TRUE)
text(c(1,1,1), c(3.25, 3, 2.75), labels = c('observed', 'default parameters', 'bias corrected'), pos = 4, cex = 1)
dev.off()













