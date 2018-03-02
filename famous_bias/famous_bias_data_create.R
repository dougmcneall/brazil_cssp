# famous_bias_data_create.R
# R script to illustrate fixing a bias in famous by using
# temperature and precipitation as inputs to the emulator.
# dougmcneall@gmail.com

library(DiceKriging)
library(RColorBrewer)
library(MASS)
library(fields)
library(parallel)

load('famous_forest_fraction.RData')

source('https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R')
source('https://raw.githubusercontent.com/dougmcneall/packages-git/master/imptools.R')
source('https://raw.githubusercontent.com/dougmcneall/packages-git/master/vistools.R')

load('famous_params.RData')
load('famous_params_beta.RData')

# Normalise the input space
X <- params_beta[ ,4:10]
X.norm <- normalize(X)
# standard set of parameters
X.standard <- c(0.875, 3, 0.03, 0.25, 36, 2, 0.5)
X.stan.norm <- normalize(matrix(X.standard, nrow = 1), wrt = X)

colnames(X.stan.norm) <- colnames(X.norm)

ndims <- ncol(X)
nens <- nrow(X)

paths <- readLines('filepaths_famous.txt')
match.list <- as.character(params$FULL_ID)
modsFromFilelist <- sub("/net/home/h01/hadda/famous/data/temp_precip_famous/", "", paths)
file.ix <- pmatch(match.list, modsFromFilelist)

temps <-cbind(match.list, read.table('temps_forests.txt', header = FALSE, skip = 1)[file.ix , ])
colnames(temps) <- c('RUN', 'GLOB_MOD_TEMP', 'AMAZ_MOD_TEMP', 
                     'CONGO_MOD_TEMP', 'SEASIA_MOD_TEMP')

precips <- cbind(match.list, read.table('precips_forests.txt', header = FALSE, skip = 1)[file.ix , ])
colnames(precips) <- c('RUN', 'GLOB_MOD_PRECIP', 'AMAZ_MOD_PRECIP',
                       'CONGO_MOD_PRECIP', 'SEASIA_MOD_PRECIP')

precips_obs <- read.table('precips_obs.txt', header = FALSE, skip = 1)
colnames(precips_obs) <- c('GLOB_OBS_PRECIP', 'AMAZ_OBS_PRECIP', 'CONGO_OBS_PRECIP', 'SEASIA_OBS_PRECIP')

temps_obs <- read.table('temps_obs.txt', header = FALSE, skip = 1)
colnames(temps_obs) <- c('GLOB_OBS_TEMP', 'AMAZ_OBS_TEMP', 'CONGO_OBS_TEMP', 'SEASIA_OBS_TEMP')

famous_agg = cbind(full_frac,temps[, 2:5], precips[,2:5])

save(file = 'famous_agg.Rdata', famous_agg, temps_obs, precips_obs)





