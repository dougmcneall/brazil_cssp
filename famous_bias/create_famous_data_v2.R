# create_famous_data_v2.R
# create the data set for the "Famous Bias" paper, 
# and check the data from the previous paper.

datadir = ('/home/h01/hadda/famous/data/')
#load('/home/h01/hadda/brazilCSSP/code/brazil_cssp/famous_bias/famous_forest_fraction.RData')
load('/home/h01/hadda/brazilCSSP/code/brazil_cssp/famous_bias/famous_agg.Rdata')


run_reorder = function(rawdat,runlist, matchruns, runidstart = 47, runidend = 53){
  # Reorder ensemble to match existing input parameter data
  runids = rep(NA, nrow(runlist))
  # Extract just the runnames
  for (i in 1:nrow(runlist)){
    runids[i] = substring(as.character(runlist[i,]), runidstart, runidend)
  }
  
  trans.ix = match(as.character(matchruns), runids)
  out = rawdat[trans.ix, ]
  
}

# Forest Fraction
famous_forest_fraction_v2_raw = read.table(paste0(datadir,"famous_forest_fraction_v2.txt"), header = TRUE)
famous_forest_fraction_runlist_v2 = read.table(paste0(datadir,"famous_forest_fraction_runlist_v2.txt"), header = FALSE)

famous_forest_fraction_v2 = run_reorder(rawdat = famous_forest_fraction_v2_raw , runlist = famous_forest_fraction_runlist_v2,
                   matchruns = famous_agg$FULL_ID)

# Surface Temperature
famous_surface_temperature_v2_raw = read.table(paste0(datadir,"famous_surface_temperature_v2.txt"), header = TRUE)
famous_surface_temperature_runlist_v2 = read.table(paste0(datadir,"famous_surface_temperature_runlist_v2.txt"), header = FALSE)

famous_surface_temperature_v2 = run_reorder(rawdat = famous_surface_temperature_v2_raw , runlist = famous_surface_temperature_runlist_v2,
                                            runidstart = 52, runidend = 58,
                                        matchruns = famous_agg$FULL_ID)

famous_precipitation_flux_v2_raw = read.table(paste0(datadir, "famous_precipitation_flux_v2.txt"), header = TRUE)
famous_precipitation_flux_runlist_v2 = read.table(paste0(datadir,"famous_precipitation_flux_runlist_v2.txt"), header = FALSE)

famous_precipitation_flux_v2 = run_reorder(rawdat = famous_precipitation_flux_v2_raw , 
                                           runlist = famous_precipitation_flux_runlist_v2,
                                           runidstart = 52, runidend = 58,
                                            matchruns = famous_agg$FULL_ID)




# dev.new(width = 8, height = 10)
# par(mfrow = c(3,2))
# plot(famous_agg$AMAZ_MOD_FRAC, famous_forest_fraction_v2$AMAZ_MOD_FRAC)
# abline(0,1)
# plot(famous_agg$SEASIA_MOD_FRAC, famous_forest_fraction_v2$SEASIA_MOD_FRAC)
# abline(0,1)
# plot(famous_agg$CONGO_MOD_FRAC, famous_forest_fraction_v2$CONGO_MOD_FRAC)
# abline(0,1)
# plot(famous_agg$NAMERICA_MOD_FRAC, famous_forest_fraction_v2$NAMERICA_MOD_FRAC)
# abline(0,1)
# plot(famous_agg$GLOB_MOD_FRAC, famous_forest_fraction_v2$GLOB_MOD_FRAC)
# abline(0,1)

# dev.new(width = 8, height = 10)
# par(mfrow = c(3,2))
# plot(famous_agg$AMAZ_MOD_TEMP, famous_surface_temperature_v2$AMAZ_MOD_TEMP)
# abline(0,1)
# plot(famous_agg$SEASIA_MOD_TEMP, famous_surface_temperature_v2$SEASIA_MOD_TEMP)
# abline(0,1)
# plot(famous_agg$CONGO_MOD_TEMP, famous_surface_temperature_v2$CONGO_MOD_TEMP)
# abline(0,1)
# plot(famous_agg$GLOB_MOD_TEMP, famous_surface_temperature_v2$GLOB_MOD_TEMP)
# abline(0,1)
# 
# dev.new(width = 8, height = 10)
# par(mfrow = c(3,2))
# plot(famous_agg$AMAZ_MOD_PRECIP, famous_precipitation_flux_v2$AMAZ_MOD_PRECIP)
# abline(0,1)
# plot(famous_agg$SEASIA_MOD_PRECIP, famous_precipitation_flux_v2$SEASIA_MOD_PRECIP)
# abline(0,1)
# plot(famous_agg$CONGO_MOD_PRECIP, famous_precipitation_flux_v2$CONGO_MOD_PRECIP)
# abline(0,1)
# plot(famous_agg$GLOB_MOD_PRECIP, famous_precipitation_flux_v2$GLOB_MOD_PRECIP)
# abline(0,1)

famous_agg_v2 = data.frame(famous_agg[, 1:12], 
                           FMEC_MOD_FRAC = famous_forest_fraction_v2$FMEC_MOD_FRAC,
                           FMECNL_MOD_FRAC = famous_forest_fraction_v2$FMECNL_MOD_FRAC,
                           FMECBL_MOD_FRAC = famous_forest_fraction_v2$FMECBL_MOD_FRAC,
                           EURASIA_MOD_FRAC = famous_forest_fraction_v2$EURASIA_MOD_FRAC,
                           GLOB_MOD_FRAC = famous_forest_fraction_v2$GLOB_MOD_FRAC,
                           AMAZ_MOD_TEMP = famous_surface_temperature_v2$AMAZ_MOD_TEMP,
                           SEASIA_MOD_TEMP = famous_surface_temperature_v2$SEASIA_MOD_TEMP,
                           CONGO_MOD_TEMP = famous_surface_temperature_v2$CONGO_MOD_TEMP,
                           NAMERICA_MOD_TEMP =famous_surface_temperature_v2$NAMERICA_MOD_TEMP,
                           FMEC_MOD_TEMP = famous_surface_temperature_v2$FMEC_MOD_TEMP,
                           EURASIA_MOD_TEMP = famous_surface_temperature_v2$EURASIA_MOD_TEMP,
                           GLOB_MOD_TEMP = famous_surface_temperature_v2$GLOB_MOD_TEMP,
                           AMAZ_MOD_PRECIP = famous_precipitation_flux_v2$AMAZ_MOD_PRECIP,
                           SEASIA_MOD_PRECIP = famous_precipitation_flux_v2$SEASIA_MOD_PRECIP,
                           CONGO_MOD_PRECIP = famous_precipitation_flux_v2$CONGO_MOD_PRECIP,
                           NAMERICA_MOD_PRECIP = famous_precipitation_flux_v2$NAMERICA_MOD_PRECIP,
                           FMEC_MOD_PRECIP = famous_precipitation_flux_v2$FMEC_MOD_PRECIP,
                           EURASIA_MOD_PRECIP = famous_precipitation_flux_v2$EURASIA_MOD_PRECIP,
                           GLOB_MOD_PRECIP = famous_precipitation_flux_v2$GLOB_MOD_PRECIP
                      )

temps_obs_v2 = read.table(paste0(datadir,"temps_obs_v2.txt"), header = TRUE)
precips_obs_v2 = read.table(paste0(datadir,"precips_obs_v2.txt"), header = TRUE)
forest_fraction_obs_v2 = read.table(paste0(datadir,"forest_fraction_obs_v2.txt"), header = TRUE)


#dev.new(width = 12, height = 12)
#pairs(famous_agg_v2, pch = 20, gap = 0)

# Add Broadleaf and needleleaf fraction, and temperature and precipitation maps.
climatedir = '/home/h01/hadda/famous/data/temp_precip_famous/'

# ----------------------------------------------------
library(ncdf4)

# Helper functions for working with FAMOUS data
f <- function(s){
  strsplit(s, split = "a.pt")[[1]][1]
}

open.frac.field = function(fn, var, lev){
  # helper function to load a map of var from nc file
  nc = nc_open(fn)
  nc.var = ncvar_get(nc, var)
  nc.var.lev = nc.var[ , ,lev]
  nc_close(nc)
  nc.var.lev
}


open.field = function(fn, var){
  # helper function to load a map of var from nc file
  nc = nc_open(fn)
  nc.var = ncvar_get(nc, var)
  nc_close(nc)
  nc.var
}

load.frac.ens = function(fn.list, var, lev){
  # open all nc files in a list, vectorise, and concatenate to
  # an ensemble matrix, each row is a map
  field.list = lapply(fn.list, FUN=open.frac.field, var=var, lev=lev)
  out = t(sapply(field.list,cbind)) # should do by columns
  out
}



load.ens = function(fn.list, var){
  # open all nc files in a list, vectorise, and concatenate to
  # an ensemble matrix, each row is a map
  field.list = lapply(fn.list, FUN=open.field, var=var)
  out = t(sapply(field.list,cbind)) # should do by columns
  out
}


remap.famous = function(dat,longs,lats, shift = FALSE){
  # reshape a map in vector form so that image() like functions
  # will plot it correctly
  mat = matrix(dat, nrow=length(longs), ncol=length(lats))[ ,length(lats):1]
  if(shift){
    block1.ix = which(longs < shift)
    block2.ix = which(longs > shift)
    mat.shift = rbind(mat[ block2.ix, ], mat[block1.ix, ]) 
    out = mat.shift
  }
  else{
    out = mat
  }
  out
}

match.list <- as.character(famous_agg_v2$FULL_ID)
file.ix <- pmatch(match.list,dir(climatedir))
# filename list
fn.list <- as.list(paste(climatedir, dir(climatedir), sep = ""))

# test open a file
nc <- nc_open(fn.list[[1]])
temp.nc <- ncvar_get(nc, "temp_mm_srf")
lats <- ncvar_get(nc,"latitude")
longs <- ncvar_get(nc,"longitude")

# lev = 1 is broadleaf forest (2 would be needleleaf)
temp.ens.v2 = load.ens(fn.list=fn.list[file.ix], var='temp_mm_srf')
precip.ens.v2 = load.ens(fn.list=fn.list[file.ix], var='precip_mm_srf')

datadir = ('/home/h01/hadda/famous/data/famous_rename/')

# broadleaf and needleleaf forest fraction
file.ix <- pmatch(match.list,dir(datadir))
fn.list <- as.list(paste(datadir, dir(datadir), sep = ""))

bl.frac.ens.v2 = load.frac.ens(fn.list=fn.list[file.ix], var='fracPFTs_mm_srf', lev=1)
nl.frac.ens.v2 = load.frac.ens(fn.list=fn.list[file.ix], var='fracPFTs_mm_srf', lev=2)


# find range of data for plotting
#zl <- c(0,1)

#library(fields)
# plot first ensemble member
#image.plot(longs, rev(lats), remap.famous(bl.frac.ens[1,], longs, lats), col=yg)

#save(file='famous_forest_fraction.RData',full_frac,obs,lats,longs,bl.frac.ens, nl.frac.ens, X.standard)
savedir = '/net/home/h01/hadda/brazilCSSP/code/brazil_cssp/famous_bias/'
save(famous_agg_v2, temps_obs_v2,precips_obs_v2,forest_fraction_obs_v2,temp.ens.v2, precip.ens.v2,bl.frac.ens.v2,nl.frac.ens.v2, file = paste0(savedir, 'famous_forest_fraction_v2.RData'))


















