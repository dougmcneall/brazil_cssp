# create_famous_data_v2.R
# create the data set for the "Famous Bias" paper, 
# and check the data from the previous paper.

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
famous_forest_fraction_v2_raw = read.table("famous_forest_fraction_v2.txt", header = TRUE)
famous_forest_fraction_runlist_v2 = read.table("famous_forest_fraction_runlist_v2.txt", header = FALSE)

famous_forest_fraction_v2 = run_reorder(rawdat = famous_forest_fraction_v2_raw , runlist = famous_forest_fraction_runlist_v2,
                   matchruns = famous_agg$FULL_ID)

# Surface Temperature
famous_surface_temperature_v2_raw = read.table("famous_surface_temperature_v2.txt", header = TRUE)
famous_surface_temperature_runlist_v2 = read.table("famous_surface_temperature_runlist_v2.txt", header = FALSE)

famous_surface_temperature_v2 = run_reorder(rawdat = famous_surface_temperature_v2_raw , runlist = famous_surface_temperature_runlist_v2,
                                            runidstart = 52, runidend = 58,
                                        matchruns = famous_agg$FULL_ID)

famous_precipitation_flux_v2_raw = read.table("famous_precipitation_flux_v2.txt", header = TRUE)
famous_precipitation_flux_runlist_v2 = read.table("famous_precipitation_flux_runlist_v2.txt", header = FALSE)

famous_precipitation_flux_v2 = run_reorder(rawdat = famous_precipitation_flux_v2_raw , 
                                           runlist = famous_precipitation_flux_runlist_v2,
                                           runidstart = 52, runidend = 58,
                                            matchruns = famous_agg$FULL_ID)




dev.new(width = 8, height = 10)
par(mfrow = c(3,2))
plot(famous_agg$AMAZ_MOD_FRAC, famous_forest_fraction_v2$AMAZ_MOD_FRAC)
abline(0,1)
plot(famous_agg$SEASIA_MOD_FRAC, famous_forest_fraction_v2$SEASIA_MOD_FRAC)
abline(0,1)
plot(famous_agg$CONGO_MOD_FRAC, famous_forest_fraction_v2$CONGO_MOD_FRAC)
abline(0,1)
plot(famous_agg$NAMERICA_MOD_FRAC, famous_forest_fraction_v2$NAMERICA_MOD_FRAC)
abline(0,1)
plot(famous_agg$GLOB_MOD_FRAC, famous_forest_fraction_v2$GLOB_MOD_FRAC)
abline(0,1)

dev.new(width = 8, height = 10)
par(mfrow = c(3,2))
plot(famous_agg$AMAZ_MOD_TEMP, famous_surface_temperature_v2$AMAZ_MOD_TEMP)
abline(0,1)
plot(famous_agg$SEASIA_MOD_TEMP, famous_surface_temperature_v2$SEASIA_MOD_TEMP)
abline(0,1)
plot(famous_agg$CONGO_MOD_TEMP, famous_surface_temperature_v2$CONGO_MOD_TEMP)
abline(0,1)
plot(famous_agg$GLOB_MOD_TEMP, famous_surface_temperature_v2$GLOB_MOD_TEMP)
abline(0,1)

dev.new(width = 8, height = 10)
par(mfrow = c(3,2))
plot(famous_agg$AMAZ_MOD_PRECIP, famous_precipitation_flux_v2$AMAZ_MOD_PRECIP)
abline(0,1)
plot(famous_agg$SEASIA_MOD_PRECIP, famous_precipitation_flux_v2$SEASIA_MOD_PRECIP)
abline(0,1)
plot(famous_agg$CONGO_MOD_PRECIP, famous_precipitation_flux_v2$CONGO_MOD_PRECIP)
abline(0,1)
plot(famous_agg$GLOB_MOD_PRECIP, famous_precipitation_flux_v2$GLOB_MOD_PRECIP)
abline(0,1)

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
                           GLOB_MOD_PRECIP = famous_precipitation_flux_v2$GLOB_MOD_PRECIP
                      )

dev.new(width = 12, height = 12)
pairs(famous_agg_v2, pch = 20, gap = 0)


















