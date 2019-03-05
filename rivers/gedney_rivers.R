# gedney_rivers.R
# Plotting the Obidos Amazon river guage data
# obidos.dat and coastal-stns-Vol-monthly.updated-oct2007.nc are
# very nearly the same data, covering slightly different time frames.
# Here, we append the data and plot it monthly and annually, as 
# well as looking at the rolling 20 year means.


library(zoo)
library(ncdf4)
# This is the .dat file from Nic Gedney
#obidos = read.table('obidos.dat', skip = 6, header = TRUE, na.strings = '-999.0')
#plot(obidos$datem, obidos$rflowobs, type = 'l', xlim = c(1960, 2020))

obidos = read.table("data/obidos.dat", 
                    skip = 6, header = TRUE, na.strings = "-999.0")

obs.ts = ts(obidos$rflowobs, start = 1901, end = 2015, frequency = 12)
obs.end.ts = window(obs.ts, start = c(2006,12), end = end(obs.ts))

# this is the .netcdf file - not sure where it's from.
nc = nc_open(file = "data/coastal-stns-Vol-monthly.updated-oct2007.nc")
stn_name = ncvar_get(nc, "stn_name") #Obidos is number 1
flow = ncvar_get(nc, 'FLOW')
tme = ncvar_get(nc, 'time')

startyear = ncvar_get(nc, 'yrb')
endyear = ncvar_get(nc, 'yre')

obidos.flow = flow[1,]
obidos.flow.ts = ts(obidos.flow, start = c(1900,1), frequency = 12)

# Join the two time series into a single timeseries
obidos.merged = c(obidos.flow.ts[-length(obidos.flow.ts)],obs.end.ts )
obidos.merged.ts = ts(obidos.merged, start = start(obidos.flow.ts), frequency = 12)
# Aggregate to annual
obidos.merged.agg = aggregate(obidos.merged.ts , freq = 1, FUN = mean, na.rm = TRUE)

flow.mean = mean(obidos.merged.ts , na.rm = TRUE)
flow.mean.all.basin.sv = (flow.mean/fac) * (1/0.77)

ma <- function(x,n=20){filter(x,rep(1/n,n), sides=2)}
yr = 1900:2014
fac = 1e6

pdf(file = 'graphics/obidos_river_flow.pdf', width = 9, height = 5)
par(las = 1)
plot(obidos.merged.ts/fac, col = 'grey', axes = FALSE, xlim = c(1925, 2015),
     lwd = 1.5,
     xlab = 'year', ylab = 'River flow (Sv)',
     main = 'Streamflow from Obidos station')
lines(obidos.merged.agg/fac, col = 'black', lwd = 1.5)
lines(yr, ma(obidos.merged.agg)/fac, col = 'red', lwd = 1.5)
abline(h = flow.mean/fac, col = 'blue', lty = 'dashed')
axis(1)
axis(2)
legend('topleft', legend = c('Monthly', 'Annual', '20 year rolling mean'),
       col = c('grey', 'black', 'red'), bty = 'n', lty = 'solid')
dev.off()

# How much does the streamflow from Obidos change?
start.flow = mean(obidos.merged.agg[28:47])
end.flow = mean(obidos.merged.agg[94:114])


runoff.obs.abs.change.sv = (end.flow - start.flow) / fac
runoff.obs.perc.change = ((end.flow - start.flow) / start.flow) * 100

save(runoff.obs.perc.change,
	runoff.obs.abs.change.sv,
	flow.mean.all.basin.sv,
 	file = 'data/perc_change.Rdata')

# NetCDF stuff.
# Useful way of finding the names of the variables
attributes(nc$var)$names
# And the dimensions
attributes(nc$dim)$names
ncatt_get(nc, attributes(nc$var)$names[1])
test  <- ncvar_get(nc, attributes(nc$var)$names[1])

