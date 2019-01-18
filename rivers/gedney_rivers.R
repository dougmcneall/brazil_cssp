# gedney_rivers.R
# Plotting the Obidos Amazon river guage data

# This is the .dat file from Nic Gedney
obidos = read.table('obidos.dat', skip = 6, header = TRUE, na.strings = '-999.0')
plot(obidos$datem, obidos$rflowobs, type = 'l', xlim = c(1960, 2020))


# this is the .netcdf file - not sure where it's from.
library(ncdf4)
library(zoo)

nc = nc_open(file = "data/coastal-stns-Vol-monthly.updated-oct2007.nc")

stn_name = ncvar_get(nc, "stn_name") #Obidos is number 1
flow = ncvar_get(nc, 'FLOW')
tme = ncvar_get(nc, 'time')

startyear = ncvar_get(nc, 'yrb')
endyear = ncvar_get(nc, 'yre')

obidos.flow = flow[1,] / 1e3

obidos.flow.ts = ts(obidos.flow, start = c(1900,1), frequency = 12)
obidos.flow.ann = aggregate(obidos.flow.ts , freq = 1, FUN = mean)


yr = seq(from = 1900+(1/12), to = 2007, by = 1/12) 
plot(obidos.flow.ts)

obidos.flow.rm = rollmean(obidos.flow, k = 30)

ma <- function(x,n=(30*12)){filter(x,rep(1/n,n), sides=2)}


dev.new(width = 10, height = 5)
plot(obidos.flow.ts, col = 'grey')
lines(obidos.flow.ann, type = 'l', col = 'black')

lines(yr, ma(obidos.flow), col = 'red')



# Useful way of finding the names of the variables
attributes(nc$var)$names
# And the dimensions
attributes(nc$dim)$names

ncatt_get(nc, attributes(nc$var)$time)


ncatt_get(nc, attributes(nc$var)$names[1])


test  <- ncvar_get(nc, attributes(nc$var)$names[1])

