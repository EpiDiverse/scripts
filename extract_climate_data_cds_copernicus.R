#Author: Iris Sammarco
#Date: 24/03/2020
#Description: extract climatic variables from cds copernicus https://cds.climate.copernicus.eu/cdsapp#!/home from the GPS location of interest

library(ncdf4)
library(chron)

ncin <- nc_open("tg_0.1deg_day_2019_grid_ensmean.nc") #load the climatic dataset in nc format downloaded by cds copernicus
ncin
lat <- ncvar_get(ncin, "latitude")
lon <- ncvar_get(ncin, "longitude")
time <- ncvar_get(ncin, "time")
tunits <- ncatt_get(ncin, "time", "units")
tunits
nt <- dim(time)
nt
tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr) [3], "-")
tmonth = as.integer(unlist(tdstr) [2])
tday = as.integer(unlist(tdstr) [3])
tyear = as.integer(unlist(tdstr) [1])
tnew <- chron(time, origin=c(tmonth,tday,tyear))
tnew

## Extract the climatic variable of interest for your GPS location
tmp.array <- ncvar_get(ncin, "tg") #tg is the name of the variable, change it according to the variable downloaded
tmp.array
dim(tmp.array)
nc_close(ncin) #I can close the file and work with the variables
tmp.array[1,1,1] #it's a 3 variables array
lon[1]
lat[1]
time[1]
tnew[1]
which.min(abs(lat-46.47120)) #I need to find the index in the latitude variable which is closest to the latitude of interest. Replace 46.47120 with the latitude of interest
lat[215] #The funciton above returns 215 as index, I can check if the latitude at this index is similar to the one I need
which.min(abs(lon-11.34305)) #Same as above, find the the index in the longitude variable closest to the longitude of interest
lon[364]
tmp.slice <- tmp.array[364,215,] # Extract from the dataset longitute, latitude and the time period of interest [lon,lat,all the time period]
tmp.slice #all the temperatures for your location
plot(tnew, tmp.slice)
data_table <- data.frame(tnew, tmp.slice)
write.table(data_table, "climate.csv", row.names = FALSE, col.names = TRUE, sep = ",")