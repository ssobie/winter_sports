##Script to merge the files for ERA Interim from 8 times daily to once daily

library(ncdf4)
library(PCICt)
library(abind)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R')

data.dir <- '/storage/data/projects/rci/data/winter_sports/ERA_INTERIM/'
var.name <- 'mx2t'

nc1 <- nc_open(paste0(data.dir,'era_interim_tasmax_0_3_6.nc'))
time1 <- netcdf.calendar(nc1)

nc2 <- nc_open(paste0(data.dir,'era_interim_tasmax_0_9_12.nc'))
time2 <- netcdf.calendar(nc2)

nc3 <- nc_open(paste0(data.dir,'era_interim_tasmax_12_3_6.nc'))
time3 <- netcdf.calendar(nc3)

nc4 <- nc_open(paste0(data.dir,'era_interim_tasmax_12_9_12.nc'))
time4 <- netcdf.calendar(nc4)

times <- list(time1,time2,time3,time4)
ncs <- list(nc1,nc2,nc3,nc4)

yrs <- 1979:2016
yr.daily <- vector(mode='list',length=4)
tasmax.daily <- c()
for (i in seq_along(yrs)) {
  print(yrs[i])
  for (n in 1:4) {    
    yr.sub <- grep(yrs[i],times[[n]])
    yr.data <- ncvar_get(ncs[[n]],var.name,start=c(1,1,yr.sub[1]),count=c(-1,-1,length(yr.sub)))-273
    dy.fac <- as.factor(format(times[[n]][yr.sub],'%Y-%m-%d'))
    yr.daily[[n]] <- aperm(apply(yr.data,c(1,2),function(x,fac){tapply(x,fac,max,na.rm=T)},dy.fac),c(2,3,1))   
  }
  dy.len <- length(levels(dy.fac))
  
  for (d in 1:dy.len) {

   yr.sub <- abind(yr.daily[[1]][,,d],
                   yr.daily[[2]][,,d],
                   yr.daily[[3]][,,d],
                   yr.daily[[4]][,,d],along=3)
   tasmax.daily <- abind(tasmax.daily,apply(yr.sub,c(1,2),max),along=3)
  }
}

