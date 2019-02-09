library(ncdf4)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

var.name <- 'pr'

##TPS scale
tps.enc <- nc_open(paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/ERA/',var.name,'_day_QDM_ERA_19790101-20181031.nc'))
tps.etx <- ncvar_get(tps.enc,var.name,start=c(11,11,1),count=c(1,1,-1))

tps.nnc <- nc_open(paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/NCEP2/',var.name,'_day_QDM_NCEP2_19790101-20181031.nc'))
tps.ntx <- ncvar_get(tps.nnc,var.name,start=c(11,11,1),count=c(1,1,-1))


anoms.enc <- nc_open(paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/ERA/',var.name,'_anoms_QDM_ERA_19790101-20181031.nc'))
anoms.etx <- ncvar_get(anoms.enc,var.name,start=c(11,11,1),count=c(1,1,-1))

anoms.nnc <- nc_open(paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/NCEP2/',var.name,'_anoms_QDM_NCEP2_19790101-20181031.nc'))
anoms.ntx <- ncvar_get(anoms.nnc,var.name,start=c(11,11,1),count=c(1,1,-1))

interp.enc <- nc_open(paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/ERA/',var.name,'_anoms_interp_ERA_19790101-20181031.nc'))
interp.etx <- ncvar_get(interp.enc,var.name,start=c(11,11,1),count=c(1,1,-1))

interp.nnc <- nc_open(paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/NCEP2/',var.name,'_anoms_interp_NCEP2_19790101-20181031.nc'))
interp.ntx <- ncvar_get(interp.nnc,var.name,start=c(11,11,1),count=c(1,1,-1))

enc <- nc_open(paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/ERA/',var.name,'_gcm_prism_ERA_19790101-20181031.nc'))
etx <- ncvar_get(enc,var.name,start=c(40,40,1),count=c(1,1,-1))

nnc <- nc_open(paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/NCEP2/',var.name,'_gcm_prism_NCEP2_19790101-20181031.nc'))
ntx <- ncvar_get(nnc,var.name,start=c(40,40,1),count=c(1,1,-1))


par(mfrow=c(4,1))
ix <- 10001:11000
plot(tps.ntx[ix],type='l')
lines(tps.etx[ix],col='red')

plot(anoms.ntx[ix],type='l')
lines(anoms.etx[ix],col='red')

plot(interp.ntx[ix],type='l')
lines(interp.etx[ix],col='red')

plot(ntx[ix],type='l')
lines(etx[ix],col='red')


nc_close(tps.enc)
nc_close(tps.nnc)

nc_close(anoms.enc)
nc_close(anoms.nnc)

nc_close(interp.enc)
nc_close(interp.nnc)
 
nc_close(enc)
nc_close(nnc)

