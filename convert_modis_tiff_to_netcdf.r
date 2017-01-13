##Script to convert the MODIS tiff files into netcdf for easier use

library(raster)
library(ncdf4)


##Function to create month file

make.month.file <- function(prefix,
                            proj.dir,write.dir,
                            start.date,end.date) {

  write.file <- paste0(write.dir,'snc.modis.terra.',gsub('-','',start.date),'-',gsub('-','',end.date),'.nc')

  dates <- seq(from=as.Date(start.date),by='day',to=as.Date(end.date))
  dates.under <- gsub('-','_',dates)
  files <- paste0(prefix,dates.under,'.clipped.tif')
  all.files <- list.files(path=proj.dir,pattern='.tif')
  yr <- format(as.Date(start.date),'%Y')
  mn <- format(as.Date(end.date),'%m')

  year.files <- all.files[grep(paste0(prefix,yr,'_',mn,'*_'),all.files)]

  r.list <- vector(mode='list',length=length(files))

  flags <- files %in% year.files
  flags.ix <- which(flags)
  miss.flag <- !flags
  mlen <- sum(miss.flag)

  for (i in seq_along(flags.ix)) {
    ix <- flags.ix[i]
    r.list[[ix]] <- raster(paste0(proj.dir,year.files[i]))
  }

  filling.file <- paste0(proj.dir,prefix,'2001_01_01.clipped.tif')
  if (mlen > 0) {
    print('Some missing days')
    miss.ix <- which(miss.flag)
    miss.raster <- raster(filling.file)*NA  
    for (j in miss.ix) {      
       print(files[j])              
       r.list[[j]] <- miss.raster
    }
  }
  ##if (mn=='12') {
  ##  browser()
  ##}
  r.stack <- stack(r.list)
  r.array <- as.array(r.stack)
  r.perm <- aperm(r.array,c(2,1,3))
  d2 <- dim(r.perm)[2]
  r.flip <- r.perm[,d2:1,]

  ##Coordinates
  r1 <- coordinates(r.list[[1]])
  lon <- sort(unique(r1[,1]))
  lat <- sort(unique(r1[,2]))

  time.units <- 'days since 1950-01-01'
  time.vals <- dates - as.Date('1950-01-01')
  time.calendar <- 'gregorian'
  ##--------------------------------------------------------------
  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', time.units, as.numeric(time.vals),
                        unlim=TRUE, calendar=time.calendar)

  var.geog <- ncvar_def('snc', units='bands', dim=list(x.geog, y.geog, t.geog),
                        missval=-32768)
  file.nc <- nc_create(write.file, var.geog)

  ncvar_put(file.nc,varid='snc',vals=r.flip,
                    start=c(1,1,1),count=c(-1,-1,-1))

  nc_close(file.nc)

}

proj.dir <- '/storage/data/projects/rci/data/winter_sports/MODIS_AQUA/'
write.dir <- '/storage/data/projects/rci/data/winter_sports/MODIS_AQUA_NETCDF/'
terra.prefix <- 'MOD10A1.'
aqua.prefix <-  'MYD10A1.'

prefix <- aqua.prefix
yr <- '2015'
year.dates <- seq(from=as.Date(paste0(yr,'-01-01')),by='day',to=as.Date(paste0(yr,'-12-31')))
months <- sprintf('%02d',1:12)

for (m in seq_along(months)) {
  print(m)
  m.ix <- grep(paste0('-',months[m],'-'),year.dates)
  start.date <- head(year.dates[m.ix],1)
  end.date   <- tail(year.dates[m.ix],1)

  make.month.file(prefix,
                  proj.dir,write.dir,
                  start.date,end.date)
}