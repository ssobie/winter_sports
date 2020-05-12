##Script to calculate the ARI (antecedent rainfall index)

ptm <- proc.time()

library(PCICt)
library(ncdf4)

##------------------------------------------------------------
##Pull the time series
time_series_from_nc <- function(nc) {

  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units

  time.start <- as.Date(strsplit(time.units, ' ')[[1]][3])
  origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                     cal=time.calendar)
  time.values <- ncvar_get(nc,'time')
  time.series <- origin + time.values*86400

  time.min <- as.character(format(head(time.series,1),'%Y%m%d'))
  time.max <- as.character(format(tail(time.series,1),'%Y%m%d'))

  rv <- list(units=time.units,
             values=time.values,
             series=time.series,
             calendar=time.calendar,
             tmax=time.max,
             tmin=time.min)
  return(rv)
}

##------------------------------------------------------------

get_var_units <- function(var.name) {
   rv <- switch(var.name,
                pr='kg m-2 d-1',
                tasmax='degC',
                tasmin='degC',
                swe='mm',
                snowdepth='m')
   return(rv)
}


##------------------------------------------------------------

insert_split_file_data <- function(var.name,write.nc,tmp.dir,split.dir) {

   splits <- c(rep(1,13),rep(2:32,each=10))
   
   for (i in 1:32) {
      ix <- which(splits %in% i)
      lat.ix <- which(splits %in% i)
      lat.st <- lat.ix[1]
      lat.en <- tail(lat.ix,1)
      lat.cnt <- diff(range(lat.ix))+1

      ix.st <- ix[1]
      ix.cnt <- diff(range(ix))+1
     
      split.file <- list.files(path=split.dir, pattern=paste0(var.name,'_L',sprintf('%02d',i),'_',sprintf('%02d',lat.st)))
      print(split.file)
      
      Sys.sleep(5)

      split.nc <- nc_open(paste0(split.dir,split.file))
      ntime <- length(ncvar_get(split.nc,'time'))
      nlon <- length(ncvar_get(split.nc,'lon'))
      nlat <- length(ncvar_get(split.nc,'lat'))
      split.insert <- array(NA,c(nlon,nlat,ntime))
      for (j in 1:nlat) {      
         split.sub <- ncvar_get(split.nc,var.name,start=c(1,j,1),count=c(-1,1,-1))
         split.insert[,j,] <- split.sub
         rm(split.sub)
      }
      ncvar_put(write.nc,var.name,split.insert,start=c(1,ix.st,1),count=c(nlon,ix.cnt,ntime))
      rm(split.insert)
      nc_close(split.nc)
      file.remove(paste0(split.dir,split.file))
      print('Copied back and cleaned up') 
      gc()
   }
}

##------------------------------------------------------------

##************************************************************

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

##varname <- 'snowdepth'
##gcm <- 'ERA5'
##reanalysis <- 'PNWNAmet'

##


tmpdir <- paste0(tmpdir,'/') ##'/local_temp/ssobie/snow_assembly/' ##
model <- gcm
var.name <- varname

if (model=='PNWNAmet' | model=='ERA5') {
   empty.dir <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/',model,'/')        
   split.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/snow_model/',        
                        'calibrated_',model,'_',reanalysis,'_prism_tps/',varname,'_',model,'_split10/')
   write.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/snow_model/',        
                        'calibrated_',model,'_',reanalysis,'_prism_tps/')
   if (model=='ERA5') {interval <- '1980-2018'}
   if (model=='PNWNAmet') {interval <- '1945-2012'}

} else {
   interval <- '19500101-21001231'
   empty.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/snow_model/templates/',model,'/')        
   split.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/snow_model/',        
                        'calibrated_',model,'_',reanalysis,'_prism_tps/',varname,'_',model,'_split10/')
   write.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/snow_model/',        
                        'calibrated_',model,'_',reanalysis,'_prism_tps/')
}


tmp.dir <- paste0(tmpdir,'combine_',model,'_',varname,'_tmp/')
if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=T)
}
 
file.copy(from=split.dir,to=tmp.dir,recursive=TRUE)
print('Done copying DS directory')
print('Done Copying')

split.dir <- paste0(tmp.dir,var.name,'_',model,'_split10/')

##-----------------------------------------------
##Copy empty snow file
empty.file <- paste0(var.name,'_BCCAQ2-PRISM_',model,'_',interval,'.nc')
write.file <- paste0(var.name,'_BCCAQ2-PRISM_',model,'_',reanalysis,'_',interval,'.nc')
file.copy(from=paste0(empty.dir,empty.file),to=paste0(tmp.dir,write.file))

write.nc <- nc_open(paste0(tmp.dir,write.file),write=TRUE)

print('Split file data insertion') 
insert_split_file_data(var.name,write.nc,tmp.dir,split.dir)

print('Copying to write dir')
file.copy(from=paste0(tmp.dir,write.file),to=write.dir,overwrite=TRUE)

file.remove(paste0(tmp.dir,write.file))
rm.files <- list.files(path=tmp.dir,recursive=T,full.name=T)
file.remove(rm.files)
print('Copied back and cleaned up')

print('Elapsed time')
print(proc.time()-ptm)

