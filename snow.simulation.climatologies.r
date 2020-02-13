##Extract time series of ARI and Hazards for a shapefile region and
##save the output series as RData

library(ncdf4)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)


##--------------------------------------------------------------------

make_annual_total_file <- function(file,tmp.dir) {
   ann.file <- gsub(pattern='day',
                    replacement='annual_total',
                    file)
   work <- paste0('cdo -O yearsum ',tmp.dir,file,' ',tmp.dir,ann.file)
   system(work)
   return(ann.file)
}

make_annual_max_file <- function(file,tmp.dir) {
   ann.file <- gsub(pattern='day',
                    replacement='annual_maximum',
                    file)
   work <- paste0('cdo -O yearmax ',tmp.dir,file,' ',tmp.dir,ann.file)
   system(work)
   return(ann.file)
}

make_monthly_mean_file <- function(file,tmp.dir) {
   mon.file <- gsub(pattern='day',
                    replacement='monthly_mean',
                    file)
   work <- paste0('cdo -O monmean ',tmp.dir,file,' ',tmp.dir,mon.file)
   system(work)
   return(mon.file)
}

make_seasonal_mean_file <- function(file,tmp.dir) {
   seas.file <- gsub(pattern='day',
                    replacement='seasonal_mean',
                    file)
   work <- paste0('cdo -O seasmean ',tmp.dir,file,' ',tmp.dir,seas.file)
   system(work)
   return(seas.file)
}


make_monthly_total_file <- function(file,tmp.dir) {
   mon.file <- gsub(pattern='day',
                    replacement='monthly_total',
                    file)
   work <- paste0('cdo -O monsum ',tmp.dir,file,' ',tmp.dir,mon.file)
   system(work)
   return(mon.file)
}

make_seasonal_total_file <- function(file,tmp.dir) {
   seas.file <- gsub(pattern='day',
                    replacement='seasonal_total',
                    file)
   work <- paste0('cdo -O seassum ',tmp.dir,file,' ',tmp.dir,seas.file)
   system(work)
   return(seas.file)
}


time_subset <- function(file,bnds,tmp.dir) {

  time.write <- gsub(pattern='[0-9]{8}-[0-9]{8}',replacement=paste0(bnds,collapse='-'),file)
  work <- paste0('cdo -O seldate,',bnds[1],'-01-01T00:00:00,',bnds[2],'-12-31T23:59:59 ',
                  tmp.dir,file,' ',tmp.dir,time.write)
  print(work)
  system(work)

  return(time.write)
}

climatology <- function(file,fxn,tmp.dir) {
  ##clim.file <- gsub(pattern='vw',replacement="vw_climatology",
  ##                  file)
  clim.file <- gsub(pattern='BCCAQ2',replacement="climatology_BCCAQ2",
                    file)

  work <- paste0('cdo -O ',fxn,' ',tmp.dir,file,' ',tmp.dir,clim.file)
  print(work)
  system(work)
  return(clim.file)
}
 
make_ann_seas_mon_totals <- function(gcm.list,read.dir) {

   for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      print(gcm)

      ###gcm.dir <- paste0(read.dir,gcm,'/')
      ###file <- paste0(var.name,'_BCCAQ2-PRISM_',gcm,'_19500101-21001231.nc')
      gcm.dir <- read.dir
      file <- list.files(path=gcm.dir,pattern=paste0('swe_day_vw_',gcm,'_A2_run1_'))
      
      file.copy(from=paste0(gcm.dir,file),to=tmp.dir,overwrite=TRUE)

      ann.file <- make_annual_max_file(file,tmp.dir)    
      file.copy(from=paste0(tmp.dir,ann.file),to=gcm.dir,overwrite=TRUE)

       seas.file <- make_seasonal_mean_file(file,tmp.dir)     
       file.copy(from=paste0(tmp.dir,seas.file),to=gcm.dir,overwrite=TRUE)
       mon.file <- make_monthly_mean_file(file,tmp.dir)     
       file.copy(from=paste0(tmp.dir,mon.file),to=gcm.dir,overwrite=TRUE)

      ann.file <- make_annual_total_file(file,tmp.dir)    
      file.copy(from=paste0(tmp.dir,ann.file),to=gcm.dir,overwrite=TRUE)
      mon.file <- make_monthly_total_file(file,tmp.dir)    
      file.copy(from=paste0(tmp.dir,mon.file),to=gcm.dir,overwrite=TRUE)
      seas.file <- make_seasonal_total_file(file,tmp.dir)    
      file.copy(from=paste0(tmp.dir,seas.file),to=gcm.dir,overwrite=TRUE)
   }

}


##--------------------------------------------------------------------

##args <- commandArgs(trailingOnly=TRUE)
##for(i in 1:length(args)){
##    eval(parse(text=args[[i]]))
##}


tmp.dir <- paste0('/local_temp/ssobie/vic/')
###tmp.dir <- paste0('/local_temp/ssobie/snow/',gcm,'/')
if (!dir.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
}

var.name <- 'swe' ##varname ##'swe'
###gcm.list <- gcm


read.dir <- "/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/"
##
gcm.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')
##gcm.list <- 'GFDL-ESM2G'
intervals <- c('1971-2000','1981-2010','1980-2018','2011-2041','2041-2070','2071-2100')

##gcm.list <- c('CCSM3','CGCM3','CSIRO35','ECHAM5','GFDL2.1','HadCM','HadGEM1','MIROC3.2')
##intervals <- c('1971-2000','1981-2010','2041-2070')
##read.dir <- '/storage/data/projects/rci/data/winter_sports/VIC/'

file.type <- c('annual','monthly','seasonal')
file.fxn <- c('timmean','ymonmean','yseasmean')

file.type <- 'seasonal_mean'
file.fxn <- 'yseasmean'

##file.type <- 'annual_maximum'
##file.fxn <- 'timmean'

##make_ann_seas_mon_totals(gcm.list,read.dir)

if (1==1) {

for (g in seq_along(gcm.list)) {
   gcm <- gcm.list[g]
   print(gcm)
   gcm.dir <- paste0(read.dir,gcm,'/')

   for (f in seq_along(file.type)) {
      file <- paste0(var.name,'_',file.type[f],'_BCCAQ2-PRISM_',gcm,'_19500101-21001231.nc')
      ###file <- list.files(path=read.dir,pattern=paste0(var.name,'_',file.type[f],'_vw_',gcm))      
      file.copy(from=paste0(gcm.dir,file),to=tmp.dir,overwrite=TRUE)  

      for (i in seq_along(intervals)) {
         interval <- intervals[i]
         print(interval)
         bnds <- strsplit(interval,'-')[[1]]
         time.file <- time_subset(file,bnds,tmp.dir)
         clim.file <- climatology(time.file,file.fxn[f],tmp.dir)

         file.copy(from=paste0(tmp.dir,clim.file),to=gcm.dir,overwrite=TRUE)
         file.remove(paste0(tmp.dir,time.file))
         file.remove(paste0(tmp.dir,clim.file))

      }
   }
   file.remove(paste0(tmp.dir,file))

}

}