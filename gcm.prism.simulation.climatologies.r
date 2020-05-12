##Extract time series of ARI and Hazards for a shapefile region and
##save the output series as RData

library(ncdf4)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)


##--------------------------------------------------------------------
make_annual_mean_file <- function(var.name,file,tmp.dir) {
   ann.file <- gsub(pattern=var.name,
                    replacement=paste0(var.name,'_annual_mean'),
                    file)
   work <- paste0('cdo -O yearmean ',tmp.dir,file,' ',tmp.dir,ann.file)
   system(work)
   return(ann.file)
}

make_annual_total_file <- function(var.name,file,tmp.dir) {
   ann.file <- gsub(pattern=var.name,
                    replacement=paste0(var.name,'_annual_total'),
                    file)
   work <- paste0('cdo -O yearsum ',tmp.dir,file,' ',tmp.dir,ann.file)
   system(work)
   return(ann.file)
}

make_monthly_mean_file <- function(var.name,file,tmp.dir) {
   mon.file <- gsub(pattern=var.name,
                    replacement=paste0(var.name,'_monthly_mean'),
                    file)
   work <- paste0('cdo -O monmean ',tmp.dir,file,' ',tmp.dir,mon.file)
   system(work)
   return(mon.file)
}

make_seasonal_mean_file <- function(var.name,file,tmp.dir) {
   seas.file <- gsub(pattern=var.name,
                    replacement=paste0(var.name,'_seasonal_mean'),
                    file)
   work <- paste0('cdo -O seasmean ',tmp.dir,file,' ',tmp.dir,seas.file)
   system(work)
   return(seas.file)
}


time_subset <- function(file,bnds,tmp.dir) {

  ###time.write <- gsub(pattern='[0-9]{8}-[0-9]{8}',replacement=paste0(bnds,collapse='-'),file)
  time.write <- gsub(pattern='[0-9]{4}-[0-9]{4}',replacement=paste0(bnds,collapse='-'),file)
  work <- paste0('cdo -O seldate,',bnds[1],'-01-01T00:00:00,',bnds[2],'-12-31T23:59:59 ',
                  tmp.dir,file,' ',tmp.dir,time.write)
  print(work)
  system(work)

  return(time.write)
}

climatology <- function(file,fxn,tmp.dir) {
  ##clim.file <- gsub(pattern='vw',replacement="vw_climatology",
  ##                  file)
  clim.file <- gsub(pattern='gcm_prism',replacement="gcm_prism_climatology",
                    file)

  work <- paste0('cdo -O ',fxn,' ',tmp.dir,file,' ',tmp.dir,clim.file)
  print(work)
  system(work)
  return(clim.file)
}
 
make_ann_seas_mon_totals <- function(var.name,gcm,reanalysis,read.dir,write.dir) {

      file <- list.files(read.dir,pattern=paste0(var.name,'_gcm_prism_',gcm,'_VW'))
      file.copy(from=paste0(read.dir,file),to=tmp.dir,overwrite=TRUE)

      ann.file <- make_annual_mean_file(var.name,file,tmp.dir)    
      file.copy(from=paste0(tmp.dir,ann.file),to=write.dir,overwrite=TRUE)

      seas.file <- make_seasonal_mean_file(var.name,file,tmp.dir)     
      file.copy(from=paste0(tmp.dir,seas.file),to=write.dir,overwrite=TRUE)

      mon.file <- make_monthly_mean_file(var.name,file,tmp.dir)     
      file.copy(from=paste0(tmp.dir,mon.file),to=write.dir,overwrite=TRUE)

   ##rv <- list(apr=apr.file,ann=ann.file,seas=seas.file,mon=mon.file)
   rv <- list(ann=ann.file,seas=seas.file,mon=mon.file)
   return(rv)
}


##--------------------------------------------------------------------

##args <- commandArgs(trailingOnly=TRUE)
##for(i in 1:length(args)){
##    eval(parse(text=args[[i]]))
##}


var.name <- 'tasmax' ##varname ##'swe'
##gcm <- 'ACCESS1-0'

intervals <- c('1950-2012','1971-2000','1981-2010','1981-2018','2011-2040','2041-2070','2071-2100')

tmp.dir <- paste0('/local_temp/ssobie/gcm_prism_tasmax_clims/')
###tmp.dir <- paste0('/local_temp/ssobie/snow/',gcm,'/')
if (!dir.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
}
##'ACCESS1-0',
gcm.list <- c('CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')

file.type <- list(ann='annual_mean',seas='seasonal_mean') ##mon='monthly_mean',
file.fxn <- list(ann='timmean',seas='yseasmean') ##mon='ymonmean',

for (gcm in gcm.list) {
   proj.dir <- paste0("/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/",gcm,"/")
   write.dir <- paste0("/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/",gcm,"/tas_climatology/")
   if (!dir.exists(write.dir)) {
      dir.create(write.dir,recursive=TRUE)
   }

   avg.files <- make_ann_seas_mon_totals(var.name,gcm,reanalysis,proj.dir,write.dir)

   print(avg.files)
if (1==1) {

   print(gcm)
   gcm.dir <- proj.dir ##paste0(read.dir,gcm,'/')
   
   for (f in seq_along(names(file.type))) {
      type <- names(file.type)[f]
      file <- avg.files[[type]]
      print(file)

      for (i in seq_along(intervals)) {
         interval <- intervals[i]
         print(interval)
         bnds <- strsplit(interval,'-')[[1]]
         time.file <- time_subset(file,bnds,tmp.dir)
         clim.file <- climatology(time.file,file.fxn[[type]],tmp.dir)

         file.copy(from=paste0(tmp.dir,clim.file),to=write.dir,overwrite=TRUE)
         file.remove(paste0(tmp.dir,time.file))
         file.remove(paste0(tmp.dir,clim.file))
      }
   }
   file <- list.files(tmp.dir,pattern=paste0(var.name,'_gcm_prism_',gcm,'_VW'))
   file.remove(paste0(tmp.dir,file))
   file.remove(paste0(tmp.dir,avg.files))
}

}