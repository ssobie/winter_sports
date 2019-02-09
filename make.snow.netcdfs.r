##Script to create SWE and Snowdepth netcdf files for the Walter Snow Model output

library(ncdf4)
library(PCICt)

##----------------------------------------------------------------------
##Global Attributes
get_global_atts <- function(drive.centre,drive.centre.name,gcm,rcp,run) {

  global.atts <- list(institution="Pacific Climate Impacts Consortium (PCIC), Victoria, BC, www.pacificclimate.org",
                   contact="Pacific Climate Impacts Consortium",
                   Conventions="CF-1.4",
                   institute_id ="PCIC",
                   domain='VanWhistler',
                   creation_date=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%Z'),
                   frequency="day",
                   product="snow model output",
                   modeling_realm="sfc",
                   project_id='ERA',
                   references=paste0("Walter et al., 2005 M.T. Walter, E.S. Brooks, D.K. McCool, L.G. King, M. Molnau, J. Boll, Process-based snowmelt modeling: does it require more input data than temperature-index modeling J. Hydrol., 300 (2005), pp. 65-75"),
                   downscaling_method="Quantile Delta Mapping + Climate Imprint",
                   downscaling_method_id='BCCAQv2+PRISM',
                   downscaling_package_id='github.com/pacificclimate/ClimDown',
                   driving_experiment=paste("historical,",rcp,sep=''),
                   driving_experiment_id=paste("historical,",rcp,sep=''),
                   driving_institution = drive.centre.name, ##Full name
                   driving_institute_id = drive.centre, ##Acronym
                   driving_model_id = gcm,
                   driving_realization = substr(run,2,2), ##These integers are from the 'r1i1p1' code
                   driving_initialization_method='1',
                   driving_physics_version = '1',
                   target_institution = "Pacific Climate Impacts Consortium",
                   target_institute_id = "PCIC",
                   target_dataset = "PCIC meteorology for NWNA",
                   target_dataset_id = "PNWNAMet",
                   target_references = "Werner, A.T., R.R. Shrestha, A.J. Cannon, M.S. Schnorbus, F.W. Zwiers, G. Dayon and F. Anslow, 2019: A long-term, temporally consistent, gridded daily meteorological dataset for northwestern North America. Scientific Data, 6, 180299, doi:10.1038/sdata.2018.299.",
                   target_version = "obtained: 14 Sept 2018",
                   target_contact = "Arelia Werner (wernera@uvic.ca)",
                   title = "Walter Snow Model Output for Vancouver Whistler")
  return(global.atts)
}


##----------------------------------------------------------------------
get_standard_atts <- function(var.name) {
  lon.atts <- list(standard_name="longitude",
                   long_name = "longitude",
                   units = "degrees_east",
                   axis = "X")
  lat.atts <- list(standard_name="latitude",
                   long_name = "latitude",
                   units = "degrees_north",
                   axis = "Y")
  snowdepth.atts <- list(standard_name = "snowdepth",
                    long_name = "Snow Depth",
                    missing_value = -9999.0,
                    cell_methods = "time: sum",
                    units = "m")
  swe.atts <- list(standard_name = "swe",
                  long_name = "Snow Water Equivalent",
                  missing_value = -9999.0,
                  cell_methods = "time: sum",
                  units = 'm')

  var.atts <- switch(var.name,
                     snowdepth=snowdepth.atts,
                     swe=swe.atts)

  rv <- list(lon=lon.atts,
             lat=lat.atts,
             var=var.atts)
  return(rv)
}



create.base.files <- function(var.name,gcm,
                              interval,
                              data.dir,write.dir) {
  
  past.file <- list.files(path=data.dir,pattern='pr_gcm_prism',full.name=TRUE)
  write.clim.name <- paste0(var.name,'_BCCAQ2-PRISM_',gcm,'_',interval,'.nc') 

  nc <- nc_open(past.file,write=FALSE)
  
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units

  time.start <- as.Date(strsplit(time.units, ' ')[[1]][3])
  past.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)

  past.values <- ncvar_get(nc,'time')

  ##Original time method
  full.values <- seq(past.values[1],tail(past.values,1),by=1)
  full.series <- format(past.origin + (past.values)*86400,'%Y-%m-%d')

  dates <- as.numeric(past.values) 

  ##Attributes to retain
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')  
  atts <- get_standard_atts(var.name)

  n.lon <- length(lon)
  n.lat <- length(lat)

  ##--------------------------------------------------------------
  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', time.units, dates,
                      unlim=FALSE, calendar=time.calendar)

  var.geog <- ncvar_def(var.name, units='m', dim=list(x.geog, y.geog, t.geog),
                        missval=atts[['missing_value']])

  file.nc <- nc_create(paste(write.dir,write.clim.name,sep=''), var.geog,h_minfree=102400)

  ##Time Attributes
  ncatt_put(file.nc,'time','standard_name','Time')
  ncatt_put(file.nc,'time','long_name','Time')

  ##Standard Attributes

  print('Lon names')
  lon.names <- names(atts$lon)
  for (j in 1:length(atts$lon))
    ncatt_put(file.nc,varid='lon',attname=lon.names[j],attval=atts$lon[[j]])
  print('Lat names')
  lat.names <- names(atts$lat)
  for (j in 1:length(atts$lat))
    ncatt_put(file.nc,varid='lat',attname=lat.names[j],attval=atts$lat[[j]])
  print('Var names')
  var.names <- names(atts$var)
  for (j in 1:length(atts$var))
    ncatt_put(file.nc,varid=var.name,attname=var.names[j],attval=atts$var[[j]])
  ##Global Attributes
  global.atts <- get_global_atts('ECMWF','ECMWF','ERA-Interim','','r1i1p1')
  global.names <- names(global.atts)
  for (g in 1:length(global.atts))
    ncatt_put(file.nc,varid=0,attname=global.names[g],attval=global.atts[[g]])

  ncvar_put(file.nc,'lon',lon)
  ncvar_put(file.nc,'lat',lat)

  nc_close(file.nc)
}


##**************************************************************************************
##-----------------------------------------------------------------------------
##Climdex from the 800 BCCAQ-PRISM output
run.bccaq.prism <- function() {
  
  interval <- '19790101-20181031'
  
  data.dir <-  '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/NCEP2/'  
  write.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/NCEP2/'  

  gcm.list <- 'NCEP2'
  for (model in gcm.list) {
    gcm <- model
    print(gcm)

    first <- create.base.files(var.name='snowdepth',gcm,
                               interval,
                               data.dir,write.dir)
    first <- create.base.files(var.name='swe',gcm,
                               interval,
                               data.dir,write.dir)
  }  
}

##**************************************************************************************
run.bccaq.prism()
