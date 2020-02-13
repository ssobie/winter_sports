##Script to pull out the time series of pr, tasmax, tamsin from the driving models 
library(ncdf4)
library(raster)


get.coordinates <- function(site) {

  coordinates <- list(callaghan=c(-123.1036,50.1383278,1009),
                      orchid_lake=c(-123.0519638,49.53678,1178),
                      palisade_lake=c(-123.0321944,49.454433,898),
                      grouse_mountain=c(-123.0774472,49.383655,1126),
                      dog_mountain=c(-122.96255,49.37251944,1007),
                      dickson_lake=c(-122.06984166,49.3168194,1147),
                      stave_lake=c(-122.315805,49.58030277,1211),
                      nahatlatch=c(-122.059261,49.825866,1530),
                      wahleach=c(-121.57945,49.2298694,1395),
                      klesilkwa=c(-121.3086527,49.129438,610),
                      lightning_lake=c(-120.850205,49.044788,1254),
                      brookmere=c(-120.87397,49.815027,994),
                      shovelnose_mountain=c(-120.864175,49.8546305,1456),
                      hamilton_hill=c(-120.7955805,49.4988027,1477),
                      spuzzum_creek=c(-121.686,49.74,1197),
                      chilliwack_river=c(-121.71667,49.0333,1600),
                      upper_squamish=c(-123.4333,50.1500,1340),
                      wahleach_lake=c(-121.5833,49.2333,1400),
                      tenquille_lake=c(-122.9333,50.5333,1680))

  rv <- coordinates[[site]]
  return(rv)
}

get_PRISM_data <- function(site, dem) {

  coords <- get.coordinates(site)
  plot.title <- site

  dem.coords <- coordinates(dem)

  lon <- unique(dem.coords[,1])
  lat <- unique(dem.coords[,2])

  lon.bnds <- coords[1]
  lat.bnds <- coords[2]
  elev <- coords[3]

  lon.ix <- which.min(abs(lon-lon.bnds))
  lat.ix <- which.min(abs(lat-lat.bnds))

  dem.elev <- as.numeric(dem[lat.ix,lon.ix,1])

  return(list(site=elev,dem=dem.elev))
}

  sites <- c('callaghan',
             'orchid_lake',
             'palisade_lake',
             'grouse_mountain',
             'dog_mountain',
             'stave_lake',
             'nahatlatch',
             'wahleach',
             'klesilkwa',
             'lightning_lake',
             'brookmere',
             'shovelnose_mountain',
             'hamilton_hill',
             'spuzzum_creek',
             'chilliwack_river',
             'upper_squamish',
             'tenquille_lake')

base.dir <- "/storage/data/projects/PRISM/bc_climate/grids/"

dem.tas <- brick(paste0(base.dir,"PRISM_BC_Domain_30s_dem.grass.tiff"))
dem.pr <-  brick(paste0(base.dir,"PRISM_BC_Domain_30s375m_dem.grass.test.tiff"))

site.elevs <- matrix(NA,nrow=length(sites),ncol=5)

  for (s in seq_along(sites)) {
    site <- sites[s]
    print(site)
    tas.elev <- get_PRISM_data(site, dem.tas)
    pr.elev <- get_PRISM_data(site, dem.pr)
    site.elevs[s,] <- c(tas.elev$site,
                        tas.elev$dem,tas.elev$dem-tas.elev$site,
                        pr.elev$dem,pr.elev$dem-pr.elev$site)
    
  }

site.info <- cbind(sites,site.elevs)
colnames(site.info) <- c('Site','Site Elev','TAS DEM','TAS DEM-Site','PR DEM','PR DEM-Site')
print(site.info)

write.file <- '/storage/data/projects/rci/data/winter_sports/course_and_pillow_prism_elevation_biases.csv'

write.table(site.info,file=write.file,sep=',',col.name=T,row.name=F,quote=F)