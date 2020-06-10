##Script to pull out the time series of pr, tasmax, tamsin from the driving models 
library(ncdf4)
library(raster)

source('/storage/home/ssobie/code/repos/winter_sports/site.coordinates.r',chdir=T)

get_PRISM_data <- function(site, dem) {

  coords <- get_coordinates(site)
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

sites <- c('blackwall_peak_course','boston_bar_lower','boston_bar_upper','burwell_lake',
           'chapman_creek','cornwall_hills','diamond_head','edwards_lake',
           'hollyburn','hope',
           'loch_lomond','lytton','mount_seymour','new_tashme',
           'ottomite','pavilion_mountain',
           'sumallo_river','tenquille_course','whistler_mountain','wolverine_creek' )



sites <- c('shovelnose_mountain','brookmere','lightning_lake','callaghan','orchid_lake',
           'palisade_lake','grouse_mountain','dog_mountain','stave_lake','nahatlatch',
           'wahleach','klesilkwa','hamilton_hill','dickson_lake','disappointment_lake',
           'duffey_lake','gnawed_mountain','highland_valley','mcgillivray_pass',
           'sumallo_river_west','great_bear','upper_squamish','spuzzum_creek','chilliwack_river','tenquille_lake',
           'wahleach_lake','blackwall_peak_pillow')


sites <- sort(sites)


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

write.file <- '/storage/data/projects/rci/data/winter_sports/course_and_pillow_prism_elevation_biases_calibration.csv'

##write.table(site.info,file=write.file,sep=',',col.name=T,row.name=F,quote=F)

##---
##Compare site elevation differences with SWE biases

cal.sites <- c('shovelnose_mountain','brookmere','lightning_lake','callaghan','orchid_lake',
           'palisade_lake','grouse_mountain','dog_mountain','stave_lake','nahatlatch',
           'wahleach','klesilkwa','hamilton_hill','dickson_lake','disappointment_lake',
           'duffey_lake','gnawed_mountain','highland_valley','mcgillivray_pass',
           'sumallo_river_west','great_bear','upper_squamish','spuzzum_creek','chilliwack_river','tenquille_lake',
           'wahleach_lake','blackwall_peak_pillow')

cal.elevs <- matrix(NA,nrow=length(cal.sites),ncol=5)
  for (s in seq_along(cal.sites)) {
    site <- cal.sites[s]
    print(site)
    tas.elev <- get_PRISM_data(site, dem.tas)
    pr.elev <- get_PRISM_data(site, dem.pr)
    cal.elevs[s,] <- c(tas.elev$site,
                        tas.elev$dem,tas.elev$dem-tas.elev$site,
                        pr.elev$dem,pr.elev$dem-pr.elev$site)
    
  }


val.sites <- c('blackwall_peak_course','boston_bar_lower','boston_bar_upper','burwell_lake',
           'chapman_creek','cornwall_hills','diamond_head','edwards_lake',
           'hollyburn','hope',
           'loch_lomond','lytton','mount_seymour','new_tashme',
           'ottomite','pavilion_mountain',
           'sumallo_river','tenquille_course','whistler_mountain','wolverine_creek' )

val.elevs <- matrix(NA,nrow=length(val.sites),ncol=5)
  for (s in seq_along(val.sites)) {
    site <- val.sites[s]
    print(site)
    tas.elev <- get_PRISM_data(site, dem.tas)
    pr.elev <- get_PRISM_data(site, dem.pr)
    val.elevs[s,] <- c(tas.elev$site,
                        tas.elev$dem,tas.elev$dem-tas.elev$site,
                        pr.elev$dem,pr.elev$dem-pr.elev$site)
    
  }

