##Script to convert the MODIS tiff files into netcdf for easier use

library(ncdf4)

##Function to create month file

modis.dir <- '/storage/data/projects/rci/data/winter_sports/MODIS_MERGED'

year.dates <- seq(from=as.Date('2001-01-01'),by='month',to=as.Date('2015-12-31'))

yrs <- format(year.dates,'%Y')
mns <- format(year.dates,'%m')

dates <- paste0(yrs,mns)

prefix <- 'snc.modis.merged.'

for (i in seq_along(dates)) {
  print(paste0(i,' in ',length(dates)))
  merged.file <- list.files(path=modis.dir,pattern=dates[i])
  subset.file <- gsub(pattern='merged',replacement='van_whistler',merged.file)
  merged.subset <- paste0('ncks -d lat,48.,51. -d lon,-124.,-120. ',modis.dir,'/',merged.file,' /storage/data/projects/rci/data/winter_sports/MODIS_VAN_WHISTLER/',subset.file)
  print(merged.subset)
  system(merged.subset)
}