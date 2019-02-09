library(R.utils)
library(curl)

ftp.dir <- 'ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/unmasked/'

tar.dir <- '/storage/data/projects/rci/data/winter_sports/obs/SNODAS/tar_files/'
gz.dir <- '/storage/data/projects/rci/data/winter_sports/obs/SNODAS/gz_files/'
tif.dir <- '/storage/data/projects/rci/data/winter_sports/obs/SNODAS/tif_files/'

dates <- seq(from=as.Date('2010-03-01'),by='month',to=as.Date('2018-12-31'))
daily <- seq(from=as.Date('2010-03-01'),by='day',to=as.Date('2018-12-31'))
missing.dates <- c()

for (i in seq_along(dates)) {
  date <- dates[i]
  print(date)
  year <- format(date,'%Y')
  month <- format(date,'%m')
  date.ix <- grep(format(date,'%Y%m'),format(daily,'%Y%m'))

  url <- paste0(ftp.dir,year,'/',month,'_',month.abb[as.numeric(month)],'/')
  h <- new_handle(dirlistonly=TRUE)
  con <- curl(url, "r", h)
  tbl <- read.table(con, stringsAsFactors=TRUE, fill=TRUE)
  close(con)
  tar.files <- as.vector(tbl[-c(1,2),])
  tar.dates <- sapply(tar.files,substr,17,24)
  mon.days <- format(daily[date.ix],'%Y%m%d')
  
  miss.ix <- mon.days %in% tar.dates
  missing.dates <- c(missing.dates,as.character(daily[date.ix][!miss.ix]))
  print(paste0('Missing Date: ',daily[date.ix][!miss.ix]))


  for (j in seq_along(tar.files)) {

    tar.file <- tar.files[j]
    day <- substr(tar.file,23,24)
    print(tar.file)  
    ##Download the unmasked files
    ##wget ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/unmasked/2015/04_Apr/SNODAS_unmasked_20150401.tar
    ##tar.file <- paste0('SNODAS_unmasked_',year,month,day,'.tar')

    wget.file <- paste0(ftp.dir,year,'/',month,'_',month.abb[as.numeric(month)],'/',tar.file)

    download <- paste0('wget ',wget.file,' -P ',tar.dir)
    system(download)
  
    ##First extraction
    ##tar -zxvf SNODAS_unmasked_20150401.tar
    ##untar <- paste0('tar -xvf ',tar.dir,tar.file,' ',tar.dir)
    ##system(untar)
    untar(paste0(tar.dir,tar.file),exdir=tar.dir)

    ##11034 is SWE
    ##11036 is Snowdepth
    gz.files <- list.files(pattern='.gz',path=tar.dir,full.name=T)
    keep.files <- grep('11034|11036',gz.files)
    file.remove(gz.files[-keep.files])
    gz.files <- list.files(pattern='.gz',path=tar.dir)
    dat.files <- gz.files[grep('.dat.gz',gz.files)]

    depth.file <- dat.files[grep('11036',dat.files)]
    gunzip(paste0(tar.dir,depth.file),overwrite=TRUE,destname=paste0(gz.dir,gsub('dat.gz','dat',depth.file)))
    swe.file <- dat.files[grep('11034',dat.files)]
    gunzip(paste0(tar.dir,swe.file),overwrite=TRUE,destname=paste0(gz.dir,gsub('dat.gz','dat',swe.file)))

    ##Create header file
    template.file <- paste0(gz.dir,'template.hdr')  

    swe.file <- list.files(path=gz.dir,pattern='11034')
    swe.hdr <- gsub('.dat','.hdr',swe.file)

    file.copy(from=template.file,to=paste0(gz.dir,swe.hdr))
  
    depth.file <- list.files(path=gz.dir,pattern='11036')
    depth.hdr <- gsub('.dat','.hdr',depth.file)
    file.copy(from=template.file,to=paste0(gz.dir,depth.hdr))

    ##Convert dat file to TIf. This must have a header file with the same name, but file extension .hdr
    swe.tif <- paste0(tif.dir,'SWE_SNODAS_UNMASKED_DOMAIN_',year,month,day,'.tif')
    make.swe.tif <- paste0("gdal_translate -of GTiff -a_srs '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' -a_ullr -130.516666666661 58.2333333333310 -62.2499999999975 24.0999999999990 -a_nodata -9999 ",gz.dir,swe.file," ", swe.tif)
    system(make.swe.tif)

  depth.tif <- paste0(tif.dir,'SNOWDEPTH_SNODAS_UNMASKED_DOMAIN_',year,month,day,'.tif')
  make.depth.tif <- paste0("gdal_translate -of GTiff -a_srs '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' -a_ullr -130.516666666661 58.2333333333310 -62.2499999999975 24.0999999999990 -a_nodata -9999 ",gz.dir,depth.file," ", depth.tif)
  system(make.depth.tif)

 
    file.remove(paste0(tar.dir,tar.file))
    ##file.remove(paste0(tar.dir,gz.files))
    dat.files <- list.files(path=gz.dir,pattern='.dat')
    file.remove(paste0(gz.dir,dat.files))
    hdr.files <- list.files(path=gz.dir,pattern='1.hdr')
    file.remove(paste0(gz.dir,hdr.files))
 }
   
}