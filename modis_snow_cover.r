##Script to download MODIS Data and produce snow cover maps

#tiles of interest
tiles <- c('h10v03','h10v04', 'h09v03', 'h09v04', 'h11v03')
#collapse for use in grep
tilepattern <- paste0(tiles, collapse = '|')

gettile <- function(tile, path, ftp){
  #command to dl the file
  args = paste0('curl ', ftp, tile, ' -o ', path, tile)
  system(args)
  
  #command to change to a geotiff from hdf.  only getting binary data here.
  args = paste0('gdal_translate -of GTiff HDF4_EOS:EOS_GRID:', path_to_working_directory, tile, ':MOD_Grid_Snow_500m:NDSI_Snow_Cover ', path, strsplit(tile, '.hdf')[[1]], '.tif')
  system(args)
  browser()  
#clean up the junk in the folder, note that if you are on windows, you'll need to use 'DEL' instead of 'RM
##  args = paste('rm ',path,'*.hdf',sep='')
##  system(args)
}

getdata <- function(yr,m,d,sat,path,shp, tilepattern){
  
  #format posix day
  date <- as.POSIXct(paste0(yr,'-',m,'-',d), tz = 'Etc/GMT+8')
  
  #formatted date for file folder on ftp
  folderdate <- format(date, format = '%Y.%m.%d')
  
  #days in a month, have to adjust for leap years.
  if (yr == 2004 | yr == 2008 | yr == 2012){
    mdays = c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)#leap year
  }else{
    mdays = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) #non leap year
  }
  
  #format DOY for actual mod10 file name.
  DOY <- ifelse(m==1, d, sum(mdays[1:(m-1)])+d)
  jd <- sprintf('%03d',DOY)
  
  #ftp folder location on the NSIDC server.  Ifelse is used to decide if its Aqua or Terra that you want.
  ftp <- ifelse(sat == 'MYD10A1', paste0('ftp://n5eil01u.ecs.nsidc.org/SAN/MOSA/', sat, '.006/', folderdate, '/'), paste0('ftp://n5eil01u.ecs.nsidc.org/SAN/MOST/', sat, '.006/', folderdate, '/'))
  
  #first need to list the contents of a directory and pipe to a file
  args <- paste0('curl -l ', ftp, ' > ', path, 'contents.txt')
  system(args)
    
  #read in the contents of the file and get the names of the tiles that we want for the day.
  list <- as.character(read.table(paste(path,'contents.txt',sep=''), header = F)$V1)
  #get only the hdfs
  list <- list[grep(list,pattern = '.hdf$')]
  #get only your tiles
  hdfs <- list[grep(list,pattern = tilepattern)]
  
  #apply the gettile function for all tiles (no output in R, only in your folder)
  lapply(hdfs,gettile, path = path, ftp = ftp)
  browser()  
  #merge into one file and clip to cheakenv, have to change this slightly if you are using a different number of tiles.
  #get all the tifs for your date.
  tifs <- list.files(path, pattern = paste0(sat,'.A', yr, jd))
  
  #mergegrid name
  mergegrid <- paste0(path,sat,'.A', yr, jd, '.merge.tif')
  
  args <- paste('gdal_merge.py -of GTiff -o' , paste0(path , mergegrid), paste0(path , tifs[1]) , paste0(path , tifs[2]) , paste0(path , tifs[3]) , paste0(path , tifs[4]), paste0(path , tifs[5]))
  
  system(args)
  browser()  
  #projected grid name
  projgrid <- paste0(sat, '.A', yr, jd, '.merge.bcalb.tif')
  
  #project from MODIS sinuisoidal projection to bcalbers and force to 500*500 using gdalwarp command line tool.  If a different output projection is needed, change the EPSG code.
  args <- paste0("gdalwarp -overwrite -tr 500 500 -s_srs '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_def' -t_srs EPSG:3005 -of GTiff " , path , mergegrid , ' ' , path , projgrid)
  system(args)
  
  #clipped grid name (dropped the .A part for easier deleting of unneeded files)
  clipgrid <- paste0(sat, '.', yr,'_',sprintf('%02d',m),'_', sprintf('%02d',d) , '.clipped.tif')
  
  #clip grid to Province of BC shapefile
  #clip to BC grid.
  args <- paste0('gdalwarp -overwrite -dstnodata 255 -q -cutline ' , path, shp , ' -crop_to_cutline -of GTiff ' , path , projgrid , ' ' , path , clipgrid)
  
  system(args)
  
  #clean up the junk, same as before, if you are on windows, you want DEL, not rm
  args <- paste0('rm ', sat, '.A*')
  ##system(args)
  
  
  #return the file name, in case you want to read in to R
  return(clipgrid)
  
}

path_to_working_directory <- '/storage/data/projects/rci/data/winter_sports/MODIS/'

snow <- getdata(2014,4,1,'MOD10A1', path_to_working_directory,'province_bcalbers.shp', tilepattern)

library(raster)
library(ggplot2)
library(dplyr)

snowimage <- raster(snow) %>%
  rasterToPoints %>% 
  as.data.frame %>% 
  set_names(c('Easting', 'Northing', 'SCA')) %>%
  mutate(SCA=ifelse(SCA==25 |SCA==50|SCA==200,SCA,NA),
         SCA=factor(SCA,levels=c(25,200,50),labels=c('No Snow','Snow','Cloud')),
         Src='Terra')


snowimage %>% ggplot(aes(x=Easting,y=Northing, fill=SCA)) + 
  geom_raster() +
  scale_fill_manual(values=c('#4daf4a', '#377eb8','#E41C1C')) + 
  theme_minimal()




