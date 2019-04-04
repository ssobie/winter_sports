##Script to plot the fraction of useable MODIS data
library(raster)
library(rgdal)
library(TeachingDemos)
library(maps)

library(ncdf4)
        
##region.shp <- readOGR(shape.dir, region, stringsAsFactors=F, verbose=F)

site.dir <- '/storage/data/projects/rci/data/winter_sports/study_map/'
course.left <- read.csv(paste0(site.dir,'snow_course_left.csv'),header=T,as.is=T)
course.right <- read.csv(paste0(site.dir,'snow_course_right.csv'),header=T,as.is=T)
snow.pillows <- read.csv(paste0(site.dir,'snow_pillow_locations.csv'),header=T,as.is=T)

nc <- nc_open('/storage/data/projects/rci/data/prism/bc_prism_dem_elevations.nc')
lon <- ncvar_get(nc,'lon')
lat <- ncvar_get(nc,'lat')


lon.st <- which.min(abs(-123.6-lon))
lon.en <- which.min(abs(-120.75-lon))
lat.st <- which.min(abs(48.85-lat))
lat.en <- which.min(abs(50.75-lat))

dem <- ncvar_get(nc,'elevation',start=c(lon.st,lat.st,1),
                                count=c(lon.en-lon.st+1,lat.en-lat.st+1,1))
nc_close(nc)
dem.range <- range(dem)

dem.breaks <- seq(0,2750,25)
dem.raster <- raster('/storage/data/projects/rci/data/prism/bc_prism_dem_elevations.nc')

colour.map <- function() {
  green.to.tan <- colorRampPalette(colors=c('#99FF99','#FFDD97'),bias=1, space = "Lab", interpolate = "linear")
  tan.to.brown <- colorRampPalette(colors=c('#FFDD97','#FFFF99','#C1934D'),bias=1, space = "Lab", interpolate = "linear")
  brown.to.purple <- colorRampPalette(colors=c('#C1934D','#CB9ED0'),bias=1, space = "Lab", interpolate = "linear")
  purple.to.white <- colorRampPalette(colors=c('#CB9ED0','#F7F5F7'),bias=1, space = "Lab", interpolate = "linear")

  rv <- c(green.to.tan(21),tan.to.brown(30),brown.to.purple(30),purple.to.white(29))
  return(rv)
}


##---------------------------------------------------------------------------------
##Make Plot

plot.file <- '/storage/data/projects/rci/data/winter_sports/plots/study.elev.map.2018.png'
shape.dir <- '/storage/data/projects/rci/data/winter_sports/study_map/'

  plot.window.xlim <- c(-123.6,-120.75)
  plot.window.ylim <- c(48.85,50.75)

  width <- 1100 ##3926
  height <- 500 ##2383
  png(file=plot.file,width=10,height=6,units='in',res=600,pointsize=6,bg='white')
  layout(matrix(c(1,1,1,1,1,1,1,1,2,3),nrow=2))
  par(mar=c(6,6,6,1))
  plot(c(),xlim=plot.window.xlim,ylim=plot.window.ylim,xaxs='i',yaxs='i',
     bg='white',# 'lightgray',
     xlab='Longitude (\u00B0E)',ylab='Latitude (\u00B0N)',main='',
     cex.axis=2.5,cex.lab=2.5,cex.main=2.4,mgp=c(3.5,2,0))
     rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='white')

     bc.dir <- '/storage/data/gis/basedata/base_layers/'
     bc.overlay <- 'bc_province_wgs84'
     bc.shp <- readOGR(bc.dir, bc.overlay, stringsAsFactors=F, verbose=F)
     us.shp <- readOGR(bc.dir, 'united_states', stringsAsFactors=F, verbose=F)

     ##image(dem.raster, col=colour.map(),breaks=dem.breaks,ribbon=FALSE,xlim=plot.window.xlim, ylim=plot.window.ylim, add=TRUE)
     shade <- raster(paste0(shape.dir,'van_whistler.tif'))
     ##img <- readTIFF(paste0(shape.dir,'van_whistler_hillshade.tif'), native=TRUE,convert=T)
     ##values(shade) <- img

     ##test <- readTIFF(paste0(shape.dir,'van_whistler_hillshade.tif'))
     image(shade,add=T,col = colour.map())    

     rivers.shp <- spTransform(readOGR(shape.dir, 'van_whistler_rivers', stringsAsFactors=F, verbose=F),CRS("+init=epsg:4326"))  
     lakes.shp <- spTransform(readOGR(shape.dir, 'van_whistler_lakes', stringsAsFactors=F, verbose=F),CRS("+init=epsg:4326"))  

     plot(lakes.shp,add=TRUE,col='lightsteelblue',border='lightblue',xlim=plot.window.xlim,ylim=plot.window.ylim)
     plot(rivers.shp,add=TRUE,col='lightsteelblue',border='lightblue',xlim=plot.window.xlim,ylim=plot.window.ylim)

     ocean.shp <- readOGR('/storage/data/projects/rci/data/assessments/crd/shapefiles/','west_coast_ocean',stringsAsFactors=F, verbose=F)
     plot(spTransform(ocean.shp,CRS("+init=epsg:4326")),add=TRUE,col='lightsteelblue',border='lightsteelblue')
     plot(spTransform(ocean.shp,CRS("+init=epsg:4326")),add=TRUE,border='black')

     plot(spTransform(us.shp,CRS("+init=epsg:4326")),add=TRUE,border='black',lwd=2)

     points(course.left$Lon,course.left$Lat,pch=18,cex=4)
     points(course.right$Lon,course.right$Lat,pch=18,cex=4)
     points(snow.pillows$Lon,snow.pillows$Lat,pch=17,cex=3)
 
     rs <- 0.25
     font <- 'Helvetica'
     tex <- 2.5
     text(-123.05196,49.63678,labels='Orchid Lake',cex=tex,col='black',family=font)
     text(-122.75219,49.47443,labels='Palisade Lake',cex=tex,col='black',family=font)
     text(-123.27744,49.44365,labels='Grouse\nMountain',cex=tex,col='black',family=font)
     text(-122.31581,49.68030,labels='Stave Lake',cex=tex,col='black',family=font)
     text(-122.05926,49.92587,labels='Nahatlatch',cex=tex,col='black',family=font)
     text(-121.58945,49.32987,labels='Wahleach',cex=tex,col='black',family=font)
     text(-121.30865,49.22944,labels='Klesilkwa',cex=tex,col='black',family=font)
     text(-120.91021,49.18479,labels='Lightning\nLake',cex=tex,col='black',family=font)
     text(-120.96418,49.98463,labels='Shovelnose\nMountain',cex=tex,col='black',family=font)
     text(-120.90558,49.54080,labels='Hamilton\nHill',cex=tex,col='black',family=font)

     text(-123.10360,50.038228,labels='Callaghan',cex=tex,col='black',family=font)
     text(-122.72256,49.29252,labels='Dog Mountain',cex=tex,col='black',family=font)
     text(-120.97397,49.72503,labels='Brookmere',cex=tex,col='black',family=font)

     text(-121.686,49.574,labels='Spuzzum Creek',cex=tex,col='black',family=font)
     text(-121.83667,49.1033,labels='Chilliwack\nRiver',cex=tex,col='black',family=font)
     text(-123.4033,50.32,labels='Upper\nSquamish',cex=tex,col='black',family=font)
     text(-122.9333,50.6333,labels='Tenquille Lake',cex=tex,col='black',family=font)

     ##text(course.right$Lon,course.right$Lat-0.06,labels=course.right$Name,cex=tex,col='black',bg='white',r=rs)
     ##text(course.left$Lon-0.125,course.left$Lat+0.05,labels=course.left$Name,cex=tex,col='black',bg='white',r=rs)
     ##text(snow.pillows$Lon,snow.pillows$Lat+0.1,labels=snow.pillows$Name,cex=tex,col='black',bg='white',r=rs)

     points(-123.117,49.264,pch=16,cex=4)
     text(-123.117,49.2,labels='Vancouver',cex=tex,col='black',family=font)

     text(-121.75,48.9,labels='USA',cex=tex,col='black',family=font)
 ##    map.scale(x=-123.5, y=50.7, ratio=TRUE, relwidth=0.1)
     box(which='plot',lwd=2)
 
u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)     

##    legend('bottomright',legend=c(col=colour.map(),cex=2,pch=15,title='Elev. (m)')


v <- c( 0.8*v[2], v[2], 0.75*v[4], v[4] )
##v <- c( v[2], 1.2*v[2], 0.75*v[4], v[4] )
     bc.proj <- spTransform(bc.shp,CRS("+init=epsg:3005"))  
     snow.ex <- spTransform(readOGR(shape.dir, 'snow_extent', stringsAsFactors=F, verbose=F),CRS("+init=epsg:3005"))  
     pnw.bnds <- spTransform(readOGR(shape.dir, 'north_america_state_provincial_boundaries', stringsAsFactors=F, verbose=F),CRS("+init=epsg:3005"))  
     bounds <- extent(bc.proj)
     xlim <- c(bounds@xmin,bounds@xmax)
     ylim <- c(bounds@ymin,bounds@ymax)

##par( fig=v, new=TRUE, mar=c(0,0,0,0) )
  par(mar=c(10,0,6,3))
  plot(c(),xlim=xlim,ylim=ylim,xaxs='i',yaxs='i',
     bg='gray94',# 'lightgray',
     xlab='',ylab='',main='',axes=F)
     rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='lightblue')
     plot(pnw.bnds,add=T,col='darkgray')
     plot(snow.ex,add=T,border='black',lwd=2)
     text(0.6*xlim[2],0.6*ylim[2],'British\nColumbia',cex=1.7,family=font)
     box(which='plot',lwd=2)

legend_image <- as.raster(matrix(rev(colour.map()), ncol=1))
  par(mar=c(1,0,1,1))
  plot(c(0,3),c(0,3),type = 'n', axes = F,xlab = '', ylab = '', main='',yaxs='i',xaxs='i')
  rect(0.01,0.35,1,2.76,border='black',lwd=2)
  text(x=0.5, y = 2.9, labels ='Elev. (m)',cex=2.75)
  text(x=1.3, y = c(0.4,1.5,2.7), labels =c('0   ','1375','2750'),cex=2.5)
  rasterImage(legend_image, 0.03, 0.36, 0.98,2.75)
##  box(which='plot')
  dev.off()