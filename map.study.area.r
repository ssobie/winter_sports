##Script to plot the fraction of useable MODIS data
library(raster)
library(rgdal)
library(TeachingDemos)
library(maps)
        
##region.shp <- readOGR(shape.dir, region, stringsAsFactors=F, verbose=F)

site.dir <- '/storage/data/projects/rci/data/winter_sports/study_map/'
course.left <- read.csv(paste0(site.dir,'snow_course_left.csv'),header=T,as.is=T)
course.right <- read.csv(paste0(site.dir,'snow_course_right.csv'),header=T,as.is=T)
snow.pillows <- read.csv(paste0(site.dir,'snow_pillow_locations.csv'),header=T,as.is=T)

##---------------------------------------------------------------------------------
##Make Plot

plot.file <- '/storage/data/projects/rci/data/winter_sports/plots/study.area.map.2018.png'
shape.dir <- '/storage/data/projects/rci/data/winter_sports/study_map/'

  plot.window.xlim <- c(-123.6,-120.75)
  plot.window.ylim <- c(48.85,50.75)

  width <- 1000 ##3926
  height <- 800 ##2383
  png(file=plot.file,width=10,height=7,units='in',res=600,pointsize=6,bg='white')
  par(mar=c(6,6,6,5))
  plot(c(),xlim=plot.window.xlim,ylim=plot.window.ylim,xaxs='i',yaxs='i',
     bg='white',# 'lightgray',
     xlab='Longitude (\u00B0E)',ylab='Latitude (\u00B0N)',main='',
     cex.axis=2.2,cex.lab=2.2,cex.main=2.4,mgp=c(3.5,2,0))
     rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='white')

     bc.dir <- '/storage/data/gis/basedata/base_layers/'
     bc.overlay <- 'bc_province_wgs84'
     bc.shp <- readOGR(bc.dir, bc.overlay, stringsAsFactors=F, verbose=F)
     us.shp <- readOGR(bc.dir, 'united_states', stringsAsFactors=F, verbose=F)

     shade <- raster(paste0(shape.dir,'van_whistler_hillshade.tif'))
     ##img <- readTIFF(paste0(shape.dir,'van_whistler_hillshade.tif'), native=TRUE,convert=T)
     ##values(shade) <- img

     ##test <- readTIFF(paste0(shape.dir,'van_whistler_hillshade.tif'))
     image(shade,add=T,col = grey(1:100/100))

     rivers.shp <- spTransform(readOGR(shape.dir, 'van_whistler_rivers', stringsAsFactors=F, verbose=F),CRS("+init=epsg:4326"))  
     lakes.shp <- spTransform(readOGR(shape.dir, 'van_whistler_lakes', stringsAsFactors=F, verbose=F),CRS("+init=epsg:4326"))  

     plot(lakes.shp,add=TRUE,col='lightsteelblue',border='lightblue',xlim=plot.window.xlim,ylim=plot.window.ylim)
     plot(rivers.shp,add=TRUE,col='lightsteelblue',border='lightblue',xlim=plot.window.xlim,ylim=plot.window.ylim)

     ocean.shp <- readOGR('/storage/data/projects/rci/data/assessments/crd/shapefiles/','west_coast_ocean',stringsAsFactors=F, verbose=F)
     plot(spTransform(ocean.shp,CRS("+init=epsg:4326")),add=TRUE,col='lightsteelblue',border='lightsteelblue')
     plot(spTransform(ocean.shp,CRS("+init=epsg:4326")),add=TRUE,border='black')

     plot(spTransform(us.shp,CRS("+init=epsg:4326")),add=TRUE,border='black',lwd=2)

     points(course.left$Lon,course.left$Lat,pch=18,cex=5)
     points(course.right$Lon,course.right$Lat,pch=18,cex=5)
     points(snow.pillows$Lon,snow.pillows$Lat,pch=17,cex=4)
 
     rs <- 0.25
     shadowtext(-123.05196,49.63678,labels='Orchid Lake',cex=2.75,col='black',bg='white',r=rs)
     shadowtext(-122.75219,49.47443,labels='Palisade Lake',cex=2.75,col='black',bg='white',r=rs)
     shadowtext(-123.27744,49.44365,labels='Grouse\nMountain',cex=2.75,col='black',bg='white',r=rs)
     shadowtext(-122.31581,49.68030,labels='Stave Lake',cex=2.75,col='black',bg='white',r=rs)
     shadowtext(-122.05926,49.92587,labels='Nahatlatch',cex=2.75,col='black',bg='white',r=rs)
     shadowtext(-121.58945,49.32987,labels='Wahleach',cex=2.75,col='black',bg='white',r=rs)
     shadowtext(-121.30865,49.22944,labels='Klesilkwa',cex=2.75,col='black',bg='white',r=rs)
     shadowtext(-120.91021,49.18479,labels='Lightning\nLake',cex=2.75,col='black',bg='white',r=rs)
     shadowtext(-120.96418,49.98463,labels='Shovelnose\nMountain',cex=2.75,col='black',bg='white',r=rs)
     shadowtext(-120.90558,49.54080,labels='Hamilton\nHill',cex=2.75,col='black',bg='white',r=rs)

     shadowtext(-123.10360,50.038228,labels='Callaghan',cex=2.75,col='black',bg='white',r=rs)
     shadowtext(-122.72256,49.29252,labels='Dog Mountain',cex=2.75,col='black',bg='white',r=rs)
     shadowtext(-120.97397,49.72503,labels='Brookmere',cex=2.75,col='black',bg='white',r=rs)

     shadowtext(-121.686,49.574,labels='Spuzzum Creek',cex=2.75,col='black',bg='white',r=rs)
     shadowtext(-121.83667,49.1033,labels='Chilliwack\nRiver',cex=2.75,col='black',bg='white',r=rs)
     shadowtext(-123.4033,50.32,labels='Upper\nSquamish',cex=2.75,col='black',bg='white',r=rs)
     shadowtext(-122.9333,50.6333,labels='Tenquille Lake',cex=2.75,col='black',bg='white',r=rs)

     ##shadowtext(course.right$Lon,course.right$Lat-0.06,labels=course.right$Name,cex=2.75,col='black',bg='white',r=rs)
     ##shadowtext(course.left$Lon-0.125,course.left$Lat+0.05,labels=course.left$Name,cex=2.75,col='black',bg='white',r=rs)
     ##shadowtext(snow.pillows$Lon,snow.pillows$Lat+0.1,labels=snow.pillows$Name,cex=2.75,col='black',bg='white',r=rs)

     points(-123.117,49.264,pch=16,cex=4)
     shadowtext(-123.117,49.2,labels='Vancouver',cex=2.75,col='black',bg='white',r=rs)

     shadowtext(-121.75,48.9,labels='USA',cex=2.75,col='black',bg='white',r=rs)
 ##    map.scale(x=-123.5, y=50.7, ratio=TRUE, relwidth=0.1)
     box(which='plot',lwd=3)
 
u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)     
v <- c( 0.8*v[2], v[2], 0.75*v[4], v[4] )
     bc.proj <- spTransform(bc.shp,CRS("+init=epsg:3005"))  
     snow.ex <- spTransform(readOGR(shape.dir, 'snow_extent', stringsAsFactors=F, verbose=F),CRS("+init=epsg:3005"))  
     bounds <- extent(bc.proj)
     xlim <- c(bounds@xmin,bounds@xmax)
     ylim <- c(bounds@ymin,bounds@ymax)

par( fig=v, new=TRUE, mar=c(0,0,0,0) )

  plot(c(),xlim=xlim,ylim=ylim,xaxs='i',yaxs='i',
     bg='gray94',# 'lightgray',
     xlab='',ylab='',main='',axes=F)
     rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='gray94')
     plot(bc.proj,add=T,col='darkgray')
     plot(snow.ex,add=T,border='black',lwd=2)
     text(0.6*xlim[2],0.6*ylim[2],'British\nColumbia',cex=1.7)
     box(which='plot')
    dev.off()


