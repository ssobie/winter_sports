##Script to plot the fraction of useable MODIS data
library(raster)
library(rgdal)
library(TeachingDemos)
library(maps)
library(scales)
library(ncdf4)
        
##region.shp <- readOGR(shape.dir, region, stringsAsFactors=F, verbose=F)

site.dir <- '/storage/data/projects/rci/data/winter_sports/study_map/'
active.courses <- read.csv(paste0(site.dir,'active_snow_courses.csv'),header=T,as.is=T)
inactive.courses <- read.csv(paste0(site.dir,'inactive_snow_courses.csv'),header=T,as.is=T)
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

plot.file <- '/storage/data/projects/rci/data/winter_sports/plots/study.elev.map.2020.png'
shape.dir <- '/storage/data/projects/rci/data/winter_sports/study_map/'

  plot.window.xlim <- c(-123.7,-120.65)
  plot.window.ylim <- c(48.85,50.8)

  width <- 1100 ##3926
  height <- 500 ##2383
  png(file=plot.file,width=10,height=6,units='in',res=600,pointsize=6,bg='white')
  ##layout(matrix(c(1,1,1,1,1,1,1,1,2,3),nrow=2))
  layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,2,3,4),nrow=3))
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
     hillshade <- raster(paste0(shape.dir,'van_whistler_hillshade.tif'))

     ##img <- readTIFF(paste0(shape.dir,'van_whistler_hillshade.tif'), native=TRUE,convert=T)
     ##values(shade) <- img

     ##test <- readTIFF(paste0(shape.dir,'van_whistler_hillshade.tif'))
     image(hillshade,add=T,col = grey(1:100/100))
     image(shade,add=T,col = alpha(colour.map(),0.8))    

     rivers.shp <- spTransform(readOGR(shape.dir, 'van_whistler_rivers', stringsAsFactors=F, verbose=F),CRS("+init=epsg:4326"))  
     ##lakes.shp <- spTransform(readOGR(shape.dir, 'van_whistler_lakes', stringsAsFactors=F, verbose=F),CRS("+init=epsg:4326"))  
     lakes.shp <- spTransform(readOGR(shape.dir, 'van_whistler_lakes_simple_0.5%', stringsAsFactors=F, verbose=F),CRS("+init=epsg:4326"))

     plot(lakes.shp,add=TRUE,col='lightsteelblue',border='lightblue',xlim=plot.window.xlim,ylim=plot.window.ylim)
     plot(rivers.shp,add=TRUE,col='lightsteelblue',border='lightblue',xlim=plot.window.xlim,ylim=plot.window.ylim)

     ocean.shp <- readOGR('/storage/data/projects/rci/data/assessments/shapefiles/bc_common','west_coast_ocean',stringsAsFactors=F, verbose=F)
     plot(spTransform(ocean.shp,CRS("+init=epsg:4326")),add=TRUE,col='lightsteelblue',border='lightsteelblue')
     plot(spTransform(ocean.shp,CRS("+init=epsg:4326")),add=TRUE,border='black')

     plot(spTransform(us.shp,CRS("+init=epsg:4326")),add=TRUE,border='black',lwd=2)

     points(active.courses$Lon,active.courses$Lat,pch=24,col='black',bg='green',cex=2.5)
     points(inactive.courses$Lon,inactive.courses$Lat,pch=25,col='black',bg='red',cex=2.5)
     points(snow.pillows$Lon,snow.pillows$Lat,pch=23,col='black',bg='orange',cex=2.5)
 
     rs <- 0.25
     font <- 'Helvetica'
     tex <- 2.5

     shadowtext(-123.10,50.22,labels='CG',cex=tex,col='black',bg='white',r=0.2,family=font) ##Callaghan
     shadowtext(-122.98,49.525,labels='OL',cex=tex,col='black',bg='white',r=0.2,family=font) ##Orchid Lake
     shadowtext(-123.14,49.6,labels='LO',cex=tex,col='black',bg='white',r=0.2,family=font) ##Loch Lomond
     shadowtext(-123.1,49.35,labels='GR',cex=tex,col='black',bg='white',r=0.2,family=font) ##Grouse Mountain
     shadowtext(-123.2,49.35,labels='HB',cex=tex,col='black',bg='white',r=0.2,family=font) ##Hollyburn
     shadowtext(-123.0,49.325,labels='DM',cex=tex,col='black',bg='white',r=0.2,family=font) ##Dog Mountain
     shadowtext(-122.87,49.37,labels='MS',cex=tex,col='black',bg='white',r=0.2,family=font) ##Mount Seymour
     shadowtext(-123.1,49.474,labels='PL',cex=tex,col='black',bg='white',r=0.2,family=font) ##Palisade Lake
     shadowtext(-122.7,49.55,labels='DI',cex=tex,col='black',bg='white',r=0.2,family=font) ##Disappointment Lake
     shadowtext(-122.94,49.46,labels='BW',cex=tex,col='black',bg='white',r=0.2,family=font) ##Burwell Lake

     shadowtext(-122.97,49.75,labels='DH',cex=tex,col='black',bg='white',r=0.2,family=font) ##Diamond Head
     shadowtext(-123.65,49.7,labels='EL',cex=tex,col='black',bg='white',r=0.2,family=font) ##Edwards Lake
     shadowtext(-123.55,49.52,labels='CC',cex=tex,col='black',bg='white',r=0.2,family=font) ##Chapman Creek
##     shadowtext(-123.15,49.95,labels='GL',cex=tex,col='black',bg='white',r=0.2,family=font) ##Garibaldi Lake
     shadowtext(-122.9,50.1,labels='WM',cex=tex,col='black',bg='white',r=0.2,family=font) ##Whistler Mountain
     shadowtext(-123.05,50.45,labels='WC',cex=tex,col='black',bg='white',r=0.2,family=font) ##Wolverine Creek
     shadowtext(-122.85,50.53,labels='TQ',cex=tex,col='black',bg='white',r=0.2,family=font) ##Tenquille Snow Course
     shadowtext(-122.55,50.44,labels='DU',cex=tex,col='black',bg='white',r=0.2,family=font) ##Duffey Lake
     shadowtext(-122.65,50.7,labels='MP',cex=tex,col='black',bg='white',r=0.2,family=font) ##McGillivray Pass
##     shadowtext(-122.2,50.66,labels='SH',cex=tex,col='black',bg='white',r=0.2,family=font) ##Shalalth
     shadowtext(-122.1,49.76,labels='NA',cex=tex,col='black',bg='white',r=0.2,family=font) ##Nahatlatch
     shadowtext(-122.3,49.51,labels='SL',cex=tex,col='black',bg='white',r=0.2,family=font) ##Stave Lake
     shadowtext(-122.11,49.26,labels='DL',cex=tex,col='black',bg='white',r=0.2,family=font) ##Dickson Lake
     shadowtext(-121.49,50.2,labels='LY',cex=tex,col='black',bg='white',r=0.2,family=font) ##Lytton
##     shadowtext(-121.4,50.65,labels='PM',cex=tex,col='black',bg='white',r=0.2,family=font) ##Pavilion Mountain
     shadowtext(-120.9,50.55,labels='HV',cex=tex,col='black',bg='white',r=0.2,family=font) ##Highland Valley
     shadowtext(-121.15,50.49,labels='GM',cex=tex,col='black',bg='white',r=0.2,family=font) ##Gnawed Mountain
     shadowtext(-120.85,49.95,labels='SM',cex=tex,col='black',bg='white',r=0.2,family=font) ##Shovelnose Mountain
     shadowtext(-120.92,49.75,labels='BR',cex=tex,col='black',bg='white',r=0.2,family=font) ##Brookmere
     shadowtext(-120.85,49.45,labels='HH',cex=tex,col='black',bg='white',r=0.2,family=font) ##Hamilton Hill 
     shadowtext(-120.75,49.15,labels='BP',cex=tex,col='black',bg='white',r=0.2,family=font) ##Blackwall Peak 
     shadowtext(-120.92,49.05,labels='LL',cex=tex,col='black',bg='white',r=0.2,family=font) ##Lightning Lake 
     shadowtext(-121.37,49.09,labels='KL',cex=tex,col='black',bg='white',r=0.2,family=font) ##Klesilkwa
     shadowtext(-121.46,49.42,labels='HO',cex=tex,col='black',bg='white',r=0.2,family=font) ##Hope
     shadowtext(-121.32,49.25,labels='SW',cex=tex,col='black',bg='white',r=0.2,family=font) ##Sumallo River West
     shadowtext(-121.2,49.17,labels='SR',cex=tex,col='black',bg='white',r=0.2,family=font) ##Sumallo River
     shadowtext(-121.2,49.33,labels='NT',cex=tex,col='black',bg='white',r=0.2,family=font) ##New Tashme
     shadowtext(-121.32,49.55,labels='BU',cex=tex,col='black',bg='white',r=0.2,family=font) ##Boston Bar Upper
     shadowtext(-121.24,49.64,labels='GB',cex=tex,col='black',bg='white',r=0.2,family=font) ##Great Bear
     shadowtext(-121.14,49.66,labels='OT',cex=tex,col='black',bg='white',r=0.2,family=font) ##Ottomite
     shadowtext(-121.01,49.61,labels='BL',cex=tex,col='black',bg='white',r=0.2,family=font) ##Boston Bar Lower


     shadowtext(-121.7,49.61,labels='SC',cex=tex,col='black',bg='white',r=0.2,family=font) ##Spuzzum Creek
     ##shadowtext(-122.98,50.6,labels='TL',cex=tex,col='black',bg='white',r=0.2,family=font) ##Tenquille Lake Pillow
     shadowtext(-123.45,50.22,labels='US',cex=tex,col='black',bg='white',r=0.2,family=font) ##Upper Squamish
     shadowtext(-121.8,49.05,labels='CR',cex=tex,col='black',bg='white',r=0.2,family=font) ##Chilliwack River
     shadowtext(-121.67,49.27,labels='WL',cex=tex,col='black',bg='white',r=0.2,family=font) ##Wahleach Lake

     #text(-122.31581,49.68030,labels='Stave Lake',cex=tex,col='black',family=font)
     #text(-122.05926,49.92587,labels='Nahatlatch',cex=tex,col='black',family=font)
     #text(-121.58945,49.32987,labels='Wahleach',cex=tex,col='black',family=font)
     #text(-121.30865,49.22944,labels='Klesilkwa',cex=tex,col='black',family=font)
     #text(-120.91021,49.18479,labels='Lightning\nLake',cex=tex,col='black',family=font)
     #text(-120.96418,49.98463,labels='Shovelnose\nMountain',cex=tex,col='black',family=font)
     #text(-120.90558,49.54080,labels='Hamilton\nHill',cex=tex,col='black',family=font)



     #text(-120.97397,49.72503,labels='Brookmere',cex=tex,col='black',family=font)

     #shadowtext(-121.686,49.574,labels='Spuzzum Creek',cex=tex,col='black',bg='white',r=0.1,family=font)
     #shadowtext(-121.83667,49.1033,labels='Chilliwack\nRiver',cex=tex,col='black',bg='white',r=0.1,family=font)
     #shadowtext(-123.41,50.25,labels='Upper\nSquamish',cex=tex,col='black',bg='white',r=0.1,family=font)
     #shadowtext(-122.9333,50.6333,labels='Tenquille\nLake',cex=tex,col='black',bg='white',r=0.1,family=font)

     ##text(course.right$Lon,course.right$Lat-0.06,labels=course.right$Name,cex=tex,col='black',bg='white',r=rs)
     ##text(course.left$Lon-0.125,course.left$Lat+0.05,labels=course.left$Name,cex=tex,col='black',bg='white',r=rs)
     ##text(snow.pillows$Lon,snow.pillows$Lat+0.1,labels=snow.pillows$Name,cex=tex,col='black',bg='white',r=rs)

     ##points(-123.117,49.264,pch=16,cex=4)
     ##text(-123.117,49.2,labels='Vancouver',cex=tex,col='black',family=font)

     #text(-121.75,48.9,labels='USA',cex=tex,col='black',family=font)
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
  par(mar=c(1,0,6,3))
  plot(c(),xlim=xlim,ylim=ylim,xaxs='i',yaxs='i',
     bg='gray94',# 'lightgray',
     xlab='',ylab='',main='',axes=F)
     rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='lightblue')
     plot(pnw.bnds,add=T,col='darkgray')
     plot(snow.ex,add=T,border='black',lwd=2)
     text(0.6*xlim[2],0.6*ylim[2],'British\nColumbia',cex=1.7,family=font)
     box(which='plot',lwd=2)

  par(mar=c(3,0,3,3))
  plot(c(),xlim=xlim,ylim=ylim,xaxs='i',yaxs='i',
     bg='gray94',# 'lightgray',
     xlab='',ylab='',main='',axes=F)
  text(x=0.5*xlim[2], y = 0.9*ylim[2], labels ='Calibration\nSnow Courses',cex=2.75)
  points(x=0.9*xlim[2], y = 0.9*ylim[2], pch=24,col='black',bg='green',cex=3.25)
  text(x=0.5*xlim[2], y = 0.6*ylim[2], labels ='Evaluation\n Snow Courses',cex=2.75)
  points(x=0.9*xlim[2], y = 0.6*ylim[2], pch=25,col='black',bg='red',cex=3.25)
  text(x=0.5*xlim[2], y = 0.3*ylim[2], labels ='Automated\nSnow Pillows',cex=2.75)
  points(x=0.9*xlim[2], y = 0.3*ylim[2], pch=23,col='black',bg='orange',cex=3.25)

  box(which='plot',lwd=2)

legend_image <- as.raster(matrix(rev(colour.map()), ncol=1))
  par(mar=c(2.95,0,1.75,1))
  plot(c(0,3),c(0,3),type = 'n', axes = F,xlab = '', ylab = '', main='',yaxs='i',xaxs='i')
  rect(0.01,0.35,1,2.71,border='black',lwd=2)
  text(x=0.5, y = 2.85, labels ='Elev. (m)',cex=2.75)
  text(x=1.3, y = c(0.35,1.45,2.65), labels =c('0   ','1375','2750'),cex=2.5)
  rasterImage(legend_image, 0.03, 0.36, 0.98,2.7)
##  box(which='plot')
  dev.off()