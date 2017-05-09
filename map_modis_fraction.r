##Script to plot the fraction of useable MODIS data
source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)
source('/storage/data/projects/rci/bcgov/moti/nrcan-precip_case_studies/code/moti.climdex.robjects.r',chdir=TRUE)

get.region.shape <- function(region,shape.dir) {
  region.shp <- readOGR(shape.dir, region, stringsAsFactors=F, verbose=F)
  return(region.shp)
}


shape.dir <- paste0('/storage/data/projects/rci/data/assessments/metro_van/shapefiles/')
region <- 'metro_van'
region.shp <- spTransform(get.region.shape(region,shape.dir),CRS("+init=epsg:4326")) 


modis.dir <- '/storage/data/projects/rci/data/winter_sports/MODIS_VAN_WHISTLER/'
snc.file <- paste0(modis.dir,'snc.modis.van_whistler.20010101-20151231.nc')
snc.nc <- nc_open(snc.file)
snc.data <- ncvar_get(snc.nc,'snc')

snc.useful <- apply(snc.data < 100,c(1,2),sum,na.rm=T)
snc.fraction <- snc.useful/(dim(snc.data)[3])*100

##image(snc.fraction)

plot.file <- '/storage/data/projects/rci/data/winter_sports/plots/modis.useable.data.png'
plot.title <- 'Fraction of Useable Days from MODIS'
plot.window.xlim <- c(-123.5,-121.5)
plot.window.ylim <- c(49.0,50.0)


width <- 1200
height <- 900
png(file=plot.file,width=width,height=height,bg='gray94')
par(mar=c(6,6,7,6))
plot(c(),xlim=plot.window.xlim,ylim=plot.window.ylim,xaxs='i',yaxs='i',
   bg='white',# 'lightgray',
     xlab='Longitude (\u00B0E)',ylab='Latitude (\u00B0N)',main=plot.title,
     cex.axis=2,cex.lab=2,cex.main=2.2,mgp=c(3.5,2,0))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='lightsteelblue') ##'lightgray')

shape.dir <- '/storage/data/gis/basedata/base_layers/'
bc.overlay <- 'bc_province_wgs84'
bc.shp <- readOGR(shape.dir, bc.overlay, stringsAsFactors=F, verbose=F)
plot(spTransform(bc.shp,CRS("+init=epsg:4326")),add=TRUE,xlim=plot.window.xlim,ylim=plot.window.ylim,col='lightgray')

map.range <- range(snc.fraction,na.rm=T)
class.breaks <- get.class.breaks('snc',type='past',map.range,manual.breaks='')

colour.ramp <- get.legend.colourbar(var.name='snc',map.range=map.range,
                                    my.bp=0,class.breaks=class.breaks,
                                    type)

image(snc.fraction, col=colour.ramp,breaks=class.breaks,ribbon=FALSE,xlim=plot.window.xlim, ylim=plot.window.ylim, add=TRUE,alpha=0.3)


ocean.shp <- readOGR('/storage/data/projects/rci/data/assessments/crd/shapefiles/','west_coast_ocean',stringsAsFactors=F, verbose=F)

##plot(spTransform(ocean.shp,CRS("+init=epsg:4326")),add=TRUE,col='lightgray',border='lightgray')

box(which='plot',lwd=3)

dev.off()




