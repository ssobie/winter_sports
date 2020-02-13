##Script to plot the fraction of useable MODIS data
source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)
library(scales)
library(RColorBrewer) # to get the sweet color ramps


##---------------------------------------------------------------------------------
##Snow Colour ramps based on the National Weather Service

snow_colour_ramp <- function(len,type) {

  if (type == 'past') {
    color.brewer.swe <- c("#CCFFE5","#99FFFF","#66B2FF","#3333FF","#7F00FF","#CC00CC")
    swe.ramp <- colorRampPalette(colors=color.brewer.swe, bias=1, space = "Lab", interpolate = "linear")
    rv <- swe.ramp(len-1)
  }
  if (type == 'anomaly' | type =='percent') {
    color.brewer.swe <- c("#FFFF99","#FFFF66","#FF9933","#FF8000","#CC6600","#990000")
    swe.ramp <- colorRampPalette(colors=color.brewer.swe, bias=1, space = "Lab", interpolate = "linear")
    rv <- swe.ramp(len-1)
  }
  if (type =='percent') {  
    rv <- rev(swe.ramp(len-1))
  }
  return(rv)
}

##---------------------------------------------------------------------------------
##Mapping Component

make_van_whistler_plot <- function(var.name,plot.type,var.title,plot.data,plot.file,class.breaks=NULL,mark,leg.title) {
   map.range <- range(as.matrix(plot.data),na.rm=T)
   if (!is.null(class.breaks)) {
      class.breaks <- class.breaks
   } else {
      class.breaks <- get.class.breaks(var.name,plot.type,map.range,manual.breaks='')      
   }

   
   if (var.name == 'swe' | var.name=='snowdepth') {
      colour.ramp <- snow_colour_ramp(length(class.breaks),plot.type)      
   } else {
      colour.ramp <- get.legend.colourbar(var.name=var.name,map.range=map.range,
                                          my.bp=0,class.breaks=class.breaks,
                                         type=plot.type)
   }
   map.class.breaks.labels <- get.class.break.labels(class.breaks)

   if ((var.name=='swe'|var.name=='snowdepth') & plot.type=='past') {
      white.class.breaks <- c(class.breaks,100000)                                       
      white.colour.ramp <- c(colour.ramp,'#FFFFFF')
   } else if ((var.name=='swe'|var.name=='snowdepth') & plot.type=='anomaly') {
      white.class.breaks <- class.breaks                                     
      white.colour.ramp <- colour.ramp
   } else {
      white.class.breaks <- class.breaks                                       
      white.colour.ramp <- colour.ramp
   }

   alb.crs <- "+init=epsg:4326"
   lons <- c(-124.0,-123.5, -123.0,-122.5,-122.0,-121.5,-121.0,-120.5,-120.0)
   lats <- c(48.75, 49.0, 49.25, 49.5, 49.75, 50.0, 50.25, 50.5, 50.75)

   shape.dir <- '/storage/data/projects/rci/data/assessments/shapefiles/bc'
   bc.shp <- get.region.shape('bc',shape.dir)

   map.extent <- extent(c(-123.55,-120.7,48.9,50.65))
   crop.extent <- map.extent

   plot.window.xlim <- c(map.extent@xmin,map.extent@xmax)
   plot.window.ylim <- c(map.extent@ymin,map.extent@ymax)

   par(mar=c(5,5.5,2,2))
   ##par(mar=mark)
   plot(c(),xlim=plot.window.xlim,ylim=plot.window.ylim,xaxs='i',yaxs='i',
   bg='white',# 'gray94',
   xlab='Longitude (\u00B0E)',ylab='Latitude (\u00B0N)',main='',##plot.title,
   cex.axis=2.5,cex.lab=2.5,cex.main=2.5,mgp=c(3.5,2,0),axes=F)
   rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='lightgray')

   axis(2,at=lats,label=lats,cex.axis=2.25)
   axis(1,at=lons,label=lons,cex.axis=2.25)

   bc.overlay <- 'north_america_state_provincial_boundaries'
   borders.shp <- readOGR(shape.dir, bc.overlay, stringsAsFactors=F, verbose=F)
   wco.shp <- readOGR('/storage/data/projects/rci/data/assessments/shapefiles/bc_common',
                           'ocean_mask', stringsAsFactors=F, verbose=F)
   shape.dir <- paste0('/storage/data/projects/rci/data/assessments/metro_van/shapefiles/')
   region <- 'metro_van'
   region.shp <- spTransform(get.region.shape(region,shape.dir),CRS("+init=epsg:4326")) 
   whistler.dir <- '/storage/data/projects/rci/data/assessments/whistler/shapefiles/'
   whistler.shp  <- spTransform(get.region.shape('WhistlerLandscapeUnit',whistler.dir),CRS("+init=epsg:4326")) 
   shade.dir <- '/storage/data/projects/rci/data/winter_sports/study_map/'
   shade <- raster(paste0(shade.dir,'van_whistler_hillshade.tif'))

   glacier.dir <- '/storage/data/gis/basedata/randolph_glacier_inventory/v32/02_rgi32_WesternCanadaUS'
   glacier.name <- '02_rgi32_WesternCanadaUS'
   glacier.shp <- spTransform(get.region.shape(glacier.name,glacier.dir),CRS("+init=epsg:4326"))

   wash.shp <- spTransform(readOGR(shade.dir, 'washington', stringsAsFactors=F, verbose=F),CRS("+init=epsg:4326"))
   rivers.shp <- spTransform(readOGR(shade.dir, 'van_whistler_rivers', stringsAsFactors=F, verbose=F),CRS("+init=epsg:4326"))
   lakes.shp <- spTransform(readOGR(shade.dir, 'van_whistler_lakes_simple_0.5%', stringsAsFactors=F, verbose=F),CRS("+init=epsg:4326"))

   image(shade,add=T,col = grey(1:100/100))

   image(plot.data,add=T,col=alpha(white.colour.ramp,0.8),breaks=white.class.breaks)

   ##plot(spTransform(borders.shp,CRS(alb.crs)),add=TRUE,border='black',cex=0.5)
   if (var.name=='swe'|var.name=='snowdepth') {
      plot(spTransform(glacier.shp,CRS(alb.crs)),add=TRUE,col=alpha('white',0.85),border=alpha('white',0.85),cex=0.5)
   }

   plot(lakes.shp,add=TRUE,col='lightgray',border='lightgray',xlim=plot.window.xlim,ylim=plot.window.ylim)
   plot(rivers.shp,add=TRUE,col='lightgray',border='lightgray',xlim=plot.window.xlim,ylim=plot.window.ylim)
 
   plot(spTransform(wco.shp,CRS(alb.crs)),add=TRUE,col='lightgray',border='darkgray',lwd=0.5)
   plot(spTransform(region.shp,CRS(alb.crs)),add=TRUE,border='black',cex=0.5)
   plot(spTransform(whistler.shp,CRS(alb.crs)),add=TRUE,border='black',cex=0.5)
   plot(spTransform(wash.shp,CRS(alb.crs)),add=TRUE,col='gray',border='black',cex=0.5)
##   abline(v=lons,lty=3,col='gray',lwd=0.7)
##   abline(h=lats,lty=3,col='gray',lwd=0.7)

   ##par(xpd=NA)
   
   legend('topright', col = "black", legend=rev(map.class.breaks.labels), pch=22, pt.bg = rev(alpha(colour.ramp,0.8)),
         pt.cex=1.55, y.intersp=0.8, title.adj=0.2, title=leg.title, xjust=0, cex=1.8)

   box(which='plot',lwd=2)

##   dev.off()
##browser()
}




##-------------------------------------------------------------------------

get.region.shape <- function(region,shape.dir) {
  region.shp <- readOGR(shape.dir, region, stringsAsFactors=F, verbose=F)
  return(region.shp)
}

plot.dir <- '/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/'
##  png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
##dev.off()