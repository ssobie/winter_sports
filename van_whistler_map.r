##Script to plot the fraction of useable MODIS data
source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)
source('/storage/data/projects/rci/bcgov/moti/nrcan-precip_case_studies/code/moti.climdex.robjects.r',chdir=TRUE)
source('/storage/home/ssobie/code/hg/pievc/spatial.r',chdir=T)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
        
get.region.shape <- function(region,shape.dir) {
  region.shp <- readOGR(shape.dir, region, stringsAsFactors=F, verbose=F)
  return(region.shp)
}




##---------------------------------------------------------------------------------
##Make Plot

##plot.file <- '/storage/data/projects/rci/data/winter_sports/plots/era.snc.success.rate.data.png'
##plot.title <- 'ERA Snow Model-MODIS Snow Cover Comparison'
vw.plot <- function(raster.object,
                    class.breaks,map.class.breaks.labels,colour.ramp,
                    plot.file,plot.title,leg.title,
                    glaciers=FALSE,bias=FALSE) {

###  plot.window.xlim <- c(-123.75,-120.5)
###  plot.window.ylim <- c(48.75,51.0)

   map.extent <- extent(c(-123.55,-120.7,48.9,50.65))
   plot.window.xlim <- c(map.extent@xmin,map.extent@xmax)
   plot.window.ylim <- c(map.extent@ymin,map.extent@ymax)


  width <- 1400
  height <- 900
  ##png(file=plot.file,width=width,height=height,bg='white') 
##  png(file=plot.file,width=9,height=6,units='in',res=600,pointsize=6,bg='white')
##  par(mar=c(6,6,6,5))
##   par(mar=c(5,5.5,2,2))

  plot(c(),xlim=plot.window.xlim,ylim=plot.window.ylim,xaxs='i',yaxs='i',
     bg='white',# 'lightgray',
     xlab='Longitude (\u00B0E)',ylab='Latitude (\u00B0N)',main=plot.title,
     cex.axis=2.2,cex.lab=2.2,cex.main=2.4,mgp=c(3.5,2,0))
     rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='lightsteelblue')

     shape.dir <- '/storage/data/gis/basedata/base_layers/'
     bc.overlay <- 'bc_province_wgs84'
     bc.shp <- readOGR(shape.dir, bc.overlay, stringsAsFactors=F, verbose=F)
     us.shp <- readOGR(shape.dir, 'united_states', stringsAsFactors=F, verbose=F)
     plot(spTransform(bc.shp,CRS("+init=epsg:4326")),add=TRUE,col='gray')
     image(raster.object, col=colour.ramp,breaks=class.breaks,ribbon=FALSE,xlim=plot.window.xlim, ylim=plot.window.ylim, add=TRUE,alpha=0.3)

     shape.dir <- paste0('/storage/data/projects/rci/data/assessments/metro_van/shapefiles/')
     region <- 'metro_van'
     region.shp <- spTransform(get.region.shape(region,shape.dir),CRS("+init=epsg:4326")) 
     whistler.dir <- '/storage/data/projects/rci/data/assessments/whistler/shapefiles/'
     whistler.shp  <- spTransform(get.region.shape('WhistlerLandscapeUnit',whistler.dir),CRS("+init=epsg:4326"))
 
     rivers.shp <- spTransform(get.region.shape('van_whistler_rivers',shape.dir),CRS("+init=epsg:4326"))  

     glacier.dir <- '/storage/data/gis/basedata/randolph_glacier_inventory/v32/02_rgi32_WesternCanadaUS'
     glacier.name <- '02_rgi32_WesternCanadaUS'
     glacier.shp <- spTransform(get.region.shape(glacier.name,glacier.dir),CRS("+init=epsg:4326"))
     
     lakes.shp <- readOGR('/storage/data/projects/rci/data/assessments/metro_van/shapefiles','van_whistler_lakes',stringsAsFactors=F, verbose=F)
     plot(spTransform(lakes.shp,CRS("+init=epsg:4326")),add=TRUE,col='lightsteelblue',border='lightsteelblue',xlim=plot.window.xlim,ylim=plot.window.ylim)
     plot(rivers.shp,add=TRUE,col='lightsteelblue',border='lightsteelblue',xlim=plot.window.xlim,ylim=plot.window.ylim)
     plot(spTransform(lakes.shp,CRS("+init=epsg:4326")),add=TRUE,border='gray',xlim=plot.window.xlim,ylim=plot.window.ylim)

     if (glaciers) {     
        plot(glacier.shp,add=TRUE,xlim=plot.window.xlim,ylim=plot.window.ylim,col='white',border='white')
     }
     plot(spTransform(us.shp,CRS("+init=epsg:4326")),add=TRUE,col='gray',xlim=plot.window.xlim,ylim=plot.window.ylim)

     ocean.shp <- readOGR('/storage/data/projects/rci/data/assessments/shapefiles/bc_common','west_coast_ocean',stringsAsFactors=F, verbose=F)
     plot(spTransform(ocean.shp,CRS("+init=epsg:4326")),add=TRUE,col='lightsteelblue',border='lightsteelblue')
     plot(spTransform(ocean.shp,CRS("+init=epsg:4326")),add=TRUE,border='black')
     ###plot(spTransform(region.shp,CRS("+init=epsg:4326")),add=TRUE,border='black',cex=0.5)
     ###plot(spTransform(whistler.shp,CRS("+init=epsg:4326")),add=TRUE,border='black',cex=0.5)


     if (bias) {
        load('/storage/data/projects/rci/data/winter_sports/snodas.bias.RData')
        col.ix <- unlist(sapply(snodas.bias$bias,function(x,y){tail(which((y - x) < 0),1)},class.breaks))
        points(snodas.bias$lons[1:13],snodas.bias$lats[1:13],col='black',pch=18,cex=2.5)                
        points(snodas.bias$lons[14:17],snodas.bias$lats[14:17],col='black',pch=17,cex=2.5)                
##        points(snodas.bias$lons,snodas.bias$lats,col=colour.ramp[col.ix],pch=18,cex=1.6)                
     }

     legend('topright', col = "black", legend=rev(map.class.breaks.labels), pch=22, pt.bg = rev(colour.ramp), 
         pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=leg.title, xjust=0, cex=1.9)
     box(which='plot',lwd=2)
##     dev.off()
}

