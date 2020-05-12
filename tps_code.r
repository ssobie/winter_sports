#!/usr/bin/Rscript

rm(list=ls())

library(fields)
library(raster)
library(session)
options(warn=-1)

################################################################################

sfc <- '_tx' # Surface variable to be gridded
method <- 'tps' # ptps or tps
frac.df <- 1 # maximum fraction of no. of observations allowed as model df
frac.df.fixed <- NA # NA; !is.na(frac.df.fixed) sets a fixed model df


z.extent <- c(-169, 40, -101, 72) # Full domain c(-179, 24.5, -100, 80)
z.var <- 'tmax' # z coordinate name in station metadata
z.dir <- 'CWNA_5.10.tif/' # directory with raster files for z coordinate
z.ext <- '.tif' # raster file extension
z.div <- 20 # Divisor for z (_p = prec 300; _tx = tmax 20; _tn = tmin 20)
fact <- 1 # if > 1, coarsen raster surfaces prior to modelling
plot.sfc <- FALSE # Plot modelled surfaces
save.sfc <- TRUE # Save modelled surfaces
NAflag <- -1e34 # missing value
################################################################################

# Interpolation extents
wna.extent <- matrix(z.extent, ncol=2, byrow=TRUE)
colnames(wna.extent) <- c('lon', 'lat')

# Station data
wna <- as.matrix(read.csv('subset.cwna-20CR2.1945-2012.gmted.csv', header=TRUE,
                          as.is=TRUE, check.names=FALSE))[,-c(1:3)]
wna.dates <- seq.Date(as.Date('1945-1-1'), as.Date('2012-12-31'), by='day')

# Station metadata
wna.meta <- read.csv('subset.cwna-20CR2.meta.gmted.climatewna.txt', head=TRUE,
                     as.is=TRUE, check.names=FALSE)
wna.stations <- grepl(sfc, wna.meta[,'station'], fixed=TRUE)
wna.stations <- wna.stations & ((wna.meta[,'lon'] >= wna.extent[1,'lon']) &
                                (wna.meta[,'lon'] <= wna.extent[2,'lon']))
wna.stations <- wna.stations & ((wna.meta[,'lat'] >= wna.extent[1,'lat']) &
                                (wna.meta[,'lat'] <= wna.extent[2,'lat']))
wna.x <- as.numeric(wna.meta[wna.stations,'lon'])
wna.y <- as.numeric(wna.meta[wna.stations,'lat'])

# Raster and station z coordinates
ras.mn <- wna.xyz.mn <- list()
for(mn in 1:12){
    # Raster grids
    if(z.var=='alt')
        ras <- raster(paste0(z.dir, z.var, z.ext))
    else{
        ras <- raster(paste0(z.dir, z.var, mn, z.ext))
    }
    if(fact > 1) ras <- aggregate(ras, fact=fact)
    if(mn==1){
        wna.extent <- wna.extent + matrix(c(-xres(ras), xres(ras),
                                            -yres(ras), yres(ras)), ncol=2)
    }
    ras <- crop(ras, extent(c(wna.extent[,1], wna.extent[,2])))
    ras.x <- xFromCol(ras, 1:ncol(ras))
    ras.y <- yFromRow(ras, 1:nrow(ras))
    ras.y <- rev(ras.y)
    ras.z <- as.matrix(ras[,,1,drop=FALSE])/z.div
    ras.z <- t(ras.z)[,length(ras.y):1]
    ras.mask <- ras.z*0+1
    ras.z[is.na(ras.z)] <- mean(ras.z, na.rm=TRUE)
    ras.mn[[mn]] <- ras.z
    # Station metadata
    if(z.var=='alt')
        wna.z <- as.numeric(wna.meta[wna.stations,z.var])/z.div
    else{
        wna.z <- as.numeric(wna.meta[wna.stations,paste0(z.var, mn)])/z.div
    }
    wna.xyz <- cbind(wna.x, wna.y, wna.z)
    colnames(wna.xyz) <- c('lon', 'lat', z.var)
    rownames(wna.xyz) <- wna.meta[wna.stations,'station']
    wna.xyz.mn[[mn]] <- wna.xyz
}
ras.xy <- cbind(c(matrix(ras.x, nrow=nrow(ras.z), ncol=ncol(ras.z))),
                c(matrix(ras.y, nrow=nrow(ras.z), ncol=ncol(ras.z),
                         byrow=TRUE)))

################################################################################

# Data for modelling
date.i <- date.seq[i]
mn.i <- as.integer(format(date.i, '%m'))
wna.var <- wna[wna.dates==date.i,wna.stations]
wna.xyz <- wna.xyz.mn[[mn.i]]
ras.z <- ras.mn[[mn.i]]

# Screen missing values
wna.var.cases <- !is.na(wna.var)
wna.var <- wna.var[wna.var.cases]
wna.xyz <- wna.xyz[wna.var.cases,]

# Trivariate thin plate spline

    # Minimum and maximum temperature
    tps.df <- NA
    if(!is.na(frac.df.fixed))
        tps.df <- frac.df.fixed*length(wna.var)
    fit.tps <- Tps(x=wna.xyz, Y=wna.var, scale.type='unscaled', df=tps.df)
    tps.df <- fit.tps$eff.df
    if(tps.df > frac.df*length(wna.var)) tps.df <- frac.df*length(wna.var)
    tps.grid <- predict(fit.tps, x=cbind(ras.xy, c(ras.z)),
                        df=tps.df)

tps.z <- matrix(tps.grid, nrow=nrow(ras.z), ncol=ncol(ras.z))*ras.mask
tps.z[tps.z < quantile(tps.z, 0.0001, na.rm=TRUE)] <-
    quantile(tps.z, 0.0001, na.rm=TRUE)-0.01
tps.xyz <- list(x=ras.x, y=ras.y, z=tps.z)


# Save raster surface
if(save.sfc){
    writeRaster(raster(tps.xyz), file=paste0('tps.sfc/tx/', date.i, sfc, '.nc'),
                overwrite=TRUE, NAflag=NAflag)
}

################################################################################

# Plots and diagnostics
if(plot.sfc){
    plot(raster(tps.xyz), col=tim.colors(256), useRaster=TRUE,
         main=paste0(date.i, sfc), xlab='Lon', ylab='Lat')
    points(wna.x, wna.y, pch=20, cex=0.2)
    grid()
}
print(date.i)
print(round(quantile(wna.var), 1))
print(round(quantile(tps.z, na.rm=TRUE), 1))



################################################################################
