##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(plotrix)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/site.coordinates.r',chdir=T)



##Modified version of Taylor diagram to enable better size control

taylor2.diagram <- function(ref, model, add = FALSE, col = "red", pch = 19, pos.cor = TRUE, 
    xlab = "", ylab = "", main = "Taylor Diagram", show.gamma = FALSE, 
    ngamma = 3, gamma.col = 8, sd.arcs = 0, ref.sd = FALSE, sd.method = "sample", 
    grad.corr.lines = c(0.2, 0.4, 0.6, 0.8, 0.9), pcex = 1, cex.axis = 1, 
    normalize = FALSE, mar = c(5, 5, 6, 6), ...) {

    grad.corr.full <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 
        1)
    R <- cor(ref, model, use = "pairwise")
    if (is.list(ref)) 
        ref <- unlist(ref)
    if (is.list(model)) 
        ref <- unlist(model)
    SD <- function(x, subn) {
        meanx <- mean(x, na.rm = TRUE)
        devx <- x - meanx
        ssd <- sqrt(sum(devx * devx, na.rm = TRUE)/(length(x[!is.na(x)]) - 
            subn))
        return(ssd)
    }
    subn <- sd.method != "sample"
    sd.r <- SD(ref, subn)
    sd.f <- SD(model, subn)
    if (normalize) {
        sd.f <- sd.f/sd.r
        sd.r <- 1
    }
    maxsd <- 1.5 * max(sd.f, sd.r)
    oldpar <- par("mar", "xpd", "xaxs", "yaxs")
    if (!add) {
        if (pos.cor) {
            if (nchar(ylab) == 0) 
                ylab = "Standard deviation"
            par(mar = mar)
            plot(0, xlim = c(0, maxsd), ylim = c(0, maxsd), xaxs = "i", 
                yaxs = "i", axes = FALSE, main = main, xlab = xlab, 
                ylab = ylab, type = "n", cex = cex.axis, ...)
            if (grad.corr.lines[1]) {
                for (gcl in grad.corr.lines) lines(c(0, maxsd * 
                  gcl), c(0, maxsd * sqrt(1 - gcl^2)), lty = 3)
            }
            segments(c(0, 0), c(0, 0), c(0, maxsd), c(maxsd, 
                0))
            axis.ticks <- pretty(c(0, maxsd))
            axis.ticks <- axis.ticks[axis.ticks <= maxsd]
            axis(1, at = axis.ticks, cex.axis = cex.axis)
            axis(2, at = axis.ticks, cex.axis = cex.axis)
            if (sd.arcs[1]) {
                if (length(sd.arcs) == 1) 
                  sd.arcs <- axis.ticks
                for (sdarc in sd.arcs) {
                  xcurve <- cos(seq(0, pi/2, by = 0.03)) * sdarc
                  ycurve <- sin(seq(0, pi/2, by = 0.03)) * sdarc
                  lines(xcurve, ycurve, col = "blue", lty = 3)
                }
            }
           if (show.gamma[1]) {
                if (length(show.gamma) > 1) 
                  gamma <- show.gamma
                else gamma <- pretty(c(0, maxsd), n = ngamma)[-1]
                if (gamma[length(gamma)] > maxsd) 
                  gamma <- gamma[-length(gamma)]
                labelpos <- seq(45, 70, length.out = length(gamma))
                for (gindex in 1:length(gamma)) {
                  xcurve <- cos(seq(0, pi, by = 0.03)) * gamma[gindex] + 
                    sd.r
                  endcurve <- which(xcurve < 0)
                  endcurve <- ifelse(length(endcurve), min(endcurve) - 
                    1, 105)
                  ycurve <- sin(seq(0, pi, by = 0.03)) * gamma[gindex]
                  maxcurve <- xcurve * xcurve + ycurve * ycurve
                  startcurve <- which(maxcurve > maxsd * maxsd)
                  startcurve <- ifelse(length(startcurve), max(startcurve) + 
                    1, 0)
                  lines(xcurve[startcurve:endcurve], ycurve[startcurve:endcurve], 
                    col = gamma.col)
                  if (xcurve[labelpos[gindex]] > 0) 
                    boxed.labels(xcurve[labelpos[gindex]], ycurve[labelpos[gindex]], 
                      gamma[gindex], border = FALSE)
                }
            }
            xcurve <- cos(seq(0, pi/2, by = 0.01)) * maxsd
            ycurve <- sin(seq(0, pi/2, by = 0.01)) * maxsd
            lines(xcurve, ycurve)
            bigtickangles <- acos(seq(0.1, 0.9, by = 0.1))
            medtickangles <- acos(seq(0.05, 0.95, by = 0.1))
            smltickangles <- acos(seq(0.91, 0.99, by = 0.01))
            segments(cos(bigtickangles) * maxsd, sin(bigtickangles) * 
                maxsd, cos(bigtickangles) * 0.97 * maxsd, sin(bigtickangles) * 
                0.97 * maxsd)
            par(xpd = TRUE)
            if (ref.sd) {
                xcurve <- cos(seq(0, pi/2, by = 0.01)) * sd.r
                ycurve <- sin(seq(0, pi/2, by = 0.01)) * sd.r
                lines(xcurve, ycurve)
            }
            points(sd.r, 0, cex = pcex)
            text(cos(c(bigtickangles, acos(c(0.95, 0.99)))) * 
                1.05 * maxsd, sin(c(bigtickangles, acos(c(0.95, 
                0.99)))) * 1.05 * maxsd, c(seq(0.1, 0.9, by = 0.1), 
                0.95, 0.99),cex=1.5)
            text(maxsd * 0.8, maxsd * 0.8, "Correlation", srt = 315,cex=1.5)
            segments(cos(medtickangles) * maxsd, sin(medtickangles) * 
                maxsd, cos(medtickangles) * 0.98 * maxsd, sin(medtickangles) * 
                0.98 * maxsd)
            segments(cos(smltickangles) * maxsd, sin(smltickangles) * 
                maxsd, cos(smltickangles) * 0.99 * maxsd, sin(smltickangles) * 
                0.99 * maxsd)
        }
        else {
            x <- ref
            y <- model
            R <- cor(x, y, use = "pairwise.complete.obs")
            E <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
            xprime <- x - mean(x, na.rm = TRUE)
            yprime <- y - mean(y, na.rm = TRUE)
            sumofsquares <- (xprime - yprime)^2
            Eprime <- sqrt(sum(sumofsquares)/length(complete.cases(x)))
            E2 <- E^2 + Eprime^2
            if (add == FALSE) {
                maxray <- 1.5 * max(sd.f, sd.r)
                plot(c(-maxray, maxray), c(0, maxray), type = "n", 
                  asp = 1, bty = "n", xaxt = "n", yaxt = "n", 
                  xlab = xlab, ylab = ylab, main = main, cex = cex.axis)
                discrete <- seq(180, 0, by = -1)
                listepoints <- NULL
                for (i in discrete) {
                  listepoints <- cbind(listepoints, maxray * 
                    cos(i * pi/180), maxray * sin(i * pi/180))
                }
                listepoints <- matrix(listepoints, 2, length(listepoints)/2)
                listepoints <- t(listepoints)
                lines(listepoints[, 1], listepoints[, 2])
                lines(c(-maxray, maxray), c(0, 0))
                lines(c(0, 0), c(0, maxray))
                for (i in grad.corr.lines) {
                  lines(c(0, maxray * i), c(0, maxray * sqrt(1 - 
                    i^2)), lty = 3)
                  lines(c(0, -maxray * i), c(0, maxray * sqrt(1 - 
                    i^2)), lty = 3)
                }
                for (i in grad.corr.full) {
                  text(1.05 * maxray * i, 1.05 * maxray * sqrt(1 - 
                    i^2), i, cex = 0.6)
                  text(-1.05 * maxray * i, 1.05 * maxray * sqrt(1 - 
                    i^2), -i, cex = 0.6)
                }
                seq.sd <- seq.int(0, 2 * maxray, by = (maxray/10))[-1]
                for (i in seq.sd) {
                  xcircle <- sd.r + (cos(discrete * pi/180) * 
                    i)
                  ycircle <- sin(discrete * pi/180) * i
                  for (j in 1:length(xcircle)) {
                    if ((xcircle[j]^2 + ycircle[j]^2) < (maxray^2)) {
                      points(xcircle[j], ycircle[j], col = "darkgreen", 
                        pch = ".")
                      if (j == 10) 
                        text(xcircle[j], ycircle[j], signif(i, 
                          2), cex = 1.5, col = "darkgreen")
                    }
                  }
                }
                seq.sd <- seq.int(0, maxray, length.out = 5)
                for (i in seq.sd) {
                  xcircle <- (cos(discrete * pi/180) * i)
                  ycircle <- sin(discrete * pi/180) * i
                  if (i) 
                    lines(xcircle, ycircle, lty = 3, col = "blue")
                  text(min(xcircle), -0.03 * maxray, signif(i, 
                    2), cex = 1.5, col = "blue")
                  text(max(xcircle), -0.03 * maxray, signif(i, 
                    2), cex = 1.5, col = "blue")
                }
                text(0, -0.08 * maxray, "Standard Deviation", 
                  cex = 0.7, col = "blue")
                text(0, -0.12 * maxray, "Centered RMS Difference", 
                  cex = 0.7, col = "darkgreen")
                points(sd.r, 0, pch = 22, bg = "darkgreen", cex = 1.1)
                text(0, 1.1 * maxray, "Correlation Coefficient", 
                  cex = 0.7)
            }
            S <- (2 * (1 + R))/(sd.f + (1/sd.f))^2
        }
    }
    
    print(paste0('TD SD: ',round(sd.f,2)))
    print(paste0('TD Cor: ',round(R,2)))
    points(sd.f * R, sd.f * sin(acos(R)), pch = pch, col = col, 
        cex = pcex)

    invisible(oldpar)
}

##-----------------------------------------------------------
##Snow Courses

if (1==1) {

##Old Sites 
if (1==0) {
sites <- c('shovelnose_mountain','brookmere','lightning_lake','callaghan','orchid_lake','palisade_lake',
           'grouse_mountain','dog_mountain','stave_lake','nahatlatch','wahleach','klesilkwa','hamilton_hill',
           'chilliwack_river','spuzzum_creek','tenquille_lake','upper_squamish')
obs.type <- c(rep('course',13),rep('asp',4))
site.letter <- c('M','B','L','C','O','P','G','D','S','N','W','K','H','R','Z','T','U')
}

'blackwall_peak_pillow',

  sites <- c('blackwall_peak_course','boston_bar_lower','boston_bar_upper','brookmere','burwell_lake',
             'callaghan','chapman_creek','chilliwack_river','cornwall_hills','diamond_head','dickson_lake',
             'disappointment_lake','dog_mountain','duffey_lake','edwards_lake','garibaldi_lake','gnawed_mountain',
             'grouse_mountain','hamilton_hill','highland_valley','hollyburn','hope','klesilkwa','lightning_lake',
             'loch_lomond','lytton','mcgillivray_pass','mount_seymour','nahatlatch','new_tashme','orchid_lake',
             'ottomite','palisade_lake','pavilion_mountain','shalalth','shovelnose_mountain','spuzzum_creek',
             'stave_lake','sumallo_river','sumallo_river_west','tenquille_lake','tenquille_course','upper_squamish',
             'wahleach','wahleach_lake','whistler_mountain','wolverine_creek')






model.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/calibrated/series/'

ncep2.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/grouse_mountain_NCEP2_800m_data.csv'

course.site.swe <- vector(mode='list',length=length(sites))
ncep2.site.swe <- vector(mode='list',length=length(sites))
era.site.swe <- vector(mode='list',length=length(sites))

course.site.pack <- vector(mode='list',length=length(sites))
ncep2.site.pack <- vector(mode='list',length=length(sites))
era.site.pack <- vector(mode='list',length=length(sites))

ncep2.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/grouse_mountain_NCEP2_800m_data.csv'
ncep2.data <- read.csv(ncep2.file,header=T,as.is=T)

ncep2.swe.sims <- matrix(0,nrow=10,ncol=dim(ncep2.data)[1])
ncep2.snow.sims <- matrix(0,nrow=10,ncol=dim(ncep2.data)[1])

era.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/grouse_mountain_ERA_800m_data.csv'
era.data <- read.csv(era.file,header=T,as.is=T)

era.swe.sims <- matrix(0,nrow=10,ncol=dim(era.data)[1])
era.snow.sims <- matrix(0,nrow=10,ncol=dim(era.data)[1])

snodas.file <- "/storage/data/projects/rci/data/winter_sports/obs/SNODAS/ncdf4_files/swe_snodas_modis_grid_van_whistler_20100101-20181231.nc"
snc <- nc_open(snodas.file)
lon <- ncvar_get(snc,'lon')
lat <- ncvar_get(snc,'lat')
snodas.dates <- as.character(netcdf.calendar(snc))

##Loop over sites
plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'
type <- 'SWE'
png(file=paste0(plot.dir,'ncep2.era.swe.',type,'.taylor.diagram.2020.png'),width=800,height=800)

for (i in seq_along(sites)) {
    site <- sites[i]
    print(site)
    ##Reanalysis 800m data
    ncep2.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_NCEP2_800m_data.csv')
    ncep2.data <- read.csv(ncep2.file,header=T,as.is=T)
    ncep2.pr <- ncep2.data$Pr
    ncep2.tasmax <- ncep2.data$Tasmax
    ncep2.tasmin <- ncep2.data$Tasmin
    ncep2.tas <- ncep2.data$Tas
    ncep2.dates <- ncep2.data$Dates

    era.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_ERA_800m_data.csv')
    era.data <- read.csv(era.file,header=T,as.is=T)
    era.pr <- era.data$Pr
    era.tasmax <- era.data$Tasmax
    era.tasmin <- era.data$Tasmin
    era.tas <- era.data$Tas
    era.dates <- era.data$Dates

    coords <- get.coordinates(site)
    lat.bnds <- coords[2]
    elev <- coords[3]


    ##Observation data
    if (obs.type[i] == 'course') {
      obs.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_courses/',site,'_snow_course.csv',sep='')
      obs.data <- read.csv(obs.file,header=T,as.is=T)
      obs.dates <- format(as.Date(obs.data[,1]),'%Y-%m-%d')
      obs.swe <- obs.data[,3] ##mm
      obs.na <- is.na(obs.swe)
      obs.pack <- obs.data[,2] ##cm
      obs.dense <-  obs.data[,4]
      obs.swe <- obs.swe[!obs.na]
      obs.dates <- obs.dates[!obs.na]

    } else {
       obs.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_pillow/',site,'_asp.csv',sep='')
       obs.data <- read.csv(obs.file,header=T,as.is=T)
       obs.dates <- format(as.Date(obs.data[,2]),'%Y-%m-%d')
       obs.tasmax <- obs.data[,3]
       obs.tasmin <- obs.data[,5]
       obs.tas <- (obs.tasmax + obs.tasmin)/2
       obs.precip <- obs.data[,7]##mm
       obs.swe <- obs.data[,11] ##mm
       obs.na <- is.na(obs.swe)
       obs.pack <- obs.data[,13] ##cm
       obs.swe <- obs.swe[!obs.na]
       obs.dates <- obs.dates[!obs.na]
    }        

    ##SNODAS Data at Courses
    lon.ix <- which.min(abs(coords[1]-lon))
    lat.ix <- which.min(abs(coords[2]-lat))
    snodas.swe <- ncvar_get(snc,'swe',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
    snodas.date.subset <- format(as.Date(snodas.dates),'%Y-%m-%d') %in% obs.dates
    snodas.obs.subset <- obs.dates %in% format(as.Date(snodas.dates),'%Y-%m-%d')
  

##    print(order(obs.dates) - 1:length(obs.dates))

    ncep2.date.subset <- format(as.Date(ncep2.dates),'%Y-%m-%d') %in% obs.dates
    ncep2.obs.subset <- obs.dates %in% format(as.Date(ncep2.dates),'%Y-%m-%d')
    snodas.ncep2.subset <- format(as.Date(snodas.dates),'%Y-%m-%d') %in% format(as.Date(ncep2.dates),'%Y-%m-%d')
    ncep2.snodas.subset <- format(as.Date(ncep2.dates),'%Y-%m-%d') %in% format(as.Date(snodas.dates),'%Y-%m-%d')

    era.date.subset <- format(as.Date(era.dates),'%Y-%m-%d') %in% obs.dates
    era.obs.subset <- obs.dates %in% format(as.Date(era.dates),'%Y-%m-%d')   
    snodas.era.subset <- format(as.Date(snodas.dates),'%Y-%m-%d') %in% format(as.Date(era.dates),'%Y-%m-%d')
    era.snodas.subset <- format(as.Date(era.dates),'%Y-%m-%d') %in% format(as.Date(snodas.dates),'%Y-%m-%d')

    model.match <- format(as.Date(ncep2.dates),'%Y-%m-%d') %in% format(as.Date(era.dates),'%Y-%m-%d')

    era.swe.sims <- read.csv(paste0(model.dir,site,'_ERA_snow_model_data.csv'),header=T,as.is=T)
    era.swe.mean <- era.swe.sims$SWE*1000 ##apply(era.swe.sims,1,mean,na.rm=T)

    ncep2.swe.sims <- read.csv(paste0(model.dir,site,'_NCEP2_snow_model_data.csv'),header=T,as.is=T) ##[model.match,]
    ncep2.swe.mean <- ncep2.swe.sims$SWE*1000 ###apply(ncep2.swe.sims,1,mean,na.rm=T)

    print('Inputs')
    print(paste0('NCEP2 SD: ',round(sd(ncep2.swe.mean[ncep2.date.subset] / sd(obs.swe[ncep2.obs.subset])),2)))
    print(paste0('NCEP2 Cor: ',round(cor(obs.swe[ncep2.obs.subset],ncep2.swe.mean[ncep2.date.subset]),2)))

    print(paste0('ERA SD: ',round(sd(era.swe.mean[era.date.subset] / sd(obs.swe[era.obs.subset])),2)))
    print(paste0('ERA Cor: ',round(cor(obs.swe[era.obs.subset],era.swe.mean[era.date.subset]),2)))


    if (type=='SNODAS') {
      if (i==1) {
        taylor2.diagram(snodas.swe[snodas.ncep2.subset],ncep2.swe.mean[ncep2.snodas.subset],sd.arcs=TRUE,normalize=T,
                   main='SNODAS SWE Comparison',pch=site.letter[i],col='green',pcex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.75)                  
        taylor2.diagram(snodas.swe[snodas.era.subset],era.swe.mean[era.snodas.subset],sd.arcs=TRUE,normalize=T,
                   pch=site.letter[i],col='blue',add=TRUE,pcex=1.5,cex=1.5)                        
      } else {
        taylor2.diagram(snodas.swe[snodas.ncep2.subset],ncep2.swe.mean[ncep2.snodas.subset],sd.arcs=TRUE,normalize=T,
                   pch=site.letter[i],add=T,col='green',pcex=1.5)
        taylor2.diagram(snodas.swe[snodas.era.subset],era.swe.mean[era.snodas.subset],sd.arcs=TRUE,normalize=T,
                   pch=site.letter[i],col='blue',add=TRUE,pcex=1.5)                  
      }
    } else {
      if (i==1) {
        taylor2.diagram(obs.swe[ncep2.obs.subset],ncep2.swe.mean[ncep2.date.subset],sd.arcs=TRUE,normalize=T,
                   main='',pch=site.letter[i],col='green',pcex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.75)                  
        taylor2.diagram(obs.swe[era.obs.subset],era.swe.mean[era.date.subset],sd.arcs=TRUE,normalize=T,
                   pch=site.letter[i],col='blue',add=TRUE,pcex=1.5,cex=1.5)                     
      } else {
        taylor2.diagram(obs.swe[ncep2.obs.subset],ncep2.swe.mean[ncep2.date.subset],sd.arcs=TRUE,normalize=T,
                   pch=site.letter[i],add=T,col='green',pcex=1.5)
        taylor2.diagram(obs.swe[era.obs.subset],era.swe.mean[era.date.subset],sd.arcs=TRUE,normalize=T,
                   pch=site.letter[i],col='blue',add=TRUE,pcex=1.5)                  
      }
    }
}        

##legend('topright',legend=c('ERA','NCEP2','SNODAS'),col=c('blue','green','red'),pch=16,cex=1.5)
legend('topright',legend=c('ERA','NCEP2'),col=c('blue','green'),pch=16,cex=1.5)

dev.off()

##


}

