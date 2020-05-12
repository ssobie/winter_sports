##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(plotrix)
library(TeachingDemos)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/site.coordinates.r',chdir=T)

##Modified version of Taylor diagram to enable better size control

taylor2_diagram <- function(ref, model, add = FALSE, col = "red", pch = 19, pos.cor = TRUE, 
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
    maxsd <- 2.0 ##1.5 * max(sd.f, sd.r)
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
                0),lwd=1.25)
            print('At segments')
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
    segments(c(0, 0), c(0, 0), c(0, maxsd), c(maxsd, 
    0),lwd=1.25)
    xcurve <- cos(seq(0, pi/2, by = 0.01)) * maxsd
    ycurve <- sin(seq(0, pi/2, by = 0.01)) * maxsd
    lines(xcurve, ycurve,lwd=1.25)

    
    #print(paste0('TD SD: ',round(sd.f,2)))
    #print(paste0('TD Cor: ',round(R,2)))
###   points(sd.f * R, sd.f * sin(acos(R)), pch = pch, col = 'black',bg=col, 
###        cex = pcex)
    shadowtext(sd.f * R, sd.f * sin(acos(R)), pch, col = col,bg='black', 
               cex = 1.2,r=0.15)

    invisible(oldpar)
    return(R)
}

##--------------------------------------

read_snow_sim <- function(site,model,reanalysis,type,model.dir) {

    ##PNWNAmet PRISM calibration
    swe.file <- paste0(model.dir,site,'_',model,'_',reanalysis,'_',type,'_snow_model_data.csv')
    swe.data <- read.csv(swe.file,header=T,as.is=T)
    swe.values <- swe.data$SWE*1000
    swe.dates <- as.Date(swe.data$Dates)
    rv <- list(dates=swe.dates,swe=swe.values)
    return(rv)
}

##---------------------------------------
##Observation data
read_course_obs <- function(site) {

   obs.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_courses/',site,'_snow_course.csv',sep='')
   obs.data <- read.csv(obs.file,header=T,as.is=T)
   obs.dates <- format(as.Date(obs.data[,1]),'%Y-%m-%d')
   obs.swe <- obs.data[,3] ##mm
   obs.na <- is.na(obs.swe)
   obs.swe <- obs.swe[!obs.na]
   obs.dates <- as.Date(obs.dates[!obs.na])

   rv <- list(dates=obs.dates,swe=obs.swe)
   return(rv)
}

read_pillow_obs <- function(site) {

   obs.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_pillow/',site,'.csv',sep='')
   obs.data <- read.csv(obs.file,header=T,as.is=T)
   obs.dates <- format(as.Date(obs.data[,2]),'%Y-%m-%d')
   obs.swe <- obs.data[,11] ##mm
   obs.na <- is.na(obs.swe)
   obs.swe <- obs.swe[!obs.na]
   obs.dates <- as.Date(obs.dates[!obs.na])
   rv <- list(dates=obs.dates,swe=obs.swe)
   return(rv)
}

##--------------------------------------

make_taylor_diagram <- function(active.courses,active.letters,
                                inactive.courses,inactive.letters,
                                pillow.sites,pillow.letters,
                                model,type,colour,model.dir) {
   cor.vals <- c()
   course.nums <- c()
   for (i in seq_along(active.courses)) {
      site <- active.courses[i]
      #print(site)
      coords <- get_coordinates(site)
      course.obs <- read_course_obs(site)
      snow.sim <- read_snow_sim(site,model,reanalysis='PNWNAmet',type,model.dir)

      date.subset <- snow.sim$dates %in% course.obs$dates
      obs.subset <- course.obs$dates %in% snow.sim$dates
      
      if (i==1) { ## & model=='PNWNAmet') {
        cr <- taylor2_diagram(course.obs$swe[obs.subset],snow.sim$swe[date.subset],sd.arcs=TRUE,normalize=T,
                   main=model,col='green',pcex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=1.75,pch=active.letters[i]) ##24                 
      } else {
        cr <- taylor2_diagram(course.obs$swe[obs.subset],snow.sim$swe[date.subset],sd.arcs=TRUE,normalize=T,
                   add=TRUE,col='green',pcex=1.5,pch=active.letters[i]) ##24
      }
      cor.vals <- c(cor.vals,cr)
      course.nums <- c(course.nums,length(course.obs$swe[obs.subset]))
    }

   for (i in seq_along(inactive.courses)) {
      site <- inactive.courses[i]
      print(site)
      coords <- get_coordinates(site)
      course.obs <- read_course_obs(site)
      snow.sim <- read_snow_sim(site,model,reanalysis='PNWNAmet',type,model.dir)

      date.subset <- snow.sim$dates %in% course.obs$dates
      obs.subset <- course.obs$dates %in% snow.sim$dates
      ##if (site=='wolverine_creek') {browser()}
      cr <- taylor2_diagram(course.obs$swe[obs.subset],snow.sim$swe[date.subset],sd.arcs=TRUE,normalize=T,
                      add=TRUE,col='red',pcex=1.5,pch=inactive.letters[i]) ##25
      cor.vals <- c(cor.vals,cr)
      course.nums <- c(course.nums,length(course.obs$swe[obs.subset]))
    }

   for (j in seq_along(pillow.sites)) {
      site <- pillow.sites[j]
      print(site)
      coords <- get_coordinates(site)
      pillow.obs <- read_pillow_obs(site)
      snow.sim <- read_snow_sim(site,model,reanalysis='PNWNAmet',type,model.dir)

      date.subset <- snow.sim$dates %in% pillow.obs$dates
      obs.subset <- pillow.obs$dates %in% snow.sim$dates

      cr <- taylor2_diagram(pillow.obs$swe[obs.subset],snow.sim$swe[date.subset],sd.arcs=TRUE,normalize=T,
                      add=TRUE,col='orange',pcex=1.5,pch=pillow.letters[j]) ##23
      cor.vals <- c(cor.vals,cr)
      course.nums <- c(course.nums,length(pillow.obs$swe[obs.subset]))
    }

    ##print('Inputs')
    ##print(paste0('NCEP2 SD: ',round(sd(ncep2.swe.mean[ncep2.date.subset] / sd(obs.swe[ncep2.obs.subset])),2)))
    ##print(paste0('NCEP2 Cor: ',round(cor(obs.swe[ncep2.obs.subset],ncep2.swe.mean[ncep2.date.subset]),2)))

    ##print(paste0('ERA SD: ',round(sd(era.swe.mean[era.date.subset] / sd(obs.swe[era.obs.subset])),2)))
    ##print(paste0('ERA Cor: ',round(cor(obs.swe[era.obs.subset],era.swe.mean[era.date.subset]),2)))

    return(list(cor=cor.vals,obs=course.nums))
}

##-----------------------------------------------------------
##Snow Courses

if (1==1) {


}
 
##Old Sites 
if (1==0) {
course.sites <- c('shovelnose_mountain','brookmere','lightning_lake','callaghan','orchid_lake','palisade_lake',
                  'grouse_mountain','dog_mountain','stave_lake','nahatlatch','wahleach','klesilkwa','hamilton_hill')
course.letters <- c('M','B','L','C','O','P','G','D','S','N','W','K','H')
pillow.sites <- c('chilliwack_river','spuzzum_creek','tenquille_lake','upper_squamish')
pillow.letters <- c('R','Z','T','U')

}



model.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow_series/'

##Loop over sites
plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'
type <- 'PRISM_TPS'
plot.file <- paste0(plot.dir,'pnwnamet.era5.',type,'.swe.taylor.diagram.2020.letters.png')

png(file=plot.file,width=10,height=5,units='in',res=600,pointsize=6,bg='white')

par(mfrow=c(1,2))

##------------------------------------------

##Locations fo PNWNAmet Evaluation

cal.courses <- c('brookmere','callaghan','dickson_lake','disappointment_lake','dog_mountain','duffey_lake',
                    'gnawed_mountain','great_bear','grouse_mountain','hamilton_hill','highland_valley',
                    'klesilkwa','lightning_lake','mcgillivray_pass','nahatlatch','orchid_lake','palisade_lake',
                    'shovelnose_mountain','stave_lake','sumallo_river_west','wahleach')
cal.letters <- c('BR','CG','DL',
                 'DI','DM','DU','GM','GB',
                 'GR','HH','HV','KL','LL',
                 'MP','NA','OL',
                 'PL','SM',
                 'SL','SW','WA')             

eval.courses <- c('blackwall_peak_course','boston_bar_lower','boston_bar_upper','burwell_lake',
                  'chapman_creek','cornwall_hills','diamond_head','edwards_lake',
                  'hollyburn','hope',
                  'loch_lomond','lytton','mount_seymour','new_tashme',
                  'ottomite','pavilion_mountain',##'shalalth',
                      'sumallo_river','tenquille_course','whistler_mountain','wolverine_creek' )
eval.letters <- c('BP','BL','BU','BW','CC','CH','DH','EL',
                  'HB','HO','LO','LY','MS','NT','OT','PM',
                  'SR','TL','WM','WC') ##'SH'

pillow.sites <- c('blackwall_peak_pillow','chilliwack_river','spuzzum_creek','tenquille_lake','upper_squamish','wahleach_lake')
pillow.letters <- c('BP','CR','SC','TQ','US','WC')


pnw.cor <- make_taylor_diagram(cal.courses,cal.letters,
                    eval.courses,eval.letters,
                    pillow.sites,pillow.letters,
                    model='PNWNAmet',type=type,'blue',model.dir)

print('-----')
print('ERA5')

##Sites for ERA5 Evaluation
cal.courses <- c('brookmere','callaghan','dickson_lake','disappointment_lake','dog_mountain','duffey_lake',
                    'gnawed_mountain','great_bear','grouse_mountain','hamilton_hill','highland_valley',
                    'klesilkwa','lightning_lake','mcgillivray_pass','nahatlatch','orchid_lake','palisade_lake',
                    'shovelnose_mountain','stave_lake','sumallo_river_west','wahleach')
cal.letters <- c('BR','CG','DL',
                 'DI','DM','DU','GB','GM',
                 'GR','HH','HV','KL','LL',
                 'MP','NA','OL',
                 'PL','SM',
                 'SL','SW','WA')             

eval.courses <- c('blackwall_peak_course','boston_bar_lower','boston_bar_upper','burwell_lake',
                  'chapman_creek','cornwall_hills','diamond_head','edwards_lake',
                  'hollyburn',
                  'loch_lomond','mount_seymour','new_tashme',
                  'ottomite','pavilion_mountain',##'shalalth',
                      'sumallo_river','tenquille_course','whistler_mountain','wolverine_creek' )
eval.letters <- c('BP','BL','BU','BW','CC','CH','DH','EL',
                  'HB','LO','MS','NT','OT','PM',
                  'SR','TL','WM','WC') ##'SH'

pillow.sites <- c('blackwall_peak_pillow','chilliwack_river','spuzzum_creek','tenquille_lake','upper_squamish','wahleach_lake')
pillow.letters <- c('BP','CR','SC','TQ','US','WC')

##Excluded due to lack of dates
##Hope, Lytton, Garibaldi Lake

era.cor <- make_taylor_diagram(cal.courses,cal.letters,
                    eval.courses,eval.letters,
                    pillow.sites,pillow.letters,
                    model='ERA5',type,'darkgreen',model.dir)

##legend('topright',legend=c('Cal.','Eval.','Pillow'),col='black',pt.bg=c('green','red','orange'),pch=c(24,25,23),cex=1.5)
legend('topright',legend=c('Cal.','Eval.','Pillow'),col='black',pt.bg=c('green','red','orange'),pch=c(22,22,22),pt.cex=2.0,cex=1.5)


##make_taylor_diagram(course.sites,course.letters,
##                    pillow.sites,pillow.letters,
##                    model='PNWNAmet',type='PRISM_with_elevation','blue',model.dir)
##make_taylor_diagram(course.sites,course.letters,
##                    pillow.sites,pillow.letters,
##                    model='ERA5',type='PRISM_with_elevation','green',model.dir)
##legend('topright',legend=c('PNWNAmet','ERA5'),col=c('blue','green'),pch=16,cex=1.5)





dev.off()

##




