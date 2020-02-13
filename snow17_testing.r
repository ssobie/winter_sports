##Snow-17 accumulation and ablation model.

##library(parallel)
library(optimx)

##Based on Anderson (2006) and Mark Raleigh's matlab code.
##Primary Citations:
## 1.  Anderson, E. A. (1973), National Weather Service River Forecast System
## Snow   Accumulation   and   Ablation   Model,   NOAA   Tech.   Memo.   NWS
## HYDro-17, 217 pp., U.S. Dep. of Commer., Silver Spring, Md.
## 2.  Anderson, E. A. (1976), A point energy and mass balance model of a snow
## cover, NOAA Tech. Rep. 19, 150 pp., U.S. Dep. of Commer., Silver Spring, Md.


##-------------------------------------------------------------------------------

snow17 <- function(dates, pr, tas, 
                  lat=50, elevation=0, dt=24, scf=1.0, rvs=1,
                  uadj=0.04, mbase=1.0, mfmax=1.05, mfmin=0.6, tipm=0.1, nmf=0.15,
                  plwhc=0.04, pxtemp=1.0, pxtemp1=-1.0, pxtemp2=3.0) {

    ##Defaults held constant:
    lat <- 50
    elevation <- 1000
    dt <- 24
    rvs <- 1
    pxtemp1 <- 0.0
    pxtemp2 <- 15

    ##Optimized parameters

    scf <- 1 ##par[1]
    mfmax <- 1.05 ##par[2]
    mfmin <- 0.5 ##par[3]
    uadj <- 0.01 #par[4]
    nmf <- 0.15 #par[5]
    mbase <- 0.9 #par[6]
    pxtemp <- 2.0 #par[7]
    plwhc <- 0.04 #par[8]
    tipm <- 0.2 ##par[3] ##par[9]

    ##pr <- data$pr
    ##tas <- data$tas
    ##dates <- data$dates

    # Initialization
    # Antecedent Temperature Index, deg C
    ait <- 0.0
    # Liquid water capacity
    w_qx <- 0.0
    # Liquid water held by the snow (mm)
    w_q <- 0.0
    # accumulated water equivalent of the iceportion of the snow cover (mm)
    w_i <- 0.0
    # Heat deficit, also known as NEGHS, Negative Heat Storage
    deficit <- 0.0

    # number of time steps
    nsteps = length(dates)
    model_swe <- rep(0,nsteps)
    outflow <- rep(0,nsteps)

    fracsnows <- rep(0,length(dates))
    snowfalls <- rep(0,length(dates))
    melting <-  rep(0,length(dates))
    wice <-   rep(0,length(dates))
    wsnow <-   rep(0,length(dates))
    deficits <- rep(0,length(dates))
    meltrain <-   rep(0,length(dates))
    meltnon <- rep(0,length(dates))
    melt1 <- melt2 <- melt3 <- rep(0,length(dates))


    # Stefan-Boltzman constant (mm/K/hr)
    stefan <- 6.12 * 10E-10
    # atmospheric pressure (mb) where elevation is in HUNDREDS of meters
    # (this is incorrectly stated in the manual)
    p_atm <- 33.86 * (29.9 - (0.335 * elevation / 100) +
                     (0.00022 * ((elevation / 100) ^ 2.4)))

    transitionx <- c(pxtemp1, pxtemp2)
    transitiony <- c(1.0, 0.0)

    tipm_dt <- 1.0 - ((1.0 - tipm) ^ (dt / 6))

    # Model Execution
    for (i in seq_along(dates)) {
        date <- dates[i]
        mf <- melt_function(date, dt, lat, mfmax, mfmin)

        # air temperature at this time step (deg C)
        t_air_mean <- tas[i]
        # precipitation at this time step (mm)
        precip <- pr[i]

        # Divide rain and snow
        if (rvs == 0) {
            if (t_air_mean <= pxtemp) {
                # then the air temperature is cold enough for snow to occur
                fracsnow <- 1.0
             } else {
                # then the air temperature is warm enough for rain
                fracsnow <- 0.0
             }
        } else if (rvs == 1) {
            if (t_air_mean <= pxtemp1) {
                fracsnow <- 1.0
            } else if (t_air_mean >= pxtemp2) {
                fracsnow <- 0.0
            } else {
                ## fracsnow <- np.interp(t_air_mean, transitionx, transitiony)
                tas.seq <- seq(pxtemp1,pxtemp2,0.1)
                ix <- which.min(abs(t_air_mean - tas.seq))
                fracsnow <- approx(c(1,0),n=length(tas.seq))$y[ix]
            }
        } else if (rvs == 2) {
            fracsnow <- 1.0
        } else {
            stop('Invalid rain vs snow option')
        }
        fracsnows[i] <- fracsnow
        fracrain <- 1.0 - fracsnow

        # Snow Accumulation
        # water equivalent of new snowfall (mm)
        pn <- precip * fracsnow * scf
        snowfalls[i] <- pn
        # w_i = accumulated water equivalent of the ice portion of the snow
        # cover (mm)
        w_i <- w_i + pn
        e <- 0.0
        # amount of precip (mm) that is rain during this time step
        rain <- fracrain * precip

        # Temperature and Heat deficit from new Snow
        if (t_air_mean < 0.0) {
            t_snow_new <- t_air_mean
            # delta_hd_snow = change in the heat deficit due to snowfall (mm)
            delta_hd_snow <-  -1*(t_snow_new * pn) / (80 / 0.5)
            t_rain <- pxtemp
        } else {
            t_snow_new <- 0.0
            delta_hd_snow <- 0.0
            t_rain <- t_air_mean
        }
        # Antecedent temperature Index
        if (pn > (1.5 * dt)) {
            ait <- t_snow_new
        } else {
            # Antecedent temperature index
            ait <- ait + tipm_dt * (t_air_mean - ait)
        }
        if (ait > 0)
            ait <- 0

        # Heat Exchange when no Surface melt
        # delta_hd_t = change in heat deficit due to a temperature gradient(mm)
        delta_hd_t <- nmf * (dt / 6.0) * ((mf) / mfmax) * (ait - t_snow_new)

        # Rain-on-snow melt
        # saturated vapor pressure at t_air_mean (mb)
        e_sat <- 2.7489 * (10 ^ 8) * exp((-4278.63 / (t_air_mean + 242.792)))
        # 1.5 mm/ 6 hrs
        if (rain > (0.25 * dt)) {
            # melt (mm) during rain-on-snow periods is:
            m_ros1 <- max(c(stefan * 24 * (((t_air_mean + 273) ^ 4) - (273 ^ 4)), 0.0))
            m_ros2 <- max(c((0.0125 * rain * t_rain), 0.0))
            m_ros3 <- max(c((8.5 * uadj * (dt / 6.0) * (((0.9 * e_sat) - 6.11) +
                           (0.00057 * p_atm * t_air_mean))),0.0))
            m_ros <- m_ros1 + m_ros2 + m_ros3
        } else {
            m_ros <- 0.0
            m_ros1 <- m_ros2 <- m_ros3 <- 0
        }
        meltrain[i] <- m_ros
        melt1[i] <- m_ros1
        melt2[i] <- m_ros2
        melt3[i] <- m_ros3

        # Non-Rain melt
        if (rain <= (0.25 * dt) & (t_air_mean > mbase)) {
            # melt during non-rain periods is:
            m_nr <- (mf * (t_air_mean - mbase)) + (0.0125 * rain * t_rain)
        } else {
            m_nr <- 0.0
        }
        meltnon[i] <- m_nr
        # Ripeness of the snow cover
        melt <- m_ros + m_nr
        if (melt <= 0)
            melt <- 0.0

        if (melt < w_i) {
            w_i <- w_i - melt
        } else {
            melt <- w_i + w_q
            w_i <- 0.0 
        }

        # qw = liquid water available melted/rained at the snow surface (mm)
        qw <- melt + rain
        melting[i] <- qw
        # w_qx = liquid water capacity (mm)
        w_qx <- plwhc * w_i
        # deficit = heat deficit (mm)
        deficit <- deficit +  delta_hd_snow + delta_hd_t

        # limits of heat deficit
        if (deficit < 0) {
            deficit <- 0.0
        } else if (deficit > 0.33 * w_i) {
            deficit <-  0.33 * w_i
        }

        # Snow cover is ripe when both (deficit=0) & (w_q = w_qx)
        if (w_i > 0.0) {
            if ((qw + w_q) > ((deficit * (1 + plwhc)) + w_qx)) {
                # THEN the snow is RIPE
                # Excess liquid water (mm)
                e <- qw + w_q - w_qx - (deficit * (1 + plwhc))
                # fills liquid water capacity
                w_q <- w_qx
                # w_i increases because water refreezes as heat deficit is
                # decreased
                w_i <- w_i + deficit
                deficit <- 0.0
            } else if ((qw >= deficit) & ait * ((qw + w_q) <= ((deficit * (1 + plwhc)) + w_qx))) {
                # THEN the snow is NOT yet ripe, but ice is being melted
                e <- 0.0
                w_q <- w_q + qw - deficit
                # w_i increases because water refreezes as heat deficit is
                # decreased
                w_i <- w_i + deficit
                deficit <- 0.0
            } else {
                # (qw < deficit) %elseif ((qw + w_q) < deficit):
                # THEN the snow is NOT yet ripe
                e <- 0.0
                # w_i increases because water refreezes as heat deficit is
                # decreased
                w_i <- w_i + qw
                deficit <- deficit - qw
            }
            swe <- w_i + w_q
        } else {
            e <- qw
            swe <- 0
        }
        if (deficit == 0) 
            ait = 0
        wice[i] <- w_i
        wsnow[i] <- w_q
        deficits[i] <- deficit

        # End of model execution
        model_swe[i] <- swe  # total swe (mm) at this time step
        outflow[i] <- e
    }    
    return(list(swe=model_swe,outflow=outflow,frac=fracsnows,snowfall=snowfalls,
    melt=melting,wice=wice,wsnow=wsnow,deficits=deficits,meltrain=meltrain,meltnon=meltnon,
    melt1=melt1,melt2=melt2,melt3=melt3))

}

melt_function <- function(date, dt, lat, mfmax, mfmin) {

    jday <- as.numeric(format(date,'%j'))
    n_mar21 <- jday - 80
    days <- 365

    # seasonal variation
    sv <- (0.5 * sin((n_mar21 * 2 * pi) / days)) + 0.5
    if (lat < 54) {
        # latitude parameter, av=1.0 when lat < 54 deg N
        av <- 1.0
    } else {
        if (jday <= 77 | jday >= 267) {
            # av = 0.0 from September 24 to March 18,
            av <- 0.0
        } else if (jday >= 117 & jday <= 227) {
            # av = 1.0 from April 27 to August 15
            av <- 1.0
        } else if (jday >= 78 & jday <= 116) {
            # av varies linearly between 0.0 and 1.0 from 3/19-4/26 and
            # between 1.0 and 0.0 from 8/16-9/23.
            ##av <- np.interp(jday, [78, 116], [0, 1])
            days <- 78:116
            ix <- which(days %in% jday) 
            av <- approx(c(0,1),n=length(days))$y[ix]

        } else if (jday >= 228 & jday <= 266) {
            ##av = np.interp(jday, [228, 266], [1, 0])
            days <- 228:266
            ix <- which(days %in% jday) 
            av <- approx(c(0,1),n=length(days))$y[ix]
        }
    }
    meltf <- (dt / 6) * ((sv * av * (mfmax - mfmin)) + mfmin)

    return(meltf)
}

##-------------------------------------------------------
##Snow pillows
##site <- 'tenquille_lake'
##site <- 'upper_squamish'
site <- 'spuzzum_creek'
print(paste0('site: ',site))

pillow.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_pillow/',site,'_asp.csv',sep='')
pillow.data <- read.csv(pillow.file,header=T,as.is=T)
pillow.dates <- as.Date(pillow.data[,2])
pillow.swe <- pillow.data[,11] ##mm

##-----------------------------------------------------------------------
##Snow courses 
site <- 'grouse_mountain'
print(paste0('site: ',site))
course.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_courses/',site,'_snow_course.csv',sep='')
course.data <- read.csv(course.file,header=T,as.is=T)
course.dates <- as.Date(course.data[,1])
course.swe <- course.data[,3] ##mm
course.pack <- course.data[,2] ##cm
course.dense <-  course.data[,4]


##Test default snow17  model
clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_PNWNAmet_800m_data.csv')
clim.data <- read.csv(clim.file,header=T,as.is=T)
dates <- as.Date(clim.data$Dates)
tas <- clim.data$Tas
pr <- clim.data$Pr

###date.subset <- dates %in% pillow.dates
###pillow.subset <- pillow.dates %in% dates[date.subset]
###data <- list(dates=dates,tas=tas,pr=pr,pillow=pillow.swe,psub=pillow.subset,dsub=date.subset)

date.subset <- dates %in% course.dates
course.subset <- course.dates %in% dates[date.subset]

data <- list(dates=dates,tas=tas,pr=pr,snow=course.swe,ssub=course.subset,dsub=date.subset)

min.RSS <- function(par,data) {

   sum((data$snow[data$ssub] -snow17(data,par)[data$dsub])^2,na.rm=T)
}

par.min.RSS <- function(par,pr,snow,tas,dates,snow.sub,data.sub) {

   sum((snow[snow.sub] - snow17(par,pr,tas,dates)[data.sub])^2,na.rm=T)

}

st <- which.min(abs(as.Date('1998-08-01')-dates))
en <- which.min(abs(as.Date('1999-05-30')-dates))


sim.snow <- snow17(dates[st:en], pr[st:en], tas[st:en], 
              lat=49.5, elevation=1000, dt=24, scf=1.0, rvs=1,
              uadj=0.04, mbase=1.0, mfmax=1.05, mfmin=0.6, tipm=0.1, nmf=0.15,
              plwhc=0.04, pxtemp=1.0, pxtemp1=-1.0, pxtemp2=3.0)


##-----------------------------------------------------------------------
##Test optimization of model

##Initial parameter values
           ##SCF, mfmax, mfmin, uadj,  nmf, mbase, pxtemp, plwhc, tipm
par   <- c(  1.0,  1.05,  0.50, 0.04, 0.15,   0.9,    1.0,  0.04,  0.2)
man.par<-c( 0.75,  0.75,  0.10, 0.04, 0.15,   0.9,    1.0,  0.04,  0.02)

upper <- c(  1.4,  2.00,  0.75, 0.19, 0.50,   1.0,    2.0,  0.30,  1.0) 
lower <- c(  0.7,  0.50,  0.05, 0.01, 0.05,   0.0,   -2.0,  0.02,  0.01) 

##-----------------------------------------------------------------------
swe <- sim.snow$swe
frac <- sim.snow$frac
fall <- sim.snow$snowfall

plot(dates[st:en],swe,col='white',
        xlim=c(as.Date('1998-08-01'),as.Date('1999-05-30')),ylim=c(-3000,3500))
##points(pillow.dates[pillow.subset],pillow.swe[pillow.subset],pch=16,cex=1)
points(course.dates[course.subset],course.swe[course.subset],pch=16,cex=1)
lines(dates[st:en],swe,lwd=2,col='green')
lines(dates[st:en],cumsum(pr[st:en]),lwd=2,col='red')
lines(dates[st:en],cumsum(fall),lwd=2,col='orange')
lines(dates[st:en],sim.snow$wice,lwd=2,col='lightblue')
lines(dates[st:en],sim.snow$wsnow,lwd=2,col='blue')
lines(dates[st:en],cumsum(-1*sim.snow$melt),lwd=2,col='black')

