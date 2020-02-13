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
get_coordinates <- function(site) {

  coordinates <- list(callaghan=c(-123.1036,50.1383278,1009),
                      orchid_lake=c(-123.0519638,49.53678,1178),
                      palisade_lake=c(-123.0321944,49.454433,898),
                      grouse_mountain=c(-123.0774472,49.383655,1126),
                      dog_mountain=c(-122.96255,49.37251944,1007),
                      dickson_lake=c(-122.06984166,49.3168194,1147),
                      stave_lake=c(-122.315805,49.58030277,1211),
                      nahatlatch=c(-122.059261,49.825866,1530),
                      wahleach=c(-121.57945,49.2298694,1395),
                      klesilkwa=c(-121.3086527,49.129438,610),
                      lightning_lake=c(-120.850205,49.044788,1254),
                      brookmere=c(-120.87397,49.815027,994),
                      shovelnose_mountain=c(-120.864175,49.8546305,1456),
                      hamilton_hill=c(-120.7955805,49.4988027,1477),
                      spuzzum_creek=c(-121.686,49.674,1197),
                      chilliwack_river=c(-121.71667,49.0333,1600),
                      upper_squamish=c(-123.4333,50.1500,1340),
                      wahleach_lake=c(-121.5833,49.2333,1400),
                      tenquille_lake=c(-122.9333,50.5333,1680))

  rv <- coordinates[[site]]
  return(rv)
}

##---------------------------------------------------------------------

min.RSS <- function(par,data,loc) {

   sum((data$snow[data$ssub] -snow17(data,par,loc)[data$dsub])^2,na.rm=T)
}

snow17 <- function(data,par,loc) {

    ##Defaults held constant:
    lat <- loc[2]
    elevation <- loc[3]
    dt <- 24
    rvs <- 1
    pxtemp1 <- -1.0
    pxtemp2 <- 3.0

    ##Optimized parameters

    scf <- par[1]
    mfmax <- par[2]
    mfmin <- 0.5 ##par[3]
    uadj <- 0.04 #par[4]
    nmf <- 0.15 #par[5]
    mbase <- 0.9 #par[6]
    pxtemp <- 2.0 #par[7]
    plwhc <- 0.04 #par[8]
    tipm <- par[3] ##par[9]

    pr <- data$pr
    tas <- data$tas
    dates <- data$dates

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

        fracrain <- 1.0 - fracsnow

        # Snow Accumulation
        # water equivalent of new snowfall (mm)
        pn <- precip * fracsnow * scf
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
            m_ros1 <- max(c(stefan * dt * (((t_air_mean + 273) ^ 4) - (273 ^ 4)), 0.0))
            m_ros2 <- max(c((0.0125 * rain * t_rain), 0.0))
            m_ros3 <- max(c((8.5 * uadj * (dt / 6.0) * (((0.9 * e_sat) - 6.11) +
                           (0.00057 * p_atm * t_air_mean))),0.0))
            m_ros <- m_ros1 + m_ros2 + m_ros3
        } else {
            m_ros <- 0.0
        }

        # Non-Rain melt
        if (rain <= (0.25 * dt) & (t_air_mean > mbase)) {
            # melt during non-rain periods is:
            m_nr <- (mf * (t_air_mean - mbase)) + (0.0125 * rain * t_rain)
        } else {
            m_nr <- 0.0
        }

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

        # End of model execution
        model_swe[i] <- swe  # total swe (mm) at this time step
        outflow[i] <- e
    }    
    ##return(list(swe=model_swe,outflow=outflow))
    return(model_swe)
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

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

loc <- get_coordinates(site)

##Test default snow17  model
clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_PNWNAmet_800m_data.csv')
clim.data <- read.csv(clim.file,header=T,as.is=T)
dates <- as.Date(clim.data$Dates)
tas <- clim.data$Tas
pr <- clim.data$Pr

##Snow pillows
pillows <- c('tenquille_lake','upper_squamish','spuzzum_creek','chilliwack_river')

if (site %in% pillows) {
   pillow.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_pillow/',site,'_asp.csv',sep='')
   pillow.data <- read.csv(pillow.file,header=T,as.is=T)
   site.dates <- as.Date(pillow.data[,2])
   site.swe <- pillow.data[,11] ##mm

} else {
   ##Snow courses 
   course.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_courses/',site,'_snow_course.csv',sep='')
   course.data <- read.csv(course.file,header=T,as.is=T)
   site.dates <- as.Date(course.data[,1])
   site.swe <- course.data[,3] ##mm
}

date.subset <- dates %in% site.dates
site.subset <- site.dates %in% dates[date.subset]

data <- list(dates=dates,tas=tas,pr=pr,snow=site.swe,ssub=site.subset,dsub=date.subset)

##-----------------------------------------------------------------------
##Initial parameter values
           ##SCF, MFMax,  TIPM
sub.par <- c(1, 0.75,  0.05)
upper.par <- c(2,2,0.2)
lower.par <- c(0.6,0.2,0.01)

optim.result <- optim(par=sub.par,fn=min.RSS,data=data,loc=loc,
                      upper=upper.par,lower=lower.par,method="L-BFGS-B")
save.dir <- '/storage/data/projects/rci/data/winter_sports/obs/snow17/'
save(optim.result,file=paste0(save.dir,site,'_snow17_optim_3_parameter_fit_bounded_pnwnamet.RData'))

##nse <- 1 - (sum( (pillow.swe[pillow.subset] - sub.swe[date.subset])^2,na.rm=T) / 
##           sum( (pillow.swe[pillow.subset] - mean(pillow.swe[pillow.subset],na.rm=T))^2,na.rm=T))
