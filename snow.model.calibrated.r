source('/storage/data/projects/rci/assessments/code/eco_hydrology.R')
library(zoo)
##-------------------------------------------------------------
##Snow/Rain Phase function

hyper_snow_phase <- function(tas,precip_mm,coeffs) {
   ##D is fixed at 1.0209
   frac <- coeffs$a*(tanh(coeffs$b*(tas-coeffs$c))-1.0209)
   sample <- runif(length(tas),min=0,max=100)
   test <- sample > frac
   snow.type <- rep(TRUE,length(tas))
   snow.type[test] <- FALSE

   NewSnowWatEq <- precip_mm
   NewSnowWatEq[!snow.type] <- 0
   R_m <- precip_mm
   R_m[snow.type] <- 0

   rv <- list(swe=NewSnowWatEq,rain=R_m)
   return(rv)
}

##-------------------------------------------------------------

snow_melt <- function(precip_mm, Tmax_C, Tmin_C, Date, lat_deg, 
                      cal_scale, cal_slope, cal_freq, 
                      slope=0, aspect=0, 
                      tempHt=1, windHt=2, groundAlbedo=0.25,SurfEmissiv=0.95, windSp=1, 
                      forest=0, startingSnowDepth_m=0, startingSnowDensity_kg_m3=450) {

  ptm <- proc.time()
  coeffs <- list(a=cal_scale,b=cal_slope,c=cal_freq)
  Tmax.check <- Tmax_C <  Tmin_C
  Tmax_C[Tmax.check] <- Tmin_C[Tmax.check]+1
 
  ##-------------------------------------------------------------      
  ## Constants :
  WaterDens <- 1000			# kg/m3
  SnowDens <- 600                         # kg/m3
  lambda <- 3.35*10^5			# latent heat of fusion (kJ/m3)
  lambdaV <- 2500				# (kJ/kg) latent heat of vaporization
  SnowHeatCap <- 2.1			# kJ/kg/C
  LatHeatFreez <- 333.3		        # kJ/kg
  Cw <- 4.2*10^3				# Heat Capacity of Water (kJ/m3/C)
	
  ##	Converted Inputs :
  Tav <- (Tmax_C+Tmin_C)/2		# degrees C
  precip_m <- precip_mm*0.001	 	# precip in m 
  NewSnowDensity <- 50+3.4*(Tav+15)		# kg/m3
  NewSnowDensity[which(NewSnowDensity < 50)] <- 50

  ##--------------------
  ##Separate Precipitation by phase
  slen <- 100
  snow.series <- matrix(0,nrow=length(precip_mm),ncol=slen)
  rain.series <- matrix(0,nrow=length(precip_mm),ncol=slen)
  for (k in 1:slen) {
     snow <- hyper_snow_phase(Tav,precip_mm,coeffs)
     snow.series[,k] <- snow$swe
     rain.series[,k] <- snow$rain
  }
  NewSnowWatEq <- apply(snow.series,1,mean)/1000
  R_m <- apply(snow.series,1,mean)/1000
  R_m[NewSnowWatEq !=0] <- 0

  NewSnow <- NewSnowWatEq*WaterDens/NewSnowDensity		# m

  ##-----------------------------

  JDay <- strptime(Date, format="%Y-%m-%d")$yday+1
  lat <- lat_deg*pi/180		#	latitude in radians

  ##rh 	<- log((windHt+0.001)/0.001)*log((tempHt+0.0002)/0.0002)/(0.41*0.41*windSp*86400)	# (day/m) Thermal Resistance	 
  rh      <- log((10+0.001)/0.001)*log((10+0.0002)/0.0002)/(0.41*0.41*windSp*86400)       # (day/m) Thermal Resistance

  if (length(windSp)==1) rh <- rep(rh,length(Tmax_C))		##	creates a vector of rh values

  cloudiness <- EstCloudiness(Tmax_C,Tmin_C)
  AE <- AtmosphericEmissivity(Tav, cloudiness)	# (-) Atmospheric Emissivity

  #  New Variables	:
  SnowTemp 		<- rep(0,length(Tmax_C)) 		# Degrees C
  rhos 			<- SatVaporDensity(SnowTemp)	# 	vapor density at surface (kg/m3)
  rhoa 			<- SatVaporDensity(Tmin_C)		#	vapor density of atmoshpere (kg/m3) 
  SnowWaterEq 	<- vector(length=length(Tmax_C))		#  (m) Equiv depth of water
  TE 				<- rep(SurfEmissiv,length(Tmax_C))	#	(-) Terrestrial Emissivity
  DCoef 			<- rep(0,length(Tmax_C))				#   Density Coefficient (-) (Simplified version)
  SnowDensity 	<- rep(SnowDens,length(Tmax_C))			#  (kg/m3)  Max density is 500
  SnowDepth 		<- vector(length=length(Tmax_C))		#  (m)
  SnowMelt 		<- rep(0,length(Tmax_C))				#  (m)
  melt            <- rep(0,length(Tmax_C))                             #  (m)
  DensityPerc             <- rep(0,length(Tmax_C))
  Albedo 			<- rep(groundAlbedo,length(Tmax_C)) 	#  (-) This will change for days with snow
	
  ##	Energy Terms
  H 		<- vector(length=length(Tmax_C))	#	Sensible Heat exchanged (kJ/m2/d)
  E 		<- vector(length=length(Tmax_C))	#	Vapor Energy	(kJ/m2/d)
  S 		<- vector(length=length(Tmax_C))	#	Solar Radiation (kJ/m2/d)
  La 		<- Longwave(AE, Tav)					#	Atmospheric Longwave Radiation (kJ/m2/d)
  Lt 		<- vector(length=length(Tmax_C))	#	Terrestrial Longwave Radiation (kJ/m2/d)
  G 		<- 173								#	Ground Condution (kJ/m2/d) 
  P 		<- Cw * R_m * Tav					# 	Precipitation Heat (kJ/m2/d)
  Energy 	<- vector(length=length(Tmax_C))	# Net Energy (kJ/m2/d)

  ##  Initial Values.  
  SnowWaterEq[1] 	<- startingSnowDepth_m * startingSnowDensity_kg_m3 / WaterDens		
  SnowDepth[1] 	<- startingSnowDepth_m			
  
# If snow on the ground or new snow, assume Albedo yesterday was 0.5
  Albedo[1] <- ifelse(NewSnow[1] > 0, 0.98-(0.98-0.50)*exp(-4*NewSnow[1]*10),
               ifelse(startingSnowDepth_m == 0, groundAlbedo, max(groundAlbedo, 0.5+(groundAlbedo-0.85)/10)))  

  S[1] <- Solar(lat=lat,Jday=JDay[1], 
                Tx=Tmax_C[1], Tn=Tmin_C[1], 
                albedo=Albedo[1], forest=forest, 
                aspect=aspect, slope=slope, 
                latUnits='degrees')

  H[1] <- 1.29*(Tav[1]-SnowTemp[1])/rh[1] 
  E[1] <- lambdaV*(rhoa[1]-rhos[1])/rh[1]
  if(startingSnowDepth_m>0) TE[1] <- 0.97 
  Lt[1] <- Longwave(TE[1],SnowTemp[1])
  Energy[1] <- S[1] + La[1] - Lt[1] + H[1] + E[1] + G + P[1]

  ##Snow Density
  SnowDensity[1] <- ifelse((startingSnowDepth_m+NewSnow[1])>0, 
                          min(SnowDens, 
                          (startingSnowDensity_kg_m3*startingSnowDepth_m + NewSnowDensity[1]*NewSnow[1]) / 
                          (startingSnowDepth_m+NewSnow[1])), SnowDens)

  ## yesterday on ground + today new   
  SnowMelt[1] <- max(0,	min((startingSnowDepth_m/10+NewSnowWatEq[1]),
                (Energy[1]-SnowHeatCap*(startingSnowDepth_m/10+NewSnowWatEq[1])*
                WaterDens*(0-SnowTemp[1]))/(LatHeatFreez*WaterDens) ) )
  melt[1] <- 0
  SnowDepth[1] <- max(0,(startingSnowDepth_m/10 + NewSnowWatEq[1]-SnowMelt[1])*WaterDens/SnowDensity[1])
  SnowWaterEq[1] <- max(0,startingSnowDepth_m/10-SnowMelt[1]+NewSnowWatEq[1])	
  DensityPerc[1] <- 10
  
  ##print('Initial time')
  ##print(proc.time()-ptm)
  ##Rprof('fx.snow.out')
  ltm <- proc.time()	

##  Snow Melt Loop	
  for (i in 2:length(Tmax_C)){

      if (NewSnow[i] > 0){ 
          Albedo[i] <- 0.98-(0.98-Albedo[i-1])*exp(-4*NewSnow[i]*10)
      } else if (SnowDepth[i-1] < 0.1){ 
          Albedo[i] <- max(groundAlbedo, Albedo[i-1]+(groundAlbedo-0.85)/10)
      } else Albedo[i] <- 0.35-(0.35-0.98)*exp(-1*(0.177+(log((-0.3+0.98)/(Albedo[i-1]-0.3)))^2.16)^0.46)

       S[i] <- Solar(lat=lat,Jday=JDay[i], Tx=Tmax_C[i], Tn=Tmin_C[i], albedo=Albedo[i-1], 
                     forest=forest, aspect=aspect, slope=slope, printWarn=FALSE)

      if(SnowDepth[i-1] > 0) TE[i] <- 0.97 	#	(-) Terrestrial Emissivity
      if(SnowWaterEq[i-1] > 0 | NewSnowWatEq[i] > 0) {
          DCoef[i] <- 9.45  ##DCoef.val ##9.5
          if(SnowMelt[i-1] == 0){ 
	     SnowTemp[i] <- max(min(0,Tmin_C[i]),
                                min(0,(SnowTemp[i-1]+
                                min(-SnowTemp[i-1],Energy[i-1]/((SnowDensity[i-1]*SnowDepth[i-1]+NewSnow[i]*NewSnowDensity[i])*SnowHeatCap*1000)))))
	  }
      }

      rhos[i] <- SatVaporDensity(SnowTemp[i])
      H[i] <- 1.29*(Tav[i]-SnowTemp[i])/rh[i] 
      E[i] <- lambdaV*(rhoa[i]-rhos[i])/rh[i]
      Lt[i] <- Longwave(TE[i],SnowTemp[i])
      Energy[i] <- S[i] + La[i] - Lt[i] + H[i] + E[i] + G + P[i]
      ##print(i)
      ##print(Energy[i])

      if (is.na(Energy[i])) {
          print(i)
          browser()
      }
      if (Energy[i]>0) k <- 1 else k <- 0.5

      SnowDensity[i] <- ifelse((SnowDepth[i-1]+NewSnow[i])>0, min(SnowDens, 
                              ((SnowDensity[i-1]+k*SnowDens*(SnowDens-SnowDensity[i-1])*exp(-DCoef[i]))*SnowDepth[i-1] + 
                              NewSnowDensity[i]*NewSnow[i])/(SnowDepth[i-1]+NewSnow[i])), SnowDens)

      # yesterday on ground + today new
      SnowMelt[i] <- max(0,	min( (SnowWaterEq[i-1]+NewSnowWatEq[i]),
                        (Energy[i]-SnowHeatCap*(SnowWaterEq[i-1]+NewSnowWatEq[i])*
                        WaterDens*(0-SnowTemp[i]))/(LatHeatFreez*WaterDens) )  )

      melt[i] <- SnowHeatCap*(SnowWaterEq[i-1]+NewSnowWatEq[i])*WaterDens*(0-SnowTemp[i])

      SnowDepth[i] <- max(0,(SnowWaterEq[i-1]+NewSnowWatEq[i]-SnowMelt[i])*WaterDens/SnowDensity[i])
      SnowWaterEq[i] <- max(0,SnowWaterEq[i-1]-SnowMelt[i]+NewSnowWatEq[i])	# (m) Equiv depth of water
      DensityPerc[i] <- SnowWaterEq[i]/SnowDepth[i]*100                

  }
  ##print('Loop time')
  ##print(proc.time()-ltm)

  Results <- list(snowdepth=SnowDepth,swe=SnowWaterEq,snowfall=NewSnow,snowdense=DensityPerc,swe=SnowWaterEq)
  ##Rprof(NULL)

  return(Results)
}
