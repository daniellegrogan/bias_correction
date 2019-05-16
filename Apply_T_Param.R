# 2015-07-03
# Apply temperature correction parameters for bias correction
# based on Hempel et al, 2013

# R code written by: Danielle S. Grogan

### ARGUMENTS
## for gcm temperature data:
# use mean if mean daily temperature is available.  Set min and max variables (obs.min.path, .name, .varname) = 0
# Otherwise, use min and max temperature. Set mean variables to 0

# tmin.path = character string, path to observation (mean temp) data.  Include "/" as last character
# tmin.name = character string, name of file, up to the year part of the name
# tmin.varnm = character string, variable name to read from file
# same as above for tmax.path, tmax.name, tmax.varname
# 
# gcm.unit.conv: number, for converting gcm data to units of observation (and correction factor) data. 
#                Assumes addition/subtraction conversion  
#                For K to C conversion: -272.15.  Set = 0 for no conversion
# 
# out.path = character string.  Path to output results, include "/" as laster character
# out.name = character string. Name for the output files.  Will be attached to "_year.nc"

# start.year = beginning year of application period
# end.year   = end year of application period

library(raster)
library(rgdal)
library(rasterVis)
#########################################################################
### SUBFUNCTIONS ###
load.data.ave<-function(mean.path,  min.path,  max.path, 
                        mean.name,  min.name,  max.name,
                        mean.varnm, min.varnm, max.varnm,
                        y, network.ext){
  if(mean.path == 0){ # if no  mean temp, load min and max, get average
    # daily max temperature
    tmax.nm<-paste(max.path, max.name, y, ".nc", sep="")
    tmax.data<-brick(tmax.nm, varname= max.varnm, stopIfNotEqualSpaced = FALSE)
    
    # daily min temperature
    tmin.nm<-paste(min.path, min.name, y, ".nc", sep="")   
    tmin.data<-brick(tmin.nm, varname= min.varnm, stopIfNotEqualSpaced = FALSE)
    
    if(network.ext@xmin != -999){
      tmax.data<-crop(tmax.data, network.ext) 
      tmin.data<-crop(tmin.data, network.ext)
    }
    tave.data<-overlay(tmax.data, tmin.data, fun=function(x,y){return ((x+y)/2)})
    
  }else{
    # daily ave temperature
    tave.nm<-paste(mean.path, mean.name, y, ".nc", sep="")
    tave.data<-brick(tave.nm, varname= mean.varnm, stopIfNotEqualSpaced = FALSE)
    
    if(network.ext@xmin != -999){
      tave.data<-crop(tave.data, network.ext) 
    }
  }
  return(tave.data)
}
#########################################################################
add.layers<-function(data, add.data, criteria){
  if(criteria == T){
    if(nlayers(add.data) == 1){
      data<-brick(add.data)
    }else{
      data = add.data
    }
  }else{
    a = add.data * 1
    data<-addLayers(data, a)
  }
  return(data)
}

#########################################################################
# MAIN
#########################################################################
apply.BCTparams<-function(tmax.path, tmin.path, tmean.path,
                          tmax.name, tmin.name, tmean.name,
                          tmax.varnm, tmin.varnm, tmean.varnm,
                          B_bar.path, B_bar.varnm, C_par.path, C_par.varnm, 
                          gcm.unit.conv, 
                          start.yr, end.yr,
                          network.path,
                          out.dir, out.name, temp.dir){

  # 1. load bias correction parameter files 
  B_bar<-brick(B_bar.path, B_bar.varnm)
  C_par<-brick(C_par.path, C_par.varnm)
  
  # data for monthly mean calculations
  month.lengths<-c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  month.start.days<-c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366)
  
  if(network.path != 0){
    network<-raster(network.path)
    network.ext<-extent(network)
    # for testing: smaller extent
    #network.ext<-extent(84, 88, 24, 28)  # xmin, xmax, ymin, ymax
  }else{
    network.ext<-extent(-999, -999, -999, -999) # indicates no data, no clipping
  }
  
  for(m in 1:12){
    ind<-rep(m, month.lengths[m])
    if(m==1){
      indices = ind
    }else{
      indices = append(indices, ind)
    }
  }
  
  # 2. loop through each year of the application period
  for(y in start.yr:end.yr){
    
    
    # Read climate model (GCM) data
    ### load gcm data
    gcm.tave.25<-load.data.ave(tmean.path, tmin.path, tmax.path, 
                              tmean.name, tmin.name, tmax.name,
                              tmean.varnm, tmin.varnm, tmax.varnm,
                              y, network.ext)
   
    # convert temp units from K to deg C if needed
    if(gcm.unit.conv != 0){ 
      gcm.tave.25<-calc(gcm.tave.25, fun=function(x){return(x+gcm.unit.conv)}, 
                        filename<-paste(temp.dir, "temp_raster01", sep=""), overwrite=TRUE)
    }
    
    # resample gcm data to spatial resolution of BC params.  Method: bilinear interpolation
    gcm.tave.25<-resample(gcm.tave.25, B_bar, method='bilinear')
    
    # Calculate monthly means
    gcm.tave.MONTH<-stackApply(gcm.tave.25, indices=indices, fun=mean, na.rm=T,
                               filename<-paste(temp.dir, "temp_raster02", sep=""), overwrite=T)
    
    # test
    gcm.tave.25.1<-extract(gcm.tave.25, 5000)
    gcm.tave.MONTH.1<-extract(gcm.tave.MONTH, 5000)
    
    # Calculate daily CORRECTED temperature values (Hempel et al 2013, Eq.22)
    for(m in 1:12){
      
      # select values for this month
      C_par.m          <- subset(C_par,m)
      gcm.tave.MONTH.m <- subset(gcm.tave.MONTH, m)
      B_bar.m          <- subset(B_bar, month.start.days[m]:(month.start.days[m+1]-1)) 
      gcm.tave.m       <- subset(gcm.tave.25, month.start.days[m]:(month.start.days[m+1]-1))
      
#       # test
#       C_par.m.1<-extract(C_par.m, 5000)
#       B_bar.m.1<-extract(B_bar.m, 5000)
#       gcm.tave.m.1<-extract(gcm.tave.m, 5000)
#       gcm.tave.MONTH.m.1<-extract(gcm.tave.MONTH.m, 5000)
#       
#       test<-c()
#       for(d in 1:31){
#         test[d] = C_par.m.1 + gcm.tave.MONTH.m.1 + B_bar.m.1[d]*(gcm.tave.m.1[d] - gcm.tave.MONTH.m.1)
#       }
      
      # Eq. 22
      Tcorrect_day.m = C_par.m + gcm.tave.MONTH.m + B_bar.m*(gcm.tave.m - gcm.tave.MONTH.m)
      
      if(m==1){
        Tcorrect_day = Tcorrect_day.m
      }else{
        Tc = Tcorrect_day.m*1
        Tcorrect_day <- addLayer(Tcorrect_day, Tc)
      }
    }

     # write output
    Tcorrect.filenm<-paste(out.dir, out.name, y, ".nc", sep="")
    writeRaster(Tcorrect_day, Tcorrect.filenm, format="CDF", varname = "tave", longname="ave daily temperature, bias-corrected", 
                varunit = "C",  xname = 'latitude', yname = 'longitude', overwrite=T)
    print(y)
  }      # close year loop 
}

######################################################################### # END OF MAIN
