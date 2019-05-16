

# try to write faster version by not looping through months (use stackApply)


# 2015-07-15
# Calculate temperature correction parameters for bias correction
# based on Hempel et al, 2013

# R code written by: Danielle S. Grogan

### ARGUMENTS

## for observational and gcm temperature data:
# use mean if mean daily temperature is available.  Set min and max variables (obs.min.path, .name, .varname) = 0
# Otherwise, use min and max temperature. Set mean variables to 0

# obs.mean.path = character string, path to observation (mean temp) data.  Include "/" as last character
# obs.mean.name = character string, name of file, up to the year part of the name
# obs.mean.varnm = character string, variable name to read from file
# obs.min.path  = character string, path to observation (min temp) data.  Include "/" as last character
# obs.min.name  = character string, name of file, up to the year part of the name
# obs.min.varnm = character string, variable name to read from file
# obs.max.path  = character string, path to observation (max temp) data.  Include "/" as last character
# obs.max.name  = character string, name of file, up to the year part of the name
# obs.max.varnm = character string, variable name to read from file
# 
# same as above for:
# gcm.mean.path  
# gcm.mean.name  
# gcm.mean.varnm 
# gcm.min.path
# gcm.min.name
# gcm.min.varnm
# gcm.max.path
# gcm.max.name
# gcm.max.varnm
# 
# gcm.unit.conv: number, for converting gcm data to units of observation data. Assumes addition/subtraction conversion  
#                For K to C conversion: -272.15.  Set = 0 for no conversion
# 
# network.path = full path (with file name) to river network file.  Used to subset climate data.
#             Set = 0 for no subsetting
# 
# out.path = character string.  Path to output results, include "/" as laster character
# out.name = character string. Name for the output files.  Will be attached to "_B_bar.nc" and "_C_param.nc"
# start.year = beginning year of bias correction period
# end.year   = end year of bias correction period
# temp.dir = character, must end with "/", for writing temporary data.  


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
    data<-addLayer(data, a)
  }
  return(data)
}

#########################################################################

### MAIN ###
#########################################################################
BC_temp_params<-function(obs.mean.path, obs.min.path, obs.max.path,
                         gcm.mean.path, gcm.min.path, gcm.max.path, 
                         obs.mean.name, obs.min.name, obs.max.name,
                         gcm.mean.name, gcm.min.name, gcm.max.name,
                         obs.mean.varnm, obs.min.varnm, obs.max.varnm,
                         gcm.mean.varnm, gcm.min.varnm, gcm.max.varnm,
                         gcm.unit.conv,
                         start.year, end.year,
                         network.path, 
                         out.path, out.name, temp.dir){
  
  # use month.lengths (days) to select monthly data out of yearly data file
  month.lengths<-c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  month.start.days<-c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366)
  n.years = end.year - start.year + 1
  month.names<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  
  if(network.path != 0){
    network<-raster(network.path)
    network.ext<-extent(network)
    # for testing: smaller extent
    #network.ext<-extent(84, 88, 24, 28) # xmin, xmax, ymin, ymax
  }else{
    network.ext<-extent(-999, -999, -999, -999) # indicates no data, no clipping
  }
  
  # make monthly indices (1's for all days in Jan, 2's for all days in Feb)
  for(m in 1:12){
    ind<-rep(m, month.lengths[m])
    if(m==1){
      indices = ind
    }else{
      indices = append(indices, ind)
    }
  }
  ind.m<-rep(indices, n.years)
  
  ind.yr = indices[1:365]
  for(y in 1:(n.years-1)){
    ind.yr.add = ind.yr[((y-1)*365+1):(y*365)] + 12
    ind.yr = append(ind.yr, ind.yr.add)
  }
  
  
  for(y in start.year:end.year){
    ### load gcm data
    tave.gcm<-load.data.ave(gcm.mean.path, gcm.min.path, gcm.max.path, 
                              gcm.mean.name, gcm.min.name, gcm.max.name,
                              gcm.mean.varnm, gcm.min.varnm, gcm.max.varnm,
                              y, network.ext)
    # load obs data
    tave.obs<-load.data.ave(obs.mean.path, obs.min.path, obs.max.path, 
                              obs.mean.name, obs.min.name, obs.max.name,
                              obs.mean.varnm, obs.min.varnm, obs.max.varnm,
                              y, network.ext)
    # remove leap year extra day (quick fix for now)
    if(y %% 4 == 0){tave.obs<-subset(tave.obs, 1:365)}
    
    # convert temp units from K to deg C if needed
    if(gcm.unit.conv != 0){ 
      tave.gcm<-calc(tave.gcm, fun=function(x){return(x+gcm.unit.conv)},
                     filename<-paste(temp.dir, "temp_raster01", sep=""), overwrite=T)
    }
    
    # resample gcm data to spatial resolution of obs data.  Method: bilinear interpolation
    tave.gcm<-resample(tave.gcm, tave.obs, method='bilinear')
    
    # all daily data
    if(y==start.year){tave.obs.all = tave.gcm.all<-c()}
    tave.gcm.all<-add.layers(tave.gcm.all, tave.gcm, (y==start.year))
    tave.obs.all<-add.layers(tave.obs.all, tave.obs, (y==start.year))

    x<-paste("year 1", y)
    print(x, quote=F)
  } # end year loop
  
  # get monthly means (12 layers).  NOTE: this step takes a long time
  gcm.mave<-stackApply(tave.gcm.all, indices=ind.m, fun=mean, na.rm=T,
                       filename<-paste(temp.dir, "temp_raster02",sep=""), overwrite=T)
  obs.mave<-stackApply(tave.obs.all, indices=ind.m, fun=mean, na.rm=T,
                       filename<-paste(temp.dir,"temp_raster03", sep=""), overwrite=T)
  
  ### Correct Monthly Mean ###
  # check if this is equal to equation 1
  C_par = overlay(obs.mave, gcm.mave, fun=function(x,y){return( (x-y))},
                  filename<-paste(temp.dir,"temp_raster06",sep=""), overwrite=T)
  
  # write C_par
  C_par.filenm<-paste(out.path, out.name, "_C_param.nc", sep="")
  writeRaster(C_par, C_par.filenm, format="CDF", varname = "C_param", longname="C temperature correction parameter", 
              varunit = "C",  xname = 'latitude', yname = 'longitude', overwrite=T)
  x<-paste("c_par written")
  print(x, quote=F)
  
  # Make cover file to set NAs to 0 later, so regression will work
  cov<-raster(nrows=nrow(C_par), ncols=ncol(C_par),  xmn=xmin(C_par), xmx=xmax(C_par), 
              ymn=ymin(C_par), ymx=ymax(C_par), crs=crs(C_par) )
  cov<-setValues(cov,0)
  
  # get monthly means (12 layers * n.years), one mean per each month and year
  gcm.mave.all<-stackApply(tave.gcm.all, indices=ind.yr, fun=mean, na.rm=T,
                           filename<-paste(temp.dir,"temp_raster04", sep=""), overwrite=T)
  obs.mave.all<-stackApply(tave.obs.all, indices=ind.yr, fun=mean, na.rm=T,
                           filename<-paste(temp.dir,'temp_raster05',sep=""), overwrite=T)

  # Calculate linear transfer function between residuals
  for(m in 1:12){
    for(y in 1:n.years){
      m.select = ((y-1)*12+1) + (m-1)
      d.select.s =  (y-1)*365 + 1 
      d.select.e = d.select.s + month.lengths[m] - 1
      
      # get month ave for this year and month
      gcm.mave.ym<-subset(gcm.mave.all, m.select)
      obs.mave.ym<-subset(obs.mave.all, m.select)
      
      # daily data for this year and month
      gcm.tave.ym<-subset(tave.gcm.all, d.select.s:d.select.e)
      obs.tave.ym<-subset(tave.obs.all, d.select.s:d.select.e)

      # calculate deltas (eq. 7 Hempel et al 2013)
      Delta.temp.gcm<-overlay(gcm.tave.ym, gcm.mave.ym, fun=function(x,y){return(x-y)})
      Delta.temp.obs<-overlay(obs.tave.ym, obs.mave.ym, fun=function(x,y){return(x-y)})
      
      # stack all of this month's data
         # this month's (n.years) averages
      if(y==1){gcm.mave.all.m = obs.mave.all.m <-c()}
      gcm.mave.all.m<-add.layers(gcm.mave.all.m, gcm.mave.ym, (y==1))
      obs.mave.all.m<-add.layers(obs.mave.all.m, obs.mave.ym, (y==1))
      
        # this month's (month.lengths[m] * n.years) deltas
      if(y==1){Delta.temp.gcm.all = Delta.temp.obs.all <-c()}
      Delta.temp.gcm.all<-add.layers(Delta.temp.gcm.all, Delta.temp.gcm, (y==1))
      Delta.temp.obs.all<-add.layers(Delta.temp.obs.all, Delta.temp.obs, (y==1)) 
    } # end year loop
    
    # Set NA to 0 so sorting will work
    Delta.temp.gcm.all<-cover(Delta.temp.gcm.all, cov)
    Delta.temp.obs.all<-cover(Delta.temp.obs.all, cov)
    
    # this month's (month.lengths[m] * n.years) deltas, all months (*12)
    # used for correcting daily residuals - NOT sorted
    if(m==1){Delta.temp.gcm.all.yrs = Delta.temp.obs.all.yrs <-c()}
    Delta.temp.gcm.all.yrs<-add.layers(Delta.temp.gcm.all.yrs, Delta.temp.gcm.all, (m==1))
    #Delta.temp.obs.all.yrs<-add.layers(Delta.temp.obs.all.yrs, Delta.temp.obs.all, (m==1)) 
    
    # sort
    Delta.temp.gcm.sort<-calc(Delta.temp.gcm.all, fun = function(x){ sort(x) },
                              filename<-paste(temp.dir,'temp_raster07',sep=""), overwrite=T)
    Delta.temp.obs.sort<-calc(Delta.temp.obs.all, fun = function(x){ sort(x) },
                              filename<-paste(temp.dir,'temp_raster08',sep=""), overwrite=T)
    
    # write intermediate output (easier to re-plot or read in later)
    Delta.gcm.month.filenm<-paste(out.path, out.name, "_DeltaT_GCM_", month.names[m], ".nc", sep="")
    writeRaster(Delta.temp.gcm.sort, Delta.gcm.month.filenm, format="CDF", varname = "DeltaT", longname="sorted temperature residuals for month", 
                varunit = "C",  xname = 'latitude', yname = 'longitude', overwrite=T)
    
    Delta.obs.month.filenm<-paste(out.path, out.name, "_DeltaT_OBS_", month.names[m], ".nc", sep="")
    writeRaster(Delta.temp.obs.sort, Delta.obs.month.filenm, format="CDF", varname = "DeltaT", longname="sorted temperature residuals for month", 
                varunit = "C",  xname = 'latitude', yname = 'longitude', overwrite=T)
    
    GCM.mave.month.filenm<-paste(out.path, out.name, "_Tave_m_GCM_", month.names[m], ".nc", sep="")
    writeRaster(gcm.mave.all.m, GCM.mave.month.filenm, format="CDF", varname = "Tave_m", longname="monthly mean temperatures for this month", 
                varunit = "C",  xname = 'latitude', yname = 'longitude', overwrite=T)
    
    OBS.mave.month.filenm<-paste(out.path, out.name, "_Tave_m_OBS_", month.names[m], ".nc", sep="")
    writeRaster(obs.mave.all.m, OBS.mave.month.filenm, format="CDF", varname = "Tave_m", longname="monthly mean temperatures for this month", 
                varunit = "C",  xname = 'latitude', yname = 'longitude', overwrite=T)
    x<-paste("deltas and averages written")
    
    # correct monthly means (used in summary.figs)
    gcm.mave.corr<-overlay(gcm.mave.all.m, subset(C_par, m), fun=function(x,y){return(x+y)})
    
    # write output: corrected monthly means (one file per month)
    gcm.mave.corr.filenm<-paste(out.path, out.name, "_Tave_m_GCMCORR_", month.names[m], ".nc", sep="")
    writeRaster(gcm.mave.corr, gcm.mave.corr.filenm, format="CDF", varname = "Tave_m_corr", longname="monthly mean temperatures for this month, corrected", 
                varunit = "C",  xname = 'latitude', yname = 'longitude', overwrite=T)
    
    ### LINEAR REGRESSION OF RESIDUALS
    # get slope of linear regression (obs as function of gcm)
    s <- stack(Delta.temp.obs.sort, Delta.temp.gcm.sort)
    
    fun=function(x) { if (is.na(x[1])){ NA } else { 
      lm(x[1:(length(x)/2)] ~ x[((length(x)/2)+1):length(x)])$coefficients[2] }}
    B_slope.m <- calc(s, fun, filename<-paste(temp.dir,'temp_raster09',sep=""), overwrite=T)
    
    # r-squared value not used in bias correction; 
    # included here to show where linear regression does well or poorly
    fun=function(x) { if (is.na(x[1])){ NA } else { 
      summary(lm(x[1:(length(x)/2)] ~ x[((length(x)/2)+1):length(x)]))$r.squared }}
    r2.m <- calc(s, fun, filename<-paste(temp.dir,'temp_raster10',sep=""), overwrite=T)
    
    # stack monthly data for output
    if(m==1){B_slope = r2 <-c()}
    B_slope<-add.layers(B_slope, B_slope.m, (m==1))
    r2<-add.layers(r2, r2.m, (m==1))
    
    x<-paste("month", m)
    print(x, quote=F)
  } # end month loop
  
  # write r2 and B_slope (not used, but output for testing)
  r2.filenm<-paste(out.path, out.name, "_r2.nc", sep="")
  writeRaster(r2, r2.filenm, format="CDF", varname = "r2", longname="r2 of temperature transfer function", 
              varunit = "unitless",  xname = 'latitude', yname = 'longitude', overwrite=T)
  
  B_slope.filenm<-paste(out.path, out.name, "_B_slope.nc", sep="")
  writeRaster(B_slope, B_slope.filenm, format="CDF", varname = "B_slope", longname="slope of temperature transfer function", 
              varunit = "unitless",  xname = 'latitude', yname = 'longitude', overwrite=T)
  
  # Daily correction factors
  # loop through months again, because need values of future months
  for(m in 1:12){                    
    for(day in 1:month.lengths[m]){
      
      # weighing factors to smooth monthly borders
      d  = ((day - 1)/(month.lengths[m] - 1)) - 0.5
      dm = 0.5*(abs(d) - d)                            # previous month weight
      d0 = 1-abs(d)                                    # current month weight
      dp = 0.5*(abs(d) + d)                            # next month weight
      
      if( m == 1 ){ 
        pre.month  = 12 
        post.month = 2
      }else if( m == 11 ){ 
        pre.month  = 10
        post.month = 12 
      }else if( m == 12){
        pre.month = 11
        post.month = 1
      }else{
        pre.month = (m-1)%%12
        post.month = (m+1)%%12
      }
      B_bar.d = dm*subset(B_slope, pre.month ) + d0*subset(B_slope, m) + dp*subset(B_slope, post.month )
      
      # stack daily data for output
      if(m==1 & day==1){B_bar <-c()}
      B_bar<-add.layers(B_bar, B_bar.d, (m==1 & day==1))
    } # end day loop
    
    # correct daily residuals (not used in BC, needed for summary.figs)
    # need to use daily data and monthly ave IN ORDER of days (don't used sorted deltas)
    if(m==1){
      start=1
      end = start + n.years * month.lengths[m]-1
    }else{
      start = end+1
      end   = start + n.years * month.lengths[m]-1
    }
   
    Delta.temp.gcm.m<-subset(Delta.temp.gcm.all.yrs, start:end)  # select all data for this month
    B_bar.month<-subset(B_bar, (month.start.days[m]:(month.start.days[m+1]-1))) # selet B_bar for this month
   
    delta.gcm.corr<-overlay(Delta.temp.gcm.m, B_bar.month, fun=function(x,y){return(x*y)},
                            filename<-paste(temp.dir,'temp_raster11',sep=""), overwrite=T)

    # write output: corrected daily residuals 
    delta.gcm.corr.filenm<-paste(out.path, out.name, "_DeltaT_GCMCORR_", month.names[m], ".nc", sep="")
    writeRaster(delta.gcm.corr, delta.gcm.corr.filenm, format="CDF", varname = "DeltaT_GCMCORR", 
                longname="corrected temperature residuals for this month", 
                varunit = "C",  xname = 'latitude', yname = 'longitude', overwrite=T)
    
    x<-paste("daily correction factors,", m)
    print(x, quote=F)
  } # end month loop
  
  # write output: B_bar
  B_bar.filenm<-paste(out.path, out.name, "_B_bar.nc", sep="")
  writeRaster(B_bar, B_bar.filenm, format="CDF", varname = "B_bar", longname="B_bar temperature correction parameter", 
              varunit = "B_bar",  xname = 'latitude', yname = 'longitude', overwrite=T)

}  #################################################### # END MAIN
