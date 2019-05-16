# 2015-07-03
# Calculate temperature correction parameters for bias correction
# based on Hempel et al, 2013

# R code written by: Danielle S. Grogan

### ARGUMENTS
# obs.path = character string, path to observation (mean temp) data.  Include "/" as last character
# obs.name = character string, name of file, up to the year part of the name.  Include "_" as last character
# obs.varnm = character string, variable name to read from file
# 
# same as above for gcm data:
# gcm.path
# gcm.name
# gcm.varnm
# 
# gcm.unit.conv =  number, conversion factor for gcm data (to make same units as obs data) 
#                  assumes * conversion (not +)
# start.year = beginning year of bias correction period
# end.year   = end year of bias correction period 
# network.path = full path (with file name) to river network file.  Used to subset climate data.
#                Set = 0 for no subsetting
#
# out.path = character string.  Path to output results, include "/" as laster character
# out.name = character string. Name for the output files.  
#           Will be attached to "_P_c_corr.nc", "_P_A_int.nc" and "_P_B_slope.nc"
# summary.figs = 1 or 0.  Set to 1 to output summary figures of correction factors.

library(raster)
library(rasterVis)
library(rgdal)

### SUBFUNCTIONS ###
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
BC_precip_params<-function(obs.path, obs.name, obs.varnm,
                           gcm.path, gcm.name, gcm.varnm, 
                           gcm.unit.conv,
                           start.year, end.year,
                           network.path, 
                           out.path, out.name,
                           summary.figs=1, temp.dir){
  
  month.lengths<-c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  month.start.days<-c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366)
  month.names<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  
  n.years = end.year - start.year + 1
  
  if(network.path != 0){
    network<-raster(network.path)
    network.ext<-extent(network)
    # for testing: smaller extent
    #network.ext<-extent(84, 88, 24, 28) # xmin, xmax, ymin, ymax
  }else{
    network.ext<-extent(-999, -999, -999, -999) # indicates no data, no clipping
  }
  
  for(m in 1:12){
    for(y in start.year:end.year){
      
      # read data 
      precip.gcm.nm<-paste(gcm.path, gcm.name, y, ".nc", sep="") 
      precip.gcm<-brick(precip.gcm.nm, varname=gcm.varnm, stopIfNotEqualSpaced = FALSE)
      
      precip.obs.nm<-paste(obs.path, obs.name, y, ".nc", sep="")
      precip.obs<-brick(precip.obs.nm, varname = obs.varnm)
      
      # select days within month
      start.day = (sum(month.lengths[1:m]) - month.lengths[m] + 1)
      end.day   = (start.day + month.lengths[m] - 1)
      
      precip.gcm.m<-subset(precip.gcm,start.day:end.day)   # select this month
      precip.obs.m<-subset(precip.obs,start.day:end.day)  # select this month
      
      # cut spatial extent of obs data to network extent.  Saves computing time for sub-global
      if(network.path != 0){
        precip.obs.m<-crop(precip.obs.m, network.ext)         
      }
      
      # make NA mask for later
      if(y == start.year & m==1){ # only need to make mask once
        data.mask<-subset(precip.obs.m, 1:12)
      }
      
      # calculate month ave
      precip.obs.mave<-calc(precip.obs.m, mean, na.rm=T,
                            filename<-paste(temp.dir, "temp_raster01", sep=""), overwrite=T)
      
      # save month's data to the all-year (this month) array
      # all daily data
      if(y==start.year){precip.obs.mlist = precip.obs.dlist<-c()}
      precip.obs.mlist<-add.layers(precip.obs.mlist, precip.obs.mave, (y==start.year)) # monthly ave
      precip.obs.dlist<-add.layers(precip.obs.dlist, precip.obs.m,    (y==start.year)) # daily data
      
      # set spatial extent and resolution of gcm data
      precip.gcm.025<-resample(precip.gcm.m, precip.obs.m)
      if(gcm.unit.conv!= 0){precip.gcm.025 = precip.gcm.025*gcm.unit.conv}
      
      # cut spatial extent of obs data to network extent.  Saves computing time for sub-global
      if(network.path != 0){
        precip.gcm.025<-crop(precip.gcm.025, network.ext) 
      }
      
      # calculate month mean
      precip.gcm.mave<-calc(precip.gcm.025, mean, na.rm=T,
                            filename<-paste(temp.dir, "temp_raster03", sep=""), overwrite=T)
      
      # save month's data to the all-year (this month) array
      if(y==start.year){precip.gcm.mlist = precip.gcm.dlist<-c()}
      precip.gcm.mlist<-add.layers(precip.gcm.mlist, precip.gcm.mave, (y==start.year)) # monthly ave
      precip.gcm.dlist<-add.layers(precip.gcm.dlist, precip.gcm.025,    (y==start.year)) # daily data
      
      # remove variables no longer needed to save memory space
      #remove('precip.obs.mave', 'precip.obs.m', 'precip.gcm.mave', 'precip.gcm.025', 'precip.gcm.m', 'start.day', 'end.day')
      x<-paste("year step 1", y)
      print(x, quote=F)
    } # close year loop
    
    # write output: monthly mean precip (used in summary.figs); one file per month
    precip.obs.mlist.filenm<-paste(out.path, out.name, "_Pmean_OBS_", month.names[m], ".nc", sep="")
    writeRaster(precip.obs.mlist, precip.obs.mlist.filenm, format="CDF", varname = "Pmean", longname="mean monthly precip", 
                varunit = "mm/day",  xname = 'latitude', yname = 'longitude', overwrite=T)
    
    precip.gcm.mlist.filenm<-paste(out.path, out.name, "_Pmean_GCM_", month.names[m], ".nc", sep="")
    writeRaster(precip.gcm.mlist, precip.gcm.mlist.filenm, format="CDF", varname = "Pmean", longname="mean monthly precip", 
                varunit = "mm/day",  xname = 'latitude', yname = 'longitude', overwrite=T)
    
    ### Correct Monthly Mean (c parameter, Eq.4 in Hempel et al 2013) ###
    precip.obs.sum.months<-calc(precip.obs.mlist, sum, na.rm=T,
                                filename<-paste(temp.dir, "temp_raster04", sep=""), overwrite=T)
    precip.gcm.sum.months<-calc(precip.gcm.mlist, sum, na.rm=T,
                                filename<-paste(temp.dir, "temp_raster05", sep=""), overwrite=T)
    
    P_c_corr.m = overlay(precip.obs.sum.months, precip.gcm.sum.months, fun=function(x,y){return(x/y)},
                         filename<-paste(temp.dir, "temp_raster06", sep=""), overwrite=T)
    P_c_corr.m[P_c_corr.m > 10] = 10 # max(P_c_corr) := 10
    
    if(m==1){P_c_corr<-c()}
    P_c_corr<-add.layers(P_c_corr, P_c_corr.m, (m==1)) # monthly ave
    
    # correct monthly means (used in summary figs)
    precip.mean.corr<-overlay(precip.gcm.mlist, P_c_corr.m, fun=function(x,y){return(x*y)})
    
    # write output: corrected monthly mean precip (used in summary.figs); one file per month
    precip.mean.corr.filenm<-paste(out.path, out.name, "_Pmean_GCMCORR_", month.names[m], ".nc", sep="")
    writeRaster(precip.mean.corr, precip.mean.corr.filenm, format="CDF", varname = "Pmean", longname="mean monthly precip, corrected", 
                varunit = "mm/day",  xname = 'latitude', yname = 'longitude', overwrite=T)
    
    
    ### Calculate epsilon_m (dry month threshold) for this month
    
    # an obs month is := dry if ave monthly precip <= 0.01 mm/day
    # a  gcm month is := dry if ave monthly precip = 0
    # IF number obs dry months > number gcm 0-precip months, 
    #    THEN define driest gcm months as "dry" in order of increasing precip (starting with least precip)
    #    AND define epsilon_m as the ave monthly precip in the last "dry" gcm month
    # IF number gcm dry months > obs dry months,
    #    THEN exclude only the gcm months with precip = 0, epsilon_m = 0
    
    dry.obs = (precip.obs.mlist <= 0.01)       # define dry months in obs data
    n.dry.obs<-calc(dry.obs, sum, na.rm=F,     
                    filename<-paste(temp.dir, "temp_raster07", sep=""), overwrite=T)  # number of dry months in obs data     
    
    dry.gcm = (precip.gcm.mlist == 0)           # define dry months in gcm data  
    n.dry.gcm = calc(dry.gcm, sum, na.rm=F,    
                     filename<-paste(temp.dir, "temp_raster08", sep=""), overwrite=T)   # number of dry months in gcm data   
    
    # Set NA to 0 so sorting will work
    cov<-raster(nrows=nrow(P_c_corr), ncols=ncol(P_c_corr),  xmn=xmin(P_c_corr), xmx=xmax(P_c_corr), 
                ymn=ymin(P_c_corr), ymx=ymax(P_c_corr), crs=crs(P_c_corr) )
    cov<-setValues(cov,0)
    precip.gcm.mlist<-cover(precip.gcm.mlist, cov)
    precip.obs.mlist<-cover(precip.obs.mlist, cov)
    
    # list monthly ave precip in rank-ordered list (sort) 
    obs.sort<-calc(precip.obs.mlist, fun=function(x) { sort(x) },
                   filename<-paste(temp.dir, "temp_raster09", sep=""), overwrite=T)
    gcm.sort<-calc(precip.gcm.mlist, fun=function(x) { sort(x) },
                   filename<-paste(temp.dir, "temp_raster10", sep=""), overwrite=T)
    
    # define where eps_m > 0, get values
    # all other values are = 0
    for(i in 1:n.years){
      match = (n.dry.obs == i)
      eps.m = overlay(match, subset(gcm.sort, i), fun=function(x,y){ return(x*y) },
                      filename<-paste(temp.dir, "temp_raster11", sep=""), overwrite=T)
      if(i == 1){
        eps_m_n = eps.m
      }else{
        eps_m_n = eps_m_n + eps.m
      }
    }
    
    if(m==1){eps_m.y<-c()}
    eps_m.y<-add.layers(eps_m.y, eps_m_n, (m==1))
    
    # remove variables no longer needed to save memory space
    # remove('eps.m', 'match')
    
    ### Calculate epsilon_d (dry day threshold) for this month
    
    # an obs day is := dry if daily precip < 1.0 mm/day
    # a  gcm day is := dry if daily precip = 0
    # IF number obs dry days > number gcm 0-precip days, 
    #    THEN define driest gcm days as "dry" in order of increasing precip (starting with least precip)
    #    AND define epsilon_d as the daily precip in the last "dry" gcm day
    # IF number gcm dry days > obs dry days,
    #    THEN exclude only the gcm days with precip = 0, epsilon_d = 0
    
    # NOTE: values from dry months will not be used in estimation of transfer function.
    # the dry months may have very high eps_d, but they won't be used
    
    ### 1.0 mm/day THRESHOLD NEEDS TO BE ASSESSED (given by Hempel et al 2013)
    # set to 0.5 mm/day because 1.0 made too many dry days
    dry.obs = (precip.obs.dlist <= 0.5)             # define dry days in obs data
    n.dry.obs<-calc(dry.obs, sum, na.rm=F,          # number of dry days in obs data
                    filename<-paste(temp.dir, "temp_raster12", sep=""), overwrite=T)          
    
    dry.gcm = (precip.gcm.dlist == 0)              # define dry months in gcm data  
    n.dry.gcm = calc(dry.gcm, sum, na.rm=F,        # number of dry months in gcm data
                     filename<-paste(temp.dir, "temp_raster13", sep=""), overwrite=T)        
    
    # list daily precip in rank-ordered list 
    precip.obs.dlist<-cover(precip.obs.dlist, cov)
    precip.gcm.dlist<-cover(precip.gcm.dlist, cov)
    
    obs.sort.d<-calc(precip.obs.dlist, fun=function(x) { sort(x) },
                     filename<-paste(temp.dir, "temp_raster14", sep=""), overwrite=T)
    gcm.sort.d<-calc(precip.gcm.dlist, fun=function(x) { sort(x) },
                     filename<-paste(temp.dir, "temp_raster15", sep=""), overwrite=T)
    
    # test
    #obs.vec<-extract(obs.sort.d, 10)
    #gcm.vec<-extract(gcm.sort.d, 10)
    #plot(obs.vec~gcm.vec)
    
    # define where eps_d > 0, get values
    # all other values are = 0
    for(i in 1:nlayers(obs.sort.d)){
      match = (n.dry.obs == i)
      eps.d = (match * subset(gcm.sort.d, i))
      if(i == 1){
        eps_d = eps.d
      }else{
        eps_d = eps_d + eps.d
      }
    }
    
    if(m==1){eps_d.y<-c()}
    eps_d.y<-add.layers(eps_d.y, eps_d, (m==1))
    
    ###  Re-distribute precip from "dry" days to "wet" days (only within wet months)
    # resulting vector = P(hat)_ij, Eq.12 (Hempel et al 2013)
    
    # Define wet months
    # for observation data, dry month is defined as mean precip < 0.01, (not epsilon_m)
    obs.wet = (precip.obs.mlist > 0.01)
    obs.wet.month = obs.wet * precip.obs.mlist  # sets dry months to 0
    
    # for gcm data, dry month is defined as mean precip < epsilon_m
    gcm.wet = (precip.gcm.mlist > eps_m_n)
    gcm.wet.month = gcm.wet * precip.gcm.mlist  # sets dry months to 0
    
    # Dry days already defined: dry.obs, dry.gcm; (0,1), where 1=dry
    
    # redistribute the precip data in each month,year individually
    for(i in 1:n.years){
      # select data for this year, month
      obs.wet.y<-subset(obs.wet, i)
      gcm.wet.y<-subset(gcm.wet, i)
      
      precip.obs.y = subset(precip.obs.dlist, ((i-1)*month.lengths[m]+1):(i*month.lengths[m]))
      precip.gcm.y = subset(precip.gcm.dlist, ((i-1)*month.lengths[m]+1):(i*month.lengths[m]))
      
      # define gcm dry days based on epsilon_d (0,1 layer)
      gcm.dry.day<-overlay(precip.gcm.y, eps_d, fun=function(x,y){return(x<=y)},
                           filename<-paste(temp.dir, "temp_raster16", sep=""), overwrite=T)
      
      # defin gcm wet days based on epsilon_d (0,1 layer)
      gcm.wet.day<-overlay(precip.gcm.y, eps_d, fun=function(x,y){return(x>y)},
                           filename<-paste(temp.dir, "temp_raster17", sep=""), overwrite=T)
      
      
      ### OBSERVATION DATA: PRECIP REDISTRIBUTION
      # get values of dry day precip, in wet months only.  
      # obs.wet.y = (0,1) to define wet months
      # (precip.obs.y <= 1.0) = (0,1) defines dry days
      #                                          dry day * wet month * precip data
      precip.obs.dry.data<-overlay((precip.obs.y <= 1.0), obs.wet.y, precip.obs.y,
                                   fun=function(x,y,z){return(x*y*z)},
                                   filename<-paste(temp.dir, "temp_raster18", sep=""), overwrite=T)
      
      # sum precip in dry days (in wet months)
      precip.obs.extra<-calc(precip.obs.dry.data, fun=sum, na.rm=T,
                             filename<-paste(temp.dir, "temp_raster19", sep=""), overwrite=T)
      
      # calculate number of wet days (only in wet months)
      #                               wet day * wet month 
      n.wet.days.obs<-calc( ((precip.obs.y > 1.0)*obs.wet.y), fun=sum, na.rm=T,
                            filename<-paste(temp.dir, "temp_raster20", sep=""), overwrite=T)
      
      # divide sum of extra precip by number of wet days in the month
      # value m_ij, eq. 12 in Hempel et al 2013
      #                          extra precip / number of wet days
      precip.obs.mij<-overlay(precip.obs.extra, n.wet.days.obs, 
                              fun=function(x,y){return(x/y)},
                              filename<-paste(temp.dir, "temp_raster21", sep=""), overwrite=T)
      precip.obs.mij<-cover(precip.obs.mij, cov) # otherwise /0 causes NAs
      
      # add extra precip from dry days to the wet days (eq. 12 in Hempel et al 2013)
      # and set dry day precip = 0
      precip.obs.redist<-overlay(obs.wet.y, (precip.obs.y > 1.0), precip.obs.y, precip.obs.mij,
                                 fun=function(a,b,c,d){ return((a*b*c) + (a*b*d)) },
                                 filename<-paste(temp.dir, "temp_raster22", sep=""), overwrite=T)
      
      # add data to brick with all redistributed precip data
      if(i==1){precip.obs.redist.m<-c()}
      precip.obs.redist.m<-add.layers(precip.obs.redist.m, precip.obs.redist, (i==1))
      
      ### GCM DATA: PRECIP REDISTRIBUTION
      # get values of dry day precip, in wet months only.  
      # gcm.wet.y = (0,1) to define wet months.  based on eps_m_n
      # (precip.gcm.y <= eps_d) = (0,1) defines dry days
      #                               dry day * wet month * precip data
      precip.gcm.dry.data<-overlay(gcm.dry.day, gcm.wet.y,  precip.gcm.y,
                                   fun=function(x,y,z){return(x*y*z)},
                                   filename<-paste(temp.dir, "temp_raster23", sep=""), overwrite=T)
      
      # sum precip in dry days (in wet months)
      precip.gcm.extra<-calc(precip.gcm.dry.data, fun=sum, na.rm=T,
                             filename<-paste(temp.dir, "temp_raster24", sep=""), overwrite=T)
      
      # calculate number of wet days (only in wet months)
      #                     wet day * wet month
      n.wet.days.gcm<-calc( (gcm.wet.day*obs.wet.y), fun=sum, na.rm=T,
                            filename<-paste(temp.dir, "temp_raster25", sep=""), overwrite=T)
      
      # divide sum of extra precip by number of wet days in the month
      # value m_ij, eq. 12 in Hempel et al 2013
      #                          extra precip / number of wet days
      precip.gcm.mij<-overlay(precip.gcm.extra, n.wet.days.gcm, 
                              fun=function(x,y){return(x/y)},
                              filename<-paste(temp.dir, "temp_raster26", sep=""), overwrite=T)
      precip.gcm.mij<-cover(precip.gcm.mij, cov) # otherwise /0 causes NAs
      
      # add extra precip from dry days to the wet days (eq. 12 in Hempel et al 2013)
      # and set dry day precip = 0 (achieved by multiplying by "wet month" and "wet day")
      # wet month precip (all days) also set to 0 (not to be used in regression, assures that slope and intercept parameters (B, A) = 1, 0)
      precip.gcm.redist<-overlay(gcm.wet.y, gcm.wet.day, precip.gcm.y, precip.gcm.mij,
                                 fun=function(a,b,c,d){ return((a*b*c) + (a*b*d)) },
                                 filename<-paste(temp.dir, "temp_raster27", sep=""), overwrite=T)
      
      # add data to brick with all redistributed precip data
      if(i==1){precip.gcm.redist.m<-c()}
      precip.gcm.redist.m<-add.layers(precip.gcm.redist.m, precip.gcm.redist, (i==1))
      
      x<-paste("year step 2", i+start.year-1)
      print(x, quote=F)
    } # end year loop
    
    ###################################################################################
    # calculate this month's (all years') daily residual data from re-distributed precip values 
    # use wet months only, wet days only 
    # Eq. 13 (Hempel et al 2013)
    # These values are used to calculate the transfer function
    
    # Pi(hat) = mean over all wet days in a month (all years)
    # only use positive precip values.
    
    # remove dry months from daily observation data
    for(i in 1:n.years){
      layer.st = ((i-1)*month.lengths[m]+1)
      layer.end = i*month.lengths[m]
      precip.obs.redist.m.wet<-overlay( subset(precip.obs.redist.m, layer.st:layer.end), subset(obs.wet, i), 
                                        fun=function(x,y){return(x*y)}, 
                                        filename<-paste(temp.dir, "temp_raster28", sep=""), overwrite=T)
      precip.gcm.redist.m.wet<-overlay( subset(precip.gcm.redist.m, layer.st:layer.end), subset(gcm.wet, i), 
                                        fun=function(x,y){return(x*y)}, 
                                        filename<-paste(temp.dir, "temp_raster29", sep=""), overwrite=T)
      
      if(i==1){precip.obs.mwet = precip.gcm.mwet <-c()}
      precip.obs.mwet<-add.layers(precip.obs.mwet, precip.obs.redist.m.wet, (i==1))
      precip.gcm.mwet<-add.layers(precip.gcm.mwet, precip.gcm.redist.m.wet, (i==1))      
    }
    
    ### OBS
    # sum precip in wet days                        
    precip.obs.wet.sum<-calc( precip.obs.mwet, fun=sum, na.rm=T, 
                              filename<-paste(temp.dir, "temp_raster30", sep=""), overwrite=T)
    
    # calculate number of wet days in the wet months
    nwet.precip.obs<-calc( (precip.obs.mwet > 1.0), fun=sum, na.rm=T,
                           filename<-paste(temp.dir, "temp_raster31", sep=""), overwrite=T)
    
    # divide precip in wet days by the number of wet days in each month
    Pi_hat_obs = overlay(precip.obs.wet.sum, nwet.precip.obs, fun=function(x,y){return(x/y)},
                         filename<-paste(temp.dir, "temp_raster32", sep=""), overwrite=T)
    
    
    ### GCM
    # sum precip in wet days
    precip.gcm.wet.sum<-calc( precip.gcm.mwet, fun=sum, na.rm=T, 
                              filename<-paste(temp.dir, "temp_raster33", sep=""), overwrite=T)
    
    # calculate number of wet days in the wet months
    nwet.precip.gcm<-calc( (precip.gcm.mwet > eps_d), fun=sum, na.rm=T,
                           filename<-paste(temp.dir, "temp_raster34", sep=""), overwrite=T)
    
    # divide precip in wet days by the number of wet days in each month
    Pi_hat_gcm = overlay(precip.gcm.wet.sum, nwet.precip.gcm, fun=function(x,y){return(x/y)},
                         filename<-paste(temp.dir, "temp_raster35", sep=""), overwrite=T)
    
    ### Deltas
    delta.Pi.obs = overlay(precip.obs.mwet, Pi_hat_obs, 
                           fun=function(x,y){return(x/y)},
                           filename<-paste(temp.dir, "temp_raster36", sep=""), overwrite=T)
    
    delta.Pi.gcm = overlay(precip.gcm.mwet, Pi_hat_gcm, 
                           fun=function(x,y){return(x/y)},
                           filename<-paste(temp.dir, "temp_raster37", sep=""), overwrite=T)
    
    
    ###################################################################################
    # calculate this month's trasnfer function
    # use wet months only, wet days only (dry days & months assigned residual values of NA)
    # output (a,b: intercept, slope) used in Eq.15 and Eq.24 (Hempel et al 2013)
    
    # sort gcm and obs residuals
    # cannot sort with NAs
    delta.Pi.obs<-cover(delta.Pi.obs, cov)
    delta.Pi.gcm<-cover(delta.Pi.gcm, cov)
    
    delta.Pi.obs.sort<-calc(delta.Pi.obs, fun=function(x) { sort(x) },
                            filename<-paste(temp.dir, "temp_raster38", sep=""), overwrite=T)
    
    delta.Pi.gcm.sort<-calc(delta.Pi.gcm, fun=function(x) { sort(x) },
                            filename<-paste(temp.dir, "temp_raster39", sep=""), overwrite=T)
    
    # write output: deltas; one file per month
    delta.Pi.obs.sort.filenm<-paste(out.path, out.name, "_deltaP_OBS_", month.names[m], ".nc", sep="")
    writeRaster(delta.Pi.obs.sort, delta.Pi.obs.sort.filenm, format="CDF", varname = "deltaP", longname="delta P", 
                varunit = "unitless",  xname = 'latitude', yname = 'longitude', overwrite=T)
    
    delta.Pi.gcm.sort.filenm<-paste(out.path, out.name, "_deltaP_GCM_", month.names[m], ".nc", sep="")
    writeRaster(delta.Pi.gcm.sort, delta.Pi.gcm.sort.filenm, format="CDF", varname = "deltaP", longname="delta P", 
                varunit = "unitless",  xname = 'latitude', yname = 'longitude', overwrite=T)
    
    # get intercept(a), and slope(b) of linear regression (obs as function of gcm)
    Deltas<-stack(delta.Pi.obs.sort, delta.Pi.gcm.sort)
    Deltas<-mask(Deltas, subset(data.mask,1))
    
    fun=function(x) { if (is.na(x[1])){ NA } else { 
      lm(x[1:(length(x)/2)] ~ x[((length(x)/2)+1):length(x)])$coefficients[2] }}
    B_slope.m <- calc(Deltas, fun, filename<-paste(temp.dir, "temp_raster40", sep=""), overwrite=T)
    
    fun=function(x) { if (is.na(x[1])){ NA } else { 
      lm(x[1:(length(x)/2)] ~ x[((length(x)/2)+1):length(x)])$coefficients[1] }}
    A_int.m <- calc(Deltas, fun, filename<-paste(temp.dir, "temp_raster41", sep=""), overwrite=T)
    
    # r-squared value not used in bias correction; 
    # included here to show where linear regression does well or poorly
    fun=function(x) { if (is.na(x[1])){ NA } else { 
      summary(lm(x[1:(length(x)/2)] ~ x[((length(x)/2)+1):length(x)]))$r.squared }}
    r2.m <- calc(Deltas, fun, filename<-paste(temp.dir, "temp_raster42", sep=""), overwrite=T)
    
    # correct deltas (used in summary figures)
    Delta.gcm.corr<-overlay(delta.Pi.gcm.sort, A_int.m, B_slope.m,
                            fun=function(x,a,b){return(a + b*x)},
                            filename<-paste(temp.dir, "temp_raster43", sep=""), overwrite=T)
    
    Delta.gcm.corr.filenm<-paste(out.path, out.name, "_deltaP_GCMCORR_", month.names[m], ".nc", sep="")
    writeRaster(Delta.gcm.corr, Delta.gcm.corr.filenm, format="CDF", varname = "deltaP_corr", longname="delta P, corrected", 
                varunit = "unitless",  xname = 'latitude', yname = 'longitude', overwrite=T)
    
    # test
#         obs.vec=extract(delta.Pi.obs.sort, 10)
#         gcm.vec=extract(delta.Pi.gcm.sort, 10)
#         corr.vec=extract(Delta.gcm.corr, 10)
#         plot(obs.vec~gcm.vec)
#         plot(obs.vec~corr.vec)
#         points(obs.vec~corr.vec, col='red')
#         abline(lm(as.numeric(obs.vec)~as.numeric(gcm.vec)))
#         slope=extract(B_slope.m, 10)
#         int=extract(A_int.m, 10)
#         r2=extract(r2.m, 10)
#         res<-int+slope*(gcm.vec)
#         points(res~gcm.vec, col='red')
#         summary(lm(as.numeric(obs.vec)~as.numeric(gcm.vec)))$r.squared
#     

    if(m==1){B_slope = A_int = r2 <-c()}
    B_slope<-add.layers(B_slope, B_slope.m, (m==1))
    A_int  <-add.layers(A_int,   A_int.m, (m==1))
    r2     <-add.layers(r2,     r2.m, (m==1))
    
    x<-paste("month,", m)
    print(x, quote=F)
  } # end month loop
  
  
  # apply mask.  makes all values outside network (eg, ocean grid cells) = NA
  P_c_corr<-mask(P_c_corr, mask=data.mask)
  A_int   <-mask(A_int,    mask=data.mask)
  B_slope <-mask(B_slope,  mask=data.mask)
  eps_m.y <-mask(eps_m.y,  mask=data.mask)
  eps_d.y <-mask(eps_d.y,  mask=data.mask)
  r2.y    <-mask(r2,       mask=data.mask)

  
  # write output: P_c_corr.m, A_int, B_slope, eps_m, eps_d
  P_c_corr.filename<-paste(out.path, out.name, "_P_c_corr.nc", sep="")
  writeRaster(P_c_corr, P_c_corr.filename, format="CDF", varname = "C", longname="C, mean precip correction parameter", 
              varunit = "unitless",  xname = 'latitude', yname = 'longitude', overwrite=T)
  
  A_int.filenm<-paste(out.path, out.name, "_P_A_int.nc", sep="")
  writeRaster(A_int, A_int.filenm, format="CDF", varname = "A_int", longname="A, precip intercept correction parameter", 
              varunit = "unitless",  xname = 'latitude', yname = 'longitude', overwrite=T)
  
  B_slope.filenm<-paste(out.path, out.name, "_P_B_slope.nc", sep="")
  writeRaster(B_slope, B_slope.filenm, format="CDF", varname = "B_slope", longname="B, precip slope correction parameter", 
              varunit = "unitless",  xname = 'latitude', yname = 'longitude', overwrite=T)
  
  eps_m.filenm<-paste(out.path, out.name, "_eps_m.nc", sep="")
  writeRaster(eps_m.y, eps_m.filenm, format="CDF", varname = "eps_m", longname="epsilon_m, dry month threshold", 
              varunit = "ave mm/day",  xname = 'latitude', yname = 'longitude', overwrite=T)
  
  eps_d.filenm<-paste(out.path, out.name, "_eps_d.nc", sep="")
  writeRaster(eps_d.y, eps_d.filenm, format="CDF", varname = "eps_d", longname="epsilon_d, dry day threshold", 
              varunit = "mm/day",  xname = 'latitude', yname = 'longitude', overwrite=T)

  r2.filenm<-paste(out.path, out.name, "_r2.nc", sep="")
  writeRaster(r2.y, r2.filenm, format="CDF", varname = "r2", longname="r squared value of linear transfer function", 
              varunit = "unitless",  xname = 'latitude', yname = 'longitude', overwrite=T)
  
  if(summary.figs == 1){
    month.names<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
    names(A_int) = names(B_slope) = names(P_c_corr) = names(eps_m.y) = names(eps_d.y) = names(r2.y) <- month.names
    
    a.name<-paste(out.path, out.name, "_a_intercept_MAPS.png", sep="")
    png(a.name, height = 1500, width = 2000, res=300)
    levelplot(A_int, margin=F, main="correction factor a: intercept")
    dev.off()
    
    b.name<-paste(out.path, out.name, "_b_slope_MAPS.png", sep="")
    png(b.name, height = 1500, width = 1500, res=300)
    levelplot(B_slope, margin=F, main="correction factor b: slope")
    dev.off()
  
    c.name<-paste(out.path, out.name, "_c_corr_MAPS.png", sep="")
    png(c.name, height = 1500, width = 1500, res=300)
    levelplot(P_c_corr, margin=F, main="correction factor c: mean offset")
    dev.off()
    
    em.name<-paste(out.path, out.name, "_eps_m_MAPS.png", sep="")
    png(em.name, height = 1500, width = 1500, res=300)
    levelplot(eps_m.y, margin=F, main="epsilon_m: dry month threshold [mm/day]")
    dev.off()
    
    ed.name<-paste(out.path, out.name, "_eps_d_MAPS.png", sep="")
    png(ed.name, height = 1500, width = 1500, res=300)
    levelplot(eps_d.y, margin=F, main="epsilon_d: dry day threshold [mm/day]")
    dev.off()
    
    r2.name<-paste(out.path, out.name, "_r2_MAPS.png", sep="")
    png(r2.name, height = 1500, width = 1500, res=300)
    levelplot(r2.y, margin=F, main="r-squared value of linear transfer function")
    dev.off()

  }
} 

#################################################### # END MAIN



