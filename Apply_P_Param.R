# 2015-07-03
# Apply temperature correction parameters for bias correction
# based on Hempel et al, 2013

# R code written by: Danielle S. Grogan

### ARGUMENTS
# gcm.path
# gcm.name
# gcm.varnm
# gcm.unit.conv = 86400
# A_int.path
# A_int.name
# B_slope.path
# B_slope.name
# C_corr.path
# C_corr.name
# eps_m.path
# eps_m.name
# eps_d.path
# eps_d.name
# start.year
# end.year

library(raster)
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
apply_BCP<-function(gcm.path, gcm.name, gcm.varnm, gcm.unit.conv, 
                    A_int.path, A_int.name, B_slope.path, B_slope.name, C_corr.path, C_corr.name, 
                    eps_m.path, eps_m.name, eps_d.path, eps_d.name,
                    start.year, end.year, network.path, out.dir, out.name, temp.dir){
  
  month.start.days<-c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366)
  month.lengths<-c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  
  ### Load Bias Correction parameters (calculated using the R script "Bias_correction_Tparam.R")
  #  correction factors are monthly (12 layers)
  A_int  <-brick(paste(A_int.path,   A_int.name,   sep=""),  varname='A_int')   
  B_slope<-brick(paste(B_slope.path, B_slope.name, sep=""),  varname='B_slope')
  C_corr <-brick(paste(C_corr.path,  C_corr.name,  sep=""),  varname='C')
  eps_m  <-brick(paste(eps_m.path,   eps_m.name,   sep=""),  varname='eps_m')
  eps_d  <-brick(paste(eps_d.path,   eps_d.name,   sep=""),  varname='eps_d')
  
  # make parameters 365 days (repeat each value for all days within month)
  for(m in 1:12){
    for(d in 1:month.lengths[m]){
      if(m==1 & d==1){
        A_int.365     = brick(subset(A_int,  m))
        B_slope.365   = brick(subset(B_slope,m))
        C_corr.365    = brick(subset(C_corr, m))
        eps_m.365     = brick(subset(eps_m,  m))
        eps_d.365     = brick(subset(eps_d,  m))
      }else{
        A_int.365    <-addLayer(A_int.365,     subset(A_int,    m))
        B_slope.365  <-addLayer(B_slope.365,   subset(B_slope,  m))
        C_corr.365   <-addLayer(C_corr.365,    subset(C_corr,   m))
        eps_m.365    <-addLayer(eps_m.365,     subset(eps_m,    m))
        eps_d.365    <-addLayer(eps_d.365,     subset(eps_d,    m))
      }   
    }
  }
  
  ### Loop through each year of the application period
  for(y in start.yr:end.yr){
    
    # Read climate model (GCM) data
    precip.gcm.nm<-paste(gcm.path, gcm.name, y, ".nc", sep="")
    
    precip.gcm<-brick(precip.gcm.nm, varname=gcm.varnm, stopIfNotEqualSpaced = FALSE)
    
    # bad fix for leap year. improve later:
    if(nlayers(precip.gcm)==366){precip.gcm<-subset(precip.gcm, 1:365)}
    
    if(network.path != 0){
      network<-raster(network.path)
      network.ext<-extent(network)
      # for testing: smaller extent
      #network.ext<-extent(84, 88, 24, 28)  # xmin, xmax, ymin, ymax
      precip.gcm<-crop(precip.gcm, network.ext)
    }
    
    # Re-grid to resolution of correction parameters
    gcm.pr.25<-resample(precip.gcm, A_int, method='bilinear')
    gcm.pr.25 = gcm.pr.25*gcm.unit.conv 
    gcm.pr.25<-mask(gcm.pr.25, eps_d.365)
    
    ### Calculate monthly mean precip
    for(m in 1:12){
      ind<-rep(m, month.lengths[m])
      if(m==1){
        indices = ind
      }else{
        indices = append(indices, ind)
      }
    }
    precip.ave<-stackApply(gcm.pr.25, indices=indices, fun=mean, na.rm=T,
                           filename<-paste(temp.dir, "temp_raster01", sep=""), overwrite=T)
    
    ### Define dry and wet months
    dry.month = (precip.ave <= eps_m)
    
    # make dry and wet month definitions 365 days long
    for(m in 1:12){
      for(d in 1:month.lengths[m]){
        if(m==1 & d==1){dry.month.365 <-c()}
        dry.month.365<-add.layers(dry.month.365, subset(dry.month,m), (m==1 & d==1)) # monthly ave  
      }
    }
    wet.month.365 = (dry.month.365==0) # if a month is not wet, it is dry.
    
    ### if month is defined as dry, only apply the c parameter (dry = ave Pr < epsilon_m)
    #                     dry.month.def * precip.data * c
    P_correct_dry<-overlay(dry.month.365, gcm.pr.25, C_corr.365, 
                           fun=function(x,y,z){return(x*y*z)},
                           filename<-paste(temp.dir, "temp_raster02", sep=""), overwrite=T)
    
    # if month is defined as wet: 
    #   (1) correct frequency of dry days by setting days with less precip than eps_d to 0. Eq. 23 in Hempel et al 2013
    
    # define dry day
    dry.day<-overlay(gcm.pr.25, eps_d.365, fun=function(x,y){return(x<=y)},
                     filename<-paste(temp.dir, "temp_raster03", sep=""), overwrite=T)
    
    #   (2) redistribute precip removed from dry days by evenly adding them to wet days (m_ij)
    # select dry day (wet month) precip
    #                       wet.month.def * precip.data * dry day def
    dry.day.precip<-overlay(wet.month.365, gcm.pr.25,     dry.day, 
                            fun=function(x,y,z){return(x*y*z)}, 
                            filename<-paste(temp.dir, "temp_raster04", sep=""),overwrite=T)
    
    # get sum of dry day (wet month) precip. 12 layers
    dry.day.precip.sum<-stackApply(dry.day.precip, indices=indices,
                                   fun=sum, na.rm=F, 
                                   filename<-paste(temp.dir, "temp_raster05", sep=""), overwrite=T)
    
    # define wet day
    wet.day<-overlay(gcm.pr.25, eps_d.365, fun=function(x,y){return(x>y)},
                     filename<-paste(temp.dir, "temp_raster06", sep=""), overwrite=T)
    
    # calculate number of wet days (only in wet months). 12 layers
    #                     wet day def * wet month def
    n.wet.days<-stackApply(  (wet.day * wet.month.365), indices=indices,
                             fun=sum, na.rm=F,
                             filename<-paste(temp.dir, "temp_raster07", sep=""), overwrite=T)
    
    # divide sum of extra precip by number of wet days in the month. 12 layers
    # value m_ij, eq. 12 in Hempel et al 2013
    #                          extra precip / number of wet days
    precip.gcm.mij<-overlay(dry.day.precip.sum, n.wet.days, 
                            fun=function(x,y){return(x/y)},
                            filename<-paste(temp.dir, "temp_raster08", sep=""), overwrite=T)
    # dividing by 0 makes na values. set to 0, then mask
    precip.gcm.mij[is.na(precip.gcm.mij)]<-c(0)   ### CHECK
    precip.gcm.mij<-mask(precip.gcm.mij, eps_d)
    
    # make precip.gcm.mij 365 days long 
    for(m in 1:12){
      for(d in 1:month.lengths[m]){
        if(m==1 & d==1){precip.mij.365 <-c()}
        precip.mij.365<-add.layers(precip.mij.365, subset(precip.gcm.mij, m), (m==1 & d==1))  
      }
    }
    
    # add extra precip from dry days to the wet days (eq. 12 in Hempel et al 2013)
    # and set dry day precip = 0 (achieved by multiplying by "wet month" and "wet day")                   
    precip.gcm.redist<-overlay(wet.month.365,  wet.day,   gcm.pr.25,    precip.mij.365,
                               fun=function(a,b,c,d){ return((a*b*c) + (a*b*d)) },
                               filename<-paste(temp.dir, "temp_raster09", sep=""), overwrite=T) 
    
    #   (3) normalize redistributed precip by dividing by mean over the wet days (Eq 13)
    # Pi(hat) = mean over all wet days in a month. 
    
    # sum precip in wet days
    #                       sum(wet.day.def * wet.month.def * precip data)
    precip.wet.sum<-stackApply( (wet.day * wet.month.365 * precip.gcm.redist),
                                 indices=indices, fun=sum, na.rm=F,
                                filename<-paste(temp.dir, "temp_raster10", sep=""), overwrite=T)
    
    # divide precip in wet days by the number of wet days in each month
    Pi_hat = overlay(precip.wet.sum, n.wet.days, fun=function(x,y){return(x/y)},
                     filename<-paste(temp.dir, "temp_raster11", sep=""), overwrite=T)
    
    # make Pi_hat_gcm 365 days long 
    for(m in 1:12){
      for(d in 1:month.lengths[m]){
        if(m==1 & d==1){Pi_hat.365 <-c()}
        Pi_hat.365<-add.layers(Pi_hat.365, (subset(Pi_hat,m)), (m==1 & d==1))    
      }
    }
    
    delta.Pi.gcm = overlay(precip.gcm.redist, Pi_hat.365, 
                           fun=function(x,y){return(x/y)},
                           filename<-paste(temp.dir, "temp_raster12", sep=""), overwrite=T)
    
    #   (4) apply a,b,c correction factors
    # linear transfer function for all wet months (Eq. 24, and Eq. 15)
    # Correct precip = c * Pi_hat_gcm * (a + b * delta.Pi.gcm)
    
    P_correct_wet<-overlay(C_corr.365, Pi_hat.365, A_int.365, B_slope.365, delta.Pi.gcm,
                           fun=function(c, p, a, b, d){
                             return( c*p*(a+(b*d)) ) },
                           filename<-paste(temp.dir, "temp_raster13", sep=""), overwrite=T)
    
    P_correct = P_correct_wet + P_correct_dry
    P_correct[P_correct > 400]<-c(400)  # limit max precip to 400 mm/day 
    
    # write output
    P_correct.filenm<-paste(out.dir, out.name, y, ".nc", sep="")
    writeRaster(P_correct, P_correct.filenm, format="CDF", varname = "precip", longname="precipitation, bias-corrected", 
                varunit = "mm/day",  xname = 'latitude', yname = 'longitude', overwrite=T)
    print(y)
  }      # close year loop
  
} # END MAIN #########################################################################

