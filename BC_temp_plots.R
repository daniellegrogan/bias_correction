# 2015-07-16
# plotting: climate data, to go with Bias Correction 
library(rgdal)


month.names<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

#####################################################################################
# monthly temperature means: Box plot
mean.temp.boxplot<-function(path, obs.mave.pre, gcm.mave.pre, gcm.cor.mave.pre, n.years,
                            shape, shape.ID, out.path, out.name){
  # shape can be a shapefile or point file
  # shape.ID = characters, IDs for shapefile or point file
  for(s in 1:length(shape.ID)){
    shape.sub<-shape[s,]
    filenm<-paste(out.path, out.name, "_", shape.ID[s],"_monthly_mean_boxplot.png", sep="")
    filenm2<-sub(" ", "_", c(filenm))
    png(filenm2, height = 1800, width = 2000, res=300)
    nrow = 3
    ncol = 4
    par(mfrow=c(nrow, ncol))
    for(m in 1:12){
      obs.mave.nm<-paste(path, obs.mave.pre, month.names[m],".nc", sep="")
      gcm.mave.nm<-paste(path, gcm.mave.pre, month.names[m],".nc", sep="")
      gcm.cor.mave.nm<-paste(path, gcm.cor.mave.pre, month.names[m],".nc", sep="")
      
      obs.mave<-brick(obs.mave.nm)
      gcm.mave<-brick(gcm.mave.nm)
      gcm.cor.mave<-brick(gcm.cor.mave.nm)
      
      obs<-extract(obs.mave,           shape.sub, fun=mean, df=F, weights=T, na.rm=T)
      gcm<-extract(gcm.mave,           shape.sub, fun=mean, df=F, weights=T, na.rm=T)
      gcm.corr<-extract(gcm.cor.mave,  shape.sub, fun=mean, df=F, weights=T, na.rm=T)
      
      all.d<-data.frame(matrix(nrow=n.years, ncol=3))
      all.d[,1]<-as.numeric(obs)
      all.d[,2]<-as.numeric(gcm)
      all.d[,3]<-as.numeric(gcm.corr)
      
      colnames(all.d)<-c("obs","gcm", "gcm+c")
      boxplot(all.d, main=month.names[m], cex.main=0.8, cex.axis=0.5)
    }
    dev.off()
   
  }
}

#####################################################################################
# daily temp residuals CDF: Figures 3b and 4 in Hempel
cdf.plot<-function(path, obs.delta.pre, gcm.delta.pre, gcm.cor.delta.pre, n.years,
                   shape, shape.ID, out.path, out.name){
  
  # temperature daily residuals CDF
  for(s in 1:length(shape.ID)){
    shape.sub<-shape[s,]
    filenm<-paste(out.path, out.name, "_", shape.ID[s],"_temperature_delta_CDF.png", sep="")
    png(filenm, height = 1800, width = 2000, res=300)
    nrow = 3
    ncol = 4
    par(mfrow=c(nrow, ncol))
    for(m in 1:12){
      obs.delta.nm<-paste(path, obs.delta.pre, month.names[m],".nc", sep="")
      gcm.delta.nm<-paste(path, gcm.delta.pre, month.names[m],".nc", sep="")
      gcm.cor.delta.nm<-paste(path, gcm.cor.delta.pre, month.names[m],".nc", sep="")
      
      obs.delta<-brick(obs.delta.nm)
      gcm.delta<-brick(gcm.delta.nm)
      gcm.cor.delta<-brick(gcm.cor.delta.nm)
      
      obs<-extract(obs.delta,           shape.sub, fun=mean, df=F, weights=T, na.rm=T)
      gcm<-extract(gcm.delta,           shape.sub, fun=mean, df=F, weights=T, na.rm=T)
      gcm.corr<-extract(gcm.cor.delta,  shape.sub, fun=mean, df=F, weights=T, na.rm=T)
      
      xmin=min(obs, gcm, gcm.corr)
      xmax=max(obs, gcm, gcm.corr)
      
      plot(ecdf(obs), pch=".", main=month.names[m], xlim=c(xmin, xmax),
           xlab=c('delta T_ij'), ylab=c("cumulative probabiliyy"),  cex.main=0.5, cex.axis=0.5, cex.lab=0.5)
      lines(ecdf(gcm), pch=".",col='red')
      lines(ecdf(gcm.corr),pch=".", col='blue')
      legend('bottomright', legend=c('obs','gcm', 'gcm corr') , 
             lty=1, col=c('black', 'red', 'blue'), cex=.5)
    }
    dev.off()
  } 
  
  # temperature daily residuals: obs as a function of gcm (gcm.corr)
  for(s in 1:length(shape.ID)){
    shape.sub<-shape[s,]
    filenm<-paste(out.path, out.name, "_", shape.ID[s],"_temperature_delta_TF.png", sep="")
    png(filenm, height = 1800, width = 2000, res=300)
    nrow = 3
    ncol = 4
    par(mfrow=c(nrow, ncol))
    for(m in 1:12){
      
      obs.delta.nm<-paste(path, obs.delta.pre, month.names[m],".nc", sep="")
      gcm.delta.nm<-paste(path, gcm.delta.pre, month.names[m],".nc", sep="")
      gcm.cor.delta.nm<-paste(path, gcm.cor.delta.pre, month.names[m],".nc", sep="")
      
      obs.delta<-brick(obs.delta.nm)
      gcm.delta<-brick(gcm.delta.nm)
      gcm.cor.delta<-brick(gcm.cor.delta.nm)
      cov = gcm.cor.delta
      cov<-setValues(cov,0)
      gcm.cor.delta<-cover(gcm.cor.delta, cov)
      gcm.cor.delta.sort<-calc(gcm.cor.delta, fun = function(x){ sort(x) })
      
      obs<-extract(obs.delta,           shape.sub, fun=mean, df=F, weights=T, na.rm=T)
      gcm<-extract(gcm.delta,           shape.sub, fun=mean, df=F, weights=T, na.rm=T)
      gcm.corr<-extract(gcm.cor.delta.sort,  shape.sub, fun=mean, df=F, weights=T, na.rm=T)
      
      xmin=min(obs, gcm, gcm.corr)
      xmax=max(obs, gcm, gcm.corr)
      
      plot( obs ~ gcm, pch=16, cex=0.5, xlim=c(xmin,xmax), ylim=c(xmin,xmax),
            main=month.names[m], xlab=c('delta T gcm'), ylab=c("delta T obs"),
            cex.main=0.5, cex.axis=0.5, cex.lab=0.5, col='black')
      points(obs ~ gcm.corr, pch=16, cex=0.5, col='blue')
      abline(0,1, lty=2)
      
      legend('topleft', legend=c('obs v gcm', 'obs v gcm corr') , 
             lty=1, col=c('black', 'blue'), cex=.5)
    }
    dev.off()
  } 
}

#####################################################################################
# maps
# c.r.map<-function(C_par, r2, out.path, out.name){
#   month.names<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
#   names(C_par) = names(r2) = month.names
# 
#   C_par.name<-paste(out.path, out.name, "_C_param_MAPS.png", sep="")
#   png(C_par.name, height = 1500, width = 2000, res=300)
#   levelplot(C_par, margin=F, main="mean temperature correction factor") 
#   dev.off()
#   
#   r2.name<-paste(out.path, out.name, "_r2_MAPS.png", sep="")
#   png(r2.name, height = 1500, width = 2000, res=300)
#   levelplot(r2, margin=F, main="r2 of temperature transfer function") 
#   dev.off()
# }







