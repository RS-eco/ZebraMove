# Additional analysis

# Loop for plotting

for(j in 1:nlayers(modis_sur_refl[[1]])){
  modis <- stack(modis_sur_refl[[1]][j],modis_sur_refl[[2]][j],modis_sur_refl[[3]][j], modis_sur_refl[[4]][j], modis_sur_refl[[5]][j], modis_sur_refl[[6]][j])
  if(file.exists(paste0("Results/Modis/mod_tc_", isos[i], "_", ts_surf_refl[j], ".tif"))){modis_tc <- stack(paste0("Results/Modis/mod_tc_", isos[i], "_", ts_surf_refl[j], ".tif"))} else{
    ggRGB(modis, r = 1, g = 2, b = 3, stretch = "lin") + geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), fill="transparent", color="red") + scale_x_continuous(name="Eastings (m)", expand=c(0,0)) + scale_y_continuous(name="Northings (m)", expand=c(0,0))
    ggsave(paste0("Figures/Modis/", isos[i], "_", ts_surf_refl[j], ".pdf"), width=9, height=7)
    beginCluster(n)
    modis_tc <- tasseledCap(modis, sat="MODIS", filename=paste0("Results/Modis/mod_tc_", isos[i], "_", ts_surf_refl[j], ".tif"), format="GTiff")
    endCluster()
  }
  if(!file.exists(paste0("Figures/Modis/Brightness_", isos[i], "_", ts_surf_refl[j],".pdf")) | !file.exists(paste0("Figures/Modis/Greenness_", isos[i], "_", ts_surf_refl[j],".pdf")) | !file.exists(paste0("Figures/Modis/Wetness_", isos[i], "_", ts_surf_refl[j],".pdf"))){
    df_modis_tc <- fortify(modis_tc)
    names(df_modis_tc) <- c("x", "y", "brightness", "greenness", "wetness")
    ggplot() + geom_raster(data=df_modis_tc, aes(x,y,fill=brightness)) + scale_fill_gradientn(name="Brightness", colors = terrain.colors(255)) + scale_x_continuous(name="Eastings (m)", expand=c(0,0)) + scale_y_continuous(name="Northings (m)", expand=c(0,0)) + geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), fill="transparent", color="red")
    ggsave(paste0("Figures/Modis/Brightness_", isos[i], "_", ts_surf_refl[j],".pdf"), width=9, height=7)
    ggplot() + geom_raster(data=df_modis_tc, aes(x,y,fill=greenness)) + scale_fill_gradientn(name="Greenness", colors = terrain.colors(100)) + scale_x_continuous(name="Eastings (m)", expand=c(0,0)) + scale_y_continuous(name="Northings (m)", expand=c(0,0)) + geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), fill="transparent", color="red")
    ggsave(paste0("Figures/Modis/Greenness_", isos[i], "_", ts_surf_refl[j],".pdf"), width=9, height=7)
    ggplot() + geom_raster(data=df_modis_tc, aes(x,y,fill=wetness)) + scale_fill_gradientn(name="Wetness", colors = terrain.colors(100)) + scale_x_continuous(name="Eastings (m)", expand=c(0,0)) + scale_y_continuous(name="Northings (m)", expand=c(0,0)) + geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), fill="transparent", color="red")
    ggsave(paste0("Figures/Modis/Wetness_", isos[i], "_", ts_surf_refl[j],".pdf"), width=9, height=7)
  }
  if(!file.exists(paste0("Results/Modis/fCover_", isos[i], "_", ts_surf_refl[j],".tif"))){
    beginCluster(n)
    modis_fC <- fCover(classImage=ls_awei[[j]], predImage=modis, filename=paste0("Results/Modis/fCover_", isos[i], "_", ts_surf_refl[j], ".tif"), format="GTiff")
    endCluster()
  }
  if(!file.exists(paste0("Figures/Modis/fCover_", isos[i], "_", ts_surf_refl[j],".pdf"))){
    df_modis_fC <- fortify(modis_fC)
    names(df_modis_fC) <- c("x", "y", "fcover")
    ggplot() + geom_raster(data=df_modis_fC, aes(x,y,fill=fcover)) + scale_fill_gradientn(name="Brightness", colors = terrain.colors(255)) + scale_x_continuous(name="Eastings (m)", expand=c(0,0)) + scale_y_continuous(name="Northings (m)", expand=c(0,0)) + geom_polygon(data=studyarea_laea, aes(x=long, y=lat, group=group), fill="transparent", color="red")
    ggsave(paste0("Figures/Modis/fCover_", isos[i], "_", ts_surf_refl[j],".pdf"), width=9, height=7)
  }
}
}

# Fractional cover with AWEI

# Load Landsat AWEI from File
ls_awei <- stack(paste0(workdir, "Results/LS_awei_ns_", isos[i], ".tif"))
ls_awei <- crop(mask(ls_awei, studyarea_laa), studyarea_laea)
### Need to check if layers are order the same way!!! ###
ls_awei <- setZ(ls_awei, ls_dates) # Correct?

# Derive appropriate Landsat scene for each MODIS layers
ls_dates <- sort(ls_dates)
ls_dif <- 16
modis_dates <- data.frame(sort(modis_dates))
modis_dates$ls_dates <- modis_dates$sort.modis_dates.
for(j in 1:nrow(ls_dates)){
  modis_dates$ls_dates[modis_dates$sort.modis_dates. >= ls_dates[j,]-ls_dif/2 & modis_dates$sort.modis_dates. < ls_dates[j,]+ls_dif/2] <- ls_dates[j,]
}

# Subset & repeat stack according to time
library(rts)
ls_awei_rts <- rts(ls_awei, ls_dates)
rts::subset(x, subset=)
# ls_awei_modis <- zApply(ls_awei, by=modis_dates$ls_dates, fun=subset)

if(!file.exists(paste0("Results/fCover_", isos[i], "_brick.grd"))){
  beginCluster(n)
  modis_fC <- mapply(x=ls_awei_modis, y=stack_list, FUN=function(x) fCover(classImage=x, predImage=y))
  endCluster()
  modis_fC <- lapply(modis_fC, stack)
  writeRaster(modis_fC, filename = paste0(workdir, "Results/fCover_", isos[i], "_brick.grd"), bandorder = "BIL", options = c("COMPRESS=NONE"), dataType = "INT2U", overwrite=TRUE)
}

# Change detection

library(bfast) # change detection
library(zoo) # time series handling
library(strucchange) # break detection

# Define helper function
f_bfm <- function(x) {
  x <- ts(x, start = c(2000, 4), frequency = 23)/10000
  bfm <- bfastmonitor(data = x, start = c(2010, 1))
  return(cbind(bfm$breakpoint, bfm$magnitude))
}

for(i in 1:length(isos)){
  modis_ndvi <- stack(paste0(workdir, "Results/ndvi_", i , "_brick.grd"))
  modis_evi <- stack(paste0(workdir, "Results/evi_", i , "_brick.grd"))  
  
  beginCluster(n)
  rbfm_ndvi <- clusteR(x = modis_ndvi, fun = calc, args = list(fun = function(x) {t(apply(x, 1, f_bfm))}))
  rbfm_evi <- clusteR(x = modis_ndvi, fun = calc, args = list(fun = function(x) {t(apply(x, 1, f_bfm))}))
  endCluster()
  
  # Convert time of break and magnitude of break to raster
  timeofbreak_ndvi <- raster(rbfm_ndvi, 1)
  magnitudeofbreak_ndvi <- raster(rbfm_ndvi, 2)
  timeofbreak_evi <- raster(rbfm_evi, 1)
  magnitudeofbreak_evi <- raster(rbfm_evi, 2)
  
  # only visualise the magnitude of the detected breaks
  magnitudeofbreak_ndvi[is.na(timeofbreak_ndvi)] <- NA
  magnitudeofbreak_evi[is.na(timeofbreak_evi)] <- NA
  cairo_pdf(paste0(workdir, "Figures/Magnitude_BreakPoint_", isos[i], ".pdf"), width=8, height=8, bg="transparent")
  par(mfrow=c(2,1))
  plot(magnitudeofbreak_ndvi, zlim = c(-0.2, 0.2), col = rev(diverge_hcl(7)))
  plot(magnitudeofbreak_evi, zlim = c(-0.2, 0.2), col = rev(diverge_hcl(7)))
  dev.off()
}
