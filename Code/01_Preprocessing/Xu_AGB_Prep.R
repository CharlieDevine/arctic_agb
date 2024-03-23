# ------------------------------------------------------------------------------------------
# Xu AGB 200 - 2019
# ------------------------------------------------------------------------------------------

xu.fp = 'U:/charlie/ABoVE/Data/Xu_AGB/Original_Files'
setwd(xu.fp)

xu.10km = stack('test10a_cd_ab_pred_corr_2000_2019_v2.tif')
xu.years = seq(2000,2019,1)

xu.out.fp = 'U:/charlie/ABoVE/Data/Xu_AGB/Resampled_Files'
setwd(xu.out.fp)

xu.agb.regrid.fun = function(xu.10km){
  for (i in 1 : nlayers(xu.10km)){
    xu = xu.10km[[i]]
    xu.cropped = crop(xu, vod.ext)
    xu.agg = aggregate(xu.cropped, fact = 2, fun = mean, na.rm = TRUE)
    xu.res = resample(xu.agg, vod.mean, method = 'ngb')
    out.filename = paste(xu.years[i], 'Xu_Annual_CarbonDensity_MgCha.tif', sep = '_')
    writeRaster(xu.res, out.filename, format = 'GTiff')
  }
}

xu.agb.regrid.fun(xu.10km)

xu.agb.res = stack(list.files(xu.out.fp, pattern = '.tif'))

for (i in 1 : length(xu.years)){
  plot(xu.agb.res[[i]], zlim = c(0,250), main = paste('Xu Carbon Density',xu.years[i]))
}