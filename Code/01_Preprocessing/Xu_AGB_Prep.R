# ------------------------------------------------------------------------------------------
# Xu AGB 2000 - 2019
# ------------------------------------------------------------------------------------------

library(raster)

setwd('../../../../')
root.fp = getwd()
data.fp = paste(root.fp, 'Data', sep = '/')

vod.ext = extent(raster(paste(data.fp, 'VODCA', 'VODCA_ABoVE_Land_Mask.tif', sep = '/')))

xu.fp = paste(data.fp, 'Biomass/Xu_AGB/Original_Files', sep = '/')
setwd(xu.fp)

xu.10km = stack('test10a_cd_ab_pred_corr_2000_2019_v2.tif')
xu.years = seq(2000,2019,1)

xu.out.fp = paste(data.fp, 'Biomass/Xu_AGB/Resampled_Files', sep = '/')
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
