library(raster)

setwd('../../')
gitrepo.fp = getwd()
setwd('../../')
root.fp = getwd()
data.fp = paste(root.fp, 'Data', sep = '/')

# -----------------------
# ABoVE domain shapefile
# -----------------------
shp.fp = paste(gitrepo.fp, 'Data/Shapefiles/Domain', sep = '/')
setwd(shp.fp)

core.domain = shapefile('Core_Region_WGS.shp', stringsAsFactors = FALSE)
entire.domain = shapefile('Extended_Plus_Core_Domains_WGS.shp', stringsAsFactors = FALSE)

# ----------------------
# VODCA land mask file
# ----------------------
lm.fp = paste(data.fp, 'VODCA', sep = '/')
setwd(lm.fp)

landmask = raster('VODCA_ABoVE_Land_Mask.tif')
landmask[landmask == 0] = NA

vod.grid = landmask
vod.ext = extent(vod.grid)

# ------------------------------------------------------------------------------------------
# GlobBiomass 2010
# ------------------------------------------------------------------------------------------

glob.fp = paste(data.fp, 'Biomass', 'Santoro_ESA', '2010', sep = '/')
setwd(glob.fp)

glob.orig = raster('GlobBiomass_2010_ABOVE_nativeres_cropped.tif')
glob.agg = aggregate(glob.orig, fact = 290, fun = mean, na.rm = TRUE)
glob.res = resample(glob.agg, vod.grid, method = 'ngb')
writeRaster(glob.res, 'GlobBiomass_2010_ABOVE_025deg.tif', format = 'GTiff')

glob = raster('GlobBiomass_2010_ABOVE_025deg.tif')
glob[glob == 0] = NA

# -----------------------------------------------------------------------------------------
# Santoro AGB 
# -----------------------------------------------------------------------------------------

# ~~~~~~~~~~~~~~~~~~~~~ 2017
# 2017 AGB data
agb.fp = paste(data.fp, 'Biomass', 'Santoro_ESA', '2017', sep = '/')
setwd(agb.fp)

agb.files = list.files(pattern = 'AGB-MERGED')
n80w100 = raster(agb.files[1])
n80w140 = raster(agb.files[2])
n80w180 = raster(agb.files[3])

agb.orig.res = mosaic(n80w140, n80w180, fun = mean)
agb.orig.res = mosaic(agb.orig.res, n80w100, fun = mean)
agb.orig.res.cropped = crop(agb.orig.res, vod.ext)
agb.above.agg = aggregate(agb.orig.res.cropped, fact = 290, fun = mean, na.rm = TRUE)
agb.above.agg.res = resample(agb.above.agg, vod.grid, method = 'ngb')
writeRaster(agb.above.agg.res, 'Santoro_AGB_2017_ABoVE_025deg.tif', format = 'GTiff')

agb.2017 = raster('Santoro_AGB_2017_ABoVE_025deg.tif')

# Standard deviation of AGB data
# agb.sd.files = list.files(pattern = 'SD-MERGED')
# sd.n80w100 = raster(agb.sd.files[1])
# sd.n80w140 = raster(agb.sd.files[2])
# sd.n80w180 = raster(agb.sd.files[3])
# sd.orig.res = mosaic(sd.n80w140, sd.n80w180, fun = mean)
# sd.orig.res = mosaic(sd.orig.res, sd.n80w100, fun = mean)
# sd.orig.res.cropped = crop(sd.orig.res, vod.ext)
# sd.above.agg = aggregate(sd.orig.res.cropped, fact = 290, fun = mean, na.rm = TRUE)
# sd.above.agg.res = resample(sd.above.agg, vod.grid, method = 'ngb')
# writeRaster(sd.above.agg.res, 'Santoro_AGB_STDEV_2017_ABoVE_025deg.tif', format = 'GTiff')
# 
# agb.sd.2017 = raster('Santoro_AGB_STDEV_2017_ABoVE_025deg.tif')
