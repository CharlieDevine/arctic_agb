library(raster)

setwd('../../')
gitrepo.fp = getwd()
shp.fp = paste(gitrepo.fp, 'Data', 'Shapefiles', 'Domain', sep = '/')
setwd('../../')
root.fp = getwd()
data.fp = paste(root.fp, 'Data', sep = '/')
esalc.fp = paste(data.fp, 'ESA_CCI_Landcover', sep = '/')

# Read ABoVE domain shapefile
entire.domain = shapefile(paste(shp.fp, 'Extended_Plus_Core_Domains_WGS.shp', sep = '/'))

# Read 0.25-deg land mask grid
landmask.grid = raster(paste(gitrepo.fp, 'Data', 'VODCA_ABoVE_Land_Mask.tif', sep = '/'))

setwd(esalc.fp)

# Read global LC file
lc = raster('ESACCI-LC-L4-LCCS-Map-300m-P1Y-2010-v2.0.7.tif')

# Crop to ABoVE domain extent
lc.domain = crop(lc, entire.domain)

# Remove particular classes
lc.domain[lc.domain > 220] = NA # permanent snow and ice
lc.domain[lc.domain == 210] = NA # water bodies

# Aggregate 300m LC map to 0.25 degree pixel grid (matching VODCA VOD spatial resolution)
# Assign new LC grid value to the most common (i.e., modal) LC value within pixel footprint
lc.domain.agg = aggregate(lc.domain, fact = 90, fun = modal, na.rm = TRUE)

# Resample aggregated LC map to match the 0.25 degree pixel grid exactly
# Since the aggregated grid is already approximately 0.25 degrees, nearest neighbor matching will "snap" the grid into place
lc.domain.res = resample(lc.domain.agg, landmask.grid, method = 'ngb', na.rm = TRUE)

# Write aggregated/resampled LC file to disc
writeRaster(lc.domain.res,
            paste(esalc.fp, '2010_ESA_CCI_LC_025deg.tif', sep = '/'),
            format = 'GTiff')