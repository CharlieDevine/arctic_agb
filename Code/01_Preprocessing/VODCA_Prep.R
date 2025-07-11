# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Compute yearly and monthly statistics for VOD
# Export to geotiff
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(raster)
library(rts)
library(ncdf4)
library(xts)

setwd('../../')
gitrepo.fp = getwd()
setwd('../../')
root.fp = getwd()
data.fp = paste(root.fp, 'Data', sep = '/')
vod.fp = paste(data.fp, 'VODCA', sep = '/')

years = as.character(seq(2003,2017,1))

roi.name = 'ABoVE'
time.step = 'yearly'
vod.band = 'X_Band' # --------> set for different VOD frequencies (e.g., C_Band or X_Band)

# ~~~~~~~~~~~~~~~~~~~~~~
# ABoVE domain shapefile
# ~~~~~~~~~~~~~~~~~~~~~~
shp.fp = paste(gitrepo.fp, 'Data/Shapefiles/Domain', sep = '/')
setwd(shp.fp)

shp = shapefile('Extended_Region_Only_WGS.shp')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert VODCA data (C or X band) from netCDF to daily geoTIFF and export files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (j in 1 : length(years)){
  # X
  year = years[j]
  fp = paste0(vod.fp,'/',vod.band,'/Global/',year)
  setwd(fp)
  
  files = list.files(fp, pattern = '.nc', recursive = TRUE, full.names = FALSE)
  dates = as.Date(substr(files, start = 20, stop = 29))
  
  file = nc_open(files[1])
  file.date = substr(files[1], start = 20, stop = 29)
  vod = ncvar_get(file, 'vod')
  nc_close(file)
  lon.lat.min.max = c(-180.0,180.0,-90.0,90.0)
  vod.ras = raster(t(vod))
  vod.crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  crs(vod.ras) = vod.crs
  extent(vod.ras) = lon.lat.min.max
  names(vod.ras) = file.date
  
  out.stack = brick(vod.ras)
  out.stack = crop(out.stack, above)
  
  for (i in 2 : length(files)){
    
    file.t = nc_open(files[i])
    file.date.t = substr(files[i], start = 20, stop = 29)
    vod.t = ncvar_get(file.t, 'vod')
    nc_close(file.t)
    vod.ras.t = raster(t(vod.t))
    extent(vod.ras.t) = lon.lat.min.max
    crs(vod.ras.t) = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    names(vod.ras.t) = file.date.t
    vod.ras.t = crop(vod.ras.t, above)
    
    out.stack = stack(out.stack, vod.ras.t)
  }
  
  if (vod.band == 'C_Band') { out.names = paste(gsub('-','_',dates),'Cband_VOD_ABoVE.tif',sep = '_') }
  if (vod.band == 'X_Band') { out.names = paste(gsub('-','_',dates),'Xband_VOD_ABoVE.tif',sep = '_') }
  
  out.fp = paste(vod.fp, vod.band, 'Extended_Domain_Daily', year, sep = '/')
  dir.create(out.fp)
  setwd(out.fp)
  
  writeRaster(out.stack,
              out.names,
              format = "GTiff",
              bylayer = TRUE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in daily VOD band geoTIFF files for each year, compute annual growing season means, and export files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

out.mean.fp = paste(vod.fp, vod.band, 'Yearly_MeanGS', sep = '/')

for (i in 1 : length(years)){
  
  year = years[i]
  
  fp = paste(vod.fp, vod.band, 'Extended_Domain_Daily', year, sep = '/')
  setwd(fp)
  
  # Get daily geoTIFF filenames, dates, and month indices
  files = list.files(fp, pattern = '.tif', recursive = FALSE, full.names = FALSE)
  dates = as.Date(gsub('_','-',substr(files, start = 1, stop = 11)))
  index = as.numeric(format(dates, format = "%m"))
  
  # Subset files to growing season
  files = files[which(index %in% c(6:10))]
  
  vod = stack(files)
  
  # ~~~~~~~~~~~~~~~~~
  # Compute mean VOD
  # ~~~~~~~~~~~~~~~~~
  
  # Yearly growing season mean
  yearly.mean = stackApply(vod, indices = 1, fun = mean, na.rm = TRUE)
  
  if (vod.band == 'C_Band') {out.file.mean = paste(year, 'AnnualMeanVOD', 'VODCA', 'Cband.tif', sep = '_') }
  if (vod.band == 'X_Band') {out.file.mean = paste(year, 'AnnualMeanVOD', 'VODCA', 'Xband.tif', sep = '_') }
  
  # Yearly output
  setwd(out.mean.fp)
  writeRaster(yearly.mean,
              out.file.mean,
              format = 'GTiff')
  
  setwd(fp)
}
