library(raster)

setwd('../../')
gitrepo.fp = getwd()
setwd('../../')
root.fp = getwd()
data.fp = paste(root.fp, 'Data', sep = '/')

# ----------------------
# VODCA land mask file
# ----------------------
lm.fp = paste(data.fp, 'VODCA', sep = '/')
setwd(lm.fp)

landmask = raster('VODCA_ABoVE_Land_Mask.tif')
landmask[landmask == 0] = NA

# ---------------------------------------------------------------------------------------------
# Function for resampling MODIS MVI datasets to match VODCA grid within extent of ABoVE domain
# ---------------------------------------------------------------------------------------------

mvi.025deg.regrid.fun = function(mvi.name){
  
  modis.fp = paste(data.fp, 'MODIS', sep = '/')
  
  if (mvi.name == 'LAI'){
    mvi.fp = paste(modis.fp, 'MOD15A2H_LAI', sep = '/')
  }
  
  if (mvi.name == 'VARI'){
    mvi.fp = paste(modis.fp, 'MCD43A4_VARI', sep = '/')
    wm = raster(paste(mvi.fp, 'ABoVE_MCD43A4_Watermask.tif', sep = '/'))
    wm[wm == 1] = NA
  }
  
  if (mvi.name %in% c('NDVI','NIR')){
    mvi.fp = paste(modis.fp, 'MCD43A4_NDVI+NIRV', sep = '/')
    wm = raster(paste(mvi.fp, 'ABoVE_MCD43A4_Watermask.tif', sep = '/'))
    wm[wm == 1] = NA
  }
  
  # Set input and output filenames
  mvi.infiles = list.files(paste(mvi.fp,'500m',sep='/'), pattern = mvi.name, full.names = FALSE)
  
  if (mvi.name == 'NIR'){
    mvi.outfiles = paste0(substr(mvi.infiles, start = 1, stop = 30),'_','025deg.tif')
  } 
  
  else {mvi.outfiles = paste0(substr(mvi.infiles, start = 1, stop = 31),'_','025deg.tif')}
  
  # Function to aggregate and resample 
  resample.fun = function(mvi.infile, mvi.outfile){
    setwd(paste(mvi.fp,'500m',sep = '/'))
    mvi.orig = raster(mvi.infile)
    
    # Apply watermask for VARI MVI prior to spatial aggregation (water areas are included in MCD43A4.006)
    if (mvi.name == 'VARI'){
      # Apply water mask
      mvi.orig = mask(mvi.orig, wm)
      # Omit VARI values outside of feasible range (-1 to 1)
      mvi.orig[mvi.orig > 1] = NA
      mvi.orig[mvi.orig < -1] = NA
    }
    
    # Upper limit of LAI values is 10, set any values above this to NA
    if (mvi.name == 'LAI'){
      mvi.orig[mvi.orig > 10] = NA
    }
    
    # Apply watermask for NIR or NDVI prior to spatial aggregation (water areas are included in MCD43A4.006)
    if (mvi.name == 'NIR'){
      mvi.orig = mvi.orig * 0.0001
      mvi.orig = mask(mvi.orig, wm)
    }
    if (mvi.name == 'NDVI'){
      mvi.orig = mask(mvi.orig, wm)
    }
    
    # Spatially aggregate to ~0.25 deg.
    mvi.agg = aggregate(mvi.orig, fact = 55, fun = mean, na.rm = TRUE)
    
    # Resample to precisely match ABoVE VODCA grid area
    mvi.res = resample(mvi.agg, landmask, method = 'ngb')
    mvi.res = mask(mvi.res, landmask)
    
    # ~~~~ Write as new geotiff file
    setwd(paste0(mvi.fp,'/','025_Deg'))
    # Disable .aux file creation when writing geotiff files
    rgdal::setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE")
    writeRaster(mvi.res, mvi.outfile, format = 'GTiff', overwrite = TRUE)
  }
  
  # Loop resample function through all files for the specified MVI
  for (i in 1 : length(mvi.infiles)){
    print(paste('Processing',mvi.infiles[i]))
    resample.fun(mvi.infiles[i], mvi.outfiles[i])
  }
}

mvi.025deg.regrid.fun('NDVI')
mvi.025deg.regrid.fun('NIR')
mvi.025deg.regrid.fun('VARI')
mvi.025deg.regrid.fun('LAI')