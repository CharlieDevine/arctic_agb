library(raster)
library(stringr)

vodca.fp = 'U:/charlie/ABoVE/Data/VODCA'

landmask = raster(paste(vodca.fp,'VODCA_ABoVE_Land_Mask.tif',sep = '/'))

bands = c('C_Band','X_Band','Ku_Band')
bandlabs = c('Cband','Xband','Kuband')

years = as.character(seq(2003,2017,1))

yearly.gs.mean.fun = function(bands,years,bandlabs){
  
  # Disable .aux file creation when writing geotiff files
  rgdal::setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE")
  
  for (i in 1 : length(bands)){
    band = bands[i]
    bandlab = bandlabs[i]
    
    for (j in 1 : length(years)){
      year = years[j]
      gs.days = gsub('-','_', seq(as.Date(paste0(year,'-06-01')), as.Date(paste0(year,'-10-01')), by = 'day'))
      
      fp = paste(vodca.fp,band,'Extended_Domain_Daily',year,sep = '/')
      file.paths = vector(mode = 'character', length = length(gs.days))
      
      for (k in 1 : length(gs.days)){
        file.paths[k] = grep(gs.days[k], list.files(fp, pattern = '.tif', full.names = TRUE, recursive = FALSE), value = TRUE)
      }
      
      print(paste('Computing mean growing season VOD for',band,year))
      gs.stack = stack(file.paths)
      gs.mean = stackApply(gs.stack, indices = 1, fun = mean, na.rm = TRUE)
      
      out.filename = paste(vodca.fp, band, 'Yearly_Mean', paste0(year,'_AnnualMeanVOD_VODCA_',bandlab,'.tif'), sep = '/')
      writeRaster(gs.mean,
                  filename = out.filename,
                  format = 'GTiff')
    }
  }
}

yearly.gs.mean.fun(bands,years,bandlabs)

yearly.gs.max.fun = function(bands,years,bandlabs){
  
  # Disable .aux file creation when writing geotiff files
  rgdal::setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE")
  
  for (i in 1 : length(bands)){
    band = bands[i]
    bandlab = bandlabs[i]
    
    for (j in 1 : length(years)){
      year = years[j]
      gs.days = gsub('-','_', seq(as.Date(paste0(year,'-06-01')), as.Date(paste0(year,'-10-01')), by = 'day'))
      
      fp = paste(vodca.fp,band,'Extended_Domain_Daily',year,sep = '/')
      file.paths = vector(mode = 'character', length = length(gs.days))
      
      for (k in 1 : length(gs.days)){
        file.paths[k] = grep(gs.days[k], list.files(fp, pattern = '.tif', full.names = TRUE, recursive = FALSE), value = TRUE)
      }
      
      print(paste('Computing max growing season VOD for',band,year))
      gs.stack = stack(file.paths)
      gs.max = stackApply(gs.stack, indices = 1, fun = max, na.rm = TRUE)
      
      out.filename = paste(vodca.fp, band, 'Yearly_MaxGS', paste0(year,'_AnnualMaxVOD_VODCA_',bandlab,'.tif'), sep = '/')
      writeRaster(gs.max,
                  filename = out.filename,
                  format = 'GTiff')
    }
  }
}

yearly.gs.max.fun(bands,years,bandlabs)
