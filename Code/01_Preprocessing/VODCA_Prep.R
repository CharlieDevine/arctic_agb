# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Compute yearly and monthly statistics for VOD
# Export to geotiff
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(raster)
library(rts)
library(RColorBrewer)

years = as.character(seq(1998,1998,1))
month.num = c('01','02','03','04','05','06','07','08','09','10','11','12')

roi.name = 'ABoVE'
time.step = 'yearly'
vod.band = 'X_Band'

if (roi.name == 'ABoVE'){
  
  # ABoVE domain .shp
  shp.fp = 'U:/charlie/ABoVE/Github/ABoVE_AGB_VOD/Data/Shapefiles/Domain'
  setwd(shp.fp)
  
  shp = shapefile('Extended_Region_Only_WGS.shp.shp')
  
  # Choose hourly or monthly statistics
  if (time.step == 'yearly'){
    
    # Yearly
    out.mean.fp = paste0('C:/Users/cjdev/Desktop/ABoVE/VODCA/',vod.band,'/Yearly_Mean')
    out.max.fp = paste0('C:/Users/cjdev/Desktop/ABoVE/VODCA/',vod.band,'/Yearly_Max')
    out.std.fp = paste0('C:/Users/cjdev/Desktop/ABoVE/VODCA/',vod.band,'/Yearly_Std_Dev')
  }
  
  if (time.step == 'monthly'){
    
    # Monthly
    #out.max.fp = paste0('C:/Users/cjdev/Desktop/ABoVE/VODCA/',vod.band,'/Monthly_Max')
    out.mean.fp = paste0('C:/Users/cjdev/Desktop/ABoVE/VODCA/',vod.band,'/Monthly_Mean')
    #out.std.fp = paste0('C:/Users/cjdev/Desktop/ABoVE/VODCA/',vod.band,'/Monthly_Std_Dev')
  }
}

# ~~~~~~~~~~~~~~~~~~
# Compute statistics
# ~~~~~~~~~~~~~~~~~~
for (i in 1 : length(years)){
  
  year = years[i]
  
  if (roi.name == 'ABoVE'){
    fp = paste0('U:/charlie/ABoVE/Data/VODCA/',vod.band,'/Extended_Domain_Daily/',year)
    setwd(fp)
  }
  
  # ~~~~~~
  files = list.files(fp, pattern = '.tif', recursive = FALSE, full.names = FALSE)
  dates = as.Date(gsub('_','-',substr(files, start = 1, stop = 11)))
  index = as.numeric(format(dates, format = "%m"))
  
  vod = stack(files)
  
  # ~~~~~~~~~~~~~~~~~
  # Max and mean VOD
  # ~~~~~~~~~~~~~~~~~
  
  if (time.step == 'yearly'){
    
    # Yearly statistics
    yearly.max = stackApply(vod, indices = 1, fun = max, na.rm = TRUE)
    yearly.mean = stackApply(vod, indices = 1, fun = mean, na.rm = TRUE)
    yearly.std = stackApply(vod, indices = 1, fun = sd, na.rm = TRUE)
    
    out.file.max = paste0('Max_',vod.band,'_VOD_',year,'.tif')
    out.file.mean = paste0('Mean_',vod.band,'_VOD_',year,'.tif')
    out.file.std = paste0('Std_Dev_',vod.band,'_VOD_',year,'.tif')
    
    # Yearly output
    setwd(out.mean.fp)
    writeRaster(yearly.mean,
                out.file.mean,
                format = 'GTiff')
    
    setwd(out.max.fp)
    writeRaster(yearly.max,
                out.file.max,
                format = 'GTiff')
    
    setwd(out.std.fp)
    writeRaster(yearly.std,
                out.file.std,
                format = 'GTiff')
  }
  
  if (time.step == 'monthly'){
    
    # Monthly statistics
    #monthly.max = stackApply(vod, index, fun = max, na.rm = TRUE)
    #names(monthly.max) = paste(year, month.num, month.abb, 'Max_VOD', vod.band, sep = '_')
    monthly.mean = stackApply(vod, index, fun = mean, na.rm = TRUE)
    names(monthly.mean) = paste(year, month.num, month.abb, 'Mean_VOD', vod.band, sep = '_')
    #monthly.std = stackApply(vod, index, fun = sd, na.rm = TRUE)
    #names(monthly.std) = paste(year, month.num, month.abb, 'Std_Dev_VOD', vod.band, sep = '_')
    
    #out.file.max = paste0(substr(names(monthly.max), start=2,stop=28),'.tif')
    out.file.mean = paste0(substr(names(monthly.mean), start=2,stop=30),'.tif')
    #out.file.std = paste0(substr(names(monthly.std), start=2,stop=30),'.tif')
    
    # Monthly output
    #setwd(out.max.fp)
    #writeRaster(monthly.max,
    #            bylayer = TRUE,
    #            out.file.max,
    #            format = 'GTiff')
    
    setwd(out.mean.fp)
    writeRaster(monthly.mean,
                bylayer = TRUE,
                out.file.mean,
                format = 'GTiff')
    
    #setwd(out.std.fp)
    #writeRaster(monthly.std,
    #            bylayer = TRUE,
    #            out.file.std,
    #            format = 'GTiff')
  }
  
  # ~~~~~~
  # Plot
  #par(oma = c(0,0,0,0),
  #    mar = c(3,3,4,5.5))
  #colors = brewer.pal(9, 'OrRd')
  #col.lim = c(0,1.2)
  #image(yearly.max, main = paste('Yearly Max VOD:',year), col = colors)
  #lines(shp)
  #box()
  #plot(yearly.max, legend.only = TRUE, col = colors, legend.width = 1.6, legend.shrink = 1)
  
  setwd(fp)
}



# --------------------------------------------------------------------


library(ncdf4)
library(raster)
library(xts)

# ~~~~~~~~~~~~~~~~~~~~~~
# ABoVE domain shapefile
# ~~~~~~~~~~~~~~~~~~~~~~
shp.fp = 'U:/charlie/ABoVE/Github/ABoVE_AGB_VOD/Data/Shapefiles/Domain'
setwd(shp.fp)

above = shapefile('Extended_Region_Only_WGS.shp')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VODCA data (C, X, or Ku band)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

years = as.character(seq(1998,1998,1))

for (j in 1 : length(years)){
  # X
  year = years[j]
  fp = paste0('U:/charlie/ABoVE/Data/VODCA/X_Band/Global/',year)
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
  
  #out.names = paste(gsub('-','_',dates),'Kuband_VOD_ABoVE.tif',sep = '_')
  out.names = paste(gsub('-','_',dates),'Xband_VOD_AUS.tif',sep = '_')
  #out.names = paste(gsub('-','_',dates),'Kuband_VOD_AUS.tif',sep = '_')
  
  out.fp = paste0("U:/charlie/ABoVE/Data/VODCA/X_Band/Extended_Domain_Daily/",year)
  #out.fp = paste0("C:/Users/cjdev/Desktop/VODCA_AUS/X_Band/",year)
  #out.fp = paste0("C:/Users/cjdev/Desktop/VODCA_AUS/Ku_Band/",year)
  dir.create(out.fp)
  setwd(out.fp)
  
  writeRaster(out.stack,
              out.names,
              format = "GTiff",
              bylayer = TRUE)
}




# Ku
year = '2002'
fp = paste0('D:/Projects/Data/VODCA/Global_Data/Ku_Band/',year)
setwd(fp)

files = list.files(fp, pattern = '.nc', recursive = TRUE, full.names = FALSE)
dates = as.Date(substr(files, start = 20, stop = 29))

file = nc_open(files[1])
file.date = substr(files[1], start = 20, stop = 29)
vod = ncvar_get(file, 'vod')
nc_close(file)
#lon = ncvar_get(file, 'lon')
#lat = ncvar_get(file, 'lat')
lon.lat.min.max = c(-180.0,180.0,-90.0,90.0)
vod.ras = raster(t(vod))
#vod.crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
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
  #crs(vod.ras.t) = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
  crs(vod.ras.t) = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  names(vod.ras.t) = file.date.t
  vod.ras.t = crop(vod.ras.t, above)
  
  out.stack = stack(out.stack, vod.ras.t)
}

out.names = paste(gsub('-','_',dates),'Kuband_VOD_ABoVE.tif',sep = '_')
#out.names = paste(gsub('-','_',dates),'Xband_VOD_AUS.tif',sep = '_')
#out.names = paste(gsub('-','_',dates),'Kuband_VOD_AUS.tif',sep = '_')

out.fp = paste0("C:/Users/cjdev/Desktop/ABoVE/VODCA/Ku_Band/Extended_Domain_Daily/",year)
#out.fp = paste0("C:/Users/cjdev/Desktop/VODCA_AUS/X_Band/",year)
#out.fp = paste0("C:/Users/cjdev/Desktop/VODCA_AUS/Ku_Band/",year)
dir.create(out.fp)
setwd(out.fp)

writeRaster(out.stack,
            out.names,
            format = "GTiff",
            bylayer = TRUE)

