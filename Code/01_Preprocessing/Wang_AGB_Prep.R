library(raster)
library(gdalUtilities)
library(rgdal)
library(sf)

setwd('../../')
gitrepo.fp = getwd()
shp.fp = paste(gitrepo.fp, 'Data', 'Shapefiles', 'Domain', sep = '/')
setwd('../../')
root.fp = getwd()
data.fp = paste(root.fp, 'Data', sep = '/')

# -----------------------
# ABoVE domain shapefiles
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

# -------------------------------------------------------------------------------------------
# Wang AGB
# -------------------------------------------------------------------------------------------

# ~~~~ For reprojecting
wang.fp = paste0(data.fp, '/Biomass/Wang_AGB_1984_2014')            
agb.fp = paste0(wang.fp, '/Original/AGB')
se.fp = paste0(wang.fp, '/Original/Standard_Error')
agb.files = list.files(agb.fp, pattern = '.tif', full.names = TRUE)
agb.filenames = substr(agb.files, start = 63, stop = 79)                                   
se.files = list.files(se.fp, pattern = '.tif', full.names = TRUE)
se.filenames = substr(se.files, start = 74, stop = 93)                                   

# Scale factor
wang.scalefact = 0.01

# ~~~~ For reprojecting using GDALwarp
# Writes reprojected files to disc
wang.reproject = function(wang.files, wang.filenames) {
  
  for (i in 1 : length(wang.files)){
    wang = raster(wang.files[i])
    
    gdalwarp(wang.files[i], 
             paste0(paste0(wang.fp,'/Processed/01_Reprojected/'),wang.filenames[i],'_rep.tif'), 
             s_srs = '+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs',
             t_srs = '+proj=longlat +datum=WGS84 +no_defs',
             r = 'average')
    print(paste('gdalwarp reproj. complete for',wang.files[i]))
    
    wang.res = raster(paste0(paste0(wang.fp,'/Processed/01_Reprojected/'),wang.filenames[i],'_rep.tif'))
    wang.res.ext = extent(wang.res)
    plot(wang.res * wang.scalefact,
         main = wang.filenames[i])
  }
}

# Run for AGB files
wang.reproject(agb.files, agb.filenames)
# Run for Standard error files
wang.reproject(se.files, se.filenames)

# ~~~~ For re-gridding
# Reads reprojected files, aggrages them to coarser spatial res (0.25deg) and writes to disc
agb.files.rep = list.files(paste0(wang.fp,'/Processed/01_Reprojected'), pattern = '.tif', full.names = TRUE, recursive = FALSE)
agb.filenames.rep = substr(agb.files.rep, start = 75, stop = 91)
wang.years = as.character(seq(1984,2014,1))

wang.regrid = function(wang.files, wang.filenames) {
  
  for(j in 1 : length(wang.files)){
    
    wang.file = brick(wang.files[j])
    wang.filename = wang.filenames[j]
    tile = substr(wang.filenames[j], start = 11, stop = 17)
    reg.path = paste0(wang.fp, '/Processed/02_Regridded/', tile) 
    dir.create(reg.path)
    
    for(i in 1 : nlayers(wang.file)){
      
      temp.ras = wang.file[[i]] 
      #temp.ras[temp.ras == 65535] = 'NA'
      temp.ras = temp.ras * wang.scalefact
      temp.agg.mean = aggregate(temp.ras, fact = 475, fun = mean, na.rm = TRUE)
      temp.agg.mean.res = resample(temp.agg.mean, crop(landmask, temp.agg.mean), method = 'ngb', na.rm = TRUE)
      
      # Disable .aux file creation when writing geotiff files
      setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE")
      
      writeRaster(temp.agg.mean.res,
                  paste0(reg.path, '/', wang.filename,'_',wang.years[i],'_','025deg.tif'),
                  type = 'GTiff')
      
      plot(temp.agg.mean.res,
           main = paste(tile, wang.years[i]),
           zlim = c(0,200))
      
      print(paste('Finished', tile, wang.years[i]))
    }
  }
}


wang.regrid(agb.files.rep, agb.filenames.rep)


# ~~~~ For mosaicing
wang.reg.fp = paste0(wang.fp, '/Processed/02_Regridded')
mos.fp = paste0(wang.fp, '/Processed/03_Mosaic')

mosaic.list.fun = function(file.list){
  
  init.ras1 = raster(file.list[1])
  init.ras2 = raster(file.list[2])
  mos.ras = mosaic(init.ras1, init.ras2, fun = mean)
  
  for (i in 3 : length(file.list)){
    ras = raster(file.list[i])
    mos.ras = mosaic(mos.ras, ras, fun = mean)
  }
  return(mos.ras)
}

for (i in 1 : length(wang.years)){
  wang.year.list = list.files(wang.reg.fp, pattern = wang.years[i], recursive = TRUE, full.names = TRUE)
  wang.year.mos = mosaic.list.fun(wang.year.list)
  writeRaster(wang.year.mos,
              paste0(mos.fp, '/', wang.years[i],'_Wang_Annual_AGB_025deg.tif'),
              format = 'GTiff')
  print(paste0(wang.years[i],' mosaic finished (',i,'/31)'))
}
