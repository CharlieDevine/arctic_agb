library(ncdf4)
library(raster)

read.liu.agb.fun = function() {
  
  liu.nc = nc_open('Aboveground_Carbon_1993_2012.nc')
  
  lat = ncvar_get(liu.nc, 'latitude')
  lon = ncvar_get(liu.nc, 'longitude')
  liu.agb.arr = ncvar_get(liu.nc, 'Aboveground Biomass Carbon')
  nc_close(liu.nc)
  
  arr.to.stack = function(arr){
    ext.ras = c(min(lon), max(lon), min(lat), max(lat))
    crs.ras = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    ras.init = raster(arr[,,1])
    extent(ras.init) = ext.ras
    crs(ras.init) = crs.ras
    ras.init = crop(ras.init, entire.domain)
    out.stack = stack(ras.init)
    
    for (i in 2 : dim(arr)[3]){
      ras = raster(arr[,,i])
      extent(ras) = ext.ras
      crs(ras) = crs.ras
      ras = crop(ras, entire.domain)
      out.stack = stack(ras,out.stack)
    }
    names(out.stack) = as.character(seq(1993,2012,1))
    return(out.stack * 2.2) # ---- Convert from MgC/ha to Mg/ha
  }
  
  liu.agb = arr.to.stack(liu.agb.arr)
  names(liu.agb) = paste('Liu_AGB',as.character(seq(1993,2012,1)), sep = '_')
  
  return(liu.agb)
}