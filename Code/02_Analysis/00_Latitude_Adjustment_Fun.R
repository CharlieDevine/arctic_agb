lat.mean.stats.fun = function(area.in) {
  
  # Get global file with 0.25 degree pixel resolution
  global.grid = raster(list.files(paste0(gitrepo.fp,'/Data'), pattern = 'Global_Grid.tif', recursive = FALSE, full.names = TRUE))
  values(global.grid) = NA
  
  rows = nrow(global.grid)
  cols = ncol(global.grid)
  
  xRes = 0.25
  yRes = 0.25
  minx = -180
  miny = 90.0 + (cols*0) + (rows*yRes)
  maxx = -180 + (cols*xRes) + (rows*0)
  maxy = 90.0
  
  radians = 0.0174532925
  radius = 6378.137 #Earth's radius in km
  
  lat = seq(from = maxy, to = miny, length.out = rows)
  lon = seq(from = minx, to = maxx, length.out = cols)
  
  xGrid = matrix(lon, nrow = rows, ncol = cols, byrow = TRUE)
  yGrid = matrix(lat, nrow = rows, ncol = cols, byrow = FALSE)
  
  # Calculate latitude-adjusted per-pixel area in sq. km
  area = (sin(yGrid*radians+0.5*yRes*radians) - sin(yGrid*radians-0.5*yRes*radians)) * (xRes*radians)*radius*radius #km2
  
  # Rasterize
  area.ras = raster(area)
  extent(area.ras) = extent(global.grid)
  crs(area.ras) = crs(global.grid)
  
  # Crop to extent of ABoVE domain
  area.above = crop(area.ras, extent(landmask))
  
  return(mask(area.above, area.in))
}