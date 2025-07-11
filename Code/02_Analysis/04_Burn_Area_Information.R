library(raster)
library(sf)
library(flextable)
library(dplyr)

# This script does the following:
# 1) Reads in burn area polygon shapefiles that were sourced from CNFDB and NIFC data repositories
#    and were then a) spatially queried using the extent of the Wang dataset, b) temporally filtered 
#    to select fires from years 2006-2010, and c) fires >= 62,500 hectares
# 2) LAT/LON coordinates of polygon centroids are extracted for each burn area
# 3) Total number of 0.25-deg grid cells are summed for each burn area polygon
# 4) Burn Area Info is summarized in a data frame showing fire year, code/ID, name, area, 0.25-deg pixel sum,
#    data source, and centroid LAT/LON coordinates; table is exported to .csv
# 5) Burn Area Info flextable is generated and exported

setwd('../../')
gitrepo.fp = getwd()
setwd('../../')
root.fp = getwd()
data.fp = paste(root.fp, 'Data', sep = '/')
figs.fp = paste(gitrepo.fp, 'Figures/Analysis', sep = '/')
code.fp = paste(gitrepo.fp, 'Code', '02_Analysis', sep = '/')
tables.fp = paste(gitrepo.fp, 'CSVs_for_Tables', sep = '/')

# ----------------------
# VODCA land mask file
# ----------------------
lm.fp = paste(data.fp, 'VODCA', sep = '/')
setwd(lm.fp)

landmask = raster('VODCA_ABoVE_Land_Mask.tif')
landmask[landmask == 0] = NA

# ------------------------
# Burn Area Shapefiles
# ------------------------

ba.shp.fp = paste(gitrepo.fp, 'Data/Shapefiles/Fires', sep = '/')
setwd(ba.shp.fp)
ba = st_read('MergedBurnAreas_2006to2010.shp')
ba = ba[1:23,] # Exclude last wildfire area since it occupies <50% pixel area

# ~~~ Create burn area info table
ba.info = data.frame(ba$Year, ba$Fire_ID, ba$Name, ba$Area_HA, ba$Region, ba$Source)
colnames(ba.info) = c('Year', 'Fire ID', 'Name', 'Area [ha]', 'Region', 'Data Source')

# Create and run function to compute burn area polygon centroids
get.ba.centroids = function(){
  centroids.lat = vector(mode = 'integer', length = nrow(ba.info))
  centroids.lon = vector(mode = 'integer', length = nrow(ba.info))
  for (i in 1 : nrow(ba.info)){
    centroids.lat[i] = st_coordinates(st_centroid(st_make_valid(ba$geometry[i])))[,2]
    centroids.lon[i] = st_coordinates(st_centroid(st_make_valid(ba$geometry[i])))[,1]
  }
  centroids = data.frame(centroids.lat, centroids.lon)
  colnames(centroids) = c('Centroid (lat)','Centroid (lon)')
  return(centroids)
}

ba.centroids = get.ba.centroids()

ba.info$CentroidLAT = ba.centroids$`Centroid (lat)`
ba.info$CentroidLON = ba.centroids$`Centroid (lon)`

# Create and run function to compute the total number of 0.25 deg pixels for each burn area
get.ba.npix = function(){
  ba.info$nPixels = vector(mode = 'integer', length = nrow(ba.info))
  for (i in 1 : nrow(ba.info)){
    ba.i = ba[i,]
    # Convert fire shape to raster, get %pixel coverage (0-1)
    fire.mask = raster::rasterize(ba.i,
                                  crop(landmask,ba.i),
                                  getCover = TRUE)
    # Remove pixels with less than 50% shape coverage
    fire.mask = round(fire.mask,1)
    fire.mask[fire.mask < 0.5] = NA
    fire.mask.npix = length(fire.mask[!is.na(fire.mask),])
    ba.info$nPixels[i] = fire.mask.npix
  }
  return(ba.info)
}

ba.info = get.ba.npix()

# Reorder data frame and export to .csv
ba.info = cbind(ba.info[,1:4],ba.info[,9],ba.info[,5:8])
colnames(ba.info)[5] = 'nPixels'
ba.info = ba.info[order(ba.info$`Data Source`, decreasing = TRUE),]

# Write BA info table to .csv
write.csv(ba.info,
          file = paste(gitrepo.fp,'Data/BurnAreaInfoTable_2006to2010.csv', sep = '/'),
          row.names = FALSE)

# Plot flextable showing BA info table
ba.info.for.table = ba.info
ba.info.for.table$Year = as.character(ba.info.for.table$Year)

ba.table = flextable(ba.info.for.table) %>%
  align_text_col(align = 'center', header = TRUE) %>%
  align_nottext_col(align = 'center', header = TRUE) %>%
  footnote(i = 1, j = 7,
           value = as_paragraph(c('CNFDB: https://cwfis.cfs.nrcan.gc.ca/ha/nfdb, NIFC: https://nifc.maps.arcgis.com/home/index.html')),
           part = 'header',
           ref_symbols = c('*'),
           inline = FALSE) %>%
  valign(valign = 'bottom', part = 'header') %>%
  autofit() %>%
  bg(part = 'all', bg = 'white') %>%
  width(j = 2, width = 1.75) %>%
  width(j = 5, width = 0.5) %>%
  width(j = 6, width = 0.75)
  

# Save flextable
save_as_image(ba.table,
              path = paste(figs.fp, 'Burn_Areas_2006to2010.png', sep = '/'))
