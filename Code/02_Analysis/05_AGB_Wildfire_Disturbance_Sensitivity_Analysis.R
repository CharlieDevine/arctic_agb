# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Evaluate sensitivity of AGB_VM and other gridded AGB products (Wang, Xu, Liu) to 
# large-scale wildfire disturbances occurring years 2006-2010 within the domain.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# -----------------------
# Get libraries
# -----------------------
library(raster)
library(sf)
library(moments)
library(rstudioapi)
library(ncdf4)
library(tidyverse)
library(dplyr)
library(reshape2)
library(colorRamps)
library(ggpubr)
library(ggpointdensity)
library(ggpattern)
library(flextable)
library(magick)
library(ggthemes)
library(grid)
library(gridExtra)
library(viridis)
library(pracma)
library(rasterVis)
library(cowplot)

setwd('../../')
gitrepo.fp = getwd()
setwd('../../')
root.fp = getwd()
data.fp = paste(root.fp, 'Data', sep = '/')
figs.fp = paste(gitrepo.fp, 'Figures/Analysis', sep = '/')
code.fp = paste(gitrepo.fp, 'Code', '02_Analysis', sep = '/')

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

# ------------------------
# Wildfire Burn Areas
# ------------------------

ba.shp.fp = paste(gitrepo.fp, 'Data/Shapefiles/Fires', sep = '/')
setwd(ba.shp.fp)
ba = shapefile('MergedBurnAreas_2006to2010.shp')
ba = ba[1:23,] # Exclude last wildfire area since it occupies <50% pixel area

# Open previously created burn area info table
ba.info = read.csv(paste(gitrepo.fp,'Data/BurnAreaInfoTable_2006to2010.csv', sep = '/'))
ba.info$Year = as.character(ba.info$Year)
ba.info$Area..ha. = round(ba.info$Area..ha., 0)
ba.info$Name[is.na(ba.info$Name)] = 'Unnamed'
colnames(ba.info)[2] = 'Fire ID'
colnames(ba.info)[4] = 'Area [ha]'
colnames(ba.info)[7] = 'Data Source'

# ---------------------------------------------------------------------------------------------
# 
# Get all gridded data products
# 
# ---------------------------------------------------------------------------------------------

# ~~~~~~~~~

# ---------------------------------------------------------------------------------------------
# ESA AGB (Santoro)
# ---------------------------------------------------------------------------------------------

# 2017
esa.agb.17.fp = paste(data.fp, 'Biomass/Santoro_ESA/2017', sep = '/')
setwd(esa.agb.17.fp)
esa.agb.17 = raster('Santoro_AGB_2017_ABoVE_025deg.tif')
esa.agb.17[esa.agb.17 == 0] = NA

# 2010
esa.agb.10.fp = paste(data.fp, 'Biomass/Santoro_ESA/2010', sep = '/')
setwd(esa.agb.10.fp)
esa.agb.10 = raster('GlobBiomass_2010_ABOVE_025deg.tif')
esa.agb.10[esa.agb.10 == 0] = NA

# ---------------------------------------------------------------------------------------------
# Xu AGB
# ---------------------------------------------------------------------------------------------
xu.fp = paste(data.fp, 'Biomass/Xu_AGB/Resampled_Files', sep = '/')
setwd(xu.fp)
xu.agb = stack(list.files(xu.fp, pattern = '.tif', full.names = FALSE, recursive = FALSE)[4:18])
xu.agb = xu.agb * 2.2 # ---- Convert from MgC/ha to Mg/ha
years = seq(2003,2017,1)

# ---------------------------------------------------------------------------------------------
# Wang AGB
# ---------------------------------------------------------------------------------------------
wang.fp = paste(data.fp, 'Biomass/Wang_AGB_1984_2014/AGB/Yearly_Mosaics_025deg', sep = '/')
setwd(wang.fp)
wang.agb = stack(list.files(wang.fp, pattern = '.tif', full.names = FALSE, recursive = FALSE)[20:31])

# Extend area to cover full extent of domain
wang.agb = extend(wang.agb, extent(landmask), value = NA, snap = 'near')

# Create Wang mask
wang.area = mean(wang.agb, na.rm = TRUE)
wang.area[wang.area != 0] = 1
wang.mask = extend(wang.area, extent(landmask))
wang.mask.df = as.data.frame(as(wang.mask, 'SpatialPixelsDataFrame'))[,2:3]
wang.footprint = st_union(st_as_sf(rasterToPolygons(wang.mask, dissolve = TRUE)))

# ---------------------------------------------------------------------------------------------
# Liu X-Band AGB
# ---------------------------------------------------------------------------------------------
liu.fp = paste(data.fp, 'Biomass/Liu_AGB', sep = '/')
setwd(liu.fp)

# Get function to read Liu AGB from ncdf file and convert to raster stack
source(paste(code.fp, '00_Read_Liu_AGB.R', sep = '/'))

# Run function
liu.agb = read.liu.agb.fun()
liu.agb = liu.agb[[11:20]] # subset years 2003-2012
liu.agb.mean = cellStats(liu.agb, stat = mean, na.rm = TRUE)

# ---------------------------------------------------------------------------------------------
# VODCA
# ---------------------------------------------------------------------------------------------
# C-Band
c.vod.ym.fp = paste(data.fp, 'VODCA/C_Band/Yearly_MeanGS', sep = '/')
setwd(c.vod.ym.fp)
c.vod.ym = stack(list.files(c.vod.ym.fp, pattern = '.tif', recursive = FALSE, full.names = FALSE))
c.vod.17 = c.vod.ym[[15]]

# X-Band
x.vod.ym.fp = paste(data.fp, 'VODCA/X_Band/Yearly_MeanGS', sep = '/')
setwd(x.vod.ym.fp)
x.vod.ym = stack(list.files(x.vod.ym.fp, pattern = '.tif', recursive = FALSE, full.names = FALSE))
x.vod.17 = x.vod.ym[[15]]

# ---------------------------------------------------------------------------------------------
# VARI (MODIS MCD43A4)
# ---------------------------------------------------------------------------------------------
vari.fp = paste(data.fp, 'MODIS/MCD43A4_VARI/025_Deg', sep = '/')
setwd(vari.fp)
vari.files = list.files(vari.fp, pattern = 'GS', full.names = FALSE, recursive = FALSE)
vari = stack(vari.files)
vari.17 = vari[[15]]

# ---------------------------------------------------------------------------------------------
# NIRV (MODIS MCD43A4)
# ---------------------------------------------------------------------------------------------
nirv.fp = paste(data.fp, 'MODIS/MCD43A4_NDVI+NIRV/025_Deg', sep = '/')
setwd(nirv.fp)
ndvi.files = list.files(nirv.fp, pattern = 'NDVI', full.names = FALSE, recursive = FALSE)
nir.files = list.files(nirv.fp, pattern = 'NIR', full.names = FALSE, recursive = FALSE)
ndvi = stack(ndvi.files)
nir = stack(nir.files)
nirv = ndvi * nir
ndvi.17 = ndvi[[15]]
nirv.17 = nirv[[15]]

# ---------------------------------------------------------------------------------------------
# LAI (MODIS MOD15A2H)
# ---------------------------------------------------------------------------------------------
lai.fp = paste(data.fp, 'MODIS/MOD15A2H_LAI/025_Deg', sep = '/')
setwd(lai.fp)
lai.files = list.files(lai.fp, pattern = 'GS', full.names = FALSE)
lai = stack(lai.files)
lai.17 = lai[[15]]

# ---------------------------------------------------------------------------------------------
# AGB_VM 
# ---------------------------------------------------------------------------------------------
AGB_VM.fp = paste(gitrepo.fp, 'Data/AGB_VM', sep = '/')
setwd(AGB_VM.fp)
AGB_VM = stack(list.files(AGB_VM.fp, pattern = '.tif', recursive = FALSE, full.names = FALSE))
AGB_VM = brick(AGB_VM)


# ---------------------------------------------------------------------------------------------
# Compute annual anomalies for AGB products
# ---------------------------------------------------------------------------------------------

# Create functions
annual.anomalies.fun = function(data.in){
  years = as.character(seq(2003,2017,1))
  data.mean = stackApply(data.in, indices = 1, fun = mean, na.rm = TRUE)
  data.sd = stackApply(data.in, indices = 1, fun = sd, na.rm = TRUE)
  temp.anoms = (data.in - data.mean) / data.sd
  names(temp.anoms) = paste(years, 'Temporal_Anomalies', sep = '_')
  return(temp.anoms)
}

wang.annual.anomalies.fun = function(data.in){
  years = as.character(seq(2003,2014,1))
  data.mean = stackApply(data.in, indices = 1, fun = mean, na.rm = TRUE)
  data.sd = stackApply(data.in, indices = 1, fun = sd, na.rm = TRUE)
  temp.anoms = (data.in - data.mean) / data.sd
  names(temp.anoms) = paste(years, 'Temporal_Anomalies', sep = '_')
  return(temp.anoms)
}

liu.annual.anomalies.fun = function(data.in){
  years = as.character(seq(2003,2012,1))
  data.mean = stackApply(data.in, indices = 1, fun = mean, na.rm = TRUE)
  data.sd = stackApply(data.in, indices = 1, fun = sd, na.rm = TRUE)
  no.change = data.sd
  no.change[no.change > 0] = NA # Isolate static pixels (i.e. no change over time, SD = 0)
  temp.anoms = (data.in - data.mean) / data.sd
  temp.anoms = merge(no.change, temp.anoms) # Merge anomalies maps with "no change" map
  names(temp.anoms) = paste(years, 'Temporal_Anomalies', sep = '_')
  return(temp.anoms)
}


# Compute biomass anomalies
AGB_VM.anoms = annual.anomalies.fun(AGB_VM)
wang.anoms = wang.annual.anomalies.fun(wang.agb)
xu.anoms = annual.anomalies.fun(xu.agb)
liu.anoms = liu.annual.anomalies.fun(liu.agb)

# ---------------------------------------------------------------------------------------------
# Compute statistics for pre- and post-fire AGB for each product's standardized anomalies
# Pre-fire = mean of std. anoms. three years prior to fire
# Post-fire = mean of std. anoms. year of fire and two subsequent years
# Visualize pre- and post-fire distributions for each AGB product (AGB_VM, Wang, Xu, Liu)
# ---------------------------------------------------------------------------------------------

prefire.col = 'forestgreen' # original = yellow
postfire.col = 'yellow' # original = royalblue

# Create function to compute pre- and post-fire AGB anomalies and store within data frame
agb.pre.post.fire.fun = function(){
  
  # Create empty datas frame to store mean pre- and post-fire anomaly values for each AGB product
  agb.pre.post.fire.anoms = ba.info[,1:2]
  nrows = nrow(agb.pre.post.fire.anoms)
  empty.col = vector(mode = 'integer', length = nrows)
  agb.pre.post.fire.anoms$AGB_Anom = empty.col
  agb.pre.post.fire.anoms$Product = empty.col
  agb.pre.post.fire.anoms$PrePost = empty.col
  
  pre.fire = 'Pre fire'
  post.fire = 'Post fire'
  
  # AGB_VM (pre)
  agb.anoms.pre.AGB_VM = agb.pre.post.fire.anoms
  agb.anoms.pre.AGB_VM$Product = rep('AGB_VM', nrows)
  agb.anoms.pre.AGB_VM$PrePost = rep(pre.fire, nrows)
  # AGB_VM (post)
  agb.anoms.post.AGB_VM = agb.pre.post.fire.anoms
  agb.anoms.post.AGB_VM$Product = rep('AGB_VM', nrows)
  agb.anoms.post.AGB_VM$PrePost = rep(post.fire, nrows)
  
  # Wang (pre)
  agb.anoms.pre.wang = agb.pre.post.fire.anoms
  agb.anoms.pre.wang$Product = rep('Wang', nrows)
  agb.anoms.pre.wang$PrePost = rep(pre.fire, nrows)
  # Wang (post)
  agb.anoms.post.wang = agb.pre.post.fire.anoms
  agb.anoms.post.wang$Product = rep('Wang', nrows)
  agb.anoms.post.wang$PrePost = rep(post.fire, nrows)
  
  # Xu (pre)
  agb.anoms.pre.xu = agb.pre.post.fire.anoms
  agb.anoms.pre.xu$Product = rep('Xu', nrows)
  agb.anoms.pre.xu$PrePost = rep(pre.fire, nrows)
  # Xu (post)
  agb.anoms.post.xu = agb.pre.post.fire.anoms
  agb.anoms.post.xu$Product = rep('Xu', nrows)
  agb.anoms.post.xu$PrePost = rep(post.fire, nrows)
  
  # liu (pre)
  agb.anoms.pre.liu = agb.pre.post.fire.anoms
  agb.anoms.pre.liu$Product = rep('Liu', nrows)
  agb.anoms.pre.liu$PrePost = rep(pre.fire, nrows)
  # liu (post)
  agb.anoms.post.liu = agb.pre.post.fire.anoms
  agb.anoms.post.liu$Product = rep('Liu', nrows)
  agb.anoms.post.liu$PrePost = rep(post.fire, nrows)
  
  for (i in 1 : nrows){
    
    fire.year = agb.pre.post.fire.anoms$Year[i]
    fire.id = agb.pre.post.fire.anoms$`Fire ID`[i]
    fire.shp = ba[ba$Fire_ID == fire.id,]
    fire.mask = rasterize(fire.shp, crop(landmask,fire.shp), getCover = TRUE) # Convert fire shape to raster, get %pixel coverage (0-1)
    fire.mask = round(fire.mask, 1)
    fire.mask[fire.mask < 0.5] = NA # Remove pixels with less than 50% shape coverage
    
    print(fire.id)
    
    AGB_VM.fire = cellStats(mask(crop(AGB_VM.anoms, fire.mask), fire.mask), stat = mean)
    wang.fire = cellStats(mask(crop(wang.anoms, fire.mask), fire.mask), stat = mean)
    xu.fire = cellStats(mask(crop(xu.anoms, fire.mask), fire.mask), stat = mean)
    liu.fire = cellStats(mask(crop(liu.anoms, fire.mask), fire.mask), stat = mean)
    
    AGB_VM.index = which(seq(2003,2017,1) == fire.year)
    wang.index = which(seq(2003,2014,1) == fire.year)
    xu.index = which(seq(2003,2017,1) == fire.year)
    liu.index = which(seq(2003,2012,1) == fire.year)
    
    # Pre-fire mean anomalies
    agb.anoms.pre.AGB_VM$AGB_Anom[i] = mean(c(AGB_VM.fire[AGB_VM.index-3], AGB_VM.fire[AGB_VM.index-2], AGB_VM.fire[AGB_VM.index-1]), na.rm = TRUE)
    agb.anoms.pre.wang$AGB_Anom[i] = mean(c(wang.fire[wang.index-3], wang.fire[wang.index-2], wang.fire[wang.index-1]), na.rm = TRUE)
    agb.anoms.pre.xu$AGB_Anom[i] = mean(c(xu.fire[xu.index-3], xu.fire[xu.index-2], xu.fire[xu.index-1]), na.rm = TRUE)
    agb.anoms.pre.liu$AGB_Anom[i] = mean(c(liu.fire[liu.index-3], liu.fire[liu.index-2], liu.fire[liu.index-1]), na.rm = TRUE)
    
    # Post-fire mean anomalies
    agb.anoms.post.AGB_VM$AGB_Anom[i] = mean(c(AGB_VM.fire[AGB_VM.index], AGB_VM.fire[AGB_VM.index+1], AGB_VM.fire[AGB_VM.index+2]), na.rm = TRUE)
    agb.anoms.post.wang$AGB_Anom[i] = mean(c(wang.fire[wang.index], wang.fire[wang.index+1], wang.fire[wang.index+2]), na.rm = TRUE)
    agb.anoms.post.xu$AGB_Anom[i] = mean(c(xu.fire[xu.index], xu.fire[xu.index+1], xu.fire[xu.index+2]), na.rm = TRUE)
    agb.anoms.post.liu$AGB_Anom[i] = mean(c(liu.fire[liu.index], liu.fire[liu.index+1], liu.fire[liu.index+2]), na.rm = TRUE)
  }
  
  anoms.combined = rbind(agb.anoms.pre.AGB_VM, agb.anoms.pre.wang, agb.anoms.pre.xu, agb.anoms.pre.liu,
                         agb.anoms.post.AGB_VM, agb.anoms.post.wang, agb.anoms.post.xu, agb.anoms.post.liu)
  
  anoms.combined = anoms.combined[order(anoms.combined$`Fire ID`),]
  anoms.combined$AGB_Anom[is.nan(anoms.combined$AGB_Anom)] = NA
  
  nfires.AGB_VM = unname(nrow(ba.info) - (colSums(is.na(anoms.combined[anoms.combined$Product %in% 'AGB_VM',]))[3] / 2))
  nfires.wang = unname(nrow(ba.info) - (colSums(is.na(anoms.combined[anoms.combined$Product %in% 'Wang',]))[3] / 2))
  nfires.xu = unname(nrow(ba.info) - (colSums(is.na(anoms.combined[anoms.combined$Product %in% 'Xu',]))[3] / 2))
  nfires.liu = unname(nrow(ba.info) - (colSums(is.na(anoms.combined[anoms.combined$Product %in% 'Liu',]))[3] / 2))
  
  anoms.combined$nFires = vector('character', nrow(anoms.combined))
  
  anoms.combined[anoms.combined$Product %in% 'AGB_VM',][,6] = paste('nFires =',nfires.AGB_VM)
  anoms.combined[anoms.combined$Product %in% 'Wang',][,6] = paste('nFires =',nfires.wang)
  anoms.combined[anoms.combined$Product %in% 'Xu',][,6] = paste('nFires =',nfires.xu)
  anoms.combined[anoms.combined$Product %in% 'Liu',][,6] = paste('nFires =',nfires.liu)
  
  anoms.combined$Product = factor(anoms.combined$Product, levels = c('AGB_VM','Wang','Xu','Liu'))
  anoms.combined$PrePost = factor(anoms.combined$PrePost, levels = c(pre.fire, post.fire))
  
  return(anoms.combined)
}

# Run function
burn.area.agb.anoms = agb.pre.post.fire.fun()

# Create box and whisker plot showing statistical distribution of pre- and post-fire AGB anomalies
ba.agb.boxplot = ggplot(burn.area.agb.anoms,
                        aes(y = AGB_Anom, fill = as.factor(PrePost))) +
  scale_y_continuous(limits = c(-1.75,1.75)) +
  geom_boxplot(position = 'dodge', alpha = 0.5) +
  scale_fill_manual(labels = c('3-year pre-fire mean Z-score', '3-year post-fire mean Z-score'), values = c(prefire.col,postfire.col)) +
  ggtitle('Aggregated for 23 wildfire \nareas (years 2006-2010)') +
  #ylab('z-score') +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  #geom_text(aes(x = 1, y = -1.5, label = paste('Diff mean:',))) +
  theme(#Axes
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 16, face = 'bold'),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_blank(),
        #Main title
        plot.title = element_text(size = 20),
        #Legend
        legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-12,0,0,0),
        legend.title = element_blank(),
        legend.text = element_text(size = 19),
        #Facet grid
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # Separate into separate grid boxes for each AGB product
  facet_wrap(. ~ Product, nrow = 2)


# -----------------------------------------------------------
# Plot timeseries of AGB anomalies within four burn areas
# -----------------------------------------------------------

# Subset individual fire areas
pasfield.2010 = ba[18,]
ryan.2010 = ba[16,]
bigcreek.2009 = ba[5,]
sid.2008 = ba[11,]

# Set graphical parameters
xu.agb.col = 'black'
wang.agb.col = 'tomato'
AGB_VM.col = 'green'
liu.agb.col = 'steelblue'

xu.agb.pch = 18
wang.agb.pch = 15
AGB_VM.pch = 17
liu.agb.pch = 19
esa.agb.pch = 13

pch.cex = 2
cexaxis = 1.8
cexlegend = 2

linetypes = factor(c(1,3,2,1))
pointtypes = factor(c(AGB_VM.pch, wang.agb.pch, xu.agb.pch, liu.agb.pch))

# Create function to plot timeseries
agb.anoms.ts.fun = function(){
  
  firemask.fun = function(fire.area){
    fire.mask = rasterize(fire.area, mask(landmask,fire.area), getCover = TRUE) 
    fire.mask = round(fire.mask, 1)
    fire.mask[fire.mask < 0.5] = NA
    return(fire.mask)
  }
  
  pasfield.2010.firemask = firemask.fun(pasfield.2010)
  ryan.2010.firemask = firemask.fun(ryan.2010)
  bigcreek.2009.firemask = firemask.fun(bigcreek.2009)
  sid.2008.firemask = firemask.fun(sid.2008)
  
  # ~~~~~~~ Burn area subsets
  # Pasfield Fire
  AGB_VM.pf = cellStats(mask(AGB_VM.anoms, pasfield.2010.firemask), stat = mean)
  wang.pf = cellStats(mask(wang.anoms, pasfield.2010.firemask), stat = mean)
  xu.pf = cellStats(mask(xu.anoms, pasfield.2010.firemask), stat = mean)
  liu.pf = cellStats(mask(liu.anoms, pasfield.2010.firemask), stat = mean)
  
  pf.anoms = data.frame(Year = rep(years,4), 
                        AGB = c(AGB_VM.pf, c(wang.pf,NA,NA,NA), xu.pf, c(liu.pf,NA,NA,NA,NA,NA)),
                        Product = c(rep('AGB_VM',15), rep('Wang',15) ,rep('Xu',15), rep('Liu',15)))
  pf.anoms$Product = factor(pf.anoms$Product, levels = c('AGB_VM','Wang','Xu','Liu'))
  pf.year = 2010
  
  # Ryan Fire
  AGB_VM.rf = cellStats(mask(AGB_VM.anoms, ryan.2010.firemask), stat = mean)
  wang.rf = cellStats(mask(wang.anoms, ryan.2010.firemask), stat = mean)
  xu.rf = cellStats(mask(xu.anoms, ryan.2010.firemask), stat = mean)
  liu.rf = cellStats(mask(liu.anoms, ryan.2010.firemask), stat = mean)
  
  rf.anoms = data.frame(Year = rep(years,4), 
                        AGB = c(AGB_VM.rf, c(wang.rf,NA,NA,NA), xu.rf, c(liu.rf,NA,NA,NA,NA,NA)),
                        Product = c(rep('AGB_VM',15), rep('Wang',15) ,rep('Xu',15), rep('Liu',15)))
  rf.anoms$Product = factor(rf.anoms$Product, levels = c('AGB_VM','Wang','Xu','Liu'))
  rf.year = 2010
  
  # Big Creek Fire
  AGB_VM.bcf = cellStats(mask(AGB_VM.anoms, bigcreek.2009.firemask), stat = mean)
  wang.bcf = cellStats(mask(wang.anoms, bigcreek.2009.firemask), stat = mean)
  xu.bcf = cellStats(mask(xu.anoms, bigcreek.2009.firemask), stat = mean)
  liu.bcf = cellStats(mask(liu.anoms, bigcreek.2009.firemask), stat = mean)
  
  bcf.anoms = data.frame(Year = rep(years,4), 
                        AGB = c(AGB_VM.bcf, c(wang.bcf,NA,NA,NA), xu.bcf, c(liu.bcf,NA,NA,NA,NA,NA)),
                        Product = c(rep('AGB_VM',15), rep('Wang',15) ,rep('Xu',15), rep('Liu',15)))
  bcf.anoms$Product = factor(bcf.anoms$Product, levels = c('AGB_VM','Wang','Xu','Liu'))
  bcf.year = 2009
  
  # Sid Fire
  AGB_VM.sf = cellStats(mask(AGB_VM.anoms, sid.2008.firemask), stat = mean)
  wang.sf = cellStats(mask(wang.anoms, sid.2008.firemask), stat = mean)
  xu.sf = cellStats(mask(xu.anoms, sid.2008.firemask), stat = mean)
  liu.sf = cellStats(mask(liu.anoms, sid.2008.firemask), stat = mean)
  
  sf.anoms = data.frame(Year = rep(years,4), 
                         AGB = c(AGB_VM.sf, c(wang.sf,NA,NA,NA), xu.sf, c(liu.sf,NA,NA,NA,NA,NA)),
                         Product = c(rep('AGB_VM',15), rep('Wang',15) ,rep('Xu',15), rep('Liu',15)))
  sf.anoms$Product = factor(sf.anoms$Product, levels = c('AGB_VM','Wang','Xu','Liu'))
  sf.year = 2008
  
  # ~~~~~~~ Plot anomalies 
  pf.anoms$Fire = rep('Pasfield Fire (2010)', 60)
  rf.anoms$Fire = rep('Ryan Fire (2010)', 60)
  bcf.anoms$Fire = rep('Big Creek Fire (2009)', 60)
  sf.anoms$Fire = rep('Sid Fire (2008)', 60)
  
  agb.fire.anoms = rbind(pf.anoms, rf.anoms, bcf.anoms, sf.anoms)
  
  agb.fires.plot = ggplot(data = agb.fire.anoms, aes(x = Year, y = AGB, color = Product, linetype = Product, shape = Product)) +
    scale_y_continuous(limits = c(-2.5,2.5)) +
    scale_x_continuous(breaks = years, limits = c(2003,2017)) +
    ylab('AGB z-scores [Mg/ha]') +
    theme_bw() +
    geom_line(linewidth = 1) +
    scale_linetype_manual(values = linetypes) +
    geom_point(size = 3) +
    scale_shape_manual(values = c(AGB_VM.pch, wang.agb.pch, xu.agb.pch, liu.agb.pch)) +
    scale_color_manual(values = c(AGB_VM.col, wang.agb.col, xu.agb.col, liu.agb.col)) +
    geom_hline(yintercept = 0, color = rgb(0,0,0,0.25), linewidth = 1) +
    facet_grid(Fire ~.) +
    theme(axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 20),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 20),
          panel.grid.minor = element_blank(),
          strip.text.y = element_text(size = 14, face = 'bold'),
          legend.position = c(0.8,0.97),
          legend.direction = 'horizontal',
          legend.text = element_text(size = 20),
          legend.title = element_blank()) +
    # ~~~~ Shade pre- and post-fire years for each fire event
    geom_vline(data = data.frame(Fire = 'Pasfield Fire (2010)'), aes(xintercept = pf.year), color = rgb(1,0.5,0,0.6), linewidth = 1.5) +
    geom_rect(data = data.frame(Fire = 'Pasfield Fire (2010)'), aes(xmin = 2007, xmax = 2009, ymin = -Inf, ymax = Inf), 
                                                                    fill = prefire.col, alpha = 0.25, inherit.aes = FALSE) +
    geom_rect(data = data.frame(Fire = 'Pasfield Fire (2010)'), aes(xmin = 2010, xmax = 2012, ymin = -Inf, ymax = Inf), 
              fill = postfire.col, alpha = 0.25, inherit.aes = FALSE) +
    geom_vline(data = data.frame(Fire = 'Ryan Fire (2010)'), aes(xintercept = rf.year), color = rgb(1,0.5,0,0.6), linewidth = 1.5) +
    geom_rect(data = data.frame(Fire = 'Ryan Fire (2010)'), aes(xmin = 2007, xmax = 2009, ymin = -Inf, ymax = Inf), 
              fill = prefire.col, alpha = 0.25, inherit.aes = FALSE) +
    geom_rect(data = data.frame(Fire = 'Ryan Fire (2010)'), aes(xmin = 2010, xmax = 2012, ymin = -Inf, ymax = Inf), 
              fill = postfire.col, alpha = 0.25, inherit.aes = FALSE) +
    geom_vline(data = data.frame(Fire = 'Big Creek Fire (2009)'), aes(xintercept = bcf.year), color = rgb(1,0.5,0,0.6), linewidth = 1.5) +
    geom_rect(data = data.frame(Fire = 'Big Creek Fire (2009)'), aes(xmin = 2006, xmax = 2008, ymin = -Inf, ymax = Inf), 
              fill = prefire.col, alpha = 0.25, inherit.aes = FALSE) +
    geom_rect(data = data.frame(Fire = 'Big Creek Fire (2009)'), aes(xmin = 2009, xmax = 2011, ymin = -Inf, ymax = Inf), 
              fill = postfire.col, alpha = 0.25, inherit.aes = FALSE) +
    geom_vline(data = data.frame(Fire = 'Sid Fire (2008)'), aes(xintercept = sf.year), color = rgb(1,0.5,0,0.6), linewidth = 1.5) +
    geom_rect(data = data.frame(Fire = 'Sid Fire (2008)'), aes(xmin = 2005, xmax = 2007, ymin = -Inf, ymax = Inf), 
              fill = prefire.col, alpha = 0.25, inherit.aes = FALSE) +
    geom_rect(data = data.frame(Fire = 'Sid Fire (2008)'), aes(xmin = 2008, xmax = 2010, ymin = -Inf, ymax = Inf), 
              fill = postfire.col, alpha = 0.25, inherit.aes = FALSE)
}

# Run function
agb.anoms.ts = agb.anoms.ts.fun()

# Combine AGB anoms timeseries and boxplots into single figure
agb.anoms.ts.bp = plot_grid(plotlist = list(agb.anoms.ts, ba.agb.boxplot),
                            axis = 'lr',
                            ncol = 2,
                            align = 'h',
                            labels = c('a','b'),
                            label_size = 25,
                            hjust = c(-0.5,0),
                            rel_widths = c(1,0.3))

ggsave(filename = paste(figs.fp, 'AGB_Anoms_TS_BP.png', sep = '/'),
       plot = agb.anoms.ts.bp,
       device = 'png',
       width = 18, height = 9, units = 'in')