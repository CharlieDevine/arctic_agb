# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Calculate pixel-wise temporal trends (% change 2003 vs. 2012, and linear slope 2003-
# 2012) for each AGB product and evaluate statistical distribution of these trends 
# within spatially aggregated vegetation cover classes.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
tables.fp = paste(gitrepo.fp, 'CSVs_for_Tables', sep = '/')

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
# ESA land cover
# ---------------------------------------------------------------------------------------------
lc.fp = paste(data.fp, 'ESA_CCI_Landcover', sep = '/')
setwd(lc.fp)
esa.lc = raster('2010_ESA_CCI_LC_025deg.tif')

# Reclassify and plot land use/land cover
esa.lc.rc = esa.lc
esa.lc.rc[esa.lc.rc %in% c(30,80,160,180,190,201)] = NA
new.class.vals = c(1,2,2,3,3,4,5,6,7,8,9,10,11)
lookup = data.frame(sort(unique(esa.lc.rc[])), new.class.vals)
esa.lc.rc = reclassify(esa.lc.rc, lookup)
esa.lc.df = as.data.frame(as(esa.lc.rc, 'SpatialPixelsDataFrame'))
colnames(esa.lc.df) = c('ClassID','x','y')
esa.lc.df = esa.lc.df[order(esa.lc.df$ClassID, decreasing = FALSE),]

class.names = c('Herbaceous cropland',          #1
                'Deciduous broad leaf forest',  #2
                'Evergreen needle leaf forest', #3
                'Mixed forest',                 #4
                'Mixed forest/shrub',           #5
                'Shrub',                        #6
                'Grass',                        #7
                'Lichens and mosses',           #8
                'Sparse vegetation',            #9
                'Bare',                         #10
                'Permanent snow/ice'            #11
)

populate.class.names = function(){
  df.out = esa.lc.df
  df.out$ClassName = NA
  esa.class.ids = unique(df.out$ClassID)
  for (i in 1 : length(esa.class.ids)){
    df.out[df.out$ClassID == esa.class.ids[i],4] = class.names[i]
  }
  return(df.out)
}

esa.lc.df = populate.class.names()
class.cols = pals::tableau20(n = length(class.names))
names(class.cols) = class.names


# ---------------------------------------------------------------------------------------------
# Pixel-based percent change and slope trends of annual AGB
# ---------------------------------------------------------------------------------------------

# ---------- Percent change function
percent.change.fun = function(agb.in, agb.prod){
  
  percent.change.ras = landmask
  percent.change.ras[] = NA
  
  dims = dim(percent.change.ras)
  
  for (i in 1 : dims[1]){
    for (j in 1 : dims[2]){
      
      agb.vec = agb.in[i,j,]
      agb.pct.chng = ((agb.vec[10] - agb.vec[1]) / agb.vec[1]) * 100
      percent.change.ras[i,j] = agb.pct.chng
      
    }
  }
  percent.change.ras[is.infinite(percent.change.ras)] = NA
  #percent.change.ras[percent.change.ras > 1000] = NA
  percent.change.df = as.data.frame(as(percent.change.ras, 'SpatialPixelsDataFrame'))
  colnames(percent.change.df) = c(agb.prod,'x','y')
  percent.change.df =  melt(percent.change.df, id.vars = c('x','y'), value.name = 'Pct_Change', variable = 'Product')
  return(list(percent.change.ras, percent.change.df))
}

# Run function for each AGB product years 2003-2012
AGB_VM.pct.change = percent.change.fun(AGB_VM[[1:10]],'AGB_VM')
wang.pct.change = percent.change.fun(wang.agb[[1:10]],'Wang')
xu.pct.change = percent.change.fun(xu.agb[[1:10]],'Xu')
liu.pct.change = percent.change.fun(liu.agb,'Liu')

pct.change.df = na.omit(rbind(AGB_VM.pct.change[[2]], 
                              wang.pct.change[[2]],
                              xu.pct.change[[2]], 
                              liu.pct.change[[2]]))

pct.change.df[pct.change.df$Pct_Change == 0,4] = NA


# ---------- Slope trend function
slope.trend.fun = function(agb.in, agb.prod){
  
  trend.years = seq(2003,2012,1)
  
  slope.trend.ras = landmask
  slope.trend.ras[] = NA
  
  dims = dim(slope.trend.ras)
  
  for (i in 1 : dims[1]){
    for (j in 1 : dims[2]){
      
      agb.vec = as.vector(unname(agb.in[i,j,]))
      
      if (is.na(agb.vec[1])) {
        slope.trend.ras[i,j] = NA
      }
      
      if (!is.na(agb.vec[1])) {
        agb.vec.lm = lm(agb.vec ~ trend.years, na.action = na.omit)
        lm.slope = round(unname(agb.vec.lm$coefficients[2]),3)
        slope.trend.ras[i,j] = lm.slope
      }
    }
  }
  slope.trend.df = as.data.frame(as(slope.trend.ras, 'SpatialPixelsDataFrame'))
  colnames(slope.trend.df) = c(agb.prod,'x','y')
  slope.trend.df =  melt(slope.trend.df, id.vars = c('x','y'), value.name = 'Slope', variable = 'Product')
  return(list(slope.trend.ras, slope.trend.df))
}

# Run function for each AGB product years 2003-2012
AGB_VM.slope.trend = slope.trend.fun(AGB_VM[[1:10]], 'AGB_VM')
wang.slope.trend = slope.trend.fun(wang.agb[[1:10]], 'Wang')
xu.slope.trend = slope.trend.fun(xu.agb[[1:10]], 'Xu')
liu.slope.trend = slope.trend.fun(liu.agb, 'Liu')

slope.trend.df = na.omit(rbind(AGB_VM.slope.trend[[2]],
                               wang.slope.trend[[2]],
                               xu.slope.trend[[2]],
                               liu.slope.trend[[2]]))

# ----------------------------------------------------------------------------
# Plot maps showing pixel-wise percent change and annual (slope) trends
# ----------------------------------------------------------------------------

# Set color/graphical parameters for each AGB product
xu.agb.col = 'black'
wang.agb.col = 'tomato'
AGB_VM.col = 'green'
liu.agb.col = 'steelblue'

xu.agb.pch = 18
wang.agb.pch = 15
AGB_VM.pch = 17
liu.agb.pch = 19
esa.agb.pch = 13

# ------------- Percent change
pct.change.maps = ggplot(data = pct.change.df, aes(x = x, y = y, fill = Pct_Change)) +
  geom_raster() +
  ggtitle('Pixel-wise % AGB change 2003 vs. 2012') +
  scale_fill_stepsn(colors = c('magenta','gray96','darkgreen'),
                    # colors = c('red','orange','gray96','palegreen1','darkgreen'),
                    breaks = seq(-100,100,25),
                    na.value = 'grey70',
                    limits = c(-100,100),
                    name = '% AGB change per pixel  (grey = no change)') +
  guides(fill = guide_colorbar(title.position = 'bottom', title.hjust = 0.75, frame.colour = 'black', frame.linewidth = 0.5)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_test() +
  geom_sf(data = st_as_sf(entire.domain), color = 'black', fill = NA, inherit.aes = FALSE) +
  # geom_sf(data = st_as_sf(core.domain), color = 'tomato', fill = NA, inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(wang.footprint), color = wang.agb.col, fill = NA, inherit.aes = FALSE, linewidth = 0.75) +
  geom_point(data = ba.info, aes(x = CentroidLON, y = CentroidLAT), pch = 1, size = 1.5, inherit.aes = FALSE) +
  facet_wrap(.~Product, ncol = 2) +
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 17, face = 'bold'),
        legend.text = element_text(size = 15),
        legend.key.width = unit(3.5, 'cm'),
        legend.key.height = unit(0.5, 'cm'),
        legend.position = 'bottom',
        strip.background = element_rect(fill = 'gray27'),
        strip.text = element_text(size = 17, color = 'white'))

pct.change.maps

# Violin plot showing the distribution shape of pixel-wise percent change values over the extent of the Wang dataset
pct.change.violin = ggplot(data = pct.change.df[pct.change.df$x %in% wang.mask.df$x & pct.change.df$y %in% wang.mask.df$y,],
                           aes(x = 1, y = Pct_Change, color = Product)) +
  geom_hline(yintercept = 0, linewidth = 1.25, color = 'grey70') +
  geom_violin(trim = TRUE,
              linewidth = 1.25,
              show.legend = FALSE,
              na.rm = TRUE,
              fill = 'transparent') +
  scale_color_manual(labels = c('AGB_VM','Wang','Xu','Liu','ESA_CCI'),
                     values = c(AGB_VM.col,wang.agb.col,xu.agb.col,liu.agb.col)) +
  scale_y_continuous(limits = c(-100,100), position = 'right') +
  geom_boxplot(width = 0.1, linewidth = 0.25, na.rm = TRUE, show.legend = FALSE) +
  ylab('Distribution of % AGB change') +
  xlab('') +
  facet_wrap(.~Product, ncol = 1) +
  theme_bw() +
  labs(caption = 'All values masked to \nextent of Wang product \n(panel a, red polygon)') +
  theme(plot.margin = unit(c(0.1,0.1,0.1,-0.5), 'cm'),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.caption = element_text(hjust = 0, size = 15, vjust = 4, face = 'italic'),
        plot.caption.position = 'plot',
        strip.text = element_text(size = 17),
        strip.background = element_rect(fill = 'transparent', color = 'white'))


pct.change.violin

# Combine map and violin plots on single figure
pct.change.maps.and.violin = plot_grid(plotlist = list(pct.change.maps,pct.change.violin),
                                       axis = 'lr',
                                       ncol = 2,
                                       align = 'h',
                                       labels = c('a','b'),
                                       label_size = 20,
                                       hjust = c(-2,1.25),
                                       rel_widths = c(0.7,0.175))


# Save map and combined map-violin figures
ggsave(paste(figs.fp, 'Percent_AGB_Change_2003vs2012.png', sep = '/'),
       pct.change.maps,
       width = 8, height = 8, units = 'in')

ggsave(paste(figs.fp, 'Percent_AGB_Change_2003vs2012_w_violin_plot.png', sep = '/'),
       pct.change.maps.and.violin,
       width = 13, height = 10, units = 'in',
       bg = ' white')


# ------------- Linear slope trends
slope.trend.maps = ggplot(data = slope.trend.df, aes(x = x, y = y, fill = Slope)) +
  geom_raster() +
  ggtitle('Pixel-wise trend (slope) of annual AGB [Mg/ha] 2003-2012') +
  scale_fill_stepsn(colors = c('gold','gray96','royalblue'),
                    breaks = seq(-15,15,3),
                    limits = c(-15,15),
                    name = 'Slope per pixel') +
  guides(fill = guide_colorbar(title.position = 'bottom', title.hjust = 0.5, frame.colour = 'black', frame.linewidth = 0.5)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_test() +
  geom_sf(data = st_as_sf(entire.domain), color = 'black', fill = NA, inherit.aes = FALSE) +
  # geom_sf(data = st_as_sf(core.domain), color = 'tomato', fill = NA, inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(wang.footprint), color = wang.agb.col, fill = NA, inherit.aes = FALSE, linewidth = 0.75) +
  geom_point(data = ba.info, aes(x = CentroidLON, y = CentroidLAT), pch = 1, inherit.aes = FALSE) +
  facet_wrap(.~Product, ncol = 2) +
  theme(plot.title = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12),
        legend.key.width = unit(3.5, 'cm'),
        legend.position = 'bottom',
        strip.background = element_rect(fill = 'gray27'),
        strip.text = element_text(size = 15, color = 'white'))

slope.trend.maps

# Save map figure
ggsave(paste(figs.fp, 'AGB_Pixel_Trends_2003to2012.png', sep = '/'),
       slope.trend.maps,
       width = 8, height = 8, units = 'in')


# ----------------------------------------------------------------------------
# Mask AGB products by land cover class
# ----------------------------------------------------------------------------

ras.obj = landmask

# Rasterize reclassified land cover class data frame
esa.lc.ras = rasterize(x = esa.lc.df[,2:3],
                       y = ras.obj,
                       field = esa.lc.df[,1])

# Subset land cover class masks
clf = esa.lc.ras %in% 1
clf[clf == 0] = NA

dbf = esa.lc.ras %in% 2
dbf[dbf == 0] = NA

enf = esa.lc.ras %in% 3
enf[enf == 0] = NA

mf = esa.lc.ras %in% 4
mf[mf == 0] = NA

mf.csl = esa.lc.ras %in% 5
mf.csl[mf.csl == 0] = NA

csl = esa.lc.ras %in% 6
csl[csl == 0] = NA

grl = esa.lc.ras %in% 7
grl[grl == 0] = NA

tnd = esa.lc.ras %in% 8
tnd[tnd == 0] = NA

spc = esa.lc.ras %in% 9
spc[spc == 0] = NA

bare = esa.lc.ras %in% 10
bare[bare == 0] = NA

# Creat function to mask each AGB product for each cover class and return a data frame
agb.lc.mask.fun = function(agb.in, product) {
  
  agb = mask(agb.in, wang.mask)
  
  agb.clf = mask(agb, clf)        # Herbaceous cropland
  agb.dbf = mask(agb, dbf)        # Deciduous broad leaf forest
  agb.enf = mask(agb, enf)        # Evergreen needle leaf forest
  agb.mf = mask(agb, mf)          # Mixed forest
  agb.mf.csl = mask(agb, mf.csl)  # Mixed forest/shrub
  agb.csl = mask(agb, csl)        # Shrub
  agb.grl = mask(agb, grl)        # Grass
  agb.tnd = mask(agb, tnd)        # Lichens and mosses
  agb.spc = mask(agb, spc)        # Sparse vegetation
  
  agb.lc.df = rbind(data.frame(Pct_Change_AGB = na.omit(values(agb.clf)),
                               LC = 'Herbaceous cropland'),
                    data.frame(Pct_Change_AGB = na.omit(values(agb.dbf)),
                               LC = 'Deciduous broad leaf forest'),
                    data.frame(Pct_Change_AGB = na.omit(values(agb.enf)),
                               LC = 'Evergreen needle leaf forest'),
                    data.frame(Pct_Change_AGB = na.omit(values(agb.mf)),
                               LC = 'Mixed forest'),
                    data.frame(Pct_Change_AGB = na.omit(values(agb.mf.csl)),
                               LC = 'Mixed forest/shrub'),
                    data.frame(Pct_Change_AGB = na.omit(values(agb.csl)),
                               LC = 'Shrub'),
                    data.frame(Pct_Change_AGB = na.omit(values(agb.grl)),
                               LC = 'Grass'),
                    data.frame(Pct_Change_AGB = na.omit(values(agb.tnd)),
                               LC = 'Lichens and mosses'),
                    data.frame(Pct_Change_AGB = na.omit(values(agb.spc)),
                               LC = 'Sparse'))
  
  agb.lc.df$Product = product
  agb.lc.df$Pct_Change_AGB[is.infinite(agb.lc.df$Pct_Change_AGB)] = NA
  return(agb.lc.df)
}


# ------------------ Pct. AGB change per product and land cover class
AGB_VM.pct.change.lc.df = agb.lc.mask.fun(AGB_VM.pct.change[[1]], 'AGB_VM')
wang.pct.change.lc.df = agb.lc.mask.fun(wang.pct.change[[1]], 'Wang')
xu.pct.change.lc.df = agb.lc.mask.fun(xu.pct.change[[1]], 'Xu')
liu.pct.change.lc.df = agb.lc.mask.fun(liu.pct.change[[1]], 'Liu')

# Create function to calculate pct. change quantiles per cover class
agb.lc.pct.change.median.stats = function(pct.change.lc.df, product.name) {
  
  # Get all unique class names
  lc.class.names = unique(pct.change.lc.df$LC)
  nrows.df = length(lc.class.names)
  
  # ---- Create empty data frames
  # Quantile values
  quantile.df = data.frame(LC = lc.class.names,
                           Quantile_000 = numeric(length = nrows.df),
                           Quantile_025 = numeric(length = nrows.df),
                           Quantile_050 = numeric(length = nrows.df),
                           Quantile_075 = numeric(length = nrows.df),
                           Quantile_100 = numeric(length = nrows.df),
                           Stat = 'Quantile')
  
  # pnorm
  pnorm.df = data.frame(LC = lc.class.names,
                        pnorm_000 = numeric(length = nrows.df),
                        pnorm_025 = numeric(length = nrows.df),
                        pnorm_050 = numeric(length = nrows.df),
                        pnorm_075 = numeric(length = nrows.df),
                        pnorm_100 = numeric(length = nrows.df),
                        Stat = 'pnorm')
  
  # Significance (p)
  p.df = data.frame(LC = lc.class.names,
                    p_000 = numeric(length = nrows.df),
                    p_025 = numeric(length = nrows.df),
                    p_050 = numeric(length = nrows.df),
                    p_075 = numeric(length = nrows.df),
                    p_100 = numeric(length = nrows.df),
                    Stat = 'p')
  
  # Standard Error (SE)
  se.df = data.frame(LC = lc.class.names,
                     se_000 = numeric(length = nrows.df),
                     se_025 = numeric(length = nrows.df),
                     se_050 = numeric(length = nrows.df),
                     se_075 = numeric(length = nrows.df),
                     se_100 = numeric(length = nrows.df),
                     Stat = 'SE')
  
  # Loop through LC classes
  for (i in 1 : length(lc.class.names)) {
    
    # Subset pct. change by LC class
    lc.class.name = lc.class.names[i]
    lc.pct.change = pct.change.lc.df[pct.change.lc.df$LC == lc.class.name,][['Pct_Change_AGB']]
    n = length(lc.pct.change)
    
    # Get mean and standard deviation of pct change LC subset
    lc.pct.change.mean = mean(lc.pct.change, na.rm = TRUE)
    lc.pct.change.std = sd(lc.pct.change, na.rm = TRUE)
    
    # Compute quantile values
    lc.pct.change.quantiles = quantile(lc.pct.change, na.rm = TRUE)
    
    # Assign 0%, 25%, 50%, 75%, and 100% quantile values to data frame for LC class subset
    quantile.df[i,2] = round(lc.pct.change.quantiles[1],2)
    quantile.df[i,3] = round(lc.pct.change.quantiles[2],2)
    quantile.df[i,4] = round(lc.pct.change.quantiles[3],2)
    quantile.df[i,5] = round(lc.pct.change.quantiles[4],2)
    quantile.df[i,6] = round(lc.pct.change.quantiles[5],2)
    
    # Compute pnorms for each quantile using mean and stdev of pct. change for LC class subset
    pnorm.df[i,2] = round(pnorm(q = as.numeric(quantile.df[i,2]), mean = lc.pct.change.mean, sd = lc.pct.change.std, lower.tail = FALSE),2)
    pnorm.df[i,3] = round(pnorm(q = as.numeric(quantile.df[i,3]), mean = lc.pct.change.mean, sd = lc.pct.change.std, lower.tail = FALSE),2)
    pnorm.df[i,4] = round(pnorm(q = as.numeric(quantile.df[i,4]), mean = lc.pct.change.mean, sd = lc.pct.change.std, lower.tail = FALSE),2)
    pnorm.df[i,5] = round(pnorm(q = as.numeric(quantile.df[i,5]), mean = lc.pct.change.mean, sd = lc.pct.change.std, lower.tail = FALSE),2)
    pnorm.df[i,6] = round(pnorm(q = as.numeric(quantile.df[i,6]), mean = lc.pct.change.mean, sd = lc.pct.change.std, lower.tail = FALSE),2)
    
    # Get upper and lower confidence intervals for median probabiliy value
    k = qbinom(p = (1-0.95)/2, size = n, prob = 0.5, lower.tail = TRUE) # Using a 95% two-tailed confidence level
    ci = sort(lc.pct.change)[c(k, n - k + 1)]
    attr(ci, 'conf.level') = 1 - 2 * pbinom(q = (k - 1), size = n, prob = 0.5)
    
    # Compute margin of error using confidence interval difference
    me = (ci[2] - ci[1]) / 2
    
    # Compute standard error(SE) by dividing error margin by 1.96
    se = me / 1.96
    
    # Subset pct change values using upper and lower confidence intervals
    lc.pct.change.ci.subset = lc.pct.change[lc.pct.change <= ci[2] & lc.pct.change >= ci[1]]
    
    # Compute two-sided t-test for CI-constrained pct change subset
    ttest = t.test(x = lc.pct.change.ci.subset,
                   conf.level = 1,
                   alternative = 'two.sided')
    
    # Assign output p value from t-test
    p.df[i,4] = ttest$p.value
    
    # Assign output SE value
    se.df[i,4] = se
    
  }
  
  # Combine 50% quantile (median) pct change, pnorm, SE, and p values
  out.df = rbind(data.frame(LC = lc.class.names,
                            Value = quantile.df$Quantile_050,
                            Stat = 'Median %change',
                            Product = product.name),
                 data.frame(LC = lc.class.names,
                            Value = pnorm.df$pnorm_050,
                            Stat = 'pnorm',
                            Product = product.name),
                 data.frame(LC = lc.class.names,
                            Value = se.df$se_050,
                            Stat = 'SE',
                            Product = product.name),
                 data.frame(LC = lc.class.names,
                            Value = p.df$p_050,
                            Stat = 'p',
                            Product = product.name)
  )
  
  return(out.df)
}

# ------------ Run function for each AGB product
# AGB_VM
AGB_VM.pct.change.lc.median.stats = agb.lc.pct.change.median.stats(AGB_VM.pct.change.lc.df, 'AGB_VM')

# Wang
wang.pct.change.lc.median.stats = agb.lc.pct.change.median.stats(wang.pct.change.lc.df, 'Wang')

# Xu
xu.pct.change.lc.median.stats = agb.lc.pct.change.median.stats(xu.pct.change.lc.df, 'Xu')

# Liu
liu.pct.change.lc.median.stats = agb.lc.pct.change.median.stats(liu.pct.change.lc.df, 'Liu')

# Combine quantile pct. change info
agb.pct.lc.change.median.stats = rbind(AGB_VM.pct.change.lc.median.stats,
                                       wang.pct.change.lc.median.stats,
                                       xu.pct.change.lc.median.stats,
                                       liu.pct.change.lc.median.stats)


# --- Generate table for 50% quantile (median), standard error (SE), and significance (p) of pct. AGB change per LC class and product
agb.pct.lc.change.median.stats = reshape(agb.pct.lc.change.median.stats, idvar = c('LC','Stat'), timevar = 'Product', v.names = 'Value', direction = 'wide')
colnames(agb.pct.lc.change.median.stats) = c('LC','Stat','AGB_VM','Wang','Xu','Liu')

agb.pct.lc.change.median.stats = data.frame(LandCoverClass = agb.pct.lc.change.median.stats$LC[1:9],
                                            AGBVM_median = agb.pct.lc.change.median.stats$AGB_VM[1:9],
                                            AGBVM_SE = agb.pct.lc.change.median.stats$AGB_VM[19:27],
                                            AGBVM_p = agb.pct.lc.change.median.stats$AGB_VM[28:36],
                                            Wang_median = agb.pct.lc.change.median.stats$Wang[1:9],
                                            Wang_SE = agb.pct.lc.change.median.stats$Wang[19:27],
                                            Wang_p = agb.pct.lc.change.median.stats$Wang[28:36],
                                            Xu_median = agb.pct.lc.change.median.stats$Xu[1:9],
                                            Xu_SE = agb.pct.lc.change.median.stats$Xu[19:27],
                                            Xu_p = agb.pct.lc.change.median.stats$Xu[28:36],
                                            Liu_median = agb.pct.lc.change.median.stats$Liu[1:9],
                                            Liu_SE = agb.pct.lc.change.median.stats$Liu[19:27],
                                            Liu_p = agb.pct.lc.change.median.stats$Liu[28:36])


# --- Round and relable values for plotting on table
# Median pct change, p, and SE
agb.pct.lc.change.median.stats[,c(2,3,5,6,8,9,11,12)] = round(agb.pct.lc.change.median.stats[,c(2,3,5,6,8,9,11,12)],2)
agb.pct.lc.change.median.stats[agb.pct.lc.change.median.stats$Liu_SE == 0.00,12] = NA
agb.pct.lc.change.median.stats[is.nan(agb.pct.lc.change.median.stats$Liu_p),13] = NA

# Re-label p values smaller than 0.001
agb.pct.lc.change.relabel.pvals = function() {
  
  data.out = agb.pct.lc.change.median.stats
  
  AGB_VM.pvals = data.out$AGBVM_p
  wang.pvals = data.out$Wang_p
  xu.pvals = data.out$Xu_p
  liu.pvals = data.out$Liu_p
  
  small.pval.label = '< 0.001'
  
  for (i in 1:9) {
    if (isTRUE(AGB_VM.pvals[i] < 0.001)) { data.out$AGBVM_p[i] = small.pval.label }
    if (isTRUE(AGB_VM.pvals[i] > 0.001)) { data.out$AGBVM_p[i] = round(AGB_VM.pvals[i],3) }
    if (isTRUE(wang.pvals[i] < 0.001)) { data.out$Wang_p[i] = small.pval.label }
    if (isTRUE(wang.pvals[i] > 0.001)) { data.out$Wang_p[i] = round(wang.pvals[i],3) }
    if (isTRUE(xu.pvals[i] < 0.001)) { data.out$Xu_p[i] = small.pval.label }
    if (isTRUE(xu.pvals[i] > 0.001)) { data.out$Xu_p[i] = round(xu.pvals[i],3) }
    if (isTRUE(liu.pvals[i] < 0.001)) { data.out$Liu_p[i] = small.pval.label }
  }
  return(data.out)
}

# Relabel
agb.pct.lc.change.median.stats = agb.pct.lc.change.relabel.pvals()

# Set graphical (color) parameters for table
pct.change.lc.cols = scales::col_factor(palette = unname(class.cols[1:9]),
                                        domain = agb.pct.lc.change.median.stats$LandCoverClass[1:9],
                                        levels = agb.pct.lc.change.median.stats$LandCoverClass[1:9])

pct.change.median.cols = scales::col_numeric(palette = c('magenta','gray96','darkgreen'),
                                             domain = c(-9, 0, 9))

# Create flextable showing stats for percent change per AGB product and cover class
agb.pct.lc.change.median.stats.table = flextable(agb.pct.lc.change.median.stats) %>%
  colformat_double(digits = 2) %>%
  separate_header(opts = 'span-top') %>%
  theme_vanilla() %>%
  align_text_col(align = 'center') %>%
  align_nottext_col(align = 'center') %>%
  bold(part = 'header') %>%
  vline(j = c(1,4,7,10)) %>%
  bg(i = 1:9,
     j = 1,
     bg = pct.change.lc.cols) %>%
  bg(i = 1:9,
     j = c(2,5,8,11),
     bg = pct.change.median.cols) %>%
  bg(i = 1:9,
     j = c(3,4,6,7,9,10,12,13),
     bg = 'white') %>%
  bg(i = c(1,2), 
     bg = 'white',
     part = 'header') %>%
  border_outer(part = 'all', border = fp_border_default(color = 'black', width = 2)) %>%
  fix_border_issues() %>%
  set_caption(caption = as_paragraph(as_chunk('% AGB change [Mg/ha] 2003 vs. 2012 per land cover class',
                                              props = fp_text_default(font.size = 15, bold = TRUE))),
              fp_p = officer::fp_par(text.align = 'center', 
                                     padding = 5,
                                     shading.color = 'white')) %>%
  fontsize(size = 13, part = 'body') %>%
  fontsize(size = 13, part = 'header') %>%
  width(j = c(4,7,10,13), width = 5) %>%
  width(j = 1, width = 20) %>%
  set_table_properties(layout = 'autofit',
                       width = 1) %>%
  labelizor(labels = c(LandCoverClass = 'Land Cover Class',
                       AGBVM = 'AGB_VM',
                       Wang = 'Wang',
                       Xu = 'Xu',
                       Liu = 'Liu'),
            part = 'all')


agb.pct.lc.change.median.stats.table
# Save manually

# Write data frame to .csv
write.csv(agb.pct.lc.change.median.stats,
          file = paste(tables.fp, 'Table_01_Pct_AGB_Change_2003vs2012_VegClasses.csv', sep = '/'),
          row.names = FALSE)


# Box-whisker figure showing distribution of %change 2003v2012 from each AGB product for each LC class
pct.change.boxplot.fun = function(){
  
  lc.pct.change.agb = na.omit(rbind(AGB_VM.pct.change.lc.df,
                                    wang.pct.change.lc.df,
                                    xu.pct.change.lc.df,
                                    liu.pct.change.lc.df))
  
  
  lc.pct.change.agb$Product = factor(lc.pct.change.agb$Product, levels = c('AGB_VM','Wang','Xu','Liu'))
  lc.pct.change.agb$LC = factor(lc.pct.change.agb$LC, levels = unique(lc.pct.change.agb$LC))
  
  facet.cols = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = unname(class.cols)[1:9]))
  
  # Box-whisker plot
  lc.areas.plot = ggplot(data = lc.pct.change.agb,
                         aes(x = Product, y = Pct_Change_AGB, color = Product, fill = Product)) +
    geom_boxplot() +
    stat_summary(fun = median, geom = 'crossbar', color = 'white', fatten = 0.75) +
    ylab('% change 2003 vs. 2012') +
    scale_y_continuous(limits = c(-60,60)) +
    scale_fill_manual(labels = c('AGB_VM','Wang','Xu','Liu'),
                      values = c(AGB_VM.col,wang.agb.col,xu.agb.col,liu.agb.col)) +
    scale_color_manual(labels = c('AGB_VM','Wang','Xu','Liu'),
                       values = c(AGB_VM.col,wang.agb.col,xu.agb.col,liu.agb.col)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    #facet_wrap(.~LC, ncol = 9) +
    ggh4x::facet_wrap2(.~LC, strip = facet.cols, ncol = 9, labeller = label_wrap_gen(width = 15)) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 17),
          axis.ticks.x = element_blank(),
          strip.text = element_text(size = 15, face = 'bold'),
          legend.position = 'bottom',
          legend.title = element_blank(),
          legend.text = element_text(size = 17))
  
  lc.areas.plot
}

pct.change.boxplot = pct.change.boxplot.fun()

ggsave(filename = paste(figs.fp, 'Percent_AGB_Change_LC_2003vs2012.png', sep = '/'),
       pct.change.boxplot,
       device = 'png',
       width = 13, height = 5, units = 'in')