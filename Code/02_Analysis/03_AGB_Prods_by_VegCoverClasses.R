# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Evaluate annual domain-averaged statistics for each AGB product
# Evaluate spatial distributions and temporal trends of AGB for each product
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
library(cowplot)
library(viridis)
library(pracma)
library(rasterVis)

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

# Get latitudinal pixel area adjustment function
source(paste(gitrepo.fp, 'Code', '02_Analysis', '00_Latitude_Adjustment_Fun.R', sep = '/'))

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
wang.mask = mean(wang.agb, na.rm = TRUE)
wang.mask[wang.mask != 0] = 1
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
# AGB_VM 
# ---------------------------------------------------------------------------------------------
AGB_VM.fp = paste(gitrepo.fp, 'Data/AGB_VM', sep = '/')
setwd(AGB_VM.fp)
AGB_VM = stack(list.files(AGB_VM.fp, pattern = '.tif', recursive = FALSE, full.names = FALSE)[1:15])
AGB_VM = brick(AGB_VM)

# ---------------------------------------------------------------------------------------------
# ESA land cover
# ---------------------------------------------------------------------------------------------
lc.fp = paste(data.fp, 'ESA_CCI_Landcover', sep = '/')
setwd(lc.fp)
esa.lc = raster('2010_ESA_CCI_LC_025deg.tif')

# Reclassify and plot land use/land cover
esa.lc.rc = esa.lc
# Omit classes that only have a few pixels (30,80,201), are innundated by water (160,180), or urban (190) 
esa.lc.rc[esa.lc.rc %in% c(30,80,160,180,190,201)] = NA
# Combine open/closed deciduous broadleaved (60,61) and open/closed neeleleaved evergreen (70,71) into 2 classes instead of 4
new.class.vals = c(1,2,2,3,3,4,5,6,7,8,9,10,11) 
lookup = data.frame(sort(unique(esa.lc.rc[])), new.class.vals)
esa.lc.rc = reclassify(esa.lc.rc, lookup)
esa.lc.df = as.data.frame(as(esa.lc.rc, 'SpatialPixelsDataFrame'))
colnames(esa.lc.df) = c('ClassID','x','y')
esa.lc.df = esa.lc.df[order(esa.lc.df$ClassID, decreasing = FALSE),]

class.names = c('Herbaceous cropland',          #11 to #1
                'Deciduous broad leaf forest',  #60+#61 to #2
                'Evergreen needle leaf forest', #70+#71 to #3
                'Mixed forest',                 #90 to #4
                'Mixed forest/shrub',           #100 to #5
                'Shrub',                        #120 to #6
                'Grass',                        #130 to #7
                'Lichens and mosses',           #140 to #8
                'Sparse',                       #150 to #9
                'Bare',                         #200 to #10
                'Snow/Ice')                     #220 to #11

class.names.abbr = c('CLF',
                     'DBF',
                     'ENF',
                     'MF',
                     'MF_CSL',
                     'CSL',
                     'GRL',
                     'TND',
                     'SPC',
                     'BARE',
                     'SNW_ICE')

# Create function to assign class names and abbreviations to reclassified ESA LULC
populate.class.names = function(){
  
  df.out = esa.lc.df
  df.out$ClassName = NA
  df.out$ClassABBR = NA
  
  esa.class.ids = unique(df.out$ClassID)
  
  # Assign class names
  for (i in 1 : length(esa.class.ids)){
    df.out[df.out$ClassID == esa.class.ids[i],4] = class.names[i]
  }
  
  # Assign class name abbreviations
  df.out[df.out$ClassID == 1,5] = class.names.abbr[1]
  df.out[df.out$ClassID == 2,5] = class.names.abbr[2]
  df.out[df.out$ClassID == 3,5] = class.names.abbr[3]
  df.out[df.out$ClassID == 4,5] = class.names.abbr[4]
  df.out[df.out$ClassID == 5,5] = class.names.abbr[5]
  df.out[df.out$ClassID == 6,5] = class.names.abbr[6]
  df.out[df.out$ClassID == 7,5] = class.names.abbr[7]
  df.out[df.out$ClassID == 8,5] = class.names.abbr[8]
  df.out[df.out$ClassID == 9,5] = class.names.abbr[9]
  df.out[df.out$ClassID == 10,5] = class.names.abbr[10]
  df.out[df.out$ClassID == 11,5] = class.names.abbr[11]
  
  return(df.out)
}

esa.lc.df = populate.class.names()
class.cols = pals::tableau20(n = length(class.names))
names(class.cols) = class.names.abbr

# Calculate percent coverage of each land cover class
esa.lc.df.summary.fun = function() {
  
  n = dim(esa.lc.df)[1]
  
  clf = esa.lc.df[esa.lc.df$ClassID == 1,]
  dbf = esa.lc.df[esa.lc.df$ClassID == 2,]
  enf = esa.lc.df[esa.lc.df$ClassID == 3,]
  mf = esa.lc.df[esa.lc.df$ClassID == 4,]
  mf.csl = esa.lc.df[esa.lc.df$ClassID == 5,]
  csl = esa.lc.df[esa.lc.df$ClassID == 6,]
  grl = esa.lc.df[esa.lc.df$ClassID == 7,]
  tnd = esa.lc.df[esa.lc.df$ClassID == 8,]
  spc = esa.lc.df[esa.lc.df$ClassID == 9,]
  bare = esa.lc.df[esa.lc.df$ClassID == 10,]
  snow.ice = esa.lc.df[esa.lc.df$ClassID == 11,]
  
  clf.pct = (dim(clf)[1] / n) * 100
  dbf.pct = (dim(dbf)[1] / n) * 100
  enf.pct = (dim(enf)[1] / n) * 100
  mf.pct = (dim(mf)[1] / n) * 100
  mf.csl.pct = (dim(mf.csl)[1] / n) * 100
  csl.pct = (dim(csl)[1] / n) * 100
  grl.pct = (dim(grl)[1] / n) * 100
  tnd.pct = (dim(tnd)[1] / n) * 100
  spc.pct = (dim(spc)[1] / n) * 100
  bare.pct = (dim(bare)[1] / n) * 100
  snow.ice.pct = (dim(snow.ice)[1] / n) * 100
  
  out.df = data.frame(ClassID = sort(unique(esa.lc.df$ClassID)),
                      PctArea = c(clf.pct, dbf.pct, enf.pct, mf.pct, mf.csl.pct, csl.pct, 
                                  grl.pct, tnd.pct, spc.pct, bare.pct, snow.ice.pct),
                      ClassName = unique(esa.lc.df$ClassName))
  
  out.df[,2] = round(out.df[,2],2)
  
  out.df = out.df %>%
    arrange(desc(ClassName)) %>%
    mutate(ypos = cumsum(out.df$PctArea) - (0.5 * out.df$PctArea))
  
  return(out.df)
}

esa.lc.df.pct = esa.lc.df.summary.fun()

# ------------------------------------ Land cover plots
# Map
esa.lc.df.plot = ggplot(data = esa.lc.df) +
  geom_raster(aes(x = x,
                  y = y,
                  fill = ClassABBR)) +
  scale_fill_manual(name = '',
                    values = class.cols) +
  #geom_sf(data = st_as_sf(entire.domain), color = NA, fill = NA, inherit.aes = FALSE, show.legend = FALSE) +
  geom_sf(data = wang.footprint, color = 'black', fill = NA, inherit.aes = FALSE, linewidth = 1) +
  xlab('') +
  ylab('') +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  #ylim(c(49.125,74.875)) +
  ggtitle('ESA CCI Veg. Cover Classes 2010') +
  theme_bw() +
  # annotation_custom(grobTree(textGrob('b', x = 0.01, y = 0.95, hjust = 0,
  #                                     gp = gpar(col = 'black', fontsize = 18, fontface = 'bold')))) +
  guides(fill = guide_legend(nrow = 3)) +
  theme(plot.title = element_text(size = 22),
        legend.position = 'bottom',
        legend.box.spacing = unit(-0.5, 'cm'),
        legend.text = element_text(size = 14),
        aspect.ratio = 0.6)

esa.lc.df.plot

# Pie chart 
esa.lc.df.piechart = ggplot(data = esa.lc.df.pct,
                            aes(x = '', y = PctArea, fill = ClassName)) +
  geom_col() +
  coord_polar('y') +
  theme_void() +
  geom_text(aes(x = 1.6, label = paste0(PctArea,'%')), 
            position = position_stack(vjust = 0.5), size = 5, fontface = 'bold') +
  scale_fill_manual(name = '', values = unname(class.cols)) +
  theme(legend.position = 'none')

esa.lc.df.piechart

# Annual AGB timeseries plot ---------------------------------
annual.mean.agb.fun = function(){
  
  # Apply latitudinal area adjustment function to AGB datsets
  
  pixel.area.lat.wang = lat.mean.stats.fun(wang.mask)
  
  wang = (wang.agb * pixel.area.lat.wang) / cellStats(mask(pixel.area.lat.wang, wang.agb), stat = sum)
  xu = (mask(xu.agb, wang.mask) * pixel.area.lat.wang) / cellStats(mask(pixel.area.lat.wang, xu.agb), stat = sum)
  liu = (mask(liu.agb, wang.mask) * pixel.area.lat.wang) / cellStats(mask(pixel.area.lat.wang, liu.agb), stat = sum)
  agb.vm = (mask(AGB_VM, wang.mask) * pixel.area.lat.wang) / cellStats(mask(pixel.area.lat.wang, AGB_VM), stat = sum)
  #esa.18 = (mask(esa.agb.18, wang.mask) * pixel.area.lat.wang) / cellStats(mask(pixel.area.lat.wang, esa.agb.18), stat = sum)
  esa.17 = (mask(esa.agb.17, wang.mask) * pixel.area.lat.wang) / cellStats(mask(pixel.area.lat.wang, esa.agb.17), stat = sum)
  esa.10 = (mask(esa.agb.10, wang.mask) * pixel.area.lat.wang) / cellStats(mask(pixel.area.lat.wang, esa.agb.10), stat = sum)
  
  # Sum the area-averaged annual biomass values
  wang.ym = data.frame('Year' = years,
                       'Annual_AGB' = c(unname(cellStats(wang, stat = sum)),NA,NA,NA),
                       'Product' = 'Wang')
  
  xu.ym = data.frame('Year' = years,
                     'Annual_AGB' = unname(cellStats(xu, stat = sum)),
                     'Product' = 'Xu')
  
  liu.ym = data.frame('Year' = years,
                      'Annual_AGB' = c(unname(cellStats(liu, stat = sum)),NA,NA,NA,NA,NA),
                      'Product' = 'Liu')
  
  agb.vm.ym = data.frame('Year' = years,
                         'Annual_AGB' = unname(cellStats(agb.vm, stat = sum)),
                         'Product' = 'AGB_VM')
  
  esa.ym = data.frame('Year' = years,
                      'Annual_AGB' = c(rep(NA,7),
                                       sum(esa.10[], na.rm = TRUE),
                                       rep(NA,6),
                                       sum(esa.17[], na.rm = TRUE)),
                      'Product' = 'ESA_CCI')
  
  prods.agb.ym = rbind(wang.ym, xu.ym, liu.ym, agb.vm.ym, esa.ym)
  prods.agb.ym$Product = factor(prods.agb.ym$Product, levels = c('AGB_VM','Wang','Xu','Liu','ESA_CCI'))
  
  return(prods.agb.ym)
}

agb.prods.ym = annual.mean.agb.fun()

# ---- Plot timeseries

# Set graphical parameters
xu.agb.col = 'black'
wang.agb.col = 'tomato'
AGB_VM.col = 'green'
liu.agb.col = 'steelblue'
esa.agb.col = 'magenta'

xu.agb.pch = 18
wang.agb.pch = 15
AGB_VM.pch = 17
liu.agb.pch = 16
esa.agb.pch = 19

linetypes = factor(c(1,3,2,1,0))

# Get linear models for annual domain-averaged AGB values for each product 2003-2017
AGB_VM.ts.lm = lm(agb.prods.ym[agb.prods.ym$Product == 'AGB_VM',2] ~ seq(2003,2017,1))
wang.ts.lm = lm(agb.prods.ym[agb.prods.ym$Product == 'Wang',2] ~ seq(2003,2017,1))
xu.ts.lm = lm(agb.prods.ym[agb.prods.ym$Product == 'Xu',2] ~ seq(2003,2017,1))
liu.ts.lm = lm(agb.prods.ym[agb.prods.ym$Product == 'Liu',2] ~ seq(2003,2017,1))
esa.ts.lm = lm(agb.prods.ym[agb.prods.ym$Product == 'ESA_CCI',2] ~ seq(2003,2017,1))

# Create plot
annual.mean.agb.ts = ggplot(data = agb.prods.ym,
            aes(x = Year,
                y = Annual_AGB,
                color = Product, 
                linetype = Product, 
                shape = Product)) +
  scale_y_continuous(limits = c(25,70)) +
  scale_x_continuous(breaks = years, limits = c(2003,2017)) +
  ylab('Mean annual AGB [Mg/ha]') +
  theme_bw() +
  # annotation_custom(grobTree(textGrob('c', x = 0.01, y = 0.95, hjust = 0,
  #                                     gp = gpar(col = 'black', fontsize = 18, fontface = 'bold')))) +
  geom_line(linewidth = 1) +
  geom_point(size = 4) +
  scale_linetype_manual(values = linetypes,
                        guide = guide_legend(override.aes = list(linetype = c(1,3,2,1,NA)))) + 
  scale_shape_manual(values = c(AGB_VM.pch, wang.agb.pch, xu.agb.pch, liu.agb.pch, esa.agb.pch)) +
  scale_color_manual(values = c(AGB_VM.col, wang.agb.col, xu.agb.col, liu.agb.col, esa.agb.col)) +
  geom_abline(slope = coef(AGB_VM.ts.lm)[[2]], intercept = coef(AGB_VM.ts.lm)[[1]],
              col = AGB_VM.col, alpha = 0.6) +
  geom_abline(slope = coef(wang.ts.lm)[[2]], intercept = coef(wang.ts.lm)[[1]],
              col = wang.agb.col, alpha = 0.6) +
  geom_abline(slope = coef(xu.ts.lm)[[2]], intercept = coef(xu.ts.lm)[[1]],
              col = xu.agb.col, alpha = 0.6) +
  geom_abline(slope = coef(liu.ts.lm)[[2]], intercept = coef(liu.ts.lm)[[1]],
              col = liu.agb.col, alpha = 0.6) +
  # Annotate pvals for Xu
  annotation_custom(grobTree(textGrob(paste('p =',round(summary(xu.ts.lm)$coefficients[,4][2],3)), 
                                      x = 0.8, y = 0.9, hjust = 0,
                             gp = gpar(col = xu.agb.col, fontsize = 18, fontface = 'italic')))) +
  # Annotate slope for Xu
  annotation_custom(grobTree(textGrob(paste('slope =',round(xu.ts.lm$coefficients[2],3)), 
                                      x = 0.8, y = 0.85, hjust = 0,
                                      gp = gpar(col = xu.agb.col, fontsize = 18, fontface = 'italic')))) +
  # Annotate pvals for Liu
  annotation_custom(grobTree(textGrob(paste('p =',round(summary(liu.ts.lm)$coefficients[,4][2],3)), 
                                      x = 0.02, y = 0.47, hjust = 0,
                                      gp = gpar(col = liu.agb.col, fontsize = 18, fontface = 'italic')))) +
  # Annotate slope for Liu
  annotation_custom(grobTree(textGrob(paste('slope =',round(liu.ts.lm$coefficients[2],3)), 
                                      x = 0.02, y = 0.42, hjust = 0,
                                      gp = gpar(col = liu.agb.col, fontsize = 18, fontface = 'italic')))) +
  # Annotate pvals for AGB_VM
  annotation_custom(grobTree(textGrob(paste('p =',round(summary(AGB_VM.ts.lm)$coefficients[,4][2],3)), 
                                      x = 0.8, y = 0.5, hjust = 0,
                                      gp = gpar(col = AGB_VM.col, fontsize = 18, fontface = 'italic')))) +
  # Annotate slope for AGB_VM
  annotation_custom(grobTree(textGrob(paste('slope =',round(AGB_VM.ts.lm$coefficients[2],3)), 
                                      x = 0.8, y = 0.45, hjust = 0,
                                      gp = gpar(col = AGB_VM.col, fontsize = 18, fontface = 'italic')))) +
  # Annotate pvals for Wang
  annotation_custom(grobTree(textGrob(paste('p =',round(summary(wang.ts.lm)$coefficients[,4][2],3)), 
                                      x = 0.02, y = 0.3, hjust = 0,
                                      gp = gpar(col = wang.agb.col, fontsize = 18, fontface = 'italic')))) +
  # Annotate slope for Wang
  annotation_custom(grobTree(textGrob(paste('slope =',round(wang.ts.lm$coefficients[2],3)), 
                                      x = 0.02, y = 0.25, hjust = 0,
                                      gp = gpar(col = wang.agb.col, fontsize = 18, fontface = 'italic')))) +
  
  theme(axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 17),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(size = 13, face = 'bold'),
        #legend.position = c(0.8,0.9),
        legend.position = 'right',
        legend.direction = 'vertical',
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 17, face = 'bold'))


annual.mean.agb.ts


# 2010 AGB plot ------------------------------------
AGB_VM.10.df = as.data.frame(as(AGB_VM[[8]], 'SpatialPixelsDataFrame'))
colnames(AGB_VM.10.df) = c('AGB','x','y')
# esa.agb.17.df = as.data.frame(as(esa.agb.17, 'SpatialPixelsDataFrame'))
# colnames(esa.agb.17.df) = c('AGB','x','y')

agb.map = ggplot(data = AGB_VM.10.df,
                 aes(x = x, y = y, fill = AGB)) +
  geom_raster() +
  scale_fill_stepsn(colors = c('palegoldenrod','palegreen1','palegreen4','darkgreen'),
                    breaks = seq(0,300,50),
                    limits = c(0,300),
                    guide = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
  # geom_sf(data = st_as_sf(entire.domain), color = 'black', fill = NA, inherit.aes = FALSE) +
  # geom_sf(data = st_as_sf(core.domain), color = 'tomato', fill = NA, inherit.aes = FALSE) +
  geom_sf(data = wang.footprint, color = 'black', fill = NA, inherit.aes = FALSE, linewidth = 1) +
  labs(fill = '') +
  xlab('') +
  ylab('') +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle('AGB_VM [Mg/ha] 2010') +
  theme_bw() +
  # annotation_custom(grobTree(textGrob('a', x = 0.01, y = 0.95, hjust = 0,
  #                                     gp = gpar(col = 'black', fontsize = 18, fontface = 'bold')))) +
  theme(plot.title = element_text(size = 20),
        legend.position = 'bottom',
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(2.75, 'cm'),
        legend.box.spacing = unit(0.25, 'cm'),
        legend.box.just = 'top',
        legend.text = element_text(size = 15),
        aspect.ratio = 0.6)
  
  
agb.map


# Combine plots ------------------------------------
agb.lc.fig = plot_grid(plot_grid(plotlist = list(agb.map, esa.lc.df.plot, esa.lc.df.piechart),
                                 axis = 'llt',
                                 ncol = 3,
                                 rel_widths = c(1,1,0.75),
                                 labels = c('a','b',''),
                                 label_size = 20,
                                 vjust = 1),
                       annual.mean.agb.ts,
                       nrow = 2,
                       labels = c('','c'),
                       label_size = 20,
                       vjust = 2,
                       rel_heights = c(1.5,1.25))


ggsave(paste(figs.fp,'ESA_AGB_LC_Maps_TS.png', sep = '/'),
       agb.lc.fig,
       device = 'png',
       width = 18, height = 9, units = 'in', bg = 'white')

# -------------------------------------------------------------------------------------------------
# Compute annual AGB statistics for each product for total area (Wang)
# -------------------------------------------------------------------------------------------------

agb.prods.ytot = rbind(data.frame('Year' = years,
                                  'Annual_AGB' = c(unname(cellStats(wang.agb, stat = sum)),NA,NA,NA),
                                  'Product' = 'Wang'),
                       data.frame('Year' = years,
                                  'Annual_AGB' = unname(cellStats(mask(xu.agb, wang.mask), stat = sum)),
                                  'Product' = 'Xu'),
                       data.frame('Year' = years,
                                  'Annual_AGB' = c(unname(cellStats(mask(liu.agb, wang.mask), stat = sum)),NA,NA,NA,NA,NA),
                                  'Product' = 'Liu'),
                       data.frame('Year' = years,
                                  'Annual_AGB' = unname(cellStats(mask(AGB_VM, wang.mask), stat = sum)),
                                  'Product' = 'AGB_VM'),
                       data.frame('Year' = years,
                                  'Annual_AGB' = c(NA,NA,NA,NA,NA,NA,NA,unname(cellStats(mask(esa.agb.10, wang.mask), stat = sum)),NA,NA,NA,NA,NA,NA,unname(cellStats(esa.agb.17, stat = sum))),
                                  'Product' = 'ESA_CCI'))

agb.prods.ystd = rbind(data.frame('Year' = years,
                                  'Annual_AGB' = c(unname(cellStats(wang.agb, stat = sd)),NA,NA,NA),
                                  'Product' = 'Wang'),
                       data.frame('Year' = years,
                                  'Annual_AGB' = unname(cellStats(mask(xu.agb, wang.mask), stat = sd)),
                                  'Product' = 'Xu'),
                       data.frame('Year' = years,
                                  'Annual_AGB' = c(unname(cellStats(mask(liu.agb, wang.mask), stat = sd)),NA,NA,NA,NA,NA),
                                  'Product' = 'Liu'),
                       data.frame('Year' = years,
                                  'Annual_AGB' = unname(cellStats(mask(AGB_VM, wang.mask), stat = sd)),
                                  'Product' = 'AGB_VM'),
                       data.frame('Year' = years,
                                  'Annual_AGB' = c(NA,NA,NA,NA,NA,NA,NA,unname(cellStats(mask(esa.agb.10, wang.mask), stat = sd)),NA,NA,NA,NA,NA,NA,unname(cellStats(esa.agb.17, stat = sd))),
                                  'Product' = 'ESA_CCI'))

agb.prods.ym$Stat = 'Mean'
agb.prods.ytot$Stat = 'Tot.'
agb.prods.ystd$Stat = 'STD'

agb.prods.annual.stats = rbind(agb.prods.ym, agb.prods.ytot, agb.prods.ystd)
agb.prods.annual.stats$Annual_AGB = round(agb.prods.annual.stats$Annual_AGB, 2)

agb.stats.mean = agb.prods.annual.stats[agb.prods.annual.stats$Stat == 'Mean',]
agb.stats.tot = agb.prods.annual.stats[agb.prods.annual.stats$Stat == 'Tot.',]
agb.stats.std = agb.prods.annual.stats[agb.prods.annual.stats$Stat == 'STD',]

agb.stats.mean = reshape(agb.stats.mean, idvar = 'Product', timevar = 'Year', v.names = 'Annual_AGB', direction = 'wide')
colnames(agb.stats.mean) = c('Product','Stat',years)

agb.stats.tot = reshape(agb.stats.tot, idvar = 'Product', timevar = 'Year', v.names = 'Annual_AGB', direction = 'wide')
colnames(agb.stats.tot) = c('Product','Stat',years)

agb.stats.std = reshape(agb.stats.std, idvar = 'Product', timevar = 'Year', v.names = 'Annual_AGB', direction = 'wide')
colnames(agb.stats.std) = c('Product','Stat',years)

# 2003
agb.stats.2003 = cbind(agb.stats.tot[,c(1,3)],
                       agb.stats.mean[,3],
                       agb.stats.std[,3])
colnames(agb.stats.2003) = c('Product','Tot.','Mean','STD')
agb.stats.2003$Year = '2003'
# 2004
agb.stats.2004 = cbind(agb.stats.tot[,c(1,4)],
                       agb.stats.mean[,4],
                       agb.stats.std[,4])
colnames(agb.stats.2004) = c('Product','Tot.','Mean','STD')
agb.stats.2004$Year = '2004'
# 2005
agb.stats.2005 = cbind(agb.stats.tot[,c(1,5)],
                       agb.stats.mean[,5],
                       agb.stats.std[,5])
colnames(agb.stats.2005) = c('Product','Tot.','Mean','STD')
agb.stats.2005$Year = '2005'
# 2006
agb.stats.2006 = cbind(agb.stats.tot[,c(1,6)],
                       agb.stats.mean[,6],
                       agb.stats.std[,6])
colnames(agb.stats.2006) = c('Product','Tot.','Mean','STD')
agb.stats.2006$Year = '2006'
# 2007
agb.stats.2007 = cbind(agb.stats.tot[,c(1,7)],
                       agb.stats.mean[,7],
                       agb.stats.std[,7])
colnames(agb.stats.2007) = c('Product','Tot.','Mean','STD')
agb.stats.2007$Year = '2007'
# 2008
agb.stats.2008 = cbind(agb.stats.tot[,c(1,8)],
                       agb.stats.mean[,8],
                       agb.stats.std[,8])
colnames(agb.stats.2008) = c('Product','Tot.','Mean','STD')
agb.stats.2008$Year = '2008'
# 2009
agb.stats.2009 = cbind(agb.stats.tot[,c(1,9)],
                       agb.stats.mean[,9],
                       agb.stats.std[,9])
colnames(agb.stats.2009) = c('Product','Tot.','Mean','STD')
agb.stats.2009$Year = '2009'
# 2010
agb.stats.2010 = cbind(agb.stats.tot[,c(1,10)],
                       agb.stats.mean[,10],
                       agb.stats.std[,10])
colnames(agb.stats.2010) = c('Product','Tot.','Mean','STD')
agb.stats.2010$Year = '2010'
# 2011
agb.stats.2011 = cbind(agb.stats.tot[,c(1,11)],
                       agb.stats.mean[,11],
                       agb.stats.std[,11])
colnames(agb.stats.2011) = c('Product','Tot.','Mean','STD')
agb.stats.2011$Year = '2011'
# 2012
agb.stats.2012 = cbind(agb.stats.tot[,c(1,12)],
                       agb.stats.mean[,12],
                       agb.stats.std[,12])
colnames(agb.stats.2012) = c('Product','Tot.','Mean','STD')
agb.stats.2012$Year = '2012'
# 2013
agb.stats.2013 = cbind(agb.stats.tot[,c(1,13)],
                       agb.stats.mean[,13],
                       agb.stats.std[,13])
colnames(agb.stats.2013) = c('Product','Tot.','Mean','STD')
agb.stats.2013$Year = '2013'
# 2014
agb.stats.2014 = cbind(agb.stats.tot[,c(1,14)],
                       agb.stats.mean[,14],
                       agb.stats.std[,14])
colnames(agb.stats.2014) = c('Product','Tot.','Mean','STD')
agb.stats.2014$Year = '2014'
# 2015
agb.stats.2015 = cbind(agb.stats.tot[,c(1,15)],
                       agb.stats.mean[,15],
                       agb.stats.std[,15])
colnames(agb.stats.2015) = c('Product','Tot.','Mean','STD')
agb.stats.2015$Year = '2015'
# 2016
agb.stats.2016 = cbind(agb.stats.tot[,c(1,16)],
                       agb.stats.mean[,16],
                       agb.stats.std[,16])
colnames(agb.stats.2016) = c('Product','Tot.','Mean','STD')
agb.stats.2016$Year = '2016'
# 2017
agb.stats.2017 = cbind(agb.stats.tot[,c(1,17)],
                       agb.stats.mean[,17],
                       agb.stats.std[,17])
colnames(agb.stats.2017) = c('Product','Tot.','Mean','STD')
agb.stats.2017$Year = '2017'

agb.stats = rbind(agb.stats.2003, agb.stats.2004, agb.stats.2005, agb.stats.2006, agb.stats.2007,
                  agb.stats.2008, agb.stats.2009, agb.stats.2010, agb.stats.2011, agb.stats.2012, 
                  agb.stats.2013, agb.stats.2014, agb.stats.2015, agb.stats.2016, agb.stats.2017)

agb.stats = agb.stats[,c(5,1,2,3,4)]

write.csv(agb.stats,
          file = paste(tables.fp, 'Table_S4_AGB_Product_Stats_2003to2017.csv', sep = '/'),
          row.names = FALSE)

# Create tables
agb.stats.table.03.10 = flextable(agb.stats[agb.stats$Year %in% as.character(seq(2003,2010,1)),]) %>%
  theme_alafoli() %>%
  colformat_num(na_str = '-') %>%
  bold(part = 'header') %>%
  fontsize(size = 12, part = 'body') %>%
  fontsize(size = 14, part = 'header') %>%
  merge_v(j = 'Year') %>%
  valign(j = 1, valign = 'top') %>%
  align_text_col(align = 'center') %>%
  align_nottext_col(align = 'center') %>%
  set_table_properties(layout = 'fixed') %>%
  #add_header_lines('Annual AGB statistics [Mg/ha]') %>%
  bg(bg = 'white', part = 'all') %>%
  vline(j = 1) %>%
  gen_grob()

agb.stats.table.11.17 = flextable(agb.stats[agb.stats$Year %in% as.character(seq(2011,2017,1)),]) %>%
  theme_alafoli() %>%
  colformat_num(na_str = '-') %>%
  bold(part = 'header') %>%
  fontsize(size = 12, part = 'body') %>%
  fontsize(size = 14, part = 'header') %>%
  merge_v(j = 'Year') %>%
  valign(j = 1, valign = 'top') %>%
  align_text_col(align = 'center') %>%
  align_nottext_col(align = 'center') %>%
  set_table_properties(layout = 'fixed') %>%
  #add_header_lines('Annual AGB statistics [Mg/ha]') %>%
  bg(bg = 'white', part = 'all') %>%
  vline(j = 1) %>%
  gen_grob()


grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.5, width = 0.5, just = 'right'))
grid.draw(agb.stats.table.03.10)
pushViewport(viewport(x = 1.47, y = 0.13, width = 1, height = 0.87, just = 'bottom'))
grid.draw(agb.stats.table.11.17)

# -------------------------------------------------------------------------------------------------
# Compute annual AGB statistics for each product by land cover classes for year 2010
# -------------------------------------------------------------------------------------------------

ras.obj = landmask

esa.lc.ras = rasterize(x = esa.lc.df[,2:3],
                       y = ras.obj,
                       field = esa.lc.df[,1])

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

agb.lc.mask.fun = function(agb.in, product, year) {
  
  # Mask input AGB to extent of Wang dataset
  agb = mask(agb.in, wang.mask)
  
  # Mask Wang-masked AGB to extent of each cover class
  agb.clf = mask(agb, clf)
  agb.dbf = mask(agb, dbf)
  agb.enf = mask(agb, enf)
  agb.mf = mask(agb, mf)
  agb.mf.csl = mask(agb, mf.csl)
  agb.csl = mask(agb, csl)
  agb.grl = mask(agb, grl)
  agb.tnd = mask(agb, tnd)
  agb.spc = mask(agb, spc)
  
  agb.lc.df = rbind(data.frame(AGB = na.omit(values(agb.clf)),
                               LC = 'Herbaceous cropland'),
                    data.frame(AGB = na.omit(values(agb.dbf)),
                               LC = 'Deciduous broad leaf forest'),
                    data.frame(AGB = na.omit(values(agb.enf)),
                               LC = 'Evergreen needle leaf forest'),
                    data.frame(AGB = na.omit(values(agb.mf)),
                               LC = 'Mixed forest'),
                    data.frame(AGB = na.omit(values(agb.mf.csl)),
                               LC = 'Mixed forest/shrub'),
                    data.frame(AGB = na.omit(values(agb.csl)),
                               LC = 'Shrub'),
                    data.frame(AGB = na.omit(values(agb.grl)),
                               LC = 'Grass'),
                    data.frame(AGB = na.omit(values(agb.tnd)),
                               LC = 'Lichens and mosses'),
                    data.frame(AGB = na.omit(values(agb.spc)),
                               LC = 'Sparse'))
  
  agb.lc.df$Product = product
  agb.lc.df$Year = year
  return(agb.lc.df)
}

AGB_VM.agb.lc.df = agb.lc.mask.fun(AGB_VM[[8]], 'AGB_VM', '2010')
wang.agb.lc.df = agb.lc.mask.fun(wang.agb[[8]], 'Wang', '2010')
xu.agb.lc.df = agb.lc.mask.fun(xu.agb[[8]], 'Xu', '2010')
liu.agb.lc.df = agb.lc.mask.fun(liu.agb[[8]], 'Liu', '2010')
esa.agb.lc.df = agb.lc.mask.fun(esa.agb.10, 'ESA_CCI', '2010')

agb.lc.df = rbind(AGB_VM.agb.lc.df, wang.agb.lc.df, xu.agb.lc.df, liu.agb.lc.df, esa.agb.lc.df)

lc.levels = c('Herbaceous cropland','Deciduous broad leaf forest','Evergreen needle leaf forest',
              'Mixed forest','Mixed forest/shrub','Shrub','Grass','Lichens and mosses','Sparse')

prod.levels = c('AGB_VM','Wang','Xu','Liu','ESA_CCI')

agb.lc.df$LC = factor(agb.lc.df$LC, levels = lc.levels)
agb.lc.df$Product = factor(agb.lc.df$Product, levels = prod.levels)

# Set facet colors
facet.cols = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = unname(class.cols)[1:9]))

# Create box-and-whisker boxplot for 2010 AGB product values per vegetation cover class
agb.lc.boxplot = ggplot(data = agb.lc.df, 
                        aes(x = Product, y = AGB, color = Product, fill = Product)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = 'crossbar', color = 'white', fatten = 0.5) +
  scale_color_manual(labels = c('AGB_VM','Wang','Xu','Liu','ESA_CCI'), 
                     values = c(AGB_VM.col,wang.agb.col,xu.agb.col,liu.agb.col,esa.agb.col)) +
  scale_fill_manual(labels = c('AGB_VM','Wang','Xu','Liu','ESA_CCI'), 
                    values = c(AGB_VM.col,wang.agb.col,xu.agb.col,liu.agb.col,esa.agb.col)) +
  ylab('AGB [Mg/ha]') +
  ggh4x::facet_wrap2(.~LC, strip = facet.cols, scales = 'free') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        strip.text = element_text(size = 15, face = 'bold'))

agb.lc.boxplot

ggsave(filename = paste(figs.fp, 'AGB_LULC_Boxplot_2010.png', sep = '/'),
       agb.lc.boxplot,
       width = 11, height = 9)


# Get annual AGB mean and standard deviation per vegetation cover class for each product
# AGB_VM
AGB_VM.agb.lc.df.annual = rbind(agb.lc.mask.fun(AGB_VM[[1]], 'AGB_VM', '2003'),
                                agb.lc.mask.fun(AGB_VM[[2]], 'AGB_VM', '2004'),
                                agb.lc.mask.fun(AGB_VM[[3]], 'AGB_VM', '2005'),
                                agb.lc.mask.fun(AGB_VM[[4]], 'AGB_VM', '2006'),
                                agb.lc.mask.fun(AGB_VM[[5]], 'AGB_VM', '2007'),
                                agb.lc.mask.fun(AGB_VM[[6]], 'AGB_VM', '2008'),
                                agb.lc.mask.fun(AGB_VM[[7]], 'AGB_VM', '2009'),
                                agb.lc.mask.fun(AGB_VM[[8]], 'AGB_VM', '2010'),
                                agb.lc.mask.fun(AGB_VM[[9]], 'AGB_VM', '2011'),
                                agb.lc.mask.fun(AGB_VM[[10]], 'AGB_VM', '2012')) %>%
  group_by(Year, LC) %>%
  summarize(AGB_median = median(AGB),
            AGB_sd = sd(AGB),
            Product = 'AGB_VM') %>%
  as.data.frame()

# Wang
wang.agb.lc.df.annual = rbind(agb.lc.mask.fun(wang.agb[[1]], 'Wang', '2003'),
                              agb.lc.mask.fun(wang.agb[[2]], 'Wang', '2004'),
                              agb.lc.mask.fun(wang.agb[[3]], 'Wang', '2005'),
                              agb.lc.mask.fun(wang.agb[[4]], 'Wang', '2006'),
                              agb.lc.mask.fun(wang.agb[[5]], 'Wang', '2007'),
                              agb.lc.mask.fun(wang.agb[[6]], 'Wang', '2008'),
                              agb.lc.mask.fun(wang.agb[[7]], 'Wang', '2009'),
                              agb.lc.mask.fun(wang.agb[[8]], 'Wang', '2010'),
                              agb.lc.mask.fun(wang.agb[[9]], 'Wang', '2011'),
                              agb.lc.mask.fun(wang.agb[[10]], 'Wang', '2012')) %>%
  group_by(Year, LC) %>%
  summarize(AGB_median = median(AGB),
            AGB_sd = sd(AGB),
            Product = 'Wang') %>%
  as.data.frame()

# Xu
xu.agb.lc.df.annual = rbind(agb.lc.mask.fun(xu.agb[[1]], 'Xu', '2003'),
                            agb.lc.mask.fun(xu.agb[[2]], 'Xu', '2004'),
                            agb.lc.mask.fun(xu.agb[[3]], 'Xu', '2005'),
                            agb.lc.mask.fun(xu.agb[[4]], 'Xu', '2006'),
                            agb.lc.mask.fun(xu.agb[[5]], 'Xu', '2007'),
                            agb.lc.mask.fun(xu.agb[[6]], 'Xu', '2008'),
                            agb.lc.mask.fun(xu.agb[[7]], 'Xu', '2009'),
                            agb.lc.mask.fun(xu.agb[[8]], 'Xu', '2010'),
                            agb.lc.mask.fun(xu.agb[[9]], 'Xu', '2011'),
                            agb.lc.mask.fun(xu.agb[[10]], 'Xu', '2012')) %>%
  group_by(Year, LC) %>%
  summarize(AGB_median = median(AGB),
            AGB_sd = sd(AGB),
            Product = 'Xu') %>%
  as.data.frame()

# Liu
liu.agb.lc.df.annual = rbind(agb.lc.mask.fun(liu.agb[[1]], 'Liu', '2003'),
                             agb.lc.mask.fun(liu.agb[[2]], 'Liu', '2004'),
                             agb.lc.mask.fun(liu.agb[[3]], 'Liu', '2005'),
                             agb.lc.mask.fun(liu.agb[[4]], 'Liu', '2006'),
                             agb.lc.mask.fun(liu.agb[[5]], 'Liu', '2007'),
                             agb.lc.mask.fun(liu.agb[[6]], 'Liu', '2008'),
                             agb.lc.mask.fun(liu.agb[[7]], 'Liu', '2009'),
                             agb.lc.mask.fun(liu.agb[[8]], 'Liu', '2010'),
                             agb.lc.mask.fun(liu.agb[[9]], 'Liu', '2011'),
                             agb.lc.mask.fun(liu.agb[[10]], 'Liu', '2012')) %>%
  group_by(Year, LC) %>%
  summarize(AGB_median = median(AGB),
            AGB_sd = sd(AGB),
            Product = 'Liu') %>%
  as.data.frame()

# ESA CCI
esa.agb.lc.df.annual = data.frame(Year = AGB_VM.agb.lc.df.annual$Year,
                                  LC = AGB_VM.agb.lc.df.annual$LC,
                                  AGB_median = numeric(length = dim(AGB_VM.agb.lc.df.annual)[1]),
                                  AGB_sd = numeric(length = dim(AGB_VM.agb.lc.df.annual)[1]),
                                  Product = 'ESA_CCI')

esa.agb.lc.df.summary = esa.agb.lc.df %>% 
  group_by(LC) %>% 
  summarize(AGB_median = median(AGB), 
            AGB_sd = sd(AGB)) %>% 
  as.data.frame()

esa.agb.lc.df.annual[esa.agb.lc.df.annual$Year == 2010, 3] = esa.agb.lc.df.summary$AGB_median
esa.agb.lc.df.annual[esa.agb.lc.df.annual$Year == 2010, 4] = esa.agb.lc.df.summary$AGB_sd
esa.agb.lc.df.annual[esa.agb.lc.df.annual$AGB_median == 0, 3] = NA
esa.agb.lc.df.annual[esa.agb.lc.df.annual$AGB_sd == 0, 4] = NA

# Combine into single data frame
agb.lc.df.annual = rbind(AGB_VM.agb.lc.df.annual,
                         wang.agb.lc.df.annual,
                         xu.agb.lc.df.annual,
                         liu.agb.lc.df.annual,
                         esa.agb.lc.df.annual)

# Factor
agb.lc.df.annual$Product = factor(agb.lc.df.annual$Product, levels = prod.levels)
agb.lc.df.annual$LC = factor(agb.lc.df.annual$LC, levels = lc.levels)

# Plot timeseries
agb.lc.df.annual.timeseries = ggplot(data = agb.lc.df.annual, 
                                     aes(x = Year, 
                                         y = AGB_median, 
                                         group = Product,
                                         color = Product,
                                         linetype = Product, 
                                         shape = Product)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  # geom_errorbar(aes(ymin = AGB_median - AGB_sd,
  #                   ymax = AGB_median + AGB_sd),
  #               position = position_dodge(0.05)) +
  scale_linetype_manual(values = linetypes,
                        guide = guide_legend(override.aes = list(linetype = c(1,3,2,1,NA)))) + 
  scale_shape_manual(values = c(AGB_VM.pch, wang.agb.pch, xu.agb.pch, liu.agb.pch, esa.agb.pch)) +
  scale_color_manual(labels = c('AGB_VM','Wang','Xu','Liu','ESA_CCI'), 
                     values = c(AGB_VM.col,wang.agb.col,xu.agb.col,liu.agb.col,esa.agb.col)) +
  ylab('Annual median AGB [Mg/ha]') +
  scale_x_discrete(breaks = unique(agb.lc.df.annual$Year),
                   labels = gsub('20','',unique(agb.lc.df.annual$Year))) +
  ggh4x::facet_wrap2(.~LC, strip = facet.cols, scales = 'free') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        strip.text = element_text(size = 15, face = 'bold'))

agb.lc.df.annual.timeseries

ggsave(filename = paste(figs.fp, 'AGB_LULC_Annual_2003to2012.png', sep = '/'),
       agb.lc.df.annual.timeseries,
       width = 11, height = 9)
