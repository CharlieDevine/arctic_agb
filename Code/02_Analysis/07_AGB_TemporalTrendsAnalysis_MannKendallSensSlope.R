# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Calculate pixel-wise temporal trend (Mann-Kendall Sen's Slope) for each AGB product
# and evaluate statistical distribution of slope trends within spatially aggregated
# vegetation cover classes.
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
library(cowplot)
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
library(xts)
library(Kendall)
library(wql)


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

# ----------------------
# Burn area table
# ----------------------
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


# -------------------------------------------------------------------------------------------------
# Compute pixel-wise Mann-Kendall trend statistics
# -------------------------------------------------------------------------------------------------

mannkendal.trend.fun = function(agb.in, agb.prod) {
  
  trend.years = as.Date(ISOdate(seq(2003,2012,1), 12, 31))
  
  # Blank raster for output Mann-Kendall p value
  mkp.trend.ras = landmask
  mkp.trend.ras[] = NA
  names(mkp.trend.ras) = paste(agb.prod,'MK_p_2003to2012', sep = '_')
  
  # Blank raster for output Mann-Kendall tau statistic
  mktau.trend.ras = landmask
  mktau.trend.ras[] = NA
  names(mktau.trend.ras) = paste(agb.prod,'MK_tau_2003to2012', sep = '_')
  
  # Blank raster for output Mann-Kendall Sen's slope
  mkss.trend.ras = landmask
  mkss.trend.ras[] = NA
  names(mkss.trend.ras) = paste(agb.prod,'MK_SS_2003to2012', sep = '_')
  
  dims = dim(mkp.trend.ras)
  
  for (i in 1 : dims[1]) {
    for (j in 1 : dims[2]) {
      
      agb.vec = as.vector(unname(agb.in[i,j,]))
      
      if (is.na(agb.vec[1])) {
        mkp.trend.ras[i,j] = NA
      }
      
      if (!is.na(agb.vec[1])) {
        # Get MK p and tau
        agb.vec.mk.1 = MannKendall(xts(agb.vec, trend.years))
        mkp.trend.ras[i,j] = agb.vec.mk.1$sl[1]
        mktau.trend.ras[i,j] = agb.vec.mk.1$tau[1]
        # Get MK sen's slope
        agb.vec.mk.2 = mannKen(agb.vec, plot = FALSE)
        mkss.trend.ras[i,j] = agb.vec.mk.2$sen.slope
      }
    }
  }
  # --- Convert maps to spatial pixels data frame
  # MK-p
  mkp.trend.df = as.data.frame(as(mkp.trend.ras, 'SpatialPixelsDataFrame'))
  colnames(mkp.trend.df) = c(agb.prod,'x','y')
  mkp.trend.df = melt(mkp.trend.df, id.vars = c('x','y'), value.name = 'MannKendall_p', variable = 'Product')
  # MK-tau
  mktau.trend.df = as.data.frame(as(mktau.trend.ras, 'SpatialPixelsDataFrame'))
  colnames(mktau.trend.df) = c(agb.prod,'x','y')
  mktau.trend.df = melt(mktau.trend.df, id.vars = c('x','y'), value.name = 'MannKendall_tau', variable = 'Product')
  # MK-ss
  mkss.trend.df = as.data.frame(as(mkss.trend.ras, 'SpatialPixelsDataFrame'))
  colnames(mkss.trend.df) = c(agb.prod,'x','y')
  mkss.trend.df = melt(mkss.trend.df, id.vars = c('x','y'), value.name = 'MannKendall_SS', variable = 'Product')
  
  return(list(mkp.trend.df, mktau.trend.df, mkss.trend.df, 
              mkp.trend.ras, mktau.trend.ras, mkss.trend.ras))
}

AGB_VM.mk = mannkendal.trend.fun(AGB_VM[[1:10]], 'AGB_VM')
wang.mk = mannkendal.trend.fun(wang.agb[[1:10]], 'Wang')
xu.mk = mannkendal.trend.fun(xu.agb[[1:10]], 'Xu')
liu.mk = mannkendal.trend.fun(liu.agb, 'Liu')

mkp.trend.df = na.omit(rbind(AGB_VM.mk[[1]],
                             wang.mk[[1]],
                             xu.mk[[1]],
                             liu.mk[[1]]))

mktau.trend.df = na.omit(rbind(AGB_VM.mk[[2]],
                               wang.mk[[2]],
                               xu.mk[[2]],
                               liu.mk[[2]]))

mkss.trend.df = na.omit(rbind(AGB_VM.mk[[3]],
                              wang.mk[[3]],
                              xu.mk[[3]],
                              liu.mk[[3]]))

mkss.trend.df[mkss.trend.df$MannKendall_SS == 0,4] = NA

# -----------------------------------------------------------------------------
# Plot maps of pixel-wise Mann-Kendall trend stats (tau, p, and Sen's Slope)
# -----------------------------------------------------------------------------

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

# Mann-Kendall p
mkp.maps = ggplot(mkp.trend.df, aes(x = x, y = y, fill = MannKendall_p)) +
  geom_raster() +
  ggtitle('Mann-Kendall p-value trend 2003-2012') +
  guides(fill = guide_colorbar(title.position = 'bottom', title.hjust = 0.5, frame.colour = 'black', frame.linewidth = 0.5)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_stepsn( colors = c('gray96','royalblue','tomato','gold'),
                     # colors = c('red','orange','gray96','palegreen1','darkgreen'),
                     breaks = seq(0,1,0.1),
                     limits = c(0,1),
                     name = 'Mann-Kendall p-value') +
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
        legend.key.width = unit(3, 'cm'),
        legend.position = 'bottom',
        strip.background = element_rect(fill = 'gray27'),
        strip.text = element_text(size = 15, color = 'white'))

mkp.maps

ggsave(paste(figs.fp, 'MannKendall_p_trend_AGB_2003to2012.png', sep = '/'),
       mkp.maps,
       width = 8, height = 8, units = 'in')

# Mann-Kendall tau statistic
mktau.maps = ggplot(mktau.trend.df, aes(x = x, y = y, fill = MannKendall_tau)) +
  geom_raster() +
  ggtitle('Mann-Kendall tau statistic 2003-2012') +
  guides(fill = guide_colorbar(title.position = 'bottom', title.hjust = 0.5, frame.colour = 'black', frame.linewidth = 0.5)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_stepsn( colors = c('gold','gray96','royalblue'),
                     breaks = seq(-1,1,0.25),
                     limits = c(-1,1),
                     name = 'Mann-Kendall tau statistic') +
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
        legend.key.width = unit(3, 'cm'),
        legend.position = 'bottom',
        strip.background = element_rect(fill = 'gray27'),
        strip.text = element_text(size = 15, color = 'white'))

mktau.maps

ggsave(paste(figs.fp, 'MannKendall_tau_AGB_2003to2012.png', sep = '/'),
       mktau.maps,
       width = 8, height = 8, units = 'in')

# Mann-Kendall Sen's Slope
mkss.maps = ggplot(mkss.trend.df, aes(x = x, y = y, fill = MannKendall_SS)) +
  geom_raster() +
  ggtitle(bquote("Pixel-wise Mann-Kendall Sen's Slope 2003-2012")) +
  guides(fill = guide_colorbar(title.position = 'bottom', title.hjust = 0.5, frame.colour = 'black', frame.linewidth = 0.5)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_stepsn( colors = c('gold','gray96','royalblue'),
                     na.value = 'grey70',
                     breaks = seq(-5,5,1),
                     limits = c(-5,5),
                     #labels = c(as.character(seq(-5,-1,1)),'Grey = no trend', as.character(seq(1,5,1))),
                     name = bquote("Mann-Kendall Sen's Slope (grey = no trend)")) +
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

mkss.maps

ggsave(paste(figs.fp, 'MannKendall_SS_AGB_2003to2012.png', sep = '/'),
       mkss.maps,
       width = 8, height = 8, units = 'in')

# Violin plot showing distribution of Mann-Kendall Sen's Slope values per AGB product masked to the extent of Wang 
mkss.violin = ggplot(data = mkss.trend.df[mkss.trend.df$x %in% wang.mask.df$x & mkss.trend.df$y %in% wang.mask.df$y,],
                     aes(x = 1, y = MannKendall_SS, color = Product)) +
  geom_hline(yintercept = 0, linewidth = 1.25, color = 'grey70') +
  geom_violin(trim = TRUE,
              linewidth = 0.5,
              show.legend = FALSE,
              na.rm = TRUE,
              fill = 'transparent') +
  scale_color_manual(labels = c('AGB_VM','Wang','Xu','Liu','ESA_CCI'),
                     values = c(AGB_VM.col,wang.agb.col,xu.agb.col,liu.agb.col)) +
  scale_y_continuous(limits = c(-5,5), position = 'right') +
  geom_boxplot(width = 0.1, linewidth = 0.25, na.rm = TRUE, show.legend = FALSE) +
  ylab(bquote("Distribution of Sen's Slope")) +
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

mkss.violin

# Combine Sen's Slope maps and violin plots into single figure
mkss.maps.and.violin = plot_grid(plotlist = list(mkss.maps,mkss.violin),
                                 axis = 'lr',
                                 ncol = 2,
                                 align = 'h',
                                 labels = c('a','b'),
                                 label_size = 20,
                                 hjust = c(-2,1.25),
                                 rel_widths = c(0.7,0.175))

mkss.maps.and.violin

# Save 
ggsave(paste(figs.fp, 'MannKendall_SS_2003vs2012_w_violin_plot.png', sep = '/'),
       mkss.maps.and.violin,
       width = 13, height = 10, units = 'in',
       bg = ' white')

# Mann-Kendall Sen's Slope -- 4-panel horizontal (wide) layout
mkss.maps.wide = ggplot(mkss.trend.df, aes(x = x, y = y, fill = MannKendall_SS)) +
  geom_raster() +
  ggtitle(bquote("Mann-Kendall Sen's Slope 2003-2012")) +
  guides(fill = guide_colorbar(title.position = 'bottom', title.hjust = 0.5, frame.colour = 'black', frame.linewidth = 0.5)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_stepsn( colors = c('gold','gray96','royalblue'),
                     na.value = 'grey70',
                     breaks = seq(-5,5,1),
                     limits = c(-5,5),
                     #labels = c(as.character(seq(-5,-1,1)),'Grey = no trend', as.character(seq(1,5,1))),
                     name = bquote("Mann-Kendall Sen's Slope (grey = no trend)")) +
  theme_test() +
  geom_sf(data = st_as_sf(entire.domain), color = 'black', fill = NA, inherit.aes = FALSE) +
  # geom_sf(data = st_as_sf(core.domain), color = 'tomato', fill = NA, inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(wang.footprint), color = wang.agb.col, fill = NA, inherit.aes = FALSE, linewidth = 0.75) +
  geom_point(data = ba.info, aes(x = CentroidLON, y = CentroidLAT), pch = 1, inherit.aes = FALSE) +
  facet_wrap(.~Product, ncol = 4) +
  theme(plot.title = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12),
        legend.key.width = unit(3, 'cm'),
        legend.position = 'bottom',
        strip.background = element_rect(fill = 'gray27'),
        strip.text = element_text(size = 15, color = 'white'))

mkss.maps.wide

ggsave(paste(figs.fp, 'MannKendall_SS_AGB_2003to2012_WIDE.png', sep = '/'),
       mkss.maps.wide,
       width = 12, height = 4, units = 'in')

# ----------------------------------------------------------------------------
# Mask AGB products by land cover class
# ----------------------------------------------------------------------------

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

agb.lc.mask.fun = function(agb.in, product) {
  
  agb = mask(agb.in, wang.mask)
  
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
  return(agb.lc.df)
}

AGB_VM.agb.lc.df = agb.lc.mask.fun(AGB_VM[[8]], 'AGB_VM')
wang.agb.lc.df = agb.lc.mask.fun(wang.agb[[8]], 'Wang')
xu.agb.lc.df = agb.lc.mask.fun(xu.agb[[8]], 'Xu')
liu.agb.lc.df = agb.lc.mask.fun(liu.agb[[8]], 'Liu')
esa.agb.lc.df = agb.lc.mask.fun(esa.agb.10, 'ESA_CCI')

agb.lc.df = rbind(AGB_VM.agb.lc.df, wang.agb.lc.df, xu.agb.lc.df, liu.agb.lc.df, esa.agb.lc.df)
agb.lc.df$LC = factor(agb.lc.df$LC, levels = unique(agb.lc.df$LC))


# Mask Mann-Kendall stats per vegetation cover class
AGB_VM.mkss.lc.df = agb.lc.mask.fun(AGB_VM.mk[[6]], 'AGB_VM')
colnames(AGB_VM.mkss.lc.df)[1] = 'MK_SS'

wang.mkss.lc.df = agb.lc.mask.fun(wang.mk[[6]], 'Wang')
colnames(wang.mkss.lc.df)[1] = 'MK_SS'

xu.mkss.lc.df = agb.lc.mask.fun(xu.mk[[6]], 'Xu')
colnames(xu.mkss.lc.df)[1] = 'MK_SS'

liu.mkss.lc.df = agb.lc.mask.fun(liu.mk[[6]], 'Liu')
colnames(liu.mkss.lc.df)[1] = 'MK_SS'

AGB_VM.mkss.lc.median = aggregate(AGB_VM.mkss.lc.df, by = list(AGB_VM.mkss.lc.df$LC), FUN = function(x){median(x, na.rm = TRUE)})[,1:2]
wang.mkss.lc.median = aggregate(wang.mkss.lc.df, by = list(wang.mkss.lc.df$LC), FUN = function(x){median(x, na.rm = TRUE)})[,1:2]
xu.mkss.lc.median = aggregate(xu.mkss.lc.df, by = list(xu.mkss.lc.df$LC), FUN = function(x){median(x, na.rm = TRUE)})[,1:2]
liu.mkss.lc.median = aggregate(liu.mkss.lc.df, by = list(liu.mkss.lc.df$LC), FUN = function(x){median(x, na.rm = TRUE)})[,1:2]

# Create function to calculate Mann-Kendall's Sen's Slope quantiles 
mkss.lc.median.stats = function(mkss.lc.df, product.name) {
  
  # Get all unique class names
  lc.class.names = unique(mkss.lc.df$LC)
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
    lc.mkss = mkss.lc.df[mkss.lc.df$LC == lc.class.name,][['MK_SS']]
    n = length(lc.mkss)
    
    # Get mean and standard deviation of pct change LC subset
    lc.mkss.mean = mean(lc.mkss, na.rm = TRUE)
    lc.mkss.std = sd(lc.mkss, na.rm = TRUE)
    
    # Compute quantile values
    lc.mkss.quantiles = quantile(lc.mkss, na.rm = TRUE)
    
    # Assign 0%, 25%, 50%, 75%, and 100% quantile values to data frame for LC class subset
    quantile.df[i,2] = round(lc.mkss.quantiles[1],2)
    quantile.df[i,3] = round(lc.mkss.quantiles[2],2)
    quantile.df[i,4] = round(lc.mkss.quantiles[3],2)
    quantile.df[i,5] = round(lc.mkss.quantiles[4],2)
    quantile.df[i,6] = round(lc.mkss.quantiles[5],2)
    
    # Compute pnorms for each quantile using mean and stdev of pct. change for LC class subset
    pnorm.df[i,2] = round(pnorm(q = as.numeric(quantile.df[i,2]), mean = lc.mkss.mean, sd = lc.mkss.std, lower.tail = FALSE),2)
    pnorm.df[i,3] = round(pnorm(q = as.numeric(quantile.df[i,3]), mean = lc.mkss.mean, sd = lc.mkss.std, lower.tail = FALSE),2)
    pnorm.df[i,4] = round(pnorm(q = as.numeric(quantile.df[i,4]), mean = lc.mkss.mean, sd = lc.mkss.std, lower.tail = FALSE),2)
    pnorm.df[i,5] = round(pnorm(q = as.numeric(quantile.df[i,5]), mean = lc.mkss.mean, sd = lc.mkss.std, lower.tail = FALSE),2)
    pnorm.df[i,6] = round(pnorm(q = as.numeric(quantile.df[i,6]), mean = lc.mkss.mean, sd = lc.mkss.std, lower.tail = FALSE),2)
    
    # Get upper and lower confidence intervals for median probabiliy value
    k = qbinom(p = (1-0.95)/2, size = n, prob = 0.5, lower.tail = TRUE) # Using a 95% two-tailed confidence level
    ci = sort(lc.mkss)[c(k, n - k + 1)]
    attr(ci, 'conf.level') = 1 - 2 * pbinom(q = (k - 1), size = n, prob = 0.5)
    
    # Compute margin of error using confidence interval difference
    me = (ci[2] - ci[1]) / 2
    
    # Compute standard error(SE) by dividing error margin by 1.96
    se = me / 1.96
    
    # Subset pct change values using upper and lower confidence intervals
    lc.mkss.ci.subset = lc.mkss[lc.mkss <= ci[2] & lc.mkss >= ci[1]]
    
    # Compute two-sided t-test for CI-constrained pct change subset
    ttest = t.test(x = lc.mkss.ci.subset,
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
                            Stat = 'Median_MK_SS',
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

AGB_VM.mkss.lc.median.stats = mkss.lc.median.stats(AGB_VM.mkss.lc.df, 'AGB_VM')
wang.mkss.lc.median.stats = mkss.lc.median.stats(wang.mkss.lc.df, 'Wang')
xu.mkss.lc.median.stats = mkss.lc.median.stats(xu.mkss.lc.df, 'Xu')
liu.mkss.lc.median.stats = mkss.lc.median.stats(liu.mkss.lc.df, 'Liu')

agb.mkss.median.lc.stats = rbind(AGB_VM.mkss.lc.median.stats,
                                 wang.mkss.lc.median.stats,
                                 xu.mkss.lc.median.stats,
                                 liu.mkss.lc.median.stats)

# Generate table for 50% quantile and pnorm of M-K Sen's Slope per LC class and product
agb.mkss.median.lc.stats = reshape(agb.mkss.median.lc.stats, idvar = c('LC','Stat'), timevar = 'Product', v.names = 'Value', direction = 'wide')
colnames(agb.mkss.median.lc.stats) = c('LC','Stat','AGB_VM','Wang','Xu','Liu')

agb.mkss.median.lc.stats = data.frame(LandCoverClass = agb.mkss.median.lc.stats$LC[1:9],
                                      AGBVM_median = agb.mkss.median.lc.stats$AGB_VM[1:9],
                                      AGBVM_SE = agb.mkss.median.lc.stats$AGB_VM[19:27],
                                      AGBVM_p = agb.mkss.median.lc.stats$AGB_VM[28:36],
                                      Wang_median = agb.mkss.median.lc.stats$Wang[1:9],
                                      Wang_SE = agb.mkss.median.lc.stats$Wang[19:27],
                                      Wang_p = agb.mkss.median.lc.stats$Wang[28:36],
                                      Xu_median = agb.mkss.median.lc.stats$Xu[1:9],
                                      Xu_SE = agb.mkss.median.lc.stats$Xu[19:27],
                                      Xu_p = agb.mkss.median.lc.stats$Xu[28:36],
                                      Liu_median = agb.mkss.median.lc.stats$Liu[1:9],
                                      Liu_SE = agb.mkss.median.lc.stats$Liu[19:27],
                                      Liu_p = agb.mkss.median.lc.stats$Liu[28:36])

agb.mkss.median.lc.stats[,c(2,3,5,6,8,9,11,12)] = round(agb.mkss.median.lc.stats[,c(2,3,5,6,8,9,11,12)],2)
agb.mkss.median.lc.stats[agb.mkss.median.lc.stats$Liu_SE == 0.00,12] = NA
agb.mkss.median.lc.stats[is.nan(agb.mkss.median.lc.stats$Liu_p),13] = NA

# Re-label p values smaller than 0.001
agb.mkss.stats.relabel.pvals = function() {
  
  data.out = agb.mkss.median.lc.stats
  
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

agb.mkss.median.lc.stats = agb.mkss.stats.relabel.pvals()

agb.mkss.median.lc.stats[,c(2,5,8,11) == -0.00] = 0.00

mkss.lc.cols = scales::col_factor(palette = unname(class.cols[1:9]),
                                  domain = agb.mkss.median.lc.stats$LandCoverClass[1:9],
                                  levels = agb.mkss.median.lc.stats$LandCoverClass[1:9])

mkss.median.cols = scales::col_numeric(palette = c('gold','gray96','royalblue'),
                                       domain = c(-0.5, 0, 0.5))


agb.mkss.median.lc.stats.table = flextable(agb.mkss.median.lc.stats) %>%
  colformat_double(digits = 2) %>%
  separate_header(opts = 'span-top') %>%
  theme_vanilla() %>%
  align_text_col(align = 'center') %>%
  align_nottext_col(align = 'center') %>%
  bold(part = 'header') %>%
  vline(j = c(1,4,7,10)) %>%
  bg(i = 1:9,
     j = 1,
     bg = mkss.lc.cols) %>%
  bg(i = 1:9,
     j = c(2,5,8,11),
     bg = mkss.median.cols) %>%
  bg(i = 1:9,
     j = c(3,4,6,7,9,10,12,13),
     bg = 'white') %>%
  bg(i = c(1,2), 
     bg = 'white',
     part = 'header') %>%
  border_outer(part = 'all', border = fp_border_default(color = 'black', width = 2)) %>%
  fix_border_issues() %>%
  set_caption(caption = as_paragraph(as_chunk(bquote("Mann-Kendall Sen's Slope trends 2003-2012 per land cover class"),
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


agb.mkss.median.lc.stats.table
# Manually save as .png

# Write data frame to .csv
write.csv(agb.mkss.median.lc.stats,
          file = paste(tables.fp, 'Table_02_MannKendall_SensSlope_2003to2012_VegClasses.csv', sep = '/'),
          row.names = FALSE)
