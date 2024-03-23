# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Perform Akaike Information Criterion (AIC) for different configurations of AGB retrieval model to determine which is suited best.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# -----------------------
# Get libraries
# -----------------------
library(rstudioapi)
library(raster)
library(moments)
library(tidyverse)
library(colorRamps)
library(AICcmodavg)
library(data.table)
library(cowplot)
library(grid)
library(gridExtra)
library(gridGraphics)
library(ggcorrplot)
library(RColorBrewer)
library(flextable)

setwd('../../')
git.fp = getwd()
figs.fp = paste(git.fp, 'Figures', sep = '/')
code.fp = paste(git.fp, 'Code', sep = '/')
setwd('../../')
root.fp = getwd()
data.fp = paste(root.fp, 'Data', sep = '/')
tables.fp = paste(git.fp, 'CSVs_for_Tables', sep = '/')


# -----------------------
# ABoVE domain shapefiles
# -----------------------
shp.fp = paste(git.fp, 'Data/Shapefiles/Domain', sep = '/')
setwd(shp.fp)

core.domain = shapefile('Core_Region_WGS.shp')
entire.domain = shapefile('Extended_Plus_Core_Domains_WGS.shp')

# ----------------------
# VODCA land mask file
# ----------------------
lm.fp = paste(data.fp, 'VODCA', sep = '/')
setwd(lm.fp)

landmask = raster('VODCA_ABoVE_Land_Mask.tif')
landmask[landmask == 0] = NA

# ---------------------------------------------------------------------------------------------
# ESA AGB (Santoro)
# ---------------------------------------------------------------------------------------------
# 2017
esa.agb.17.fp = paste(data.fp, 'Biomass/Santoro_ESA/2017', sep = '/')
setwd(esa.agb.17.fp)
esa.agb.17 = raster('Santoro_AGB_2017_ABoVE_025deg.tif')
esa.agb.17[esa.agb.17 == 0] = NA

log.agb.17 = log(esa.agb.17)
log.agb.17[log.agb.17 < 0] = NA

# ---------------------------------------------------------------------------------------------
# VODCA VOD
# ---------------------------------------------------------------------------------------------
# C-Band
c.vod.ym.fp = paste(data.fp, 'VODCA/C_Band/Yearly_MeanGS', sep = '/')
setwd(c.vod.ym.fp)
c.vod.ym = stack(list.files(c.vod.ym.fp, pattern = '.tif', recursive = FALSE, full.names = FALSE)[1:15])
c.vod.17 = c.vod.ym[[15]]

# X-Band
x.vod.ym.fp = paste(data.fp, 'VODCA/X_Band/Yearly_MeanGS', sep = '/')
setwd(x.vod.ym.fp)
x.vod.ym = stack(list.files(x.vod.ym.fp, pattern = '.tif', recursive = FALSE, full.names = FALSE)[1:15])
x.vod.17 = x.vod.ym[[15]]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Multispectral Vegetation Indices (MVIs)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# -----------------------
# VARI (MODIS MCD43A4)
# -----------------------

vari.fp = paste(data.fp, 'MODIS/MCD43A4_VARI/025_Deg', sep = '/')
setwd(vari.fp)
vari.files = list.files(vari.fp, pattern = 'GS', full.names = FALSE, recursive = FALSE)
vari = stack(vari.files)
vari.17 = vari[[15]]

# -----------------------
# NIRV (MODIS MCD43A4)
# -----------------------

nirv.fp = paste(data.fp, 'MODIS/MCD43A4_NDVI+NIRV/025_Deg', sep = '/')
setwd(nirv.fp)
ndvi.files = list.files(nirv.fp, pattern = 'NDVI', full.names = FALSE, recursive = FALSE)
nir.files = list.files(nirv.fp, pattern = 'NIR', full.names = FALSE, recursive = FALSE)
ndvi = stack(ndvi.files)
nir = stack(nir.files)
nirv = ndvi * nir
ndvi.17 = ndvi[[15]]
nirv.17 = nirv[[15]]

# -----------------------
# LAI (MODIS MOD15A2H)
# -----------------------

lai.fp = paste(data.fp, 'MODIS/MOD15A2H_LAI/025_Deg', sep = '/')
setwd(lai.fp)
lai.files = list.files(lai.fp, pattern = 'GS', full.names = FALSE)
lai = stack(lai.files)
lai.17 = lai[[15]]

# ---------------------------------------------------------------------------------------------
# Distribution of VOD at different AGB intervals (AGB and logAGB)
# ---------------------------------------------------------------------------------------------

vod.agb.dist.fun = function(){
  
  # AGB 0-50 Mg/ha
  agb.0.50 = esa.agb.17
  agb.0.50[agb.0.50 > 50] = NA
  # AGB 50-100 Mg/ha
  agb.50.100 = esa.agb.17
  agb.50.100[agb.50.100 < 50] = NA
  agb.50.100[agb.50.100 > 100] = NA
  # AGB 100-150 Mg/ha
  agb.100.150 = esa.agb.17
  agb.100.150[agb.100.150 < 100] = NA
  agb.100.150[agb.100.150 > 150] = NA
  # AGB 150-200 Mg/ha
  agb.150.200 = esa.agb.17
  agb.150.200[agb.150.200 < 150] = NA
  agb.150.200[agb.150.200 > 200] = NA
  # AGB 200-260 Mg/ha
  agb.200.260 = esa.agb.17
  agb.200.260[agb.200.260 < 200] = NA
  
  # --- Log transform AGB and then mask by logAGB intervals
  agb.log = log(esa.agb.17)
  agb.log[agb.log < 0] = NA
  
  # logAGB 0-1.2
  agb.log.0to1.2 = agb.log
  agb.log.0to1.2[agb.log.0to1.2 > 1.2] = NA
  # logAGB 1.2-2.4
  agb.log.1.2to2.4 = agb.log
  agb.log.1.2to2.4[agb.log.1.2to2.4 < 1.2] = NA
  agb.log.1.2to2.4[agb.log.1.2to2.4 > 2.4] = NA
  # logAGB 2.4-3.6
  agb.log.2.4to3.6 = agb.log
  agb.log.2.4to3.6[agb.log.2.4to3.6 < 2.4] = NA
  agb.log.2.4to3.6[agb.log.2.4to3.6 > 3.6] = NA
  # logAGB 3.6-4.8
  agb.log.3.6to4.8 = agb.log
  agb.log.3.6to4.8[agb.log.3.6to4.8 < 3.6] = NA
  agb.log.3.6to4.8[agb.log.3.6to4.8 > 4.8] = NA
  # log 4.8-6
  agb.log.4.8to6 = agb.log
  agb.log.4.8to6[agb.log.4.8to6 < 4.8] = NA
  
  # Mask VOD by AGB intervals
  c.vod.0.50 = na.omit(mask(c.vod.17, agb.0.50)[])
  c.vod.50.100 = na.omit(mask(c.vod.17, agb.50.100)[])
  c.vod.100.150 = na.omit(mask(c.vod.17, agb.100.150)[])
  c.vod.150.200 = na.omit(mask(c.vod.17, agb.150.200)[])
  c.vod.200.260 = na.omit(mask(c.vod.17, agb.200.260)[])
  
  l1 = c(length(c.vod.0.50),length(c.vod.50.100),length(c.vod.100.150),length(c.vod.150.200),length(c.vod.200.260))
  
  c.vod.agb.1 = rbind(data.frame(vod = c.vod.0.50, int = rep('0-50',l1[1]), scale = rep('Linear',l1[1])),
                      data.frame(vod = c.vod.50.100, int = rep('50-100',l1[2]), scale = rep('Linear',l1[2])),
                      data.frame(vod = c.vod.100.150, int = rep('100-150',l1[3]), scale = rep('Linear',l1[3])),
                      data.frame(vod = c.vod.150.200, int = rep('150-200',l1[4]), scale = rep('Linear',l1[4])),
                      data.frame(vod = c.vod.200.260, int = rep('200-260',l1[5]), scale = rep('Linear',l1[5])))
  
  c.vod.agb.1$int = factor(c.vod.agb.1$int, levels = c('0-50','50-100','100-150','150-200','200-260'))
  
  # Mask VOD by logAGB intervals
  c.vod.log1 = na.omit(mask(c.vod.17, agb.log.0to1.2)[])
  c.vod.log2 = na.omit(mask(c.vod.17, agb.log.1.2to2.4)[])
  c.vod.log3 = na.omit(mask(c.vod.17, agb.log.2.4to3.6)[])
  c.vod.log4 = na.omit(mask(c.vod.17, agb.log.3.6to4.8)[])
  c.vod.log5 = na.omit(mask(c.vod.17, agb.log.4.8to6)[])
  
  l2 = c(length(c.vod.log1),length(c.vod.log2),length(c.vod.log3),length(c.vod.log4),length(c.vod.log5))
  
  c.vod.agb.2 = rbind(data.frame(vod = c.vod.log1, int = rep('0-1.2',l2[1]), scale = rep('logAGB',l2[1])),
                      data.frame(vod = c.vod.log2, int = rep('1.2-2.4',l2[2]), scale = rep('logAGB',l2[2])),
                      data.frame(vod = c.vod.log3, int = rep('2.4-3.6',l2[3]), scale = rep('logAGB',l2[3])),
                      data.frame(vod = c.vod.log4, int = rep('3.6-4.8',l2[4]), scale = rep('logAGB',l2[4])),
                      data.frame(vod = c.vod.log5, int = rep('4.8-6',l2[5]), scale = rep('logAGB',l2[5])))
  
  c.vod.agb.2$int = factor(c.vod.agb.2$int, levels = c('0-1.2','1.2-2.4','2.4-3.6','3.6-4.8','4.8-6'))
  
  # ------------------------
  # Figures
  # ------------------------
  
  # ~~~~~~ Distributions
  agb.vals = na.omit(esa.agb.17[])
  log.agb.vals = na.omit(agb.log[])
  
  agb.df = data.frame(agb = agb.vals, scale = rep('AGB',length(agb.vals)))
  log.agb.df = data.frame(agb = log.agb.vals, scale = rep('logAGB',length(log.agb.vals)))
  
  # Original AGB
  agb.hist = agb.df %>%
    ggplot(aes(x = agb)) +
    geom_histogram(binwidth = 2, aes(y = stat(density))) +
    #ggtitle('Distiribution of original AGB values') +
    geom_density(col = 'tomato', size = 2) +
    theme_bw() +
    annotation_custom(grobTree(textGrob('a', x = 0.01, y = 0.95, hjust = 0,
                                        gp = gpar(col = 'black', fontsize = 18, fontface = 'bold')))) +
    ylab('Density') + xlab('AGB [Mg/ha]') +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_blank())
  
  # logAGB
  log.agb.hist = log.agb.df %>%
    ggplot(aes(x = agb)) +
    geom_histogram(binwidth = 0.03, aes(y = stat(density))) +
    #ggtitle('Distiribution of log-transformed AGB values') +
    geom_density(col = 'tomato', size = 2) +
    theme_bw() +
    annotation_custom(grobTree(textGrob('d', x = 0.01, y = 0.95, hjust = 0,
                                        gp = gpar(col = 'black', fontsize = 18, fontface = 'bold')))) +
    ylab('Density') + xlab('log(AGB) [Mg/ha]') +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_blank())
  
  #agb.hist
  #log.agb.hist
  
  # ~~~~~~ Scatterplots
  agb.cband = mask(esa.agb.17, c.vod.17)[]
  log.agb.cband = mask(agb.log, c.vod.17)[]
  
  agb.cband.df = na.omit(data.frame(agb = agb.cband, vod = c.vod.17[]))
  log.agb.cband.df = na.omit(data.frame(agb = log.agb.cband, vod = c.vod.17[]))
  
  agb.cband.sp = ggplot(agb.cband.df) +
    geom_point(aes(x = agb, y = vod), alpha = 0.1, pch = 16) +
    scale_y_continuous(limits = c(0,1)) +
    ylab('C-Band VOD [ ]') + 
    xlab('AGB [Mg/ha]') +
    theme_bw() +
    annotation_custom(grobTree(textGrob('b', x = 0.01, y = 0.95, hjust = 0,
                                        gp = gpar(col = 'black', fontsize = 18, fontface = 'bold')))) +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15))
  
  # logAGB
  log.agb.cband.sp = ggplot(log.agb.cband.df) +
    geom_point(aes(x = agb, y = vod), alpha = 0.1, pch = 16) +
    scale_y_continuous(limits = c(0,1)) +
    ylab('C-Band VOD [ ]') + 
    xlab('log(AGB) [Mg/ha]') +
    theme_bw() +
    annotation_custom(grobTree(textGrob('e', x = 0.01, y = 0.95, hjust = 0,
                                        gp = gpar(col = 'black', fontsize = 18, fontface = 'bold')))) +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15))
  
  #agb.cband.sp
  #log.agb.cband.sp
  
  
  # ~~~~~~ Boxplots
  # Boxplot 1
  c.vod.agb.plot = ggplot(c.vod.agb.1,
                          aes(y = vod, x = int)) +
    geom_boxplot(position = 'dodge') +
    scale_y_continuous(limits = c(0,1)) +
    ylab('C-Band VOD [ ]') + 
    xlab('AGB [Mg/ha]') +
    theme_bw() +
    annotation_custom(grobTree(textGrob('c', x = 0.01, y = 0.95, hjust = 0,
                                        gp = gpar(col = 'black', fontsize = 18, fontface = 'bold')))) +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15))
  
  # Boxplot 2
  c.vod.logagb.plot = ggplot(c.vod.agb.2,
                             aes(y = vod, x = int)) +
    geom_boxplot(position = 'dodge') +
    scale_y_continuous(limits = c(0,1)) +
    ylab('C-Band VOD [ ]') + 
    xlab('log(AGB) [Mg/ha]') +
    theme_bw() +
    annotation_custom(grobTree(textGrob('f', x = 0.01, y = 0.95, hjust = 0,
                                        gp = gpar(col = 'black', fontsize = 18, fontface = 'bold')))) +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15))
  
  #c.vod.agb.plot
  #c.vod.logagb.plot
  
  grid.arrange(agb.hist, agb.cband.sp, c.vod.agb.plot,
               log.agb.hist, log.agb.cband.sp, c.vod.logagb.plot,
               ncol = 3, nrow = 2)
}


vod.agb.dist.plots = vod.agb.dist.fun()

ggsave(paste(figs.fp, 'AIC', 'VOD_AGB_DIST.png', sep = '/'),
       vod.agb.dist.plots,
       width = 15, height = 8)

# ---------------------------------------------------------------------------------------------
# Generate list of all MLR model combinations integrating VOD, VARI, NDVI, NIRv, and LAI
# ---------------------------------------------------------------------------------------------
list.model.combos = function(agb.in, vod.in, vari.in, ndvi.in, nirv.in, lai.in) {
  
  vod = vod.in[]
  agb = mask(agb.in, vod.in)[]
  vari = mask(vari.in, vod.in)[]
  ndvi = mask(ndvi.in, vod.in)[]
  nirv = mask(nirv.in, vod.in)[]
  lai = mask(lai.in, vod.in)[]
  
  modelvars.df = data.frame(AGB = agb,
                            VOD = vod,
                            VARI = vari,
                            NDVI = ndvi,
                            NIRV = nirv,
                            LAI = lai)
  models = list()
  
  for (i in 1:5) {
    vc = combn(names(modelvars.df)[2:6], i)
    for (j in 1 : ncol(vc)) {
      model = as.formula(paste0('AGB ~', paste0(vc[,j], collapse = '+')))
      models = c(models,model)
    }
  }
  
  return(models)
}

model.list = list.model.combos(esa.agb.17, c.vod.17, vari.17, ndvi.17, nirv.17, lai.17)

# ---------------------------------------------------------------------------------------------
# Compare different model configurations using AIC test
# ---------------------------------------------------------------------------------------------

diff.pal = colorRampPalette(c('orange','gray96','forestgreen'))
breakpoints = seq(-4,4,0.5)

# RMSE function
# ------------------
rmse.fun = function(x.in, y.in){
  x = x.in[]
  y = mask(y.in, x.in)[]
  
  minus = x - y
  square = minus^2
  mean = mean(square, na.rm = TRUE)
  RMSE = sqrt(mean) 
  return(round(RMSE,0))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model AIC comparison
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#log.agb.17 = log(esa.agb.17)
#log.agb.17[log.agb.17 < 0] = NA

vod.agb.aic.test = function(agb.in, vod.in, vari.in, ndvi.in, nirv.in, lai.in, main.text, log.data){
  
  vod = vod.in[]
  agb = mask(agb.in, vod.in)[]
  vari = mask(vari.in, vod.in)[]
  ndvi = mask(ndvi.in, vod.in)[]
  nirv = mask(nirv.in, vod.in)[]
  lai = mask(lai.in, vod.in)[]
  
  # VOD only
  mod01.vod = lm(agb ~ vod, na.action = na.exclude)
  
  # VARI configurations
  mod02.vod.vari = lm(agb ~ vod + vari, na.action = na.exclude)
  mod03.vod.vari.ndvi = lm(agb ~ vod + vari + ndvi, na.action = na.exclude)
  mod04.vod.vari.ndvi.nirv = lm(agb ~ vod + vari + ndvi + nirv, na.action = na.exclude)
  mod05.vod.vari.ndvi.nirv.lai = lm(agb ~ vod + vari + ndvi + nirv + lai, na.action = na.exclude)
  mod06.vod.vari.ndvi.lai = lm(agb ~ vod + vari + ndvi + lai, na.action = na.exclude)
  mod07.vod.vari.nirv = lm(agb ~ vod + vari + nirv, na.action = na.exclude)
  mod08.vod.vari.nirv.lai = lm(agb ~ vod + vari + nirv + lai, na.action = na.exclude)
  mod09.vod.vari.lai = lm(agb ~ vod + vari + lai, na.action = na.exclude)
  
  # Configurations excluding VARI
  mod10.vod.ndvi = lm(agb ~ vod + ndvi, na.action = na.exclude)
  mod11.vod.ndvi.nirv = lm(agb ~ vod + ndvi + nirv, na.action = na.exclude)
  mod12.vod.ndvi.lai = lm(agb ~ vod + ndvi + lai, na.action = na.exclude)
  mod13.vod.ndvi.nirv.lai = lm(agb ~ vod + ndvi + nirv + lai, na.action = na.exclude)
  
  # Configurations excluding VARI and NDVI
  mod14.vod.nirv = lm(agb ~ vod + nirv, na.action = na.exclude)
  mod15.vod.nirv.lai = lm(agb ~ vod + nirv + lai, na.action = na.exclude)
  
  # Configuration excluding VARI, NDVI, and NIRv
  mod16.vod.lai = lm(agb ~ vod + lai, na.action = na.exclude)
  
  # VIs alone
  mod17.vari = lm(agb ~ vari, na.action = na.exclude)
  mod18.vari.ndvi = lm(agb ~ vari + ndvi, na.action = na.exclude)
  mod19.vari.ndvi.nirv = lm(agb ~ vari + ndvi + nirv, na.action = na.exclude)
  mod20.vari.ndvi.nirv.lai = lm(agb ~ vari + ndvi + nirv + lai, na.action = na.exclude)
  mod21.vari.nirv = lm(agb ~ vari + nirv, na.action = na.exclude)
  mod22.vari.nirv.lai = lm(agb ~ vari + nirv + lai, na.action = na.exclude)
  mod23.vari.lai = lm(agb ~ vari + lai, na.action = na.exclude)
  
  mod24.ndvi = lm(agb ~ ndvi, na.action = na.exclude)
  mod25.ndvi.nirv = lm(agb ~ ndvi + nirv, na.action = na.exclude)
  mod26.ndvi.nirv.lai = lm(agb ~ ndvi + nirv + lai, na.action = na.exclude)
  
  mod27.nirv = lm(agb ~ nirv, na.action = na.exclude)
  mod28.nirv.lai = lm(agb ~ nirv + lai, na.action = na.exclude)
  
  mod29.lai = lm(agb ~ lai, na.action = na.exclude)
  
  mod30.ndvi.lai = lm(agb ~ ndvi + lai, na.action = na.exclude)
  mod31.vari.ndvi.lai = lm(agb ~ vari + ndvi + lai, na.action = na.exclude)
  
  # Combine all models into single list
  all.mods.list = list(mod01.vod,
                       mod02.vod.vari,
                       mod03.vod.vari.ndvi,
                       mod04.vod.vari.ndvi.nirv,
                       mod05.vod.vari.ndvi.nirv.lai,
                       mod06.vod.vari.ndvi.lai,
                       mod07.vod.vari.nirv,
                       mod08.vod.vari.nirv.lai,
                       mod09.vod.vari.lai,
                       mod10.vod.ndvi,
                       mod11.vod.ndvi.nirv,
                       mod12.vod.ndvi.lai,
                       mod13.vod.ndvi.nirv.lai,
                       mod14.vod.nirv,
                       mod15.vod.nirv.lai,
                       mod16.vod.lai,
                       mod17.vari,
                       mod18.vari.ndvi,
                       mod19.vari.ndvi.nirv,
                       mod20.vari.ndvi.nirv.lai,
                       mod21.vari.nirv,
                       mod22.vari.nirv.lai,
                       mod23.vari.lai,
                       mod24.ndvi,
                       mod25.ndvi.nirv,
                       mod26.ndvi.nirv.lai,
                       mod27.nirv,
                       mod28.nirv.lai,
                       mod29.lai,
                       mod30.ndvi.lai,
                       mod31.vari.ndvi.lai)
  
  # Create vector of model names corresponding to list
  all.mods.names = c('VOD alone',  #01
                     'VOD+VARI',  #02
                     'VOD+VARI+NDVI',  #03
                     'VOD+VARI+NDVI+NIRv',  #04
                     'VOD+VARI+NDVI+NIRv+LAI',  #05
                     'VOD+VARI+NDVI+LAI',  #06
                     'VOD+VARI+NIRv',  #07
                     'VOD+VARI+NIRv+LAI',  #08
                     'VOD+VARI+LAI',  #09
                     'VOD+NDVI',  #10
                     'VOD+NDVI+NIRv',  #11
                     'VOD+NDVI+LAI',  #12
                     'VOD+NDVI+NIRv+LAI',  #13
                     'VOD+NIRv',  #14
                     'VOD+NIRv+LAI',  #15
                     'VOD+LAI',  #16
                     'VARI alone',  #17
                     'VARI+NDVI',  #18
                     'VARI+NDVI+NIRv',  #19
                     'VARI+NDVI+NIRv+LAI',  #20
                     'VARI+NIRv',  #21
                     'VARI+NIRv+LAI',  #22
                     'VARI+LAI',  #23
                     'NDVI alone',  #24
                     'NDVI+NIRv',  #25
                     'NDVI+NIRv+LAI',  #26
                     'NIRv alone',  #27
                     'NIRv+LAI',  #28
                     'LAI alone',  #29 
                     'NDVI+LAI', #30
                     'VARI+NDVI+LAI' #31
  )
  
  # Generate AIC report
  all.mods.aic = aictab(all.mods.list, all.mods.names, second.ord = FALSE, sort = TRUE)
  print(all.mods.aic)
  
  # ~~~~~~~~~ Plot AIC report as table
  # grid.newpage()
  # 
  # # Set theme/appearance parameters
  # tt3 = ttheme_minimal(
  #   core=list(bg_params = list(fill = blues9[c(3,1)], col=NA),
  #             fg_params=list()),
  #   colhead=list(fg_params=list(col="black")),
  #   rowhead=list(fg_params=list(col="black")))
  # 
  # # Plot table
  # grid.arrange(tableGrob(all.mods.aic,
  #                        theme = tt3),
  #              nrow = 1,
  #              top = main.text)
  
  # Create output data frame with AIC results
  out.df = setDT(all.mods.aic, keep.rownames = TRUE)[]
  out.df = out.df[,c(1,2,3,4,5,8)]
  out.df$rn = sprintf('%02d', as.numeric(out.df$rn))
  out.df$K = sprintf('%02d', as.numeric(out.df$K))
  out.df$mr = sprintf('%02d', seq(1,nrow(out.df),1))
  colnames(out.df)[1] = 'Model Index'
  colnames(out.df)[2] = 'Model Name'
  colnames(out.df)[3] = '# Parameters'
  colnames(out.df)[6] = 'Log-Likelihood'
  colnames(out.df)[7] = 'Model Rank'
  
  out.df = out.df[,c(7,1,2,3,4,5,6)]
  
  aic.report.tab = flextable(out.df) %>%
    theme_vanilla() %>%
    align_text_col(align = 'center') %>%
    align_nottext_col(align = 'center') %>%
    set_table_properties(layout = 'autofit') %>%
    add_header_row(values = main.text,
                   colwidths = 7, top = TRUE) %>%
    bg(part = 'all', bg = 'white')
  
  save_as_image(aic.report.tab,
                path = paste(figs.fp, 'AIC', paste0('AIC_Report_',log.data,'.png'), sep = '/'))
  
  write.csv(out.df,
            file = paste(tables.fp, paste0('AIC_Report_',log.data,'.csv'), sep = '/'),
            row.names = FALSE)
  
  return(out.df)
}

# ---------------------------------------------------------------------------------------------
# Run tests
# ---------------------------------------------------------------------------------------------

# ~~~~ c-Band VOD
# Linear test
cvod.agb.aic = vod.agb.aic.test(esa.agb.17, c.vod.17, vari.17, ndvi.17, nirv.17, lai.17,
                                'No log transformation for any model components',
                                'CbandVOD_logNone')

# Non-linear test (log transform of AGB data only)
cvod.logagb.aic = vod.agb.aic.test(log.agb.17, c.vod.17, vari.17, ndvi.17, nirv.17, lai.17,
                                   'Log transformation(s): target AGB',
                                   'CbandVOD_logAGB')

# Non-linear test (log transform of VOD data only)
logcvod.agb.aic = vod.agb.aic.test(esa.agb.17, log(c.vod.17), vari.17, ndvi.17, nirv.17, lai.17,
                                   'Log transformation(s): input C-Band VOD',
                                   'CbandVOD_logVOD')

# Non-linear test (log transform of AGB and VOD data)
logcvod.logagb.aic = vod.agb.aic.test(log.agb.17, log(c.vod.17), vari.17, ndvi.17, nirv.17, lai.17,
                                      'Log transformation(s): target AGB, input C-Band VOD',
                                      'CbandVOD_logAGBlogVOD')

# ~~~~ X-Band VOD
# Linear test
xvod.agb.aic = vod.agb.aic.test(esa.agb.17, x.vod.17, vari.17, ndvi.17, nirv.17, lai.17,
                                'No log transformation for any model components',
                                'XbandVOD_logNone')

# Non-linear test (log transform of AGB data only)
xvod.logagb.aic = vod.agb.aic.test(log.agb.17, x.vod.17, vari.17, ndvi.17, nirv.17, lai.17,
                                   'Log transformation(s): target AGB',
                                   'XbandVOD_logAGB')

# Non-linear test (log transform of VOD data only)
logxvod.agb.aic = vod.agb.aic.test(esa.agb.17, log(x.vod.17), vari.17, ndvi.17, nirv.17, lai.17,
                                   'Log transformation(s): input X-Band VOD',
                                   'XbandVOD_logVOD')

# Non-linear test (log transform of AGB and VOD data)
logxvod.logagb.aic = vod.agb.aic.test(log.agb.17, log(x.vod.17), vari.17, ndvi.17, nirv.17, lai.17,
                                      'Log transformation(s): target AGB, input X-Band VOD',
                                      'XbandVOD_logAGBlogVOD')


# Order each AIC dataframe by model index (01-31)
cvod.agb.aic = cvod.agb.aic[order(cvod.agb.aic$`Model Index`),]
cvod.logagb.aic = cvod.logagb.aic[order(cvod.logagb.aic$`Model Index`),]
logcvod.agb.aic = logcvod.agb.aic[order(logcvod.agb.aic$`Model Index`),]
logcvod.logagb.aic = logcvod.logagb.aic[order(logcvod.logagb.aic$`Model Index`),]

xvod.agb.aic = xvod.agb.aic[order(xvod.agb.aic$`Model Index`),]
xvod.logagb.aic = xvod.logagb.aic[order(xvod.agb.aic$`Model Index`),]
logxvod.agb.aic = logxvod.agb.aic[order(logxvod.agb.aic$`Model Index`),]
logxvod.logagb.aic = logxvod.logagb.aic[order(logxvod.logagb.aic$`Model Index`),]

# Combine AIC scores from all three tests into single data frame
caic.combined = data.frame(c(cvod.agb.aic$AIC,
                             cvod.logagb.aic$AIC,
                             logcvod.agb.aic$AIC,
                             logcvod.logagb.aic$AIC),
                           c(cvod.agb.aic$Delta_AIC,
                             cvod.logagb.aic$Delta_AIC,
                             logcvod.agb.aic$Delta_AIC,
                             logcvod.logagb.aic$Delta_AIC),
                           c(cvod.agb.aic$`Model Index`,
                             cvod.logagb.aic$`Model Index`,
                             logcvod.agb.aic$`Model Index`,
                             logcvod.logagb.aic$`Model Index`),
                           c(cvod.agb.aic$`Model Name`,
                             cvod.logagb.aic$`Model Name`,
                             logcvod.agb.aic$`Model Name`,
                             logcvod.logagb.aic$`Model Name`),
                           c(rep('1) Linear',nrow(cvod.agb.aic)),
                             rep('3) logAGB',nrow(cvod.agb.aic)),
                             rep('2) logVOD',nrow(cvod.agb.aic)),
                             rep('4) logAGB & logVOD',nrow(cvod.agb.aic))))

colnames(caic.combined) = c('AIC','dAIC','ModID','ModName','ModConfig')
caic.combined$VOD_Band = rep('C-Band', length(caic.combined$AIC))

xaic.combined = data.frame(c(xvod.agb.aic$AIC,
                             xvod.logagb.aic$AIC,
                             logxvod.agb.aic$AIC,
                             logxvod.logagb.aic$AIC),
                           c(xvod.agb.aic$Delta_AIC,
                             xvod.logagb.aic$Delta_AIC,
                             logxvod.agb.aic$Delta_AIC,
                             logxvod.logagb.aic$Delta_AIC),
                           c(xvod.agb.aic$`Model Index`,
                             xvod.logagb.aic$`Model Index`,
                             logxvod.agb.aic$`Model Index`,
                             logxvod.logagb.aic$`Model Index`),
                           c(xvod.agb.aic$`Model Name`,
                             xvod.logagb.aic$`Model Name`,
                             logxvod.agb.aic$`Model Name`,
                             logxvod.logagb.aic$`Model Name`),
                           c(rep('1) Linear',nrow(xvod.agb.aic)),
                             rep('3) logAGB',nrow(xvod.agb.aic)),
                             rep('2) logVOD',nrow(xvod.agb.aic)),
                             rep('4) logAGB & logVOD',nrow(xvod.agb.aic))))

colnames(xaic.combined) = c('AIC','dAIC','ModID','ModName','ModConfig')
xaic.combined$VOD_Band = rep('X-Band', length(xaic.combined$AIC))


# Combine AIC data frames 
aic.combined = rbind(caic.combined, xaic.combined)


# Create plot summarizing results of AIC tests
aic.plot = ggplot(aic.combined) +
  geom_bar(aes(x = ModName, y = AIC, fill = VOD_Band), 
           position = 'dodge', stat = 'identity') +
  scale_fill_manual(values = c('tomato','steelblue')) +
  theme_minimal() +
  theme(legend.position = 'top',
        legend.justification = 'right',
        legend.title = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1))


# ---------------------------------------------------------------------------------------------
# Subset aic.combined so it contains a fewer number of models, generate final figures
# ---------------------------------------------------------------------------------------------

aic.final = aic.combined
aic.final = aic.final[aic.final$ModName %in% c('VOD alone','LAI alone','NDVI alone','VARI alone','NIRv alone',
                                               'VOD+LAI','VOD+NDVI','VOD+NIRv','VOD+VARI',
                                               'VOD+VARI+LAI','VOD+VARI+NDVI','VOD+VARI+NIRv'),]

# Subset further to separate AIC scores for linear and log-transformed models
aic.final.linear = aic.final[aic.final$ModConfig == '1) Linear',]
aic.final.logagb = aic.final[aic.final$ModConfig == '3) logAGB',]


# Create plot showing AIC results from the subsets
linear.aic.plot = ggplot(aic.final.linear) +
  geom_bar(aes(x = reorder(ModName, -AIC), y = AIC, fill = VOD_Band), 
           position = 'dodge', stat = 'identity', 
           width = 0.80, show.legend = FALSE) +
  ggtitle(expression('AGB'['ref'] == beta[0]+beta[1]*'VOD'+beta[2]*'VI'[1]+'... ...')) +
  theme_bw() +
  scale_fill_manual(values = c('tomato','steelblue')) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 12),
        axis.title.x = element_blank())

logagb.aic.plot = ggplot(aic.final.logagb) +
  geom_bar(aes(x = reorder(ModName, -AIC), y = AIC, fill = VOD_Band), 
           position = 'dodge', stat = 'identity', 
           width = 0.80, show.legend = TRUE) +
  ggtitle(expression('logAGB'['ref'] == beta[0]+beta[1]*'VOD'+beta[2]*'VI'[1]+'... ...')) +
  theme_bw() +
  scale_fill_manual(values = c('tomato','steelblue')) +
  # annotation_custom(grobTree(textGrob('A', x = 0.01, y = 0.95, hjust = 0,
  #                                     gp = gpar(col = 'black', fontsize = 20, fontface = 'bold')))) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 12),
        axis.title.x = element_blank(),
        legend.position = c(0.98,0.83),
        legend.justification = 'right',
        legend.title = element_text(face = 'bold'))


grid.arrange(linear.aic.plot, logagb.aic.plot, ncol = 1, nrow = 2)


# ---------------------------------------------------------------------------------------------
# Create correlation matrix to visualize multiplecollinearity bewteen variables
# ---------------------------------------------------------------------------------------------

# Linear relationships
linear.df.17 = data.frame(esa.agb.17[], c.vod.17[], x.vod.17[], vari.17[], ndvi.17[], nirv.17[], lai.17[])
colnames(linear.df.17) = c('ESA_CCI_AGB','C-Band_VOD','X-Band_VOD','VARI','NDVI','NIRv','LAI')
linear.df.17.cor = cor(linear.df.17, method = 'pearson', use = 'complete.obs')

linear.corrplot = ggcorrplot(linear.df.17.cor, hc.order = FALSE, type = 'upper',
                             outline.col = 'white',
                             lab = TRUE, show.legend = FALSE,
                             title = 'Pearson R Correlation',
                             colors = c('tomato','white','seagreen')) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 12))


# Relationships with log-transformed reference AGB
logagb.df.17 = data.frame(log(esa.agb.17[]), c.vod.17[], x.vod.17[], vari.17[], ndvi.17[], nirv.17[], lai.17[])
colnames(logagb.df.17) = c('log(ESA_CCI_AGB)','C-Band_VOD', 'X-Band_VOD', 'VARI','NDVI','NIRv','LAI')
logagb.df.17.cor = cor(logagb.df.17, method = 'pearson', use = 'complete.obs')

logagb.corrplot = ggcorrplot(logagb.df.17.cor, hc.order = FALSE, type = 'upper',
                             outline.col = 'white',
                             lab = TRUE, show.legend = FALSE,
                             title = 'Pearson R Correlation',
                             colors = c('tomato','white','seagreen')) +
  # annotation_custom(grobTree(textGrob('B', x = 0.01, y = 0.95, hjust = 0,
  #                                     gp = gpar(col = 'black', fontsize = 20, fontface = 'bold')))) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 12),
        plot.background = element_rect(fill = 'white', color = 'white'))


# Plot AIC plots and correlation plots on the same figure
grid.arrange(linear.aic.plot, linear.corrplot, logagb.aic.plot, logagb.corrplot,
             ncol = 2, nrow = 2,
             widths = c(3,1.5),
             heights = c(0.5,0.5))


# ~~~~~ Final AIC and correlation plot figure
png(bg = 'white')
final.aic.corplot.fig = grid.arrange(logagb.aic.plot, logagb.corrplot,
                                     ncol = 2, nrow = 1, widths = c(3,1.5))

final.aic.corplot.fig = ggpubr::annotate_figure(final.aic.corplot.fig,
                                                fig.lab = 'a',
                                                fig.lab.pos = 'top.left',
                                                fig.lab.size = 20,
                                                fig.lab.face = 'bold')

final.aic.corplot.fig = ggpubr::annotate_figure(final.aic.corplot.fig,
                                                fig.lab = 'b',
                                                fig.lab.pos = 'top.right',
                                                fig.lab.size = 20,
                                                fig.lab.face = 'bold')

ggsave(paste(figs.fp, 'AIC','AIC_Results_Final.png', sep = '/'),
       final.aic.corplot.fig,
       width = 12, height = 4, units = 'in',
       bg = 'white')

# ------------------------------------------------------------------------------------------------------------------
# Variance inflation factor (VIF) for input variables of each model configuration containing more than one variable
# ------------------------------------------------------------------------------------------------------------------

# Read '00_Variance_Inflation_Factor_VIF_Fun.R' from Code/02_Analysis 
# Requires folliwing inputs (in order): reference AGB, VOD, VARI, NDVI, NIRv, and LAI
source(paste(code.fp, '02_Analysis', '00_Variance_Inflation_Factor_VIF_Fun.R', sep = '/'))

# Run VIF function
vif.fun(log.agb.17, c.vod.17, vari.17, ndvi.17, nirv.17, lai.17)
