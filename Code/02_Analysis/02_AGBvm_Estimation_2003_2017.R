# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Using model configuration selected through AIC tests, retrieve AGB over extended ABoVE domain by combining VOD and MVIs for years 2003-2017
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# -----------------------
# Get libraries
# -----------------------
library(raster)
library(sf)
library(colorRamps)
library(ggplot2)
library(ggpointdensity)
library(ggpattern)
library(ggpubr)
library(grid)
library(gridExtra)

setwd('../../')
gitrepo.fp = getwd()
setwd('../../')
root.fp = getwd()
data.fp = paste(root.fp, 'Data', sep = '/')
figs.fp = paste(gitrepo.fp, 'Figures/Analysis', sep = '/')

# -----------------------
# ABoVE domain shapefiles
# -----------------------
shp.fp = paste(gitrepo.fp, 'Data', 'Shapefiles', sep = '/')
setwd(shp.fp)

core.domain = shapefile(paste(shp.fp, 'Domain','Core_Region_WGS.shp', sep = '/'))
entire.domain = shapefile(paste(shp.fp, 'Domain', 'Extended_Plus_Core_Domains_WGS.shp', sep = '/'))

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

# 2010
esa.agb.10.fp = paste(data.fp, 'Biomass/Santoro_ESA/2010', sep = '/')
setwd(esa.agb.10.fp)
esa.agb.10 = raster('GlobBiomass_2010_ABOVE_025deg.tif')
esa.agb.10[esa.agb.10 == 0] = NA

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For AGB Retrieval:
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

# -----------------------
# RMSE function
# -----------------------
rmse.fun = function(x.in, y.in){
  x = x.in[]
  y = mask(y.in, x.in)[]
  
  minus = x - y
  square = minus^2
  mean = mean(square, na.rm = TRUE)
  RMSE = sqrt(mean)
  return(round(RMSE,0))
}

diff.pal = colorRampPalette(c('orange','gray96','forestgreen'))
breakpoints = seq(-4,4,0.5)


# ---------------------------------------------------------------------------------------------
# AGB estimation (calibration 2017)
# ---------------------------------------------------------------------------------------------

# Exponential MLR function fitting VOD (or logVOD), VARI, and NIRv logAGB
agb.vod.vari.nirv.exp.fun = function(agb.in, vod.in, vod.name, vari.in, nirv.in, log.agb, log.vod, vi.names, xylim){
  
  vod = vod.in
  agb = mask(agb.in, vod)
  vari = mask(vari.in, vod)
  nirv = mask(nirv.in, vod)
  
  if (log.agb == TRUE && log.vod == FALSE){
    model = lm(log(agb[]) ~ vod[] + vari[] + nirv[], na.action = na.exclude)
    # coefs = exp(model$coefficients)
    # agb.vod = coefs[1] * (coefs[2]^vod) * (coefs[3]^vari) * (coefs[4]^nirv) * (coefs[5]^lai) 
    coefs = model$coefficients
    agb.vod = exp(coefs[1] + (coefs[2]*vod) + (coefs[3]*vari) + (coefs[4]*nirv))
    plot.text = 'Log trans: target AGB (dependent)'
  }
  
  if (log.agb == TRUE && log.vod == TRUE){
    model = lm(log(agb[]) ~ log(vod[]) + vari[] + nirv[], na.action = na.exclude)
    coefs = model$coefficients
    agb.vod = exp(coefs[1] + (coefs[2]*log(vod)) + (coefs[3]*vari) + (coefs[4]*nirv))
    plot.text = 'Log trans: target AGB (dependent) and VOD (independent) '
  }
  
  rmse = rmse.fun(agb, agb.vod)
  model.rsqrd = round(summary(model)$r.squared,2)
  retrieval = lm(agb[] ~ agb.vod[], na.action = na.exclude)
  retrieval.rsqrd = round(summary(retrieval)$r.squared, 2)
  retrieval.p.cor = round(cor(agb[], agb.vod[], method = 'pearson', use = 'complete.obs'), 2)
  
  par(mar = c(5,5,4,2),
      oma = c(0,0,0,0))
  smoothScatter(x = agb[], y = agb.vod[], xlim = xylim, ylim = xylim,
                ylab = 'Estimated AGB [Mg/ha]',
                xlab = 'ESA CCI AGB [Mg/ha]',
                main = paste('MLR components:',vod.name,vi.names), 
                cex.main = 1.5, cex.lab = 1.8, 
                colramp = colorRampPalette(c('white',blues9,'tomato'), space = 'Lab'))
  abline(1,1, lwd = 1, lty = 2, col = 'grey')
  abline(retrieval, lwd = 2, lty = 2, col = 'tomato')
  legend('topleft',
         legend = plot.text, cex = 1.5, bty = 'n')
  legend('bottomright',
         legend = c(paste('RMSE =', rmse, 'Mg/ha'),
                    as.expression(bquote(R^2 == .(retrieval.rsqrd)))),
         cex = 1.8, bty = 'n')
  
  #print(coefs)
  print(summary(model))
  return(coefs)
}


# Get coefficients
b = agb.vod.vari.nirv.exp.fun(esa.agb.17, c.vod.17, 'C-Band VOD (VODCA),', vari.17, nirv.17,
                              log.agb = TRUE, log.vod = FALSE,
                              vi.names = c('VARI, NIRv'),
                              c(0,260))

# ---------------------------------------------------------------------------------------------
# Compute AGB_VM 2003-2017
# Export to disc
# ---------------------------------------------------------------------------------------------

AGB_VM = exp(b[1] + (b[2]*c.vod.ym) + (b[3]*vari) + (b[4]*nirv))

# Export AGB_VM data as .tif files
AGB_VM.fp = paste(gitrepo.fp, 'Data/AGB_VM', sep = '/')
AGB_VM.outfiles = paste(AGB_VM.fp, paste(seq(2003,2017,1),'AGB_VM.tif',sep = '_'), sep = '/')

writeRaster(AGB_VM,
            filename = AGB_VM.outfiles,
            format = 'GTiff',
            bylayer = TRUE,
            overwrite = TRUE)


# ---------------------------------------------------------------------------------------------
# Plot 2017/2010 calibration/evaluation scatterplot and map figures
# ---------------------------------------------------------------------------------------------

cal.plot.fun = function(){
  
  # Combine AGB_VM and ESA AGB maps for each year as data frames
  agb.17.df = na.omit(data.frame(AGB_VM = values(AGB_VM[[15]]),
                                 ESA = values(esa.agb.17),
                                 Year = rep('2017 (calibration)', ncell(esa.agb.17))))
  agb.17.df$RMSE = rep(rmse.fun(esa.agb.17, AGB_VM[[15]]), nrow(agb.17.df))
  
  agb.10.df = na.omit(data.frame(AGB_VM = values(AGB_VM[[8]]),
                                 ESA = values(esa.agb.10),
                                 Year = rep('2010 (evaluation)', ncell(esa.agb.10))))
  agb.10.df$RMSE = rep(rmse.fun(esa.agb.10, AGB_VM[[8]]), nrow(agb.10.df))
  
  agb.df = rbind(agb.17.df, agb.10.df)
  agb.df$Year = factor(agb.df$Year, levels = c('2017 (calibration)','2010 (evaluation)'))
  
  # Create scatterplot
  sp.labs = data.frame(Year = c('2017 (calibration)','2010 (evaluation)'), label = c('a','b'))
  sp.labs$Year = factor(sp.labs$Year, levels = c('2017 (calibration)','2010 (evaluation)'))
  
  agb.sp = ggplot(agb.df, aes(x = ESA, y = AGB_VM)) +
    geom_pointdensity(adjust = 50, size = 1) +
    scale_x_continuous(limits = c(0,260)) +
    scale_y_continuous(limits = c(0,260)) +
    scale_color_gradient2(low = 'white', mid = blues9, high = 'tomato', space = 'Lab') +
    stat_smooth(formula = y ~ x, method = 'lm', se = FALSE, color = 'tomato', linetype = 'dashed', fullrange = TRUE) +
    geom_abline(slope = 1, col = 'grey') +
    xlab(expression(AGB[ESA-CCI] ~ '[Mg/ha]')) +
    ylab(expression(AGB[VM] ~ '[Mg/ha]')) +
    theme_bw() +
    stat_regline_equation(aes(label = ..rr.label..), formula = y~x, size = 6, fontface = 'bold') +
    geom_text(aes(x = 180, y = 10, label = paste('RMSE =',RMSE,'[Mg/ha]'), family = 'sans', fontface = 'plain'), size = 5.5) +
    facet_grid(.~Year) +
    geom_text(x = 260, y = 260, data = sp.labs, aes(label = label), size = 8, fontface = 'bold') +
    theme(legend.position = 'none',
          strip.background = element_rect(fill = "grey39"),
          strip.text = element_text(color = 'white', size = 18),
          axis.text = element_text(size = 14))
  
  # Combine AGB_VM-ESA difference maps for each year as data frames
  AGB_VM.agb.17.diff = AGB_VM[[15]] - esa.agb.17
  AGB_VM.agb.17.diff = as(AGB_VM.agb.17.diff, 'SpatialPixelsDataFrame')
  AGB_VM.agb.17.diff = as.data.frame(AGB_VM.agb.17.diff)
  colnames(AGB_VM.agb.17.diff) = c('AGB_Diff','x','y')
  AGB_VM.agb.17.diff$Year = rep('2017', nrow(AGB_VM.agb.17.diff))
  
  AGB_VM.agb.10.diff = AGB_VM[[8]] - esa.agb.10
  AGB_VM.agb.10.diff = as(AGB_VM.agb.10.diff, 'SpatialPixelsDataFrame')
  AGB_VM.agb.10.diff = as.data.frame(AGB_VM.agb.10.diff)
  colnames(AGB_VM.agb.10.diff) = c('AGB_Diff','x','y')
  AGB_VM.agb.10.diff$Year = rep('2010', nrow(AGB_VM.agb.10.diff))
  
  AGB_VM.agb.diff = rbind(AGB_VM.agb.17.diff, AGB_VM.agb.10.diff)
  AGB_VM.agb.diff$Year = factor(AGB_VM.agb.diff$Year, levels = c('2017','2010'))
  
  diff.pal = colorRampPalette(c('red','orange','khaki1','gray96','palegreen','darkgreen','navy'))
  
  # Plot difference maps
  map.labs = data.frame(Year = c('2017','2010'), label = c('c','d'))
  map.labs$Year = factor(map.labs$Year, levels = c('2017','2010'))
  
  agb.map = ggplot(data = AGB_VM.agb.diff, aes(x=x, y=y, fill = AGB_Diff)) +
    geom_raster() + 
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    scale_fill_gradient2(low = 'tomato', mid = 'ivory', high = 'forestgreen', midpoint = 0, limits = c(-200,200),
                         name = 'AGB diff. [Mg/ha]') +
    ggtitle(expression(AGB[VM] ~ 'minus' ~ AGB[ESA-CCI])) +
    theme_test() + 
    geom_sf(data = st_as_sf(entire.domain), color = 'black', fill = NA, inherit.aes = FALSE, show.legend = 'line') +
    geom_sf(data = st_as_sf(core.domain), color = 'tomato', fill = NA, inherit.aes = FALSE) +
    facet_grid(.~Year) +
    geom_text(x = -165, y = 78, data = map.labs, aes(label = label), inherit.aes = FALSE, size = 8, fontface = 'bold') +
    guides(fill = guide_colorbar(title.position = 'bottom', title.hjust = 0.5, frame.colour = 'black')) +
    theme(plot.title = element_text(size = 15, hjust = 0.5),
          legend.position = 'bottom',
          legend.key.width = unit(3.5, 'cm'),
          legend.title = element_text(size = 15),
          axis.title = element_blank(),
          #strip.background = element_rect(fill = "grey39"),
          strip.background = element_blank(),
          #strip.text = element_text(color = 'white', size = 18),
          strip.text = element_blank())
  
  #final.fig = grid.arrange(agb.sp, agb.map, ncol = 1, nrow = 2)
  final.fig = ggarrange(agb.sp, agb.map, ncol = 1, nrow = 2, align = 'v', heights = c(0.8,1)) + bgcolor('white')
}

cal.plot.17.10 = cal.plot.fun()

ggsave(filename = paste(figs.fp, 'AGB_VM_CalEval.png', sep = '/'),
       cal.plot.17.10,
       device = 'png',
       width = 9, height = 9, units = 'in')
