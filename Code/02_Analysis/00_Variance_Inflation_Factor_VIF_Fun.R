vif.fun = function(agb.in, vod.in, vari.in, ndvi.in, nirv.in, lai.in) {
  
  vod = vod.in[]
  agb = mask(agb.in, vod.in)[]
  vari = mask(vari.in, vod.in)[]
  ndvi = mask(ndvi.in, vod.in)[]
  nirv = mask(nirv.in, vod.in)[]
  lai = mask(lai.in, vod.in)[]
  
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
  mod18.vari.ndvi = lm(agb ~ vari + ndvi, na.action = na.exclude)
  mod19.vari.ndvi.nirv = lm(agb ~ vari + ndvi + nirv, na.action = na.exclude)
  mod20.vari.ndvi.nirv.lai = lm(agb ~ vari + ndvi + nirv + lai, na.action = na.exclude)
  mod21.vari.nirv = lm(agb ~ vari + nirv, na.action = na.exclude)
  mod22.vari.nirv.lai = lm(agb ~ vari + nirv + lai, na.action = na.exclude)
  mod23.vari.lai = lm(agb ~ vari + lai, na.action = na.exclude)
  
  mod25.ndvi.nirv = lm(agb ~ ndvi + nirv, na.action = na.exclude)
  mod26.ndvi.nirv.lai = lm(agb ~ ndvi + nirv + lai, na.action = na.exclude)
  
  mod28.nirv.lai = lm(agb ~ nirv + lai, na.action = na.exclude)
  
  mod30.ndvi.lai = lm(agb ~ ndvi + lai, na.action = na.exclude)
  mod31.vari.ndvi.lai = lm(agb ~ vari + ndvi + lai, na.action = na.exclude)
  
  
  # Combine all models into single list
  all.mods.list = list(mod02.vod.vari,
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
                       mod18.vari.ndvi,
                       mod19.vari.ndvi.nirv,
                       mod20.vari.ndvi.nirv.lai,
                       mod21.vari.nirv,
                       mod22.vari.nirv.lai,
                       mod23.vari.lai,
                       mod25.ndvi.nirv,
                       mod26.ndvi.nirv.lai,
                       mod28.nirv.lai,
                       mod30.ndvi.lai,
                       mod31.vari.ndvi.lai)
  
  # Create vector of model names corresponding to list
  all.mods.names = c('VOD+VARI',  #02
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
                     'VARI+NDVI',  #18
                     'VARI+NDVI+NIRv',  #19
                     'VARI+NDVI+NIRv+LAI',  #20
                     'VARI+NIRv',  #21
                     'VARI+NIRv+LAI',  #22
                     'VARI+LAI',  #23
                     'NDVI+NIRv',  #25
                     'NDVI+NIRv+LAI',  #26
                     'NIRv+LAI',  #28
                     'NDVI+LAI', #30
                     'VARI+NDVI+LAI' #31
  )
  
  # Generate AIC report
  all.mods.aic = aictab(all.mods.list, all.mods.names, second.ord = FALSE)
  
  # ------ Compute VIF for each model containing multiple variables
  mod02.vod.vari.vif = car::vif(mod02.vod.vari)
  mod03.vod.vari.ndvi.vif = car::vif(mod03.vod.vari.ndvi)
  mod04.vod.vari.ndvi.nirv.vif = car::vif(mod04.vod.vari.ndvi.nirv)
  mod05.vod.vari.ndvi.nirv.lai.vif = car::vif(mod05.vod.vari.ndvi.nirv.lai)
  mod06.vod.vari.ndvi.lai.vif = car::vif(mod06.vod.vari.ndvi.lai)
  mod07.vod.vari.nirv.vif = car::vif(mod07.vod.vari.nirv)
  mod08.vod.vari.nirv.lai.vif = car::vif(mod08.vod.vari.nirv.lai)
  mod09.vod.vari.lai.vif = car::vif(mod09.vod.vari.lai)
  mod10.vod.ndvi.vif = car::vif(mod10.vod.ndvi)
  mod11.vod.ndvi.nirv.vif = car::vif(mod11.vod.ndvi.nirv)
  mod12.vod.ndvi.lai.vif = car::vif(mod12.vod.ndvi.lai)
  mod13.vod.ndvi.nirv.lai.vif = car::vif(mod13.vod.ndvi.nirv.lai)
  mod14.vod.nirv.vif = car::vif(mod14.vod.nirv)
  mod15.vod.nirv.lai.vif = car::vif(mod15.vod.nirv.lai)
  mod16.vod.lai.vif = car::vif(mod16.vod.lai)
  mod18.vari.ndvi.vif = car::vif(mod18.vari.ndvi)
  mod19.vari.ndvi.nirv.vif = car::vif(mod19.vari.ndvi.nirv)
  mod20.vari.ndvi.nirv.lai.vif = car::vif(mod20.vari.ndvi.nirv.lai)
  mod21.vari.nirv.vif = car::vif(mod21.vari.nirv)
  mod22.vari.nirv.lai.vif = car::vif(mod22.vari.nirv.lai)
  mod23.vari.lai.vif = car::vif(mod23.vari.lai)
  mod25.ndvi.nirv.vif = car::vif(mod25.ndvi.nirv)
  mod26.ndvi.nirv.lai.vif = car::vif(mod26.ndvi.nirv.lai)
  mod28.nirv.lai.vif = car::vif(mod28.nirv.lai)
  mod30.ndvi.lai.vif = car::vif(mod30.ndvi.lai)
  mod31.vari.ndvi.lai.vif = car::vif(mod31.vari.ndvi.lai)
  
  # Create data frame with AIC results
  aic.df = setDT(all.mods.aic, keep.rownames = TRUE)[]
  aic.df = aic.df[,c(1,2,3,4,5)]
  aic.df$rn = sprintf('%02d', as.numeric(aic.df$rn))
  aic.df$K = sprintf('%02d', as.numeric(aic.df$K))
  
  # Create data frame with VIF results
  vif.df = rbind(data.frame('Model' = all.mods.names[1],
                            'VOD' = mod02.vod.vari.vif[1],
                            'VARI' = mod02.vod.vari.vif[2],
                            'NDVI' = NA,
                            'NIRv' = NA,
                            'LAI' = NA),
                 data.frame('Model' = all.mods.names[2],
                            'VOD' = mod03.vod.vari.ndvi.vif[1],
                            'VARI' = mod03.vod.vari.ndvi.vif[2],
                            'NDVI' = mod03.vod.vari.ndvi.vif[3],
                            'NIRv' = NA,
                            'LAI' = NA),
                 data.frame('Model' = all.mods.names[3],
                            'VOD' = mod04.vod.vari.ndvi.nirv.vif[1],
                            'VARI' = mod04.vod.vari.ndvi.nirv.vif[2],
                            'NDVI' = mod04.vod.vari.ndvi.nirv.vif[3],
                            'NIRv' = mod04.vod.vari.ndvi.nirv.vif[4],
                            'LAI' = NA),
                 data.frame('Model' = all.mods.names[4],
                            'VOD' = mod05.vod.vari.ndvi.nirv.lai.vif[1],
                            'VARI' = mod05.vod.vari.ndvi.nirv.lai.vif[2],
                            'NDVI' = mod05.vod.vari.ndvi.nirv.lai.vif[3],
                            'NIRv' = mod05.vod.vari.ndvi.nirv.lai.vif[4],
                            'LAI' = mod05.vod.vari.ndvi.nirv.lai.vif[5]),
                 data.frame('Model' = all.mods.names[5],
                            'VOD' = mod06.vod.vari.ndvi.lai.vif[1],
                            'VARI' = mod06.vod.vari.ndvi.lai.vif[2],
                            'NDVI' = mod06.vod.vari.ndvi.lai.vif[3],
                            'NIRv' = NA,
                            'LAI' = mod06.vod.vari.ndvi.lai.vif[4]),
                 data.frame('Model' = all.mods.names[6],
                            'VOD' = mod07.vod.vari.nirv.vif[1],
                            'VARI' = mod07.vod.vari.nirv.vif[2],
                            'NDVI' = NA,
                            'NIRv' = mod07.vod.vari.nirv.vif[3],
                            'LAI' = NA),
                 data.frame('Model' = all.mods.names[7],
                            'VOD' = mod08.vod.vari.nirv.lai.vif[1],
                            'VARI' = mod08.vod.vari.nirv.lai.vif[2],
                            'NDVI' = NA,
                            'NIRv' = mod08.vod.vari.nirv.lai.vif[3],
                            'LAI' = mod08.vod.vari.nirv.lai.vif[4]),
                 data.frame('Model' = all.mods.names[8],
                            'VOD' = mod09.vod.vari.lai.vif[1],
                            'VARI' = mod08.vod.vari.nirv.lai.vif[2],
                            'NDVI' = NA,
                            'NIRv' = NA,
                            'LAI' = mod08.vod.vari.nirv.lai.vif[3]),
                 data.frame('Model' = all.mods.names[9],
                            'VOD' = mod10.vod.ndvi.vif[1],
                            'VARI' = NA,
                            'NDVI' = mod10.vod.ndvi.vif[2],
                            'NIRv' = NA,
                            'LAI' = NA),
                 data.frame('Model' = all.mods.names[10],
                            'VOD' = mod11.vod.ndvi.nirv.vif[1],
                            'VARI' = NA,
                            'NDVI' = mod11.vod.ndvi.nirv.vif[2],
                            'NIRv' = mod11.vod.ndvi.nirv.vif[3],
                            'LAI' = NA),
                 data.frame('Model' = all.mods.names[11],
                            'VOD' = mod12.vod.ndvi.lai.vif[1],
                            'VARI' = NA,
                            'NDVI' = mod12.vod.ndvi.lai.vif[2],
                            'NIRv' = NA,
                            'LAI' = mod12.vod.ndvi.lai.vif[3]),
                 data.frame('Model' = all.mods.names[12],
                            'VOD' = mod13.vod.ndvi.nirv.lai.vif[1],
                            'VARI' = NA,
                            'NDVI' = mod13.vod.ndvi.nirv.lai.vif[2],
                            'NIRv' = mod13.vod.ndvi.nirv.lai.vif[3],
                            'LAI' = mod13.vod.ndvi.nirv.lai.vif[4]),
                 data.frame('Model' = all.mods.names[13],
                            'VOD' = mod14.vod.nirv.vif[1],
                            'VARI' = NA,
                            'NDVI' = NA,
                            'NIRv' = mod14.vod.nirv.vif[2],
                            'LAI' = NA),
                 data.frame('Model' = all.mods.names[14],
                            'VOD' = mod15.vod.nirv.lai.vif[1],
                            'VARI' = NA,
                            'NDVI' = NA,
                            'NIRv' = mod15.vod.nirv.lai.vif[2],
                            'LAI' = mod15.vod.nirv.lai.vif[3]),
                 data.frame('Model' = all.mods.names[15],
                            'VOD' = mod16.vod.lai.vif[1],
                            'VARI' = NA,
                            'NDVI' = NA,
                            'NIRv' = NA,
                            'LAI' = mod16.vod.lai.vif[2]),
                 data.frame('Model' = all.mods.names[16],
                            'VOD' = NA,
                            'VARI' = mod18.vari.ndvi.vif[1],
                            'NDVI' = mod18.vari.ndvi.vif[2],
                            'NIRv' = NA,
                            'LAI' = NA),
                 data.frame('Model' = all.mods.names[17],
                            'VOD' = NA,
                            'VARI' = mod19.vari.ndvi.nirv.vif[1],
                            'NDVI' = mod19.vari.ndvi.nirv.vif[2],
                            'NIRv' = mod19.vari.ndvi.nirv.vif[3],
                            'LAI' = NA),
                 data.frame('Model' = all.mods.names[26],
                            'VOD' = NA,
                            'VARI' = mod31.vari.ndvi.lai.vif[1],
                            'NDVI' = mod31.vari.ndvi.lai.vif[2],
                            'NIRv' = NA,
                            'LAI' = mod31.vari.ndvi.lai.vif[3]),
                 data.frame('Model' = all.mods.names[18],
                            'VOD' = NA,
                            'VARI' = mod20.vari.ndvi.nirv.lai.vif[1],
                            'NDVI' = mod20.vari.ndvi.nirv.lai.vif[2],
                            'NIRv' = mod20.vari.ndvi.nirv.lai.vif[3],
                            'LAI' = mod20.vari.ndvi.nirv.lai.vif[4]),
                 data.frame('Model' = all.mods.names[19],
                            'VOD' = NA,
                            'VARI' = mod21.vari.nirv.vif[1],
                            'NDVI' = NA,
                            'NIRv' = mod21.vari.nirv.vif[2],
                            'LAI' = NA),
                 data.frame('Model' = all.mods.names[20],
                            'VOD' = NA,
                            'VARI' = mod22.vari.nirv.lai.vif[1],
                            'NDVI' = NA,
                            'NIRv' = mod22.vari.nirv.lai.vif[2],
                            'LAI' = mod22.vari.nirv.lai.vif[3]),
                 data.frame('Model' = all.mods.names[21],
                            'VOD' = NA,
                            'VARI' = mod23.vari.lai.vif[1],
                            'NDVI' = NA,
                            'NIRv' = NA,
                            'LAI' = mod23.vari.lai.vif[2]),
                 data.frame('Model' = all.mods.names[22],
                            'VOD' = NA,
                            'VARI' = NA,
                            'NDVI' = mod25.ndvi.nirv.vif[1],
                            'NIRv' = mod25.ndvi.nirv.vif[2],
                            'LAI' = NA),
                 data.frame('Model' = all.mods.names[25],
                            'VOD' = NA,
                            'VARI' = NA,
                            'NDVI' = mod30.ndvi.lai.vif[1],
                            'NIRv' = NA,
                            'LAI' = mod30.ndvi.lai.vif[2]),
                 data.frame('Model' = all.mods.names[23],
                            'VOD' = NA,
                            'VARI' = NA,
                            'NDVI' = mod26.ndvi.nirv.lai.vif[1],
                            'NIRv' = mod26.ndvi.nirv.lai.vif[2],
                            'LAI' = mod26.ndvi.nirv.lai.vif[3]),
                 data.frame('Model' = all.mods.names[24],
                            'VOD' = NA,
                            'VARI' = NA,
                            'NDVI' = NA,
                            'NIRv' = mod28.nirv.lai.vif[1],
                            'LAI' = mod28.nirv.lai.vif[2])
  )
  
  vif.df = cbind(vif.df[,1], round(vif.df[,c(2:6)], 3))
  colnames(vif.df)[1] = 'Model'
  
  vif.df$Model = factor(vif.df$Model, 
                        levels = c('VOD+VARI',
                                   'VOD+NDVI',
                                   'VOD+NIRv',
                                   'VOD+LAI',
                                   'VARI+NDVI',
                                   'VARI+NIRv',
                                   'VARI+LAI',
                                   'NDVI+NIRv',
                                   'NDVI+LAI',
                                   'NIRv+LAI',
                                   'VOD+VARI+NDVI',
                                   'VOD+VARI+NIRv',
                                   'VOD+VARI+LAI',
                                   'VOD+NDVI+NIRv',
                                   'VOD+NDVI+LAI',
                                   'VOD+NIRv+LAI',
                                   'VARI+NDVI+NIRv',
                                   'VARI+NDVI+LAI',
                                   'VARI+NIRv+LAI',
                                   'NDVI+NIRv+LAI',
                                   'VOD+VARI+NDVI+NIRv',
                                   'VOD+VARI+NDVI+LAI',
                                   'VOD+VARI+NIRv+LAI',
                                   'VOD+NDVI+NIRv+LAI',
                                   'VARI+NDVI+NIRv+LAI',
                                   'VOD+VARI+NDVI+NIRv+LAI'))
  
  vif.df = vif.df[order(vif.df$Model),]
  
  vif.colors = scales::col_numeric(palette = c('wheat','tomato'),
                                   domain = c(min(vif.df[,c(2:6)],na.rm = TRUE),max(vif.df[,c(2:6)],na.rm = TRUE)),
                                   na.color = 'gray')
  
  
  vif.table = flextable(vif.df) %>%
    theme_vanilla() %>%
    align_text_col(align = 'center') %>%
    align_nottext_col(align = 'center') %>%
    set_table_properties(layout = 'autofit') %>%
    add_header_row(values = 'Variance Inflation Factor (VIF) for MLR model input variables',
                   colwidths = 6, top = TRUE) %>%
    bg(i = 1:26,
       j = c('VOD','VARI','NDVI','NIRv','LAI'),
       bg = vif.colors) %>%
    bg(i = 1:26,
       j = 'Model',
       bg = 'white') %>%
    bg(part = 'header', bg = 'white')
  
  
  vif.table
  
  save_as_image(vif.table,
                path = paste(figs.fp,'AIC','CbandVOD_logAGB_ModInputVars_VIF_table.png', sep = '/'))
  
  write.csv(vif.df,
            file = paste(tables.fp, 'Table_S2_MLR_VIF.csv', sep = '/'),
            row.names = FALSE)
}
