library(plyr)
library(tidyverse)
library(scales) 
library(ggExtra)
library(ggpubr)
library(patchwork)
library(viridis)
library(gridExtra)
library(ggpmisc)
library(maps)
library(ggthemes)
library(mapproj)
library(ggrepel)
library(sf)
library(rnaturalearth)


library(MASS)
library(Metrics)
library(doParallel)

source(file = "SCRIPTS/Match-up_with_SOMLIT/Functions.R") 

registerDoParallel(cores = detectCores()/2)

#### Define path ####

work_dir <- "/home/terrats/Desktop/RIOMAR/"
path_to_SOMLIT_data_with_MU_results <- file.path(work_dir, "DATA/IN_SITU/")
path_to_river_data <- file.path(work_dir, "DATA/RIVERS/FOR_MAPS/FILES/")
where_to_save_MU_results <- file.path(work_dir, "RESULTS/MATCH_UP/")

#### Configuration ####

max_CV <- 30 # 30
min_n <- NA # 5
Maximum_time_difference = NA # 3 # in hour
vars_to_plot <- c("SPM-G", "CHLA") # c("SST-DAY", 'SST-NIGHT') 
# stations_to_remove <- c('Sola', 'Bouee 13', 'Eyrac', 'Comprian')
stations_to_remove <- c()
Version_best_spatio_temporal_grid <- 1 # V1 show good results
SOMLIT_QC_codes_to_keep <- c(2,6,7)

#### Load data ####

df <- load_SOMLIT_data(path_to_SOMLIT_data_with_MU_results, stations_to_remove)

color_values = c("Point L" = mako(n = 1,begin = 0.8,end = 0.8), # Manche
                 "Point C" = mako(n = 1,begin = 0.85,end = 0.85), 
                 "Luc-sur-Mer" = mako(n = 1,begin = 0.90,end = 0.9), 
                 "Smile" = mako(n = 1,begin = 0.95,end = 0.95),
                 
                 "Bizeux" = viridis(n = 1,begin = 1,end = 1), # Bretagne
                 "Le Buron" = viridis(n = 1,begin = 0.925,end = 0.925), 
                 "Cézembre" = viridis(n = 1,begin = 0.85,end = 0.85), 
                 "Estacade" = viridis(n = 1,begin = 0.775,end = 0.775), 
                 "Astan" = viridis(n = 1,begin = 0.7,end = 0.7), 
                 "Portzic" = viridis(n = 1,begin = 0.625,end = 0.625), 
                 
                 "Antioche" = plasma(n = 1,begin = 0.05,end = 0.05), # Golfe de Gascogne
                 "pk 86" = plasma(n = 1,begin = 0.1,end = 0.1), 
                 "pk 52" = plasma(n = 1,begin = 0.15,end = 0.15),
                 "pk 30" = plasma(n = 1,begin = 0.2,end = 0.2),
                 "Comprian" = plasma(n = 1,begin = 0.25,end = 0.25), 
                 "Eyrac" = plasma(n = 1,begin = 0.3,end = 0.3), 
                 "Bouee 13" = plasma(n = 1,begin = 0.35,end = 0.35), 
                 
                 "Sola" = rocket(n = 1,begin = 0.80,end = 0.80), # Golfe du Lion
                 "Sete" = rocket(n = 1,begin = 0.85,end = 0.85), 
                 "Frioul" = rocket(n = 1,begin = 0.90,end = 0.90),
                 "Point B" = rocket(n = 1,begin = 0.95,end = 0.95)
) 

SOMLIT_stations <- df %>% distinct(Latitude, Longitude, Site)

# world_map = map_data("world") %>% filter(! long > 180)
# # France_boundaries <- ne_countries(country = "France", scale = 'medium', type = 'map_units', returnclass = 'sf')  
# France_boundaries <- read_sf(file.path(work_dir, "DATA", "FRANCE_shapefile", "gadm41_FRA_0.shp"))  
# riversData <- read_sf(file.path(path_to_river_data, "HydroRIVERS_v10_eu.shp")) 
# riversData <- riversData %>% st_intersection(France_boundaries)
# # Tuto here : https://www.etiennebacher.com/posts/2021-12-27-mapping-french-rivers-network/
# # 
# SOMLIT_station_map <- ggplot() +
#   # geom_polygon(data = world_map, aes(x=long, y = lat, group = group), fill = "black") +
#   geom_sf(data=France_boundaries, color="black", inherit.aes = FALSE) +
#   geom_sf(data=riversData, color="white", inherit.aes = FALSE) +
#   # geom_sf(data=riversData, color="white", inherit.aes = FALSE) +
#   # expand_limits(x = c(-8,11), y = c(52, 40)) +
#   # coord_map("moll") +
#   coord_sf(xlim = c(-8,11), ylim = c(40,52)) +
#   geom_point(data = SOMLIT_stations, aes(x = Longitude, y = Latitude), color = "red", size = 3) +
#   geom_text_repel(data = SOMLIT_stations, aes(x = Longitude, y = Latitude, label = Site),
#                   color = "red", size = 5, max.overlaps = 20) +
#   labs(title = paste(nrow(SOMLIT_stations), "SOMLIT stations"), x = "", y = "") +
#   scale_x_continuous(label = function(x) {paste(x, "°E", sep = "")}) +
#   scale_y_continuous(label = function(x) {paste(x, "°N", sep = "")}) +
#   theme_map() +
#   ggplot_theme() +
#   theme(plot.title = element_text(color = "red", size = 30),
#         axis.text = element_text(size = 20))
# # 
# # save_plot_as_png(SOMLIT_station_map, name = "SOMLIT_station_map", width = 16, height = 15, path = where_to_save_MU_results)

sat_products_to_process <- names(df) %>% 
  str_subset("_value$") %>% 
  str_replace("_value", "") %>% 
  str_replace("CHLA_|SPM-G_|SPM-R_|SST-DAY_|SST-NIGHT_", "") %>% 
  str_replace("_[0-9]x[0-9]", "") %>% 
  unique() %>% 
  sort()

valid_MU_criteria_for_title = paste("Match-up criteria : Coefficient of Variation ≤ ", max_CV, "% / n ≥ ", min_n, ' / time difference ≤ ', Maximum_time_difference, "h", sep = "")

statistics_and_scatterplots <- sat_products_to_process %>% llply(function(sat_product_to_process) {
  
  min_n_of_the_product <- ifelse(str_detect(sat_product_to_process, pattern = "SEXTANT_merged|polymer") & is.finite(min_n), 1, min_n)
  
  SPM_CHLA_plots <- vars_to_plot %>% llply(function(variable) {
    
    insitu_variable <- case_when(variable %in% c("SPM-G", "SPM-R") ~ "SPM", .default = variable) %>% str_replace("-DAY|-NIGHT", "")
    sat_product_full_name <- paste(variable, sat_product_to_process, sep = "_")
    
    if ( any(str_detect(names(df), sat_product_full_name)) == FALSE ) { return() }
    
    df_formatted <- SOMLIT_stations$Site %>% llply(function(Site_name) {
      
      df_Site <- get_sat_data_usign_the_optimal_grid(df = df, Site_name = Site_name, 
                                                     Version_optimized_grid = Version_best_spatio_temporal_grid,
                                                     sat_product_full_name = sat_product_full_name)
      
      return(df_Site)
      
    }, .inform = TRUE) %>% rbind.fill()
    
    df_of_the_sat_product <- df_formatted %>% 
      
      dplyr::rename(insitu_value := {{insitu_variable}},
                    insitu_qc := !!sym(paste("q", insitu_variable, sep = ""))) %>% 
      
      # mutate(sat_value = case_when( grepl("polymer", sat_product_to_process) ~ sat_value, 
      #                               .default = sat_median)) %>% # 1x1 grid only for polymer atm correction
      mutate(sat_value = sat_median) %>% # 1x1 grid only for polymer atm correction
      
      filter_at(vars(all_of(c("insitu_value", "sat_value"))), ~ is.finite(.)) %>% 
      filter(insitu_qc %in% SOMLIT_QC_codes_to_keep) %>% 
      
      mutate(sat_CV = 100 * sat_sd / sat_value,
             Year = lubridate::year(DATE),
             Month = lubridate::month(DATE))  
    
    if (max_CV %>% is.finite()) {df_of_the_sat_product <- df_of_the_sat_product %>% filter(sat_CV <= max_CV) }  
    if (min_n_of_the_product %>% is.finite()) {df_of_the_sat_product <- df_of_the_sat_product %>% filter(sat_n >= min_n_of_the_product) }
    
    if (nrow(df_of_the_sat_product) == 0) { return() }
    
    df_with_only_positive_values <- df_of_the_sat_product %>% filter_at(vars(ends_with("value")), ~ . > 0)
    
    statistics_df <- compute_stats(sat_values = df_with_only_positive_values$sat_value, 
                                   insitu_values = df_with_only_positive_values$insitu_value)
    
    if (str_detect(variable, 'SST')) {
      statistics_df['bias_Pahlevan'] = statistics_df['bias_Pahlevan_linear']
      statistics_df['error_Pahlevan'] = statistics_df['error_Pahlevan_linear']
      statistics_df['slope_log'] = statistics_df['slope']
    }
    
    Figures <- make_the_figures(df_with_only_positive_values, var_to_assess = variable, statistics_outputs = statistics_df)
    
    # if (variable %>% str_detect("SST")) { 
    save_plot_as_png(Figures$scatterplot_with_side_histograms, name = paste(variable, "V", Version_best_spatio_temporal_grid, sep = "_"), 
                     path = file.path(where_to_save_MU_results, "SCATTERPLOTS", "Per_sat_product"), 
                     width = 18, height = 14)
    # }
    
    return( list("scatterplot" = Figures$scatterplot, 
                 "scatterplot_with_side_hist" = Figures$scatterplot_with_side_hist, 
                 "freq_per_year" = Figures$bar_plot_freq_per_year,
                 "freq_per_month" = Figures$bar_plot_freq_per_month,
                 "statistics" = tibble(Sat_product = sat_product_to_process,
                                       Variable = variable,
                                       statistics_df %>% purrr::discard_at(vars(starts_with("lm_line"))) %>% data.frame())) )
    
  }, .inform = TRUE) %>% setNames(vars_to_plot)
  
  if (SPM_CHLA_plots %>% llply(is.null) %>% unlist() %>% all()) {
    return()
  }
  
  return( list("statistics" = ldply(SPM_CHLA_plots, function(x) x$statistics, .id = NULL),
               "scatterplots" = SPM_CHLA_plots %>% llply(function(x) {x$scatterplot}) %>% setNames(names(SPM_CHLA_plots))) )
  
}, .inform = TRUE, .parallel = TRUE) %>% 
  set_names(sat_products_to_process)



vars_to_plot %>% l_ply(function(var) {
  
  SEXTANT_plots <- c("SEXTANT_meris", "SEXTANT_modis", "SEXTANT_seawifs", "SEXTANT_viirs", "SEXTANT_merged") %>% paste("_Standard", sep = "") %>% 
    llply(function(x) {
      if (var %in% names(statistics_and_scatterplots[[x]]$scatterplots) == FALSE) {
        return()
      }
      statistics_and_scatterplots[[x]]$scatterplots[[var]] + 
        labs(title = str_replace_all(x, "SEXTANT|_|Standard", replacement = "")) + 
        guides(color=guide_legend(ncol=1), 
               linetype = guide_legend(ncol = 1,
                                       override.aes = list(color = c("black"),
                                                           shape = c(NA),
                                                           linetype = c("dashed")))) + 
        theme(legend.box.background = element_rect(color = "white"), plot.title = element_text(size = 50))}) 
  
  SEXTANT_plots <- ggarrange(plotlist = SEXTANT_plots, ncol = 4, nrow = 2, align = 'hv', common.legend = TRUE, legend = "right")
  
  SEXTANT_plots <- annotate_figure(
    annotate_figure(SEXTANT_plots, 
                    top=text_grob(valid_MU_criteria_for_title %>% paste(" -  For the merged product : n ≥ 1 / time difference ≤ 1d"), 
                                  face = "italic", size = 35 * 1.5, color = "black")),
    top=text_grob("SEXTANT", face = "bold", size = 45*1.5, color = "black"),
  )
  
  ODATIS_plots <- c("acolite", "polymer", "nirswir") %>% 
    
    llply(function(atm_correction) {
      
      sensor_plots <- c("meris", "modis", "olcia", "olcib") %>% llply(function(sensor_name) {
        
        product_name <- paste("ODATIS", sensor_name, atm_correction, sep = "_")
        
        if (product_name %in% names(statistics_and_scatterplots) == FALSE) {
          return(NULL)
        }
        
        if (var %in% names(statistics_and_scatterplots[[product_name]]$scatterplots) == FALSE) {
          return()
        }
        
        statistics_and_scatterplots[[product_name]]$scatterplots[[var]] + 
          labs(title = str_replace_all(product_name, "ODATIS|_", replacement = " ")) + 
          guides(color=guide_legend(ncol=1), 
                 linetype = guide_legend(ncol = 1,
                                         override.aes = list(color = c("black"),
                                                             shape = c(NA),
                                                             linetype = c("dashed")))) + 
          theme(legend.box.background = element_rect(color = "white"), plot.title = element_text(size = 40))  
        
      })
      
      return(sensor_plots)
      
    }) 
  
  ODATIS_plots <- ggarrange(plotlist = ODATIS_plots %>% flatten(), 
                            ncol = 4, nrow = 3, 
                            align = 'hv', common.legend = TRUE, legend = "right")
  
  ODATIS_plots <- annotate_figure(
    annotate_figure(ODATIS_plots, 
                    top=text_grob(valid_MU_criteria_for_title, face = "italic", size = 50, color = "black")),
    top=text_grob("ODATIS", face = "bold", size = 70, color = "black"),
  )
  
  final_plots = ggarrange(ODATIS_plots, SEXTANT_plots, ncol = 1, nrow = 2, align = "hv", common.legend = TRUE, heights = c(0.15, 0.1))
  
  save_plot_as_png(final_plots, name = paste("Product_comparison", var, "V", Version_best_spatio_temporal_grid, sep = "_"), 
                   path = file.path(where_to_save_MU_results, "SCATTERPLOTS"),
                   width = 60, height = 60, res = 100)
  
}, .inform = TRUE, .parallel = TRUE)
