list_of_packages <- c("plyr", "tidyverse", "maps", "scales", "ggpubr")
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list_of_packages, require, character.only = TRUE)

#### Main ####
# where_to_save_the_figure <- "/home/terrats/Desktop/RIOMAR/TEST"

Figure_1 <- function(where_to_save_the_figure) {
  
  main_folder_of_Figure_1 <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURES", "FIGURE_1")

  SPM_map <- file.path( main_folder_of_Figure_1, "DATA", "SPM_map.csv" ) %>% read_csv()
  insitu_stations <- file.path( main_folder_of_Figure_1, "DATA", "Stations_position.csv" ) %>% read_csv()
  RIOMAR_limits <- file.path( main_folder_of_Figure_1, "DATA", "RIOMAR_limits.csv" ) %>% read_csv() %>% 
                    pivot_longer(-...1, names_to = "Zone", values_to = "Value") %>%
                    pivot_wider(names_from = ...1, values_from = Value)

  basic_map <- create_the_basic_map(map_df = SPM_map, var_name = 'SPM', in_situ_fixed_station = insitu_stations)
  
  points_for_the_legend <- data.frame(SOURCE = c('SOMLIT', 'REPHY'),
                                      longitude = c(0,0),
                                      latitude = c(0,0))
  
  Figure_1 <- basic_map + 
      geom_point(data = insitu_stations %>% filter(SOURCE == 'REPHY'), 
                 aes(x = LONGITUDE, y = LATITUDE), 
                 fill = "red", color = "black", size = 4, shape = 24, stroke = 1) +
      geom_point(data = insitu_stations %>% filter(SOURCE == 'SOMLIT'), 
                 aes(x = LONGITUDE, y = LATITUDE), 
                 fill = "red", color = "black", size = 10, shape = 21, stroke = 2) + 
      geom_rect(data = RIOMAR_limits, aes(xmin = lon_min, xmax = lon_max, ymin = lat_min, ymax = lat_max),
                fill = "transparent", color = "red", linetype = "dashed", size = 2) +
      
      geom_point(data = points_for_the_legend, aes(x = longitude, y = latitude, shape = SOURCE), size = 0.1) +
      
      scale_shape_manual(values = c('SOMLIT' = 21, "REPHY" = 24), breaks=c('SOMLIT', 'REPHY'), 
                         labels = c(paste('SOMLIT (n=', length(which(insitu_stations$SOURCE == "SOMLIT")), ")", sep = ""),
                                    paste('REPHY (n=', length(which(insitu_stations$SOURCE == "REPHY")), ")", sep = ""))) +
      guides(
        shape = guide_legend(keyheight = unit(0.3, "cm"), byrow = TRUE,
                             override.aes = list(size = c(10, 4),
                                                 shape = c(21,24),
                                                 fill = c("red", "red"),
                                                 color = c('black', 'black'),
                                                 stroke = c(2, 1)),
                             order = 1),
        fill = guide_colorbar(barwidth = 30, barheight = 2)) +
      labs(shape = "In-situ stations") + 
      theme(legend.position = "bottom",
            legend.title.position = "top",
            legend.title = element_text(angle = 0, hjust = 0.5),
            legend.spacing.x = unit(5, "cm"))
    
  save_plot_as_png(Figure_1, "Figure_1", width = 18, height = 16, path = main_folder_of_Figure_1)
  
}

# where_to_save_the_figure <- '/home/terrats/Desktop/RIOMAR/TEST/'
Figure_2 <- function(where_to_save_the_figure) {
  
  main_folder_of_Figure_2 <- file.path(where_to_save_the_figure, "ARTICLE", "FIGURES", "FIGURE_2")
  
  SPM_map_data <- file.path( main_folder_of_Figure_2, "DATA" ) %>% 
                list.files(pattern = "*.csv", full.names = TRUE) %>% 
                llply(read_csv) %>% 
                keep(~ 'analysed_spim' %in% names(.))
  
  insitu_stations <- file.path( main_folder_of_Figure_2, "DATA", "Stations_position.csv" ) %>% read_csv()

  points_for_the_legend <- data.frame(SOURCE = c('SOMLIT', 'REPHY'), longitude = c(0,0), latitude = c(0,0))
  
  SPM_maps <- SPM_map_data %>% 
                llply(function(x) {
                  
                  insitu_stations_of_the_map <- insitu_stations %>% 
                                                  filter((LATITUDE %>% between(min(x$lat), max(x$lat))) & 
                                                         (LONGITUDE %>% between(min(x$lon), max(x$lon))))
                  
                  the_map <- create_the_basic_map(x, 'SPM', legend_limits = c(3,10))
                  the_map <- the_map + 
                    
                                geom_point(data = insitu_stations_of_the_map %>% filter(SOURCE == 'REPHY'), 
                                           aes(x = LONGITUDE, y = LATITUDE), 
                                           fill = "red", color = "black", size = 6, shape = 24, stroke = 1) +
                                geom_point(data = insitu_stations_of_the_map %>% filter(SOURCE == 'SOMLIT'), 
                                           aes(x = LONGITUDE, y = LATITUDE), 
                                           fill = "red", color = "black", size = 14, shape = 21, stroke = 2) + 
                    
                                geom_point(data = points_for_the_legend, aes(x = longitude, y = latitude, shape = SOURCE), size = 0.1) +
                    
                                scale_shape_manual(values = c('SOMLIT' = 21, "REPHY" = 24), breaks=c('SOMLIT', 'REPHY'), 
                                                   labels = c('SOMLIT', 'REPHY')) +
                                guides(
                                  shape = guide_legend(keyheight = unit(0.3, "cm"), byrow = TRUE,
                                                       override.aes = list(size = c(14, 6),
                                                                           shape = c(21,24),
                                                                           fill = c("red", "red"),
                                                                           color = c('black', 'black'),
                                                                           stroke = c(2, 1)),
                                                       order = 1),
                                  fill = guide_colorbar(barwidth = 45, barheight = 2)) +
                                labs(shape = "In-situ stations") + 
                                theme(legend.position = "bottom",
                                      legend.title.position = "top",
                                      legend.title = element_text(angle = 0, hjust = 0.5),
                                      legend.spacing.x = unit(5, "cm"),
                                      axis.text = element_text(size=25, colour = "black"))

                    return(the_map)
                  }) %>% 
                ggarrange(plotlist = ., common.legend = TRUE)
  
  save_plot_as_png(SPM_maps, "Figure_2", width = 28, height = 16, path = main_folder_of_Figure_2)
  
}


# where_to_save_the_figure <- '/home/terrats/Desktop/RIOMAR/TEST/ARTICLE/FIGURES/FIGURE_4'
# name_of_the_plot <- "C"
Figure_4 <- function(where_to_save_the_figure, name_of_the_plot) {

  SPM_map_data <- read_csv(file.path( where_to_save_the_figure, "DATA", paste(name_of_the_plot, ".csv", sep = "") ))
  
  if (name_of_the_plot %in% c("A", "B")) {
    the_map <- create_the_basic_map(SPM_map_data, 'SPM', legend_limits = c(1,5))
  } else {
    the_map <- create_the_basic_map(SPM_map_data, 'plume', legend_limits = c(1,5))
  }
    
  the_map <- the_map + 
              guides(fill = guide_colorbar(barwidth = 60, barheight = 2, title.position = "top")) +
              theme(legend.position = "top",
                    legend.title = element_text(angle = 0, hjust = 0.5),
                    axis.text = element_text(size=25, colour = "black"))
  
  if (name_of_the_plot == "B") {
    points_used_for_finding_SPM_threshold <- read_csv(file.path( where_to_save_the_figure, "DATA", "B_points_used_for_finding_SPM_threshold.csv" ))
    the_map <- the_map +
      geom_point(data = points_used_for_finding_SPM_threshold, aes(x = longitude, y = latitude), color = "red")
  }
  
  save_plot_as_png(the_map, name_of_the_plot, width = 28, height = 16, path = where_to_save_the_figure)
  
}


# where_to_save_the_figure <- '/home/terrats/Desktop/RIOMAR/TEST/ARTICLE/FIGURES/FIGURE_5'
Figure_5 <- function(where_to_save_the_figure) {
  
  SPM_map_data <- where_to_save_the_figure %>% file.path('DATA') %>% 
                    list.files(pattern = "*.csv", full.names = TRUE) %>% 
                    llply(read_csv)
  
  SPM_maps <- SPM_map_data %>% llply(function(SPM_map) {
    
    create_the_basic_map(SPM_map, 'plume', legend_limits = c(1,5)) + 
      guides(fill = guide_colorbar(barwidth = 60, barheight = 2, title.position = "top")) +
      theme(legend.position = "top",
            legend.title = element_text(angle = 0, hjust = 0.5),
            axis.text = element_text(size=25, colour = "black"))
    
  })
  
  save_plot_as_png(ggarrange(plotlist = SPM_maps, common.legend = TRUE), 
                   'Figure_5', width = 28, height = 16, path = where_to_save_the_figure)
  
}



#### Utils ####

create_the_basic_map <- function(map_df, var_name, 
                           in_situ_fixed_station = NULL, cruise_stations = NULL, 
                           glider_stations = NULL,
                           legend_limits = NULL) {
  
  if (str_detect(var_name, 'chl|CHL')) {
    title = "[Chl-a]"
    unit = "mg m³"
    if (legend_limits %>% is.null()) {legend_limits <- c(1e-1, 5e0)} 
  }
  
  if (str_detect(var_name, 'tsm|SPM|TSM|plume')) {
    title = "[SPM]"
    unit = "g m³"
    # legend_limits <- map_df$analysed_spim[which(map_df$plume)] %>% quantile(probs = c(0.1, 0.9), na.rm = TRUE)
    if (legend_limits %>% is.null()) {legend_limits <- c(1e-1, 5e0)} 
  }
  
  FRANCE_shapefile <- map_data('world')[map_data('world')$region == "France",]

  the_base_map <- ggplot() + 
    geom_raster(data = map_df, aes(x = lon, y = lat, fill = analysed_spim), interpolate = FALSE) + 
    scale_fill_viridis_c(na.value = "transparent", option = "viridis", trans = "log10", 
                         limits = c(legend_limits[1], legend_limits[2]), oob = scales::squish, 
                         n.breaks = 5, name = paste(title, " (", unit, ")", sep = "")) +
    guides(fill = guide_colourbar(title.position = "right"))
    
  if (var_name == 'plume') {
    the_base_map <- the_base_map + geom_raster(data = map_df[which(map_df$plume),], aes(x = lon, y = lat), fill = "red", interpolate = FALSE) 
  }
  
  the_map <- the_base_map + 
    
    ## First layer: worldwide map
    geom_polygon(data = map_data("world"), aes(x=long, y=lat, group = group), color = 'grey60', fill = 'black') +
    ## Second layer: Country map
    geom_polygon(data = FRANCE_shapefile, aes(x=long, y=lat, group = group), color = 'grey60', fill = 'black') +
    coord_cartesian(xlim = range(map_df$lon), ylim = range(map_df$lat), expand = FALSE) +
    
    scale_x_continuous(name = "", labels = function(x) paste(x, "°E", sep = "")) +
    scale_y_continuous(name = "", labels = function(x) paste(x, "°N", sep = "")) +
    ggplot_theme() + 
    
    theme(plot.title = element_text(size = 45),
          legend.position = "right",
          legend.title = element_text(angle = -90, hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.key.height = unit(6, "lines"),
          legend.key.width = unit(3, "lines")) 
  
  return(the_map)
  
}



ggplot_theme <-   function() {
  theme(text = element_text(size=35, colour = "black"), #25
        plot.title = element_text(hjust = 0.5, size = 55),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.text = element_text(size = 35, colour = "black"),
        axis.title = element_text(size = 40, colour = "black"),
        axis.text.x=element_text(angle=0),
        axis.ticks.length=unit(.25, "cm"))}


save_plot_as_png <- function(plot, name = c(), width = 14, height = 8.27, path, res = 150) {
  
  graphics.off()
  
  if (name %>% length() == 1) {
    if (dir.exists(path) == FALSE) {dir.create(path, recursive = TRUE)}
    path <- file.path(path, paste(name, ".png", sep = ""))
  } else {
    path <- paste(path, ".png", sep = "")
  }
  
  if (grepl(pattern = ".png.png", path)) {path <- path %>% gsub(pattern = ".png.png", replacement = ".png", x = .)}
  
  png(path, width = width, height = height, units = "in", res = res)
  print(plot)
  dev.off()
  
}
