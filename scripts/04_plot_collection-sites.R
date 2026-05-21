#### plot_collection-sites ####
## Kevin Labrador
## 2024-09-07

####-----INTRODUCTION-----####
# Project collection sites from a shape file

####-----INITIALIZE-----####

#### Housekeeping ####
# Clear global environment.
rm(list = ls())

# Set-up the working directory in the source file location: 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Load Libraries ####
pacman::p_load(
  tidyverse,
  janitor,
  readxl,
  scales,
  zoo,
  ggmap,
  ggspatial,
  ggpubr,
  sf,
  cowplot
)


#### USER-DEFINED VARIABLES ####

# Assign input file paths
path_sf_rota <- 
    "../data/raw/rota_shapefile/"
  

path_sf_cnmi <- 
    "../data/raw/cnmi_map/shape files/" %>% 
  list.files(
    pattern = ".zip",
    full.names = T
  )

path_site_info <- 
    "../data/raw/rme_site-info.xlsx"

# Assign output files
path_plot_sampling_site <- 
    "../results/plot_sampling_site.svg"

path_site_info_for_supp_table <-
  "../results/rme_site-info.csv"
  
  

####----LOAD FILES----####

df_site_info <- 
  path_site_info %>% 
  read_excel() %>% 
  clean_names() %>% 
  mutate (
    
    # Reorder factors
    resiliency_category = 
      factor (
        resiliency_category,
        levels = c(
          "high", 
          "medium-high",
          "medium-low", 
          "low")
      )
  ) %>% 
  
  # Add numeric index
  mutate (index = 1:nrow(.))


sf_rota <- 
  path_sf_rota %>% 
  st_read() %>% 
  # Convert CRS from UTM Zone 55N to long/lat coords
  st_transform(crs = 4326)

sf_cnmi <- 
  map(
    path_sf_cnmi, 
    function(zip_path) {
      # Construct the path to the .shp file inside the zip
      shp_name <- tools::file_path_sans_ext(basename(zip_path)) # e.g., "agrihan_shoreline"
      shp_path <- paste0("/vsizip/", zip_path, "/", shp_name, ".shp")
      
      # Read the shapefile
      st_read(shp_path) %>% 
        
        # Transform to appropriate coordinate system
        st_transform(crs = 4326)
    }
  ) %>% 
  reduce(rbind)


####-----PLOT SAMPLING LOCATIONS-----####
# Prepare data frames for plotting

## Extract coordinates from the dataset
coords <- 
  df_site_info %>% 
  select (
    index,
    site,
    site_code,
    latitude,
    longitude,
    resiliency_category
  ) 

## Prepare bounding box
padding <- 0.2
bbox <- 
  data.frame(
    xmin = min(coords$longitude) - padding,
    xmax = max(coords$longitude) + padding,
    ymin = min(coords$latitude) - padding,
    ymax = max(coords$latitude) + padding
  )

# Plot Basemap
(basemap <- 
    ggplot () + 
    
    # Prepare basemap
    geom_sf (
      data = sf_rota,
      fill = "#EDDDDD",
      col = "black"
    ) +
    
    # Add the coordinates
    geom_point(
      data = coords,
      aes (x = longitude,
           y = latitude),
      col = "black",
      alpha = 0.50,
      size = 3.5
    ) +
    
    # Add labels
    geom_text (
      data = coords,
      aes(x=longitude,
          y=latitude,
          label = index),
      nudge_x = 0.0035,
      nudge_y = 0.0035,
      size = 5) +

    
    # Set the theme
    theme_classic() + 
    theme (
      panel.border = 
        element_rect(
          color = "black",
          fill = NA
        ),
      axis.text = 
        element_text(
          size = 10
      )
    ) + 
    
    # Add the annotations
    annotation_scale (
      location = "bl", 
      bar_cols = c(
        "black", 
        "white")
    ) +
    annotation_north_arrow (
      location = "tl",
      which_north="true",
      style = north_arrow_fancy_orienteering
    ) +
    labs (
      x = "Longitude",
      y = "Latitude"
    )
) 


(plot_inset <- 
    ggplot () + 
    
    # Prepare basemap
    geom_sf (
      data = sf_cnmi,
      fill = "#EDDDDD",
      col = "black"
    ) +
    
    # Add bounding box
    geom_rect(
      data = bbox,
      aes (
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      col = "black",
      fill = "red3",
      alpha = 0.25
    ) +
    
    scale_x_continuous(
      breaks = 145, 
      limits = c(144, 146)
    ) + 
    scale_y_continuous(
      breaks = seq(14, 22, by = 2)
    ) + 
    
    # Set the theme
    theme_pubr() + 
    theme (
      panel.border = 
        element_rect(
          color = "black",
          fill = NA
        ),
      plot.background = 
        element_rect(
          color = "black"
        ),
      axis.text = 
        element_text(
          size = 10
        )
    )
)


(p <- 
    ggdraw() +
    draw_plot(basemap) +  # Main plot
    draw_plot(plot_inset, 
              x = 0.77, y = 0.20,   
              width = 0.30, height = 0.30)
)


####-----SAVE FILES-----#### 
# Save output
ggsave(
  p, 
  file = path_plot_sampling_site,
  width=7.25, 
  height=5, 
  units = "in", 
  dpi=600
)

df_site_info %>% 
  select (
    index,
    site,
    longitude,
    latitude
  ) %>% 
  write.csv(
    file = path_site_info_for_supp_table,
    row.names = F
  )

