#### plot_timeline ####
## Kevin Labrador
## 2025-01-31

####----INTRODUCTION----####


####----INITIALIZE----####

#### HOUSEKEEPING ####
# Clear global environment.
rm(list = ls())

# Set-up the working directory in the source file location: 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Load Libraries ####
pacman::p_load(
  janitor,
  ggpubr,
  vistime,
  readxl,
  tidyverse
)

#### USER DEFINED VARIABLES ####

# Assign input file paths
path_timeline <- 
  "../data/raw/rme_timeline.xlsx"

# Assign output file paths
path_timeline_plot <- 
    "../results/plot_timeline.png"

#### Load Files ####
timeline_data <- 
  path_timeline %>% 
  read_xlsx() %>% 
  mutate (color = c(
    "#b3cde0",
    "#00a8e8",
    "#0077b6",
    "#003f5c")
  )


#### Plot Timeline ####

(plot_timeline <- 
   gg_vistime(
     timeline_data,
     col.event = "storage_condition",
     col.start = "date_start",
     col.end = "date_end",
     col.group = "storage_condition",
     show_labels = F,
     title = "Preservation Timeline"
   ) +
   theme_classic() +
   scale_x_datetime(
     date_breaks = "1 month", 
     date_labels = "%b %Y" # <-- This line adds year ticks
   )  
)



# Export plot
ggsave(
  plot_timeline,
  file = path_timeline_plot,
  width = 5,
  height = 3.5,
  units = "in",
  dpi = 330)
