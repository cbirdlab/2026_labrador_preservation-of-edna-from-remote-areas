#### wrangle_data-frame-for-statistical-tests ####
## Kevin Labrador
## 2025-06-02

####----INTRODUCTION----####
# This sets up the statistical tests for the quant data. 

####----INITIALIZE----####

#### HOUSEKEEPING ####
# Clear global environment.
rm(list = ls())

# Set-up the working directory in the source file location: 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Load Libraries ####
pacman::p_load(
  janitor,
  readxl,
  tidyverse
)

#### USER DEFINED VARIABLES ####

# Assign input file paths
path_quant_report_summarized <- 
  "../data/quant_report_compiled.csv"

path_normalized_dna_concentration_for_pcr <- 
  "../data/rme_plate_normalization_for_pcr.xlsx"

path_sampling_site_metadata <- 
  "../data/rme_site-info.xlsx"

path_pcr_scores <- 
  "../data/rme_pcr-score.xlsx"

path_pcr_quants <- 
  "../data/rme_pcr-16_sample_concentrations.csv"

path_sample_info <- 
  "../data/rme_sample-info.xlsx"

path_timeline <- 
  "../data/rme_timeline.xlsx"



# Assign output file paths
path_df_for_stats_csv <- 
  "../data/df_for_stats.csv"

path_df_for_stats_rds <- 
  "../data/df_for_stats.rds"



####----LOAD FILES----####

#### DNA concentration ####
df_dna_quants <- 
  path_quant_report_summarized %>% 
  read.csv() %>% 
  
  ## Modify names of negative extraction controls
  mutate (
    sample_id = 
      case_when (
        dna_extract_tube_id == "rme_15" ~ "nec-0124-nc01",
        dna_extract_tube_id == "rme_24" ~ "nec-0124-nc02",
        dna_extract_tube_id == "rme_35" ~ "nec-0124-nc03",
        dna_extract_tube_id == "rme_55" ~ "nec-0124-nc01",
        dna_extract_tube_id == "rme_71" ~ "nec-0124-nc02",
        dna_extract_tube_id == "rme_90" ~ "nec-0124-nc03",
        dna_extract_tube_id == "rmep_46" ~ "nec-1104-nc01",
        dna_extract_tube_id == "rmep_47" ~ "nec-1104-nc02",
        T ~ sample_id
      )
  ) %>% 
  
  ## Rename "dna_concentration_ng_per_ul"
  rename (
    dna_concentration_ng_per_ul = ng_per_ul,
    dna_concentration_ng_per_ul_upr95 = ng_per_ul_upr95,
    dna_concentration_ng_per_ul_lwr95 = ng_per_ul_lwr95,
  )


# Normalized DNA concentration for PCR
df_normalized_dna <- 
  path_normalized_dna_concentration_for_pcr %>% 
  read_xlsx() 


# PCR amplification score
df_pcr_scores <- 
  path_pcr_scores %>% 
  read_xlsx() %>%
  
  # Rename levels
  mutate (
    dna_extract_tube_id = tolower (dna_extract_tube_id),
    treatment = tolower (treatment)
  ) %>% 
  
  # Update site codes
  mutate (
    site_code = 
      case_when (
        sample_id == "C.tinkeri" ~ "cti",
        T ~ substr(site_code, 1, 3)
      )
  ) %>% 
  
  # Update factor levels
  mutate (
    treatment = 
      case_when (
        site_code == "nfc" ~ "field control",
        site_code == "nec" ~ "extraction control",
        site_code %in% c("npc", "cti") ~ "pcr control",
        T ~ treatment
      ) 
  ) 

# PCR Quants
df_pcr_quants <- 
  path_pcr_quants %>% 
  read_csv() %>% 
  ## Rename "dna_concentration_ng_per_ul"
  rename (
    pcr_concentration_ng_per_ul = ng_per_ul_mean,
    pcr_concentration_ng_per_ul_upr95 = ng_per_ul_upr95,
    pcr_concentration_ng_per_ul_lwr95 = ng_per_ul_lwr95,
  )

# Site Metadata
site_metadata <- 
  path_sampling_site_metadata %>% 
  read_xlsx() %>% 
  clean_names() %>% 
  mutate (
    site_code = as.factor(tolower(site_code)),
    resiliency_category = factor(
      resiliency_category,
      levels = c("high", 
                 "medium-high", 
                 "medium-low", 
                 "low")),
    weather_conditions = as.factor(weather_conditions),
    wave_strength = as.factor(wave_strength)
  ) 

# Sample Information
sample_info <- 
  path_sample_info %>% 
  read_xlsx() %>% 
  clean_names() %>% 
  # Change the sample_code to sample_id
  mutate (
    sample_id = 
      tolower(sample_code) 
  )

# Cold Chain Timeline 
timeline_data <- 
  path_timeline %>% 
  read_xlsx() %>% 
  clean_names()

####----PREPARE DATA FRAME FOR STATS----####
# Assign data frame
df_stats <- 
  df_dna_quants %>% 
  filter (extraction_round == "first") %>% 
  select (
    sample_id,
    site_code,
    dna_extract_tube_id,
    sample_type,
    treatment,
    source,
    dna_concentration_ng_per_ul,
    log10_dna_concentration) %>% 
  
  # Add information on normalized DNA concentrations
  left_join(
    df_normalized_dna %>% 
      select (
        dna_extract_tube_id,
        dna_concentration_in_normalized_plate_ng_per_ul
      )
  ) %>% 
  
  # Add information on PCR amplification scores. Use full join to include PCR controls
  full_join(
    df_pcr_scores %>% 
      select (
        dna_extract_tube_id,
        well_id,
        site_code,
        sample_id,
        sample_type,
        treatment,
        source,
        pcr_score
      )
  ) %>% 
  
  # Add information on PCR quants
  left_join(
    df_pcr_quants %>% 
      select (
        dna_extract_tube_id,
        pcr_concentration_ng_per_ul
      )
  ) %>% 
  
  # Add information on date and time of collection
  full_join (
    site_metadata %>% 
      select (
        site_code, 
        date_and_time_collection)
  ) %>% 
  
  # Add information on water volume filtered
  full_join (
    sample_info %>% 
      select (
        sample_id,
        volume_filtered_ml
      ) 
  )%>% 
  
  # Assign time for negative field controls
  mutate (
    date_and_time_collection = 
      case_when (
        grepl ("rme_14", dna_extract_tube_id) ~ as_datetime("2023-11-04 07:38:00"),
        grepl ("rme[p]?_(28|42)", dna_extract_tube_id) ~ as_datetime("2023-11-04 17:21:00"),
        T ~ date_and_time_collection
      )
    
  ) %>% 
  
  # Calculate duration on ice of samples during fieldwork
  mutate (
    date_and_time_end = 
      case_when(
        grepl("nec|npc|cti", site_code) ~ NA,
        T ~ timeline_data$date_end[1]),
    duration_on_ice_min = as.numeric(date_and_time_end - date_and_time_collection),
    duration_on_ice_hours = duration_on_ice_min / 60
  ) %>% 
  
  ## Change cases for treatments
  mutate (
    treatment =
      factor (
        case_when (
          treatment == "buffer tl" ~ "Buffer TL",
          treatment == "dess" ~ "DESS",
          treatment == "no preservative" ~ "No Pres",
          treatment == "extraction control" ~ "Ext Ctrl",
          treatment == "field control" ~ "Field Ctrl",
          treatment == "pcr control" ~ "PCR Ctrl",
          T ~ treatment),
        levels = c(
          "Ext Ctrl",
          "Field Ctrl",
          "PCR Ctrl",
          "No Pres",
          "DESS",
          "Buffer TL"
        )
      )
  ) %>% 
  
  ## Change cases for source
  mutate (
    source = 
      case_when(
        source == "pcr" ~ toupper(source),
        source == "preservative" ~ "Filtrate",
        T ~ str_to_sentence(source)
      )
  ) %>% 
  
  ## Change cases for sample type
  mutate (
    sample_type = 
      case_when (
        sample_type == "pcr control"  ~ "PCR Control",
        T ~ str_to_title(sample_type)
      )
  ) %>% 
  
  ## Change columns to factor
  mutate (
    site_code = as.factor (site_code),
    source = as.factor (source),
    sample_type = as.factor (sample_type)
  ) %>% 
  
  
  # Set Negative PCR Controls (millipore water) normalized DNA concentration to 0
  mutate (
    dna_concentration_in_normalized_plate_ng_per_ul = 
      case_when (
        grepl ("npc", sample_id) ~ 0,
        T ~ dna_concentration_in_normalized_plate_ng_per_ul
      )
  ) %>% 
  
  ## Rename columns
  rename (
    filter_id = sample_id,
    pcr_gel_score = pcr_score,
    pcr_plate_well_id = well_id,
    site_id = site_code,
    fraction = source,
    preservative = treatment
  )


####----SAVE FILES----####
write_csv(
  df_stats,
  file = path_df_for_stats_csv
)

write_rds(
  df_stats,
  file = path_df_for_stats_rds
)



