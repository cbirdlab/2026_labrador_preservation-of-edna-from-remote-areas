#### wrangle_quant_results ####
## Kevin Labrador
## 2025-06-02

####----INTRODUCTION----####
# Compiles all the quant reports from the [DNA Quantification Tools](https://github.com/tamucc-gcl/gcl_bioinformatic_tools/blob/main/labTools/DNA_Quant/README.md)

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
  readxl,
  tidyverse
)

#### USER DEFINED VARIABLES ####

# Assign input file paths
path_quant_filter <-
  "../data/processed/jselwyn_quant_modules/quant_module_2/rme_dna_2024-02-06_sample-concentration.csv"


path_quant_preservative<- 
  "../data/processed/jselwyn_quant_modules/quant_module_2/rme_preservative_dna_2024-11-06_sample-concentration.csv"


# Assign output file paths
path_quant_report_compiled <- 
  "../data/processed/quant_report_compiled.csv"

####----WRANGLE QUANT RESULTS----####

# Quant data
quant_data_filter <- 
  path_quant_filter %>% 
  read_csv() %>% 
  mutate (
    source = "filter"
  )

quant_data_preservative <- 
  path_quant_preservative %>%
  read_csv() %>% 
  mutate (
    source = "preservative"
  ) %>% 
  # Drop old ng columns first
  select(
    -contains("ng")
  ) %>%
  # Convert pg to ng and rename
  mutate(
    across(contains("pg"), ~ .x / 1000, 
           .names = "converted_{.col}")
  ) %>%
  rename_with(
    ~ str_replace(.x, "converted_pg", "ng"), 
    starts_with("converted_pg")
  ) %>% 
  select (- contains ("pg"))


quant_data_compiled <- 
  full_join(
    quant_data_filter,
    quant_data_preservative
  ) %>% 
  
  select (
    dna_plate_id,
    dna_extract_tube_id,
    sample_type,
    preservative,
    dna_plate_row,
    dna_plate_col,
    sample_id,
    ng_per_ul_mean,
    ng_per_ul_lwr95,
    ng_per_ul_upr95,
    source
  ) %>% 
  
  # Modify data frame
  
  ## Separate sample id to respective codes
  separate_wider_delim(
    cols = sample_id,
    delim = "-",
    names = c("site_code", NA, "preservative_code"),
    cols_remove = F
    
  ) %>% 
  
  ## Rename columns
  rename (
    
    ### "preservative" to "treatment"
    treatment = preservative,
    treatment_code = preservative_code,
    
    ### "ng_per_ul_mean" to "ng_per_ul"
    ng_per_ul = ng_per_ul_mean
    
  ) %>% 
  
  # Modify columns
  mutate (
    
    ## Add extraction_round column
    extraction_round = 
      case_when(
        str_detect(dna_extract_tube_id, "^rmep_") ~ "first",
        str_detect(dna_extract_tube_id, "^rme_") & 
          as.numeric(stringr::str_extract(dna_extract_tube_id, "\\d+")) <= 45 ~ "first",
        TRUE ~ "second"
      ),
    
    ## Change treatment code from cold storage "cs" to no preservative "np"
    treatment_code =
      case_when (
        grepl("cs", treatment_code) ~ "np",
        T ~ treatment_code
      ),
    
    ## Modify treatments based on treatment code
    treatment = 
      case_when(
        site_code == "nec" ~ "extraction control",
        site_code == "nfc" ~ "field control",
        site_code == "npc" ~ "pcr control",
        T ~ treatment
      ),
    
    ## Merge extraction round and source into a single column
    extraction_strategy = paste(source, extraction_round, sep = "_"),
    
    ## Calculate log of DNA concentration
    log10_dna_concentration = log10(ng_per_ul)
  )


####----SAVE FILES----####

# Save quant data
write.csv(
  quant_data_compiled,
  file = path_quant_report_compiled,
  row.names = F
)
