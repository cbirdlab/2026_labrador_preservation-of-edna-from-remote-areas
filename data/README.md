# DATA: Effects of preservatives on the recovery of environmental DNA (eDNA) from seawater: Implications for field sampling in remote areas

## Kevin Labrador
---

# Introduction
This README contains the description of the datasets used in the manuscript.

<details>	
<summary> cnmi_map </summary>

## Description
This directory contains the shape files used to plot the Commonwealth of the Northern Mariana Islands (CNMI), which was used as inset in Figure 1.

</details>

<details>	
<summary> rota_shapefile </summary>

## Description
This directory contains the shape files used to plot the Rota Island in CNMI. This was used to generate Figure 1.

</details>

<details>	
<summary> df_for_stats.csv </summary>

## Description
The compiled data frame used for downstream statistical analyses.

### Column Headers

- A) filter_id: ID assigned to the Sterivex filters used throughout the experiment, following the convention <site_id>-<collection_date>-<preservative_code>
- B) site_id: a 3-letter code ID for the sites
- C) dna_extract_tube_id: ID assigned to the DNA extracted from the Sterivex filters
- D) sample_type: sample classification (sample, field control, extraction control, positive control, PCR control)
- E) preservative: the preservative treatment used in the experiment
- F) fraction: source of eDNA, either from Filter or Filtrate
- G) dna_concentration_ng_per_ul: concentration of eDNA (ng/uL) calculated from the Quant Module
- H) log10_dna_concentration: log10(dna_concentration_ng_per_ul)
- I) dna_concentration_in_normalized_plate_ng_per_ul: concentration of eDNA used for PCR
- J) pcr_plate_well_id: well designation of the DNA sample in the 96-well PCR plate
- K) pcr_gel_score: a binary score indicating the presence (1) or absence (0) of the target 16S rRNA amplicon in the agarose gel
- L) pcr_concentration_ng_per_ul: concentration (ng/uL) of the PCR products
- M) date_and_time_collection: date (yyyy-mm-dd) and time (hh:mm:ss) of eDNA collection
- N) volume_filtered: volume (mL) of seawater filtered through the Sterivex filter
- O) date_and_time_collection: date (yyyy-mm-dd) and time (hh:mm:ss) when sample collection ended
- P) duration_on_ice_min: how long (min) the samples were stored on ice in the field before storage in -20C
- Q) duration_on_ice_hours: duration_on_ice_min / 60
 
</details>

<details>	
<summary> quant_report_compiled.csv </summary>

## Description
Compilation of the quant reports from the [DNA Quantification Toolkit (Selwyn et al., 2025)](https://github.com/tamucc-gcl/gcl_bioinformatic_tools/blob/main/labTools/DNA_Quant/README.md) 

### Column Headers

- A) dna_plate_id: ID assigned to the DNA plate 
- B) dna_extract_tube_id: ID assigned to the individual DNA tubes
- C) sample_type: sample classification (sample, field control, extraction control)
- D) treatment: the preservative treatment used in the experiment
- E) dna_plate_row: the row location of the DNA sample in the DNA plate
- F) dna_plate_col: the column location of the DNA sample in the DNA plate
- G) site_code: a 3-letter code ID for the sites
- H) treatment_code: a 2-letter code used for the treatments
- I) sample_id: ID assigned to the Sterivex filters used throughout the experiment, following the convention <site_id>-<collection_date>-<preservative_code>
- J) ng_per_ul: DNA concentration (ng/uL) 
- K) ng_per_ul_lwr95: lower 95% confidence interval of the estimated DNA concentration
- L) ng_per_ul_upr95: upper 95% confidence interval of the estimated DNA concentration
- M) source: source of eDNA, either from Filter or Filtrate
- N) extraction_round: indicates which round of extraction (first, second) the eDNA was obtained
- O) extraction_strategy: <source>_<extraction_round>
- P) log10_dna_concentration: log10(DNA concentration) 
 
 
</details>


<details>	
<summary> rme_dna_2024-02-06_sample_concentration.csv </summary>

## Description
Concentration of DNA extracted from Sterivex filters using the [DNA Quantification Toolkit (Selwyn et al., 2025)](https://github.com/tamucc-gcl/gcl_bioinformatic_tools/blob/main/labTools/DNA_Quant/README.md)

### Column Headers

- A) quant_stage: indicates whether the quantification was done the first time (original) or was redone due to various problems (e.g., standards were erratic)
- B) dna_plate_id: ID assigned to the DNA plate 
- C) dna_extract_tube_ID: ID assigned to the individual DNA tubes
- D) sample_type: sample classification (sample, field control, extraction control)
- E) preservative: the preservative treatment used in the experiment
- F) dna_plate_row: the row location of the DNA sample in the DNA plate
- G) dna_plate_col: the column location of the DNA sample in the DNA plate
- H) sample_id: ID assigned to the Sterivex filters used throughout the experiment, following the convention <site_id>-<collection_date>-<preservative_code>
- I) ng_per_ul_mean: mean DNA concentration (ng/uL) across all replicates
- J) ng_per_ul_lwr95: lower 95% confidence interval of the estimated DNA concentration
- K) ng_per_ul_upr95: upper 95% confidence interval of the estimated DNA concentration 
- L) flags: notes/comments on samples based on their concentration
 
 
</details>


<details>	
<summary> rme_preservative_dna_2024-11-06_sample-concentration.csv </summary>

## Description
Concentration of DNA extracted from the preservative filtrate fraction using the [DNA Quantification Toolkit (Selwyn et al., 2025)](](https://github.com/tamucc-gcl/gcl_bioinformatic_tools/blob/main/labTools/DNA_Quant/README.md)

### Column Headers

- A) quant_stage: indicates whether the quantification was done the first time (original) or was redone due to various problems (e.g., standards were erratic)
- B) dna_plate_id: ID assigned to the DNA plate 
- C) dna_extract_tube_ID: ID assigned to the individual DNA tubes
- D) sample_type: sample classification (sample, field control, extraction control)
- E) preservative: the preservative treatment used in the experiment
- F) dna_plate_row: the row location of the DNA sample in the DNA plate
- G) dna_plate_col: the column location of the DNA sample in the DNA plate
- H) sample_id: ID assigned to the Sterivex filters used throughout the experiment, following the convention <site_id>-<collection_date>-<preservative_code>
- I) ng_per_ul_mean: mean DNA concentration (ng/uL) across all replicates
- J) ng_per_ul_lwr95: lower 95% confidence interval of the estimated DNA concentration
- K) ng_per_ul_upr95: upper 95% confidence interval of the estimated DNA concentration 
- L) flags: notes/comments on samples based on their concentration
 
 
</details>

<details>	
<summary> rme_pcr-16_accuclear_models_2025-07-02.rds </summary>

## Description
Model result from calculating the concentration of the PCR products using the DNA Quantification Toolkit (Selwyn et al., 2025). This is an RDS object that can be opened using R.

 
</details>

<details>	
<summary> rme_site-info.csv </summary>

## Description
Sampling site metadata

### Column Headers

- A) site: site name
- B) site-code: 3-letter site code
- C) resiliency-category: the relative resilience category of the site based on Maynard et al., 2015
- D) date_and_time_collection: date (yyyy-mm-dd) and time (hh:mm:ss) of eDNA collection
- E) weather-conditions: weather-conditions during sample collection (sunny, partly cloudy, intermittent rain)
- F) wave-strength: wave strength during sample collection
- G) depth_ft: depth (ft) of water column where surface seawater collection took place
- H) latitude: site coordinate (latitude) in decimal degrees
- I) longitude: site coordinate (longitude) in decimal degrees
  
</details>