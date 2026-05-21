#### model_amplification-yield ####
## Kevin Labrador
## 2025-10-02

####----INTRODUCTION----####
# This sets up the statistical tests for the amplification success.
# The PCR concentration (ng/uL) were calculated using JSelwyn's quant app.
## From the Quant App's first module, we determined the limit of detection (LoD; the lowest concentration of the quant standard that was used in the assay)
## If [PCR] < LoD, then failed amplification; otherwise, successful amplification.
## I find this to be more objective than using gel scores (can be very subjective)
# Amplification success was scored as binary (0/1). 

####----INITIALIZE----####

#### HOUSEKEEPING ####
# Clear global environment.
rm(list = ls())

# Set-up the working directory in the source file location: 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Load Libraries ####
pacman::p_load(
  
  ## Used in final version of the script
  janitor,
  ggpubr,
  readxl,
  emmeans,
  multcomp,
  emmeans,
  DHARMa,
  brms,
  cmdstanr,
  ggbeeswarm,
  fitdistrplus,
  glmmTMB,
  sjPlot,
  ggeffects,
  ggpattern,
  ggheatmap,
  performance,
  tidyverse
)

#### USER DEFINED VARIABLES ####

# Assign input file paths
path_df_stats <- 
  "../data/df_for_stats.rds"

path_model_pcr_quant <- 
  "../data/rme_pcr-16_accuclear_models_2025-07-02.rds"

# Assign output file paths
path_plot_pcr_probability <- 
  "../results/plot_pcr-probability.png"

path_plot_quant_model <- 
  "../results/plot_pcr-quant-model.pdf"

path_plot_log_model <- 
  "../results/plot_logistic-regression_pcr-amplification.png"

path_plot_log_model_supplementary <- 
  "../results/plot_logistic-regression_pcr-amplification_buffer-tl-only.png"

path_model_coef_amplicon_concentration <- 
  "../results/model_coef_amplicon-concentration.csv"

# Set theme for ggplots
theme_set(theme_bw())

####----LOAD FILES----####
df_stats <- 
  path_df_stats %>% 
  readRDS() %>% 
  mutate (
    treatment_modified_controls = 
      case_when (
        grepl ("Extraction Control|Field Control", sample_type) ~ "Neg Ctrl",
        grepl ("C.tinkeri", filter_id) ~ "Pos Ctrl",
        grepl ("npc", filter_id) ~ "Neg Ctrl",
        T ~ preservative
      ) %>% 
      factor (levels = c(
        "Neg Ctrl",
        "Pos Ctrl",
        "No Pres",
        "DESS",
        "Buffer TL"
      ))
  ) %>%
  
  # Create a new column that merges treatment and fraction into a single variable
  mutate (
    treatment =
      paste (
        fraction,
        treatment_modified_controls,
        sep = "_"
      ) %>%
      factor (levels = c(
        "PCR_Neg Ctrl",
        "PCR_Pos Ctrl",
        "Filter_Neg Ctrl",
        "Filter_No Pres",
        "Filter_DESS",
        "Filter_Buffer TL",
        "Filtrate_Neg Ctrl",
        "Filtrate_DESS",
        "Filtrate_Buffer TL"
      ))
  )


# Open the regression model used in JSelwyn's Quant App
(quant_standard_model <- 
    path_model_pcr_quant %>% 
    readRDS() %>% 
    # Pull model from the list
    (function (x) x$model) %>% 
    
    # Back transform log10 values to original scale
    mutate (
      ng_per_well = 10^`log10(ng_per_well)`, # This is what was used in the model
      ng_per_ul = ng_per_well / 5, #5 = volume of standards used during quant
      rfu = 10^`log10(rfu)`
    ) 
)

# Determine the limit of detection (LoD) from the back-transformed values; the first row is the limit
limit_of_detection <- 
  quant_standard_model %>% 
  head (1) %>% 
  pull (ng_per_well)

print (paste("Limit of detection:", limit_of_detection, "ng/uL"))


####----SCORE PCR USING QUANT DATA----####
# Prepare data frame to score PCR based on quant data
df_pcr_score <-  
  df_stats %>%
  select (
    site_id,
    filter_id,
    pcr_plate_well_id,
    fraction,
    preservative,
    treatment_modified_controls,
    treatment,
    dna_concentration_in_normalized_plate_ng_per_ul,
    pcr_concentration_ng_per_ul,
    pcr_gel_score
  ) %>% 
  
  # Rename columns
  rename (
    dna_concentration = dna_concentration_in_normalized_plate_ng_per_ul,
    pcr_concentration = pcr_concentration_ng_per_ul
  ) %>% 
  
  # Reorder levels
  mutate (
    fraction =
      factor (
        fraction,
        levels = 
          c("Filter",
            "Filtrate",
            "PCR")
      )
  ) %>% 
  
  # Score amplification based on PCR quants
  mutate (
    pcr_quant_score = 
      case_when (
        pcr_concentration < limit_of_detection ~ 0,
        T ~ 1
      )
  )

# Check group sizes
df_pcr_score %>% 
  count(
    pcr_gel_score, 
    pcr_quant_score
  ) %>%   
  complete(
    pcr_gel_score, 
    pcr_quant_score, 
    fill = list(n = 0)
  ) %>% 
  ggplot(aes(
    x = pcr_gel_score %>% as.factor(),
    y = pcr_quant_score %>% as.factor(),
    fill = n
  )) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "dodgerblue") +
  geom_text(aes(label = n), color = "black") +   # <-- add values
  theme_minimal() +
  labs (
    x = "PCR Gel Score",
    y = "PCR Quant Score",
    subtitle = "Concordance of PCR scores based on gel and quant"
  )



####----VISUALIZE PCR SCORES----####

#### PCR quant score ####
df_quant_score_profile <-  
  df_pcr_score %>% 
  select (
    filter_id,
    preservative, 
    fraction,
    pcr_quant_score,
  ) %>% 
  filter (!grepl ("Ctrl", preservative)) %>% 
  pivot_wider(
    names_from = fraction,
    values_from = pcr_quant_score,
    values_fill = 0
  ) %>% 
  mutate (
    group = 
      factor(
        case_when(
          Filter == 1 & Filtrate == 1 ~ "+/+",
          Filter == 1 & Filtrate == 0 ~ "+/-",
          Filter == 0 & Filtrate == 1 ~ "-/+",
          Filter == 0 & Filtrate == 0 ~ "-/-"
        ),
        levels = c(
          "-/-",
          "+/+",
          "+/-",
          "-/+"
        )
      )
  ) 


(plot_pcr_quant_score <- 
    df_quant_score_profile %>% 
    ggplot (
      aes (
        x = preservative
        # fill = group
      )
    ) + 
    geom_bar_pattern(
      aes(
        pattern = group,
        pattern_angle = group,
        group = group
      ),
      position = "fill",
      fill = "white",
      color = "black",
      linewidth = 1,
      pattern_fill = "gray50",
      pattern_color = "white",
      pattern_density = 0.25,
      pattern_spacing = 0.05
    ) + 
    labs (
      x = NULL,
      y = "Relative Frequency (%)",
      subtitle = "Amplification Profile\nBased on Limit of Detection (LoD)"
    ) +
    theme_classic() +
    scale_y_continuous(
      labels = scales::percent_format()
    ) +
    scale_pattern_manual(
      name = "Amplification profile\n(Filter/Filtrate)",
      values = c(
        "-/-" = "none",
        "+/+" = "weave",
        "+/-" = "stripe",
        "-/+" = "stripe"
      ),
      drop = FALSE
    ) +
    scale_pattern_angle_manual(
      name = "Amplification profile\n(Filter/Filtrate)",
      values = c(
        "-/-" = 0,
        "+/+" = 45,
        "+/-" = 135,
        "-/+" = 45
      ),
      drop = FALSE
    )
)


#### PCR gel score ####
df_gel_score_profile <- 
  df_pcr_score %>% 
  select (
    filter_id,
    preservative, 
    fraction,
    pcr_gel_score,
  ) %>% 
  filter (!grepl ("Ctrl", preservative)) %>% 
  pivot_wider(
    names_from = fraction,
    values_from = pcr_gel_score,
    values_fill = 0
  ) %>% 
  mutate (
    group = 
      factor(
        case_when(
          Filter == 1 & Filtrate == 1 ~ "+/+",
          Filter == 1 & Filtrate == 0 ~ "+/-",
          Filter == 0 & Filtrate == 1 ~ "-/+",
          Filter == 0 & Filtrate == 0 ~ "-/-"
        ),
        levels = c(
          "-/-",
          "+/+",
          "+/-",
          "-/+"
        )
      )
  )

(plot_pcr_gel_score <- 
    df_gel_score_profile %>% 
    ggplot (
      aes (
        x = preservative,
        fill = group
      )
    ) + 
    geom_bar_pattern(
      aes(
        pattern = group,
        pattern_angle = group,
        group = group
      ),
      position = "fill",
      fill = "white",
      color = "black",
      linewidth = 1,
      pattern_fill = "gray50",
      pattern_color = "white",
      pattern_density = 0.25,
      pattern_spacing = 0.05
    ) + 
    labs (
      x = NULL,
      y = "Relative Frequency (%)",
      subtitle = "Amplification Profile\nBased on Presence of Target Bands on Gel"
    ) +
    theme_classic() +
    scale_y_continuous(
      labels = scales::percent_format()
    ) + 
    scale_pattern_manual(
      name = "Amplification profile\n(Filter/Filtrate)",
      values = c(
        "-/-" = "none",
        "+/+" = "weave",
        "+/-" = "stripe",
        "-/+" = "stripe"
      ),
      drop = FALSE
    ) +
    scale_pattern_angle_manual(
      name = "Amplification profile\n(Filter/Filtrate)",
      values = c(
        "-/-" = 0,
        "+/+" = 45,
        "+/-" = 135,
        "-/+" = 45
      ),
      drop = FALSE
    )
)



####---- TEST AMP FREQ BASED ON PCR AMPLIFICATION PROFILE ----####

#### Fisher's exact test of association between categorical vars: Gel Scores ####
(df_gel_amp_frequency <- 
   df_gel_score_profile %>% 
   droplevels() %>% 
   count(preservative, group) %>% 
   pivot_wider (
     names_from = group,
     values_from = n,
     values_fill = 0
   ) %>% 
   column_to_rownames("preservative") %>%
   as.matrix()
)
# Fisher's exact test of independence
fisher.test(df_gel_amp_frequency)


#### Fisher's exact test of association between categorical vars: Quant Scores ####
(df_quant_amp_frequency <- 
   df_quant_score_profile %>% 
   droplevels() %>% 
   count(preservative, group) %>% 
   pivot_wider (
     names_from = group,
     values_from = n,
     values_fill = 0
   ) %>% 
   column_to_rownames("preservative") %>%
   as.matrix()
)

# Fisher's exact test of independence
fisher.test(df_quant_amp_frequency)


# Fisher's exact test of independence for Buffer TL (TL Filter > TL Filtrate)
(df_gel_amp_frequency <- 
    df_gel_score_profile %>% 
    droplevels() %>% 
    count(preservative, group) %>% 
    pivot_wider (
      names_from = group,
      values_from = n,
      values_fill = 0
    ) %>% 
    column_to_rownames("preservative") %>%
    as.matrix()
)


#### LOGISTIC REGRESSION ON PCR QUANTS AND AMP ####
df_for_logreg <- 
  df_pcr_score %>% 
  select (treatment,
          site_id, 
          filter_id,
          preservative,
          fraction,
          pcr_concentration,
          pcr_gel_score) %>% 
  filter (!grepl ("Ctrl", treatment))



model_1 <- 
  glmmTMB (
    pcr_gel_score ~ pcr_concentration,
    data = df_for_logreg,
    family = binomial (link = "logit")
  )

model_2 <- 
  glmmTMB (
    pcr_gel_score ~ pcr_concentration + treatment,
    data = df_for_logreg,
    family = binomial (link = "logit")
  )

model_3 <- 
  glmmTMB (
    pcr_gel_score ~ pcr_concentration * treatment,
    data = df_for_logreg,
    family = binomial (link = "logit")
  )

model_4 <- 
  glmmTMB (
    pcr_gel_score ~ pcr_concentration + fraction + preservative,
    data = df_for_logreg,
    family = binomial (link = "logit")
  )

model_5 <- 
  glmmTMB (
    pcr_gel_score ~ pcr_concentration + preservative * fraction,
    data = df_for_logreg,
    family = binomial (link = "logit")
  )


model_list <- 
  list (
    model_1 = model_1,
    model_2 = model_2,
    model_3 = model_3,
    model_4 = model_4,
    model_5 = model_5
  )


(model_summary <- 
    lapply (
      model_list,
      summary
    )
)

# Test models
test_performance (model_list)

(model_selection <- 
    compare_performance(
      model_list
    )
)

# Check models 
plot(model_selection)

best_model <- model_1
summary(best_model)
check_model(best_model)

# Get p-values for best model
null_model <- 
  glmmTMB (
    pcr_gel_score ~ 1,
    data = df_for_logreg,
    family = binomial (link = "logit")
  )

(anova_out <- (anova(null_model, best_model)))

# Quick reference plot
model_prediction <- 
  predict_response(
    best_model,
    terms = c(
      "pcr_concentration [0:4.7, by=0.01]"
    )
  )


(plot_log_model <- 
    plot (
      model_prediction,
      show_data = T,
      show_ci = T
    ) + 
    labs (
      x = "[Amplicon] (ng/uL)",
      y = "Probability of Amplification",
      title = NULL,
      col = "Treatment",
      subtitle = paste0(
        "GLMM (Binomial, link = logit):\n",
        "gel_score ~ [amplicon]"
      )
    )+ 
    geom_vline(
      xintercept = limit_of_detection,
      lty = 2,
      col = "red"
    ) +
    theme_bw() 
  # scale_color_manual(
  #   values = c(
  #     #"gray50",
  #     "#D95F02",
  #     "#7570B3",
  #     "#E7298A"
  #     #"#1B9E77",
  #   )
  # )
)

#### SUPPLEMENTARY LOGISTIC REGRESSION LOOKING FOR BUFFER TL ####
supplementary_model <- 
  glmmTMB (
    pcr_gel_score ~ pcr_concentration + fraction,
    data = df_for_logreg %>%
      filter (
        preservative == "Buffer TL"
      ),
    family = binomial (link = "logit")
  )

model_prediction_supplementary <- 
  predict_response(
    supplementary_model,
    terms = c(
      "pcr_concentration [0:4.7, by=0.01]",
      "fraction"
    )
  )


(plot_log_model_supplementary <- 
    plot (
      model_prediction_supplementary,
      show_data = T,
      show_ci = T
    ) + 
    labs (
      x = "[Amplicon](ng/uL)",
      y = "Probability of Amplification",
      title = NULL,
      col = "Fraction",
      subtitle = paste0(
        "GLMM (Binomial, link = logit; Buffer TL only):\n",
        "gel_score ~ [amplicon] + fraction"
      )
    )+ 
    geom_vline(
      xintercept = limit_of_detection,
      lty = 2,
      col = "red"
    ) +
    theme_bw() +
    scale_fill_manual(
      values =
        c("dodgerblue", "darkorange")
    ) +
    scale_color_manual(
      values =
        c("dodgerblue", "darkorange")
    )
)


(df_model_coef <-
    plot_model(
      supplementary_model,
      show.intercept = F,
      show.values = T,
      transform = NULL) %>%
    
    # Modify plot
    pluck ("data") %>% 
    mutate (
      term = 
        str_remove(
          term, 
          "fraction"
        ) 
    )%>% 
    
    mutate (
      term = 
        case_when (
          grepl ("pcr_concentration", term) ~ "[Amplicon]",
          T ~ term
        )
    )
)

(plot_model_coef_supplementary <- 
    ggplot (
      data = df_model_coef,
      aes (
        x = estimate,
        y = reorder(term, estimate),
        xmin = conf.low,
        xmax = conf.high,
        label = p.label
      )
    ) + 
    geom_errorbarh(
      width = 0.1,
      col = "black"
    ) + 
    geom_point(
      size = 3,
      col = "black"
    ) +
    geom_text(
      vjust = -1,
      size = 2.5
      
    ) + 
    geom_vline(
      xintercept = 0,
      lty = 2,
      color = "gray50"
    ) +
    labs (
      x = "Model Coefficient",
      y = "Model Terms",
      subtitle = 
        "Reference Category: Filter" 
    ) + 
    guides (
      fill = "none",
      shape = "none"
    )
)


(plot_log_model_comp_supplementary <- 
    ggarrange(
      plot_log_model_supplementary +
        theme (legend.position = "top",
               plot.subtitle = element_text(size = 8)),
      plot_model_coef_supplementary +
        theme (plot.subtitle = element_text(size = 8)),
      labels = "AUTO",
      widths = c(1, 0.50)
    )
)

#### PREPARE PCR QUANT DATA ####
df_pcr_quant <- 
  df_pcr_score %>% 
  # Remove Positive Control because I do not have a DNA concentration for this
  filter (
    filter_id != "C.tinkeri"
  ) %>% 
  droplevels()


# Check group sizes
table(df_pcr_quant$treatment,
      df_pcr_quant$pcr_quant_score)
dim(df_pcr_quant)


####----ASSESS DATA DISTRIBUTION----####
# Identify the variable to assess
y_var <- 
  df_pcr_quant %>% 
  pull(pcr_concentration) 

hist(y_var)

# Describe the distribution using the Cullen and Frey Graph
descdist(
  data = y_var, 
  discrete = F,
  boot = 1000
)

# Fit the data
distributions <- c(
  "norm",
  "lnorm", # Turn this off when there are 0s in the model
  #"beta", # Turn this off when range of values go beyond [0,1]
  "gamma",
  "weibull",
  "exp",
  "unif",
  "logis"
)

fit <- lapply (
  distributions,
  function (d) fitdist(
    y_var, 
    distr = d, 
    discrete = F,
    method = "mle"
  )
); names(fit) <- distributions

# Compare the fits
lapply(fit, plot)
gofstat(fit)
denscomp(fit)


####----FIT GLMM ON PCR CONCENTRATION----####
# Check group sizes
table(df_pcr_quant$treatment)
dim(df_pcr_quant)

model_1 <- 
  glmmTMB(
    pcr_concentration ~ dna_concentration * treatment,
    data = df_pcr_quant,
    family = lognormal()
  )

model_2 <- 
  glmmTMB(
    pcr_concentration ~ dna_concentration + treatment + (1|site_id/filter_id),
    data = df_pcr_quant,
    family = lognormal()
  )

model_3 <- 
  glmmTMB(
    pcr_concentration ~ dna_concentration * treatment + (1|site_id/filter_id),
    data = df_pcr_quant,
    family = lognormal()
  )


model_list <- 
  list (
    model_1 = model_1,
    model_2 = model_2,
    model_3 = model_3
  )


(model_summary <- 
    lapply (
      model_list,
      summary
    )
)

# Test models
test_performance (model_list)

(model_selection <- 
    compare_performance(
      model_list
    )
)

# Check models 
plot(model_selection)

####----CHECK MODEL ASSUMPTIONS----####

# Assign best model
best_model <- model_2
summary(best_model)

# Get p-values for best model
null_model <- 
  glmmTMB(
    pcr_concentration ~ 1 + (1|site_id/filter_id),
    data = df_pcr_quant,
    family =
  )

(anova_out <- anova(null_model, best_model))

# Check model diagnostics
check_model(best_model)



# Check model residuals - normality, overdispersion, and outlier
model_residuals <- simulateResiduals(best_model)
testResiduals(model_residuals)

# Check model residuals - homoscedasticity
plot(model_residuals)

# Check for multicollinearity
check_collinearity(best_model)  # VIF should be < 5 for all predictors
check_autocorrelation(best_model)

# Check random effect structure
check_singularity(best_model)  # Should return F; if T, then the random effects structure is oversimplified for your data (var = 0)

# Check overdispersion
check_overdispersion(best_model)


####----PLOT MODEL COEFFICIENTS----####
reference_category <- 
  levels(df_stats$treatment)[1] %>% 
  str_replace("_", " x ")

factor_order <- 
  best_model %>% 
  summary() %>% 
  pluck ("coefficients") %>% 
  pluck ("cond") %>% 
  as.data.frame() %>% 
  clean_names() %>% 
  rownames_to_column ("term") %>% 
  mutate (
    term = 
      str_remove(
        term, 
        "treatment"
      ) 
  )%>% 
  filter (!term %in% c("(Intercept)", "dna_concentration")) %>% 
  separate (
    term, 
    into = c("fraction", "preservative"),
    remove =  F,
    sep = "_"
  ) %>% 
  mutate (
    term = 
      case_when (
        grepl ("dna_concentration", term) ~ "Template Concentration",
        T ~ str_replace(term, "_", " x ")
      )
  ) %>% 
  arrange (estimate) %>% 
  pull (term)


(df_model_coef <-
    plot_model(
      best_model,
      show.intercept = F,
      type = "est") %>%
    
    
    # Modify plot
    pluck ("data") %>% 
    mutate (
      term = 
        str_remove(
          term, 
          "treatment"
        ) 
    )%>% 
    
    separate (
      term, 
      into = c("fraction", "preservative"),
      remove =  F,
      sep = "_"
    ) %>% 
    
    mutate (
      term = 
        case_when (
          grepl ("dna_concentration", term) ~ "[Template DNA]",
          T ~ str_replace(term, "_", " x ")
        )
    ) %>% 
    
    mutate (
      fraction = 
        case_when (
          term == "[Template DNA]" ~ "X", # "X" is just a placeholder for plotting, does not mean anything
          T ~ fraction
        )
    ) %>% 
    
    mutate (
      treatment =
        case_when (
          term == "[Template DNA]" ~ NA, 
          T ~ preservative
        )
    ) %>%         
    
    mutate (
      term = factor (
        term, 
        levels = c(
          factor_order,
          "[Template DNA]"
        )
      )
    ) 
)

(plot_model_coef <- 
    ggplot (
      data = df_model_coef,
      aes (
        x = estimate,
        y = term,
        xmin = conf.low,
        xmax = conf.high,
        label = p.label,
        fill = preservative,
        pch = fraction
      )
    ) + 
    geom_errorbarh(
      width = 0.25,
      col = "black"
    ) + 
    geom_point(
      size = 3,
      col = "black",
      alpha = 0.5
    ) +
    geom_text(
      vjust = -1,
      size = 2.5
      
    ) + 
    geom_vline(
      xintercept = 0,
      lty = 2,
      color = "gray50"
    ) +
    labs (
      x = expression ("Model Coefficent ("*beta*")"),
      y = "Model Terms",
      subtitle =paste0("Reference Category:\n", reference_category)
    ) + 
    theme (
      plot.subtitle = element_text (size = 8)
    ) +
    scale_fill_manual(
      values = c(
        "#E7298A",
        "gray50",
        "#7570B3",
        "#1B9E77",
        "#D95F02"
      )
    ) +
    scale_shape_manual(
      values =
        c(21, 22, 4)
    ) +
    guides (
      fill = "none",
      shape = "none"
    )
)



####----PLOT USING GGEFFECTS----####
model <- best_model

# Quick reference plot
model_prediction <- 
  predict_response(
    model,
    terms = c(
      "dna_concentration [0.01:0.25, by=0.01]",
      "treatment"
    )
  )

# Quick reference plot
plot (
  model_prediction,
  show_data = T
)


# Fully modified plot
(df_plot <- 
    model_prediction %>% 
    as.data.frame()
)


(plot_model <- 
    df_plot %>% 
    separate (
      group,
      into = c(
        "fraction",
        "preservative"
      ),
      sep = "_",
      remove = F
    ) %>% 
    rename (
      treatment = group
    ) %>% 
    filter (
      fraction != "PCR"
    ) %>% 
    mutate (
      source =
        factor (
          fraction,
          levels = 
            c(
              "Filter",
              "Filtrate"
            )
        )
    ) %>% 
    ggplot (
      aes (
        x = x,
        y = predicted,
        group = treatment,
        lty = fraction,
        col = preservative,
        fill = preservative
      )
    ) + 
    geom_line() +
    geom_ribbon (
      aes (
        ymin = conf.low,
        ymax = conf.high
      ),
      alpha = 0.05,
      color = NA
    )  +
    # Add raw data
    geom_point(
      data = df_pcr_quant %>%
        select (-preservative) %>%
        filter (fraction != "PCR") %>%
        mutate (preservative = treatment_modified_controls),
      aes (
        x = dna_concentration,
        y = pcr_concentration,
        fill = preservative,
        shape = fraction
      ),
      size = 2.5,
      alpha = 0.50
    ) +
    
    # # Double check labels
    # geom_text(
    #   data =
    #     df_pcr_quant %>%
    #     select (-preservative) %>%
    #     filter (fraction != "PCR") %>%
    #     mutate (preservative = treatment_modified_controls,
    #             filter_id =  str_remove(filter_id, "-\\d{4}")),
    #   aes (
    #     label = filter_id,
    #     x = dna_concentration,
    #     y = pcr_concentration
    #   ),
    #   size = 3,
    #   hjust = 1,
    #   vjust =1
    # )+
    
    labs (
      x = "[Template] (ng/uL)",
      y = "[Amplicon] (ng/uL)",
      shape = "Fraction",
      lty = "Fraction",
      fill = "Preservation\nTreatment",
      color = "Preservation\nTreatment",
      subtitle =  paste(
        "GLMM (family = lognormal):\n",
        "[amplicon] ~ [template] + treatment + (1|site_id/filter_id)\n"
      )
    ) + 
    theme (
      legend.position = "top",
      legend.box = "vertical",
      legend.spacing.y = unit(0.0001, "cm"),
      plot.subtitle = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    ) +
    scale_fill_brewer(
      palette = "Dark2"
    ) + 
    scale_color_brewer(
      palette = "Dark2"
    ) +
    scale_linetype_manual(
      values = c(1, 2)
    ) + 
    scale_shape_manual(
      values = c(21, 22)
    ) +
    scale_y_log10()
)


####----COMBINE PLOTS----####
# Create a multifacet plot 
(plot_model_and_coef <- 
   ggarrange(
     plot_model,
     plot_model_coef + 
       theme (
         axis.text.y = element_text (size = 8)
       ),
     ncol = 2,
     labels = c("A", "B"),
     widths = c(1, 0.85)
   )
)


(plot_pcr_scores <- 
    ggarrange(
      plot_pcr_gel_score + 
        theme(
          plot.subtitle = element_text (size = 8)
        ),
      plot_pcr_quant_score +
        theme (axis.title.y = element_blank(),
               plot.subtitle = element_text (size = 8)
        ),
      common.legend = T,
      legend = "right",
      widths = c(1, 1),
      labels = c("C", "D")
    )
)


(plot_models <- 
    ggarrange(
      plot_model_and_coef,
      plot_pcr_scores,
      plot_log_model + theme(plot.subtitle = element_text (size = 8)) + theme_bw(),
      nrow = 3,
      labels = c(NA, NA, "E"),
      heights = c(1,0.75,0.75)
    )
)

# (plot_amp_frequency <- 
#     ggarrange(
#       plot_quant_amp_frequncy,
#       plot_gel_amp_frequncy,
#       ncol = 2,
#       common.legend = F,
#       legend = "top",
#       labels = c("C", "D")
#     )
# )


####----SAVE FILES----####
ggsave(
  plot_models,
  filename = path_plot_quant_model,
  width = 8,
  height = 10, 
  units = "in",
  dpi = 600)

ggsave(
  plot_log_model,
  filename = path_plot_log_model,
  width = 7,
  height = 5, 
  units = "in",
  dpi = 330)



ggsave(
  plot_log_model_comp_supplementary,
  filename = path_plot_log_model_supplementary,
  width = 7,
  height = 4.5, 
  units = "in",
  dpi = 330)


write.csv (
  df_model_coef,
  file = path_model_coef_amplicon_concentration
)

