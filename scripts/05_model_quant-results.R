#### model_quant-results ####
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
  
  ## Used in earlier versions of the script
  #rstatix,
  #lme4,
  #lmerTest,
  #FactoMineR,
  #factoextra,
  #PerformanceAnalytics,
  #lmtest,
  #fitur,
  #smplot2,
  #ggeffects,
  #broom.mixed,
  #bbmle,
  #multcompView,
  
  ## Used in final version of the script
  janitor,
  ggpubr,
  emmeans,
  multcomp,
  emmeans,
  DHARMa,
  ggbeeswarm,
  fitdistrplus,
  glmmTMB,
  sjPlot,
  performance,
  tidyverse
)

#### USER DEFINED VARIABLES ####

# Assign input file paths
path_df_for_stats <- 
  "../data/df_for_stats.rds"

# Assign output file paths
path_compare_control <- 
  "../results/plot_compare_control.png"

path_model_coef <- 
  "../results/plot_model_coef.pdf"

path_model_emmeans <- 
  "../results/model_emmeans.csv"

path_plot_contrasts <- 
  "../results/plot_model-contrasts.png"


####----LOAD FILES----####
df_stats <- 
  path_df_for_stats %>% 
  read_rds() 


# Set theme for ggplots
theme_set(theme_bw())

####----ASSESS DATA DISTRIBUTION----####
# Identify the variable to assess
y_var <- 
  df_stats %>% 
  select (dna_concentration_ng_per_ul) %>% 
  na.omit() %>% 
  pull()

# Describe the distribution using the Cullen and Frey Graph
descdist(
  data = y_var, 
  discrete = F,
  boot = 1000
)

# Fit the data
distributions <- c(
  "norm",
  "lnorm",
  "beta", # Turn this off when range of values go beyond [0,1]
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

####----TEST DIFFERENCES BETWEEN CONTROLS----####

# Prepare a model that compares the field and extraction controls
model_control_comparison <- 
  df_stats %>% 
  filter (
    sample_type %in% c(
      "Extraction Control",
      "Field Control"
    )
  ) %>% 
  lm (formula = log10_dna_concentration ~ sample_type * fraction)

# Check model summary
summary(model_control_comparison)

(plot_model_control_coef <- 
    plot_model(
      show.intercept = T,
      colors = "bw",
      model_control_comparison,
      show.values = T,
      show.p = T,
      vline.color="gray50",
    ) + 
    theme_bw()
)

# Perform ANOVA
(anova_out <- anova(model_control_comparison))

# Check model diagnostics
check_model(model_control_comparison)

# Check model residuals - normality, overdispersion, and outlier
model_residuals <- simulateResiduals(model_control_comparison)
testResiduals(model_control_comparison)
shapiro.test(model_residuals$fittedResiduals)

# Check model residuals - homoscedasticity
plot(model_residuals)

# Check for multicollinearity
check_collinearity(model_control_comparison)  # VIF should be < 5 for all predictors
check_autocorrelation(model_control_comparison)

# Check random effect structure
check_singularity(model_control_comparison)  # Should return F; if T, then the random effects structure is oversimplified for your data (var = 0)


####----EMMEANS CALCULATION FOR CONTROLS----####
model <- model_control_comparison

# Compute estimated marginal means for the treatment*fraction interaction
(emmeans_model <-
    emmeans(model,
            ~  sample_type * fraction,
            alpha = 0.05,
            type = "response"
    )
)

# Check contrast
contrast(
  regrid(emmeans_model), # emmeans back transformed to the original units of response var
  method = 'pairwise', 
  simple = 'each', 
  combine = F, 
  adjust = "mvt")

# group treatment combos
(groupings_model <-
    multcomp::cld(emmeans_model, 
                  alpha = 0.05,
                  Letters = letters,
                  type="response",
                  adjust = "mvt") %>%
    as.data.frame %>%
    mutate(group = str_remove_all(
      .group,
      " "),
      group = str_replace_all(
        group, 
        "(?<=.)(?=.)", 
        "\n"
      )
    )
)

####----VISUALIZE EMMEANS----####
(plot_compare_controls <- 
   
   # Plot EMMEANS
   groupings_model %>% 
   ggplot(aes(x=sample_type,
              y=10^emmean,
              fill = fraction)) +
   geom_col(
     color = "black",
     position = position_dodge(width=0.9),
     alpha = 0.50
   ) +
   geom_point(
     size = 3.5,
     position = position_dodge(width = 0.9),
     pch = 21,
     col = "black"
   ) + 
   geom_errorbar(  
     aes(ymin = 10^lower.CL,
         ymax = 10^upper.CL),
     width = 0.2,
     color = "black",
     position = position_dodge(width=0.9)
   ) +
   
   # PLOT COMP LETTERS
   geom_text(
     aes(label=group),
     position = position_dodge(width=0.9),
     color = "black",
     vjust = -1,
     hjust = -0.5,
     size = 3.5
   ) +  
   
   
   # PLOT RAW DATA
   geom_beeswarm (
     data = df_stats %>% 
       filter (
         sample_type %in% c(
           "Extraction Control",
           "Field Control"
         )
       ),
     aes (x = sample_type,
          y = dna_concentration_ng_per_ul,
          fill = fraction ),
     pch = 21,
     color = "black",
     alpha = 0.50,
     dodge.width = 0.9
     
   ) +
   
   theme_bw() +
   scale_fill_manual(
     values = c(
       "dodgerblue",
       "darkorange"
     )
   ) + 
   labs (
     x = "Sample Type",
     y = "[DNA] (ng/uL)",
     fill = "Fraction",
     subtitle = expression (
       "ANOVA:"~log[10]*"[DNA] ~ sample type * fraction")
   )+
   theme(legend.position = "right")  +
   theme(axis.title.x = element_blank(),
         plot.subtitle = element_text(size = 8)) 
)

####----MODIFY DATA FRAME FOR STATS----####
df_stats_modified <- 
  df_stats %>% 
  
  # Remove PCR controls
  filter (fraction != "PCR") %>% 
  droplevels() %>% 
  
  # Create a modified treatment column that merges field and extraction controls into a single variable
  mutate (
    treatment_merged_controls = 
      case_when (
        grepl ("Ctrl", preservative) ~ "Neg Ctrl",
        T ~ preservative
      ) %>% 
      factor (levels = c(
        "Neg Ctrl",
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
        treatment_merged_controls,
        sep = "_"
      ) %>%
      factor (levels = c(
        "Filter_Neg Ctrl",
        "Filter_No Pres",
        "Filter_DESS",
        "Filter_Buffer TL",
        "Filtrate_Neg Ctrl",
        "Filtrate_DESS",
        "Filtrate_Buffer TL"
      ))
  ) %>%
  
  # Scale continuous explanatory variables
  mutate (
    volume_filtered_ml_scaled = as.numeric(scale(volume_filtered_ml)),
    duration_on_ice_hours_scaled = as.numeric(scale(duration_on_ice_hours))
  )

# Check group sizes
table(df_stats_modified$treatment)

####----FIT GLMM----####

# Ask Jason: How do we report if a glmmTMB is significant?

# Model 1: assumes homoscedasticity
model_1 <- 
  glmmTMB(
    formula = dna_concentration_ng_per_ul ~ treatment + (1|site_id/filter_id),
    data = df_stats_modified,
    family = lognormal()
  )

# Model 2: assumes heteroscedasticity
model_2 <- 
  glmmTMB(
    formula = dna_concentration_ng_per_ul ~ treatment + (1|site_id/filter_id),
    data = df_stats_modified,
    family = lognormal(),
    dispformula = ~ treatment
  )


model_list <- 
  list (
    model_1 = model_1,
    model_2 = model_2
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
best_model<- 
  model_selection %>%
  as.data.frame() %>% 
  arrange (AIC) %>% 
  slice_head(n = 1) %>% 
  pull(Name) %>% 
  get()

summary(best_model)

# Compare best_model with the null model
null_model <- 
  glmmTMB (
    dna_concentration_ng_per_ul ~ 1 + (1|site_id/filter_id),
    data = df_stats_modified,
    family = lognormal())

anova(null_model, best_model)

# Check model diagnostics
check_model(best_model)

# Check model residuals - normality, overdispersion, and outlier
model_residuals <- simulateResiduals(best_model)
testResiduals(model_residuals)
shapiro.test(model_residuals$fittedResiduals)

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
  best_model$frame %>% 
  (function (x) levels (x$treatment)[1]) %>% 
  str_replace("_", " x ")



(model_coef <-
    plot_model(
      best_model,
      show.intercept = F,
      type = "est") %>%
    
    
    # Modify plot
    pluck ("data") %>% 
    filter (wrap.facet == "conditional") %>% 
    mutate (
      term = str_remove(term, "treatment"),
      term = str_replace(term, "_", " x "), 
      term = factor (
        term, 
        levels = c(
          "Filtrate x Neg Ctrl",
          "Filtrate x DESS",
          "Filtrate x Buffer TL",
          "Filter x Buffer TL",
          "Filter x DESS",
          "Filter x No Pres"
        )
      )
    ) %>% 
    separate_wider_delim (
      cols = term, 
      names = c("fraction", "preservative"),
      delim = " x ",
      cols_remove =  F
    )
)

write.csv (
  model_coef,
  file = "../results/model-coef_dna-yield.csv"
)
    
    
(plot_model_coef <- 
    model_coef %>% 
    ggplot (
      aes (
        x = estimate,
        y = term,
        xmin = conf.low,
        xmax = conf.high,
        label = p.label,
        fill = fraction
      )
    ) + 
    geom_errorbarh(
      width = 0.25,
      col = "black"
    ) + 
    geom_point(
      size = 2,
      pch = 21,
      col = "black"
    ) +
    geom_text(
      vjust = -0.75,
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
        "dodgerblue",
        "darkorange"
      )
    ) +
    guides (
      fill = "none"
    )
)


####----EMMEANS CALCULATION - CATEGORICAL PREDICTORS ONLY----####
model <- best_model

# Compute estimated marginal means for the treatment*fraction interaction
(emmeans_model <-
    emmeans(model,
            ~  treatment,
            alpha = 0.05
    ) %>% 
    regrid (type = "response")
)

####----ASSIGN GROUPING USING CLD----####
# group treatment combos
(groupings_model <-
   multcomp::cld(
     emmeans_model, 
     alpha = 0.05,
     Letters = letters,
     type="response",
     adjust = "mvt"
   ) %>%
   as.data.frame %>%
   mutate(
     group = str_remove_all(
       .group,
       " "),
     group = str_replace_all(
       group, 
       "(?<=.)(?=.)", 
       "\n"
     )
   ) %>% 
   
   # separate "treatment" into "fraction" and "treatment" columns
   separate(
     treatment, 
     into = c(
       "fraction", 
       "preservative"), 
     sep = "_") %>% 
   
   # Turn implicit missing values to explicit missing values
   complete(preservative, fraction) %>% 
   mutate (
     preservative = 
       factor (
         preservative,
         levels = c(
           "Neg Ctrl",
           "No Pres",
           "DESS",
           "Buffer TL"
         ))
   )
)

####----VISUALIZE EMMEANS----####
(plot_model <- 
   
   # Plot EMMEANS
   groupings_model %>% 
   ggplot(aes(x=preservative,
              y=response,
              fill = fraction)) +
   geom_col(
     color = "black",
     position = position_dodge(width=0.9),
     alpha = 0.50
   ) +
   geom_point(
     size = 3.5,
     position = position_dodge(width = 0.9),
     pch = 21,
     col = "black"
   ) + 
   geom_errorbar(  
     aes(ymin = asymp.LCL,
         ymax = asymp.UCL),
     width = 0.2,
     color = "black",
     position = position_dodge(width=0.9)
   ) +
   
   # PLOT COMP LETTERS
   geom_text(
     aes(label=group),
     position = position_dodge(width=0.9),
     color = "black",
     vjust = -0.25,
     hjust = -0.75,
     size = 3
   ) +  
   
   
   # PLOT RAW DATA
   geom_beeswarm (
     data = df_stats_modified %>% 
       # Turn implicit missing values to explicit missing values
       complete(treatment_merged_controls, fraction),
     aes (x = treatment_merged_controls,
          y = dna_concentration_ng_per_ul,
          fill = fraction ),
     pch = 21,
     color = "black",
     alpha = 0.50,
     dodge.width = 0.9
     
   ) +
   
   theme_bw() +
   scale_fill_manual(
     values = c(
       "dodgerblue",
       "darkorange"
     )
   ) + 
   labs (
     x = NULL,
     y = "[DNA] (ng/uL)",
     fill = "Fraction",
     subtitle = paste("GLMM (family = lognormal):\n",
     "[DNA] ~ treatment + (1|site_id/filter_id)\ndispersal = ~treatment")
   )+
   theme(legend.position = "top",
         plot.subtitle = element_text(size = 8,
                                      margin = margin(t = 5))) +
   
   
   # Annotate (no data, "ND") for no preservative x preservative combination
   annotate(
     geom = "text",
     x = 2.25,
     y = 0.025,
     label = "ND"
   )
) 

####----VISUALIZE RELATIVE DNA CONCENTRATION----####
adj_emmeans <-
  
  # Perform contrast on the emmeans model by subtracting the Neg Ctrl for each fraction
  contrast(
    emmeans_model,
    list(
      "Filter_No Pres" =
        c(-1, 1, 0, 0, 0, 0, 0),
      
      "Filter_DESS" =
        c(-1, 0, 1, 0, 0, 0, 0),
      
      "Filter_Buffer TL" =
        c(-1, 0, 0, 1, 0, 0, 0),
      
      "Filtrate_DESS" =
        c( 0, 0, 0, 0, -1, 1, 0),
      
      "Filtrate_Buffer TL" =
        c( 0, 0, 0, 0, -1, 0, 1)
    )
  ) %>%   
  as.data.frame() %>% 
  separate_wider_delim (
    cols = contrast,
    names = c("fraction", "preservative"),
    delim = "_",
    cols_remove = F
  ) %>% 
  mutate (estimate = round(estimate, 2))
  

## 3. Sum DNA yield per preservative
df_sum <- 
  adj_emmeans %>% 
  summarize (
    sum_emmean = sum(estimate),
    .by = preservative
  ) 

# 4. Get the max value among preservative: max(no preservative, dess, buffer tl) and set as denominator for model
denominator <- 
  df_sum %>% 
  pull (sum_emmean) %>% 
  max() 

# 5. Calculate % recovery
df_pct_recovery <- 
  adj_emmeans %>%
  
  # Calculate pct_recovery based on max concentration
  mutate(
    pct_recovered = (estimate / denominator) * 100
  ) %>%
  select(
    preservative, 
    fraction, 
    pct_recovered
  ) 

# Create unrecovered rows per treatment
df_pct_unrecovered <- 
  df_pct_recovery %>%
  summarize (
    sum_pct_recovered = sum(pct_recovered),
    pct_unrecovered = 100 - sum_pct_recovered,
    .by = preservative
  ) %>% 
  mutate (fraction = "Unrecovered DNA") %>% 
  rename (pct_recovered = pct_unrecovered) %>% 
  select (-sum_pct_recovered)


# Bind the two together
(plot_relative_concentration <- 
    bind_rows(
      df_pct_recovery, 
      df_pct_unrecovered
    ) %>% 
    filter (preservative != "Neg Ctrl") %>% 
    mutate (
      preservative = 
        factor(preservative,
               levels = c(
                 #"Neg Ctrl",
                 "No Pres",
                 "DESS",
                 "Buffer TL"
               )),
      fraction = 
        factor (fraction,
                levels = c(
                  "Unrecovered DNA",
                  "Filter",
                  "Filtrate"
                ))
    ) %>% 
    ggplot (
      aes (x = preservative,
           y = pct_recovered,
           fill = fraction
      )
    ) + 
    geom_col(
      position = "fill",
      col = "black",
      alpha = 0.50,
      width = 0.75
    ) +
    scale_fill_manual(
      values = c(
        "white",
        "dodgerblue",
        "darkorange"
      )
    ) + 
    labs (
      x = NULL,
      y = "Relative Mean\nDNA Recovery (%)",
      fill = "Fraction"
    ) +
    scale_y_continuous(
      labels = scales::percent
    ) +
    coord_cartesian(
      expand = F
    ) +
    theme_classic()
)



####----CONTRASTS FOR FILTER + FILTRATE EMM----####
pooled_contrasts <-
  contrast(
    emmeans_model,
    list(
      "No Pres" =
      
        # Filter: NoPres - Neg Ctrl
        c(-1, 1, 0, 0, 0, 0, 0),
      
      "DESS" =
        # Filter: DESS - Neg Ctrl
        c(-1, 0, 1, 0, 0, 0, 0) +
        # Filtrate: DESS - Neg Ctrl
        c(0, 0, 0, 0, -1, 1, 0),
      
      "Buffer TL" =
        # Filter: Buffer TL - Neg Ctrl
        c(-1, 0, 0, 1, 0, 0, 0) +
        # Filtrate: Buffer TL - Neg Ctrl
        c( 0, 0, 0, 0, -1, 0, 1)
    )
  )

summary(pooled_contrasts)

contrast_df_preservative <-   
    pairs(
      pooled_contrasts,
      adjust = "mvt"
    ) %>% 
    as.data.frame() %>% 
    arrange (estimate) 


factor_level <- contrast_df_preservative$contrast

(plot_contrast <- 
    contrast_df_preservative %>% 
    mutate (
      sig_symbol =
        case_when (
          p.value < 0.001 ~ "***",
          p.value < 0.01 ~ "**",
          p.value < 0.05 ~ "*",
          T ~ "ns"
        )
    ) %>% 
    mutate (
      label = paste(
        round (estimate, 3),
        sig_symbol,
        sep = " "
      )
    ) %>% 
    mutate (
      contrast = 
        fct_relevel(
          contrast,
          factor_level
        )
    )%>% 
    clean_names() %>%  
    ggplot (
      aes (
        x = estimate,
        y = contrast,
        xmin = estimate - se,
        xmax = estimate + se,
        label = label      
      )
    ) + 
    geom_errorbarh(
      width = 0.25,
      col = "black"
    ) + 
    geom_point(
      size = 3,
      pch = 19,
      col = "black"
    ) +
    geom_text(
      vjust = -0.75,
      size = 2.5
      
    ) + 
    labs (
      x = expression ("Difference ("*Delta*") in EMM [DNA] (ng/uL)"),
      y = "Contrasts",
    ) 
)



####----COMBINE PLOTS----####
# Create a multifacet plot 
(plot_model_and_coef <- 
   ggarrange (
     ggarrange(
       plot_model,
       plot_model_coef + 
         theme (
           axis.text.y = element_text (size = 8)
         ),
       ncol = 2,
       labels = "AUTO",
       widths = c(1, 0.85)
       #align = "h"
     ),
     plot_relative_concentration,
     nrow = 2,
     heights = c(1, 0.65),
     labels = c(NA, "C")
   )
)



####----SAVE FILES----####
ggsave(
  plot_model_and_coef,
  filename = path_model_coef,
  width = 7,
  height = 7, 
  units = "in",
  dpi = 600)

ggsave(
  plot_compare_controls,
  filename = path_compare_control,
  width = 5,
  height = 5, 
  units = "in",
  dpi = 330)


ggsave(
  plot_contrast,
  filename = path_plot_contrasts,
  width = 5,
  height = 3.5, 
  units = "in",
  dpi = 330)



write.csv(
  groupings_model %>% 
    mutate (dna_concentration = response),
  path_model_emmeans 
)

# Project table
tab_model(best_model)



