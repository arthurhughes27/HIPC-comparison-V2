suppressPackageStartupMessages({
  library(shiny)
  library(dplyr)
  library(tidyr)
  library(GSA)
  library(stringr)
  library(reshape2)
  library(foreach)
  library(doParallel)
  library(purrr)
  library(circlize)
  library(tidyverse)
  library(ggh4x)
  library(ggradar)
  library(ggnewscale)
  library(gridExtra)
  library(ggpubr)
  library(grid)
  library(data.table)
  library(ggtext)
  library(shinyBS)
  library(ggforce)
  library(scales)
  library(glue)
})

p_dearseq_dgsa_results_processed = fs::path("output", "results", "dearseq", "dearseq_dgsa_results_processed.rds")
p_qusage_dgsa_results_processed = fs::path("output", "results", "qusage", "qusage_dgsa_results_processed.rds")

# I am including results from dearseq and a different DGSA method (QuSAGE), so I bind these dataframes together.
# You can just load one, but call it must be called "results_df".

results_dearseq = readRDS(p_dearseq_dgsa_results_processed)
results_qusage = readRDS(p_qusage_dgsa_results_processed)

results_df = bind_rows(results_dearseq, results_qusage)

rm(results_dearseq, results_qusage)

# Define some meta options for colours

# 1. define ordering & colours
condition_order <- c(
  "Tuberculosis (RVV)",
  "Varicella Zoster (LV)",
  "Yellow Fever (LV)",
  "Ebola (RVV)",
  "Hepatitis A/B (IN/RP)",
  "HIV (RVV)",
  "Influenza (IN)",
  "Influenza (LV)",
  "Malaria (RP)",
  "Meningococcus (CJ)",
  "Meningococcus (PS)",
  "Pneumococcus (PS)",
  "Smallpox (LV)"
)
condition_colors <- c(
  "#b94a73",
  "#c6aa3c",
  "#6f71d9",
  "#64c46a",
  "#be62c2",
  "#7d973c",
  "#563382",
  "#4ea76e",
  "#bc69b0",
  "#33d4d1",
  "#bb4c41",
  "#6a87d3",
  "#b57736"
)
names(condition_colors) <- condition_order

aggregate_order <- c(
  "Antigen Presentation",
  "Inflammatory/TLR/Chemokines",
  "Interferon/Antiviral Sensing",
  "Monocytes",
  "DC Activation",
  "Neutrophils",
  "NK Cells",
  "Signal Transduction",
  "ECM And Migration",
  "Energy Metabolism",
  "Cell Cycle",
  "Platelets",
  "T Cells",
  "B Cells",
  "Plasma Cells",
  "NA"
)
aggregate_colors <- c(
  "#7c5fcd",
  "#57c39d",
  "#c1121f",
  "#55c463",
  "#7082ca",
  "#64a332",
  "#45aecf",
  "#df9545",
  "#b7b238",
  "#a6b36c",
  "#667328",
  "#662d2e",
  "#ff8fa3",
  "#c05299",
  "#8f2d56",
  "#adb5bd"
)
names(aggregate_colors) <- aggregate_order

### DEFINE USER PARAMETERS ###

### COMMON PARAMETERS ###
# Timepoints
timepoints_set = levels(results_df$time)

# Conditions
conditions_set = levels(results_df$condition)

# Aggregates
aggregates_set = levels(results_df$gs.aggregate)

# Score options
scores_set = c("Activation scores" = "activation.score",
               "Mean of fold changes" = "fc.score")

# Correlation options
correlation_set = c("Mean of gene correlations" = "mean.corr",
                    "Correlation of mean gene expression" = "corr.mean")

# Effect size filtration
filtration_set = c("No Filtration on effect size" = "none",
                   "Filter by absolute activation score" = "activation.score",
                   "Filter by absolute mean fold change" = "fc.score")

# How to filter by effect size
filtration_mode_set = c("Data-driven" = "data",
                        "User threshold" = "user")

# P-value multiplicity correction method
p_correction_set = c("No correction" = "none",
                     "Bonferroni" = "bonferroni",
                     "Benjamini-Hochberg/FDR" = "BH",
                     "Benjamini-Yekutieli" = "BY",
                     "Hochberg" = "hochberg",
                     "Holm" ="holm",
                     "Hommel" = "hommel")

# P-value multiplicity correction
sensitivity_heatmap_p_correction_set = c(
  "Bonferroni" = "bonferroni",
  "Benjamini-Hochberg/FDR" = "BH",
  "Benjamini-Yekutieli" = "BY",
  "Holm" ="holm",
  "Hommel" = "hommel",
  "Hochberg" = "hochberg")

# P-value correction/effect size filtration scale (global, within-time, within-comparison)
p_approach_set = c("Global" = "global",
                   "Within-time" = "withinTime",
                   "Within-comparison" = "withinComparison")


### CIRCOS PLOT OPTIONS ###
# What is the rule to plot a link?
arc_set = c("Common differential regulation" = "any",
            "Common differential regulation in same direction" = "direction",
            "Positive only" = "positive",
            "Negative only" = "negative")

# Which conditions should be displayed and how?
circos_order_set = c("Default fixed order and display all" = "set_all",
                     "Default fixed order and only display available" = "set_available",
                     "Hierarchical Clustering" = "cluster")

circos_ring_set = c("All inner rings" = "all",
                    "Expression profiles only" = "expression",
                    "No inner rings" = "none")

### HEATMAP OPTIONS ###
# How should the y-axis labels (gene sets) be ordered?
y_order_set = c("Hierarchical clustering" = "cluster",
                "By aggregate" = "aggregate")

# How should the x-axis labels (conditions) be ordered?
x_order_set = c("Default fixed order" = "set",
                "Cluster within-time" = "cluster-time")

# Filter genesets by common differential expression?
filter_commonDE_set = c("No filtration" = "none",
                        "Filter on across-time common DE" = "global",
                        "Filter on within-time common DE" = "withinTime",
                        "Filter by sharing score" = "score")

### LONGITUDINAL EXPRESSION PROFILES OPTIONS ###
# How should the expression profile bars by coloured?
exprprofiles_barcolours_set = c("Direction of regulation" = "direction",
                                "Aggregates" = "aggregate")

### SPIDER PLOTS OPTIONS ###
# Which plots should be shown?
spider_grouping_set = c("Selected conditions within one time" = "withinTime",
                        "One condition across selected times" = "withinCondition")

# Should grid positions be conserved?
spider_grid_set = c("Fix positions on grid" = "set",
                    "Only display available vaccines" = "available")

### SENSITIVITY HEATMAP OPTIONS ###
# Define set of reasonable adjusted p-value thresholds
p_thresholds_set =  c(0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1)

# Define set of reasonable FC thresholds
fc_thresholds_set =  seq(0,2,0.1)

# Methods to use
methods_set = results_df$method %>% unique()

### INDIVIDUAL PARAMETER SENSITIVITY OPTIONS ###
indiv_sensitivity_parameter_set = c("Differential gene-set analysis method" = "method",
                                    "P-value correction method" = "p_correction",
                                    "P-value correction subset" = "p_approach",
                                    "Adjusted p-value threshold" = "p_threshold",
                                    "Absolute log2-fold-change threshold" = "fc_threshold")

indiv_sensitivity_groupby_set = c("By parameter to vary" = "byparameter", 
                                  "By vaccine" = "byvaccine")

method_unique_filtration_set = c("No filtration" = "none",
                                 "Positively regulated only" = "positive",
                                 "Negatively regulated only" = "negative")

### DEFINE PAGE LAYOUT ###
ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
      /* whole page background */
      body {
        background-color: #ffffff;
      }
      /* tab panel content background */
      .tab-content {
        background-color: #ffffff;
        padding: 20px;
        border-radius: 5px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
      }
      /* optional: change the tab nav background too */
      .nav.nav-tabs {
        background-color: #ffffff;
      }
      .nav-tabs > li > a {
        color: #023047;
      }
      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:focus,
      .nav-tabs > li.active > a:hover {
        background-color: #ffffff;
        color: #023047;
      }
    "))
  ),
  
  titlePanel("Vaccine Comparison Application"),
  
  tabsetPanel(
    id = "main_tabs",
    
    ## 0) About tab
    tabPanel(
      title = "About",
      value = "about_tab",
      fluidRow(
        column(
          width = 12,
          h2("Application Background"),
          p("This RShiny app allows you to explore differential gene‐set activation across various conditions and times."),
          p("The 'Circos plots' tab allows for the exploration of commonly differentially expressed gene sets across conditions within a single time."),
          p("The 'Heatmaps' tab generates a heatmap to represent differential gene-set expression profiles of conditions across time."),
          p("The 'Longitudinal expression profiles' tab represents these profiles in a different way."),
          p("The 'Spider plots' tab represents average differential activation of gene-set aggregates across conditions or time."),
          p("The 'Sensitivity heatmap' tab represents geneset-aggregate-level results as a function of user differential expression specifications (e.g. significance level, effect size cutoff) in order to explore the sensitivity of the results."),
          p("The 'Individual parameter sensitivity' tab represents vaccine-level results as a function of one user specification, with the other parameters fixed, to explore the parameters which have the most impact on results."),
        )
      )
    ),
    
    ### TAB 1 : CIRCOS PLOTS ###
    tabPanel(
      title = "Circos plots", value = "circos_tab",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          
          # Update plot button
          actionButton(
            inputId = "circos_Update",
            label   = "Update Circos Plot",
            class   = "btn-primary",
            width   = '100%'
          ),
          
          # Data selection group
          wellPanel(
            style = "background-color: #c6d2ed;",
            h4("Data Selection"),
            
            # Results from which method?
            tags$label(
              "DGSA Method:",
              tags$span(
                icon("info-circle"),
                id    = "circos_method_info",
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            selectInput(
              inputId = "circos_method",
              label   = NULL,
              choices = methods_set,
              selected = unique(results_df$method)[1],
              multiple = FALSE
            ),
            bsTooltip(
              id      = "circos_method_info",
              title   = "Choose the method of differential gene set expression analysis whose results are to be viewed.",
              placement = "right",
              trigger   = "hover"
            ),
            
            # Conditions to visualise
            tags$label(
              "Conditions:",
              tags$span(
                icon("info-circle"),
                id    = "circos_conditions_info",
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            selectInput(
              inputId = "circos_conditions",
              label   = NULL,
              choices = conditions_set,
              selected = conditions_set,
              multiple = TRUE
            ),
            bsTooltip(
              id      = "circos_conditions_info",
              title   = "Choose one or more experimental conditions to include in the circos plot.",
              placement = "right",
              trigger   = "hover"
            ),
            
            # Timepoint to visualise
            tags$label(
              "Timepoint (Days):",
              tags$span(
                icon("info-circle"),
                id    = "circos_time_info",       # target for tooltip
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            selectInput(
              inputId = "circos_time",
              label   = NULL,
              choices = timepoints_set,
              selected = 1
            ),
            bsTooltip(
              id      = "circos_time_info",
              title   = "Choose the timepoint to visualise in the circos plot.",
              placement = "right",
              trigger   = "hover"
            ),
            
            # Aggregates to visualise
            tags$label(
              "Gene-set Aggregates:",
              tags$span(
                icon("info-circle"),
                id    = "circos_aggregates_info",       # target for tooltip
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            selectInput(
              inputId = "circos_aggregates",
              label   = NULL,
              choices = aggregates_set,
              selected = aggregates_set,
              multiple = TRUE
            )
          ),
          bsTooltip(
            id      = "circos_aggregates_info",
            title   = "Choose the gene-set aggregates to visualise in the circos plot.",
            placement = "right",
            trigger   = "hover"
          ),
          
          # Determining significance
          wellPanel(
            style = "background-color: #cce3de;",
            h4("Statistical Significance"),
            
            # p-value correction method
            tags$label(
              "P-value correction method",
              tags$span(
                icon("info-circle"),
                id    = "circos_p_correction_info",       # target for tooltip
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            selectInput(
              inputId = "circos_p_correction",
              label   = NULL,
              choices = p_correction_set,
              selected = "BH"
            ),
            bsTooltip(
              id      = "circos_p_correction_info",
              title   = "Choose the p-value multiplicity correction method.",
              placement = "right",
              trigger   = "hover"
            ),
            
            # Subset of p-values to correct
            tags$label(
              "Correction & Filtering Strategy:",
              tags$span(
                icon("info-circle"),
                id    = "circos_p_approach_info",       # target for tooltip
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            selectInput(
              inputId = "circos_p_approach",
              label   = NULL,
              choices = p_approach_set,
              selected = "global"
            ),
            bsTooltip(
              id      = "circos_p_approach_info",
              title   = "Choose the strategy for correcting p-values and for filtration on effect size. 'Global' refers to correcting all p-values across conditions and timepoints at once. 'Within-timepoint' refers to correcting p-values within each distinct timepoint. 'Within-comparison' refers to correcting p-values within each condition-timepoint combination.",
              placement = "right",
              trigger   = "hover"
            ),
            
            # p-value threshold
            tags$label(
              "P-value Threshold:",
              tags$span(
                icon("info-circle"),
                id    = "circos_p_threshold_info",       # target for tooltip
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            numericInput(
              inputId = "circos_p_threshold",
              label   = NULL,
              value   = 0.05,
              min     = 0,
              max     = 1,
              step    = 0.01
            ),
            bsTooltip(
              id      = "circos_p_threshold_info",
              title   = "Threshold for determining statistical significance based on adjusted p-values.",
              placement = "right",
              trigger   = "hover"
            ),
            
            # Filter on effect size? If so, which effect size?
            tags$label(
              "Effect size to filter on:",
              tags$span(
                icon("info-circle"),
                id    = "circos_filter_variable_info",       # target for tooltip
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            selectInput(
              inputId = "circos_filter_variable",
              label   = NULL,
              choices = filtration_set,
              selected = "none"
            ),
            bsTooltip(
              id      = "circos_filter_variable_info",
              title   = "If the user desires to add a condition on the effect size to define statistical significance, on which variable should this be?",
              placement = "right",
              trigger   = "hover"
            ),
            
            # If filtering on effect size, how? user-threshold or data-driven approach
            conditionalPanel(
              condition = "input.circos_filter_variable != 'none'",
              tags$label(
                "How to filter on effect size?",
                tags$span(
                  icon("info-circle"),
                  id    = "circos_filter_mode_info",       # target for tooltip
                  style = "cursor: help; color: #337ab7; margin-left: 4px;"
                )
              ),
              selectInput(
                inputId = "circos_filter_mode",
                label   = NULL,
                choices = filtration_mode_set,
                selected = "user"
              ),
              bsTooltip(
                id      = "circos_filter_mode_info",
                title   = "How should filtration on effect size be performed? This can either be a threshold determined by the user, or data-driven by filtering out all results whose scores are below a given quantile. If the data-driven option is chosen, the strategy for determining the data-driven thresholds are given by the user-defined strategy above (i.e. globally, within-time, or within-comparison).",
                placement = "right",
                trigger   = "hover"
              )
            ),
            
            # If user threshold, specify it
            conditionalPanel(
              condition = "input.circos_filter_variable != 'none' &&
                input.circos_filter_mode == 'user'",
              tags$label(
                "Threshold for effect size filtering:",
                tags$span(
                  icon("info-circle"),
                  id    = "circos_user_threshold_info",       # target for tooltip
                  style = "cursor: help; color: #337ab7; margin-left: 4px;"
                )
              ),
              numericInput(
                inputId = "circos_user_threshold",
                label   = NULL,
                value   = 0.5,
                min     = 0,
                max     = 3,
                step    = 0.1
              ),
              bsTooltip(
                id      = "circos_user_threshold_info",
                title   = "Results whose chosen scores absolute values are below this threshold will not be defined as statistically significant.",
                placement = "right",
                trigger   = "hover"
              )
            ),
            
            # If data-driven, which quantile to filter by
            conditionalPanel(
              condition = "input.circos_filter_variable != 'none' &&
                input.circos_filter_mode == 'data'",
              tags$label(
                "Quantile for effect size filtering:",
                tags$span(
                  icon("info-circle"),
                  id    = "circos_quantile_threshold_info",       # target for tooltip
                  style = "cursor: help; color: #337ab7; margin-left: 4px;"
                )
              ),
              numericInput(
                inputId = "circos_quantile_threshold",
                label   = NULL,
                value   = 0.5,
                min     = 0,
                max     = 1,
                step    = 0.1
              ),
              bsTooltip(
                id      = "circos_quantile_threshold_info",
                title   = "Results whose chosen scores absolute values are below the threshold(s) given by this quantile will not be defined as statistically significant.",
                placement = "right",
                trigger   = "hover"
              )
            )
          ),
          
          # Visualization settings group
          wellPanel(
            style = "background-color: #ffcad4;",
            h4("Visualisation Settings"),
            
            # Which scores to display in the second ring
            tags$label(
              "Which scores to display?",
              tags$span(
                icon("info-circle"),
                id    = "circos_scores_info",       # target for tooltip
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            selectInput(
              inputId = "circos_scores",
              label   = NULL,
              choices = scores_set,
              selected = "fc.score"
            ),
            bsTooltip(
              id      = "circos_scores_info",
              title   = "Which scores to display on the second ring of the circos plot. Note that this choice is purely visual and is independent from the determination of statistical significance.",
              placement = "right",
              trigger   = "hover"
            ),
            
            # Which correlation scores to display in the outer ring
            tags$label(
              "Correlation metric to Display:",
              tags$span(
                icon("info-circle"),
                id    = "circos_correlation_info",       # target for tooltip
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            selectInput(
              inputId = "circos_correlation",
              label   = NULL,
              choices = correlation_set,
              selected = "corr.mean"
            ),
            bsTooltip(
              id      = "circos_correlation_info",
              title   = "Which correlation metric with a given response to display on the first ring of the circos plot. The options are the mean of genewise correlations within each set, or the correlation of the mean expression of the genes in each set.",
              placement = "right",
              trigger   = "hover"
            ),
            
            # What is the rule for plotting an arc between segments?
            tags$label(
              "What is the rule for plotting arcs?",
              tags$span(
                icon("info-circle"),
                id    = "circos_arc_info",       # target for tooltip
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            selectInput(
              inputId = "circos_arc",
              label   = NULL,
              choices = arc_set,
              selected = "direction"
            ),
            bsTooltip(
              id      = "circos_arc_info",
              title   = "Rule for plotting arcs between segments.",
              placement = "right",
              trigger   = "hover"
            ),
            
            # Which inner rings should be shown?
            tags$label(
              "Which inner rings should be shown?",
              tags$span(
                icon("info-circle"),
                id    = "circos_ring_info",       # target for tooltip
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            selectInput(
              inputId = "circos_ring",
              label   = NULL,
              choices = circos_ring_set,
              selected = "expression"
            ),
            bsTooltip(
              id      = "circos_ring_info",
              title   = "User can specify which inner rings of the circos plot should be shown. The options are correlation + expression, expression only, or none.",
              placement = "right",
              trigger   = "hover"
            ),
            
            # What is the order of the conditions? Fixed or adapted to data-availability
            tags$label(
              "Order Segments By:",
              tags$span(
                icon("info-circle"),
                id    = "circos_order_info",       # target for tooltip
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            selectInput(
              inputId = "circos_order",
              label   = NULL,
              choices = circos_order_set,
              selected = "set_all"
            ),
            bsTooltip(
              id      = "circos_order_info",
              title   = "Which conditions to show in the circos plot? This can either be a fixed order, or adapted based on the conditions who have available data at the selected timepoint.",
              placement = "right",
              trigger   = "hover"
            ),
            
            # Should the scores be clipped to avoid distortion from extreme values?
            tags$label(
              "Clip Score At Quantile:",
              tags$span(
                icon("info-circle"),
                id    = "circos_quantile_scoreclip_info",       # target for tooltip
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            numericInput(
              inputId = "circos_quantile_scoreclip",
              label   = NULL,
              value   = 0.99,
              min     = 0,
              max     = 1,
              step    = 0.01
            ),
            bsTooltip(
              id      = "circos_quantile_scoreclip_info",
              title   = "Should the scores in the second layer be clipped to avoid extreme values distorting the scale? If so, at which quantile should the scores be clipped? For example, if the user selects 0.95, every score above the 95th percentile will be limited to that value.",
              placement = "right",
              trigger   = "hover"
            ),
          ),
          
          # Download section
          wellPanel(
            style = "background-color: #ced4da;",h4("Download Options"),
            # choose format
            selectInput(
              inputId = "circos_download_format",
              label   = "Format of download:",
              choices = c("PDF" = "pdf", "PNG" = "png", "JPEG" = "jpeg"),
              selected = "pdf"
            ),
            
            # size controls (in inches)
            numericInput(
              inputId = "circos_download_width",
              label   = "Width of download (inches):",
              value   = 35,
              min     = 1
            ),
            numericInput(
              inputId = "circos_download_height",
              label   = "Height of download (inches):",
              value   = 15,
              min     = 1
            ),
            
            numericInput(
              inputId = "circos_download_dpi",
              label   = "Resolution (DPI):",
              value   = 300,
              min     = 72
            ),
            
            
            
            ## Download plot button
            div(
              style = "width: 100%;",
              downloadButton(
                outputId = "download_circos_plot",
                label    = "Download Circos Plot",
                class    = "btn-success"
              )
            )
          )
        ), # sidebarPanel end
        # Plot settings
        mainPanel(
          width = 9,
          div(
            style = "text-align: center;",
            plotOutput(
              outputId = "circosPlot",
              height = "1000px",
              width  = "160%"
            )
          )
        )
      )
    )
    ,
    
    ### TAB 2 : HEATMAPS ###
    
    tabPanel(
      title = "Heatmaps", value = "heatmap_tab",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          # Action button
          actionButton(
            inputId = "heatmap_Update",
            label   = "Update Heatmap",
            class   = "btn-primary",
            width   = '100%'
          ),
          # Data selection group
          wellPanel(
            style = "background-color: #c6d2ed;",
            h4("Data Selection"),
            
            # Results from which method?
            tags$label(
              "DGSA Method:",
              tags$span(
                icon("info-circle"),
                id    = "heatmap_method_info",
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            selectInput(
              inputId = "heatmap_method",
              label   = NULL,
              choices = methods_set,
              selected = unique(results_df$method)[1],
              multiple = FALSE
            ),
            bsTooltip(
              id      = "heatmap_method_info",
              title   = "Choose the method of differential gene set expression analysis whose results are to be viewed.",
              placement = "right",
              trigger   = "hover"
            ),
            
            selectInput(
              inputId = "heatmap_conditions",
              label   = "Conditions:",
              choices = conditions_set,
              selected = conditions_set,
              multiple = TRUE
            ),
            selectInput(
              inputId = "heatmap_times",
              label   = "Time (days):",
              choices = timepoints_set,
              selected = timepoints_set,
              multiple = TRUE
            ),
            selectInput(
              inputId = "heatmap_aggregates",
              label   = "Aggregates:",
              choices = aggregates_set,
              selected = aggregates_set,
              multiple = TRUE
            )
          ),
          
          # Statistical settings group
          wellPanel(
            style = "background-color: #cce3de;",
            h4("Statistical Significance"),
            selectInput(
              inputId = "heatmap_p_correction",
              label   = "P-value Correction:",
              choices = p_correction_set,
              selected = "BH"
            ),
            selectInput(
              inputId = "heatmap_p_approach",
              label   = "Correction & Filtering Strategy:",
              choices = p_approach_set,
              selected = "global"
            ),
            numericInput(
              inputId = "heatmap_p_threshold",
              label   = "P-value Threshold:",
              value   = 0.05,
              min     = 0,
              max     = 1,
              step    = 0.01
            ),
            selectInput(
              inputId = "heatmap_filter_variable",
              label   = "Filter Variable:",
              choices = filtration_set,
              selected = "none"
            ),
            conditionalPanel(
              condition = "input.heatmap_filter_variable != 'none'",
              numericInput(
                inputId = "heatmap_quantile_threshold",
                label   = "Filter Below Quantile:",
                value   = 0.5,
                min     = 0,
                max     = 1,
                step    = 0.1
              )
            )
          ),
          
          # Visualization settings group
          wellPanel(
            style = "background-color: #ffcad4;",
            h4("Visualization Settings"),
            selectInput(
              inputId = "heatmap_scores",
              label   = "Scores to Display:",
              choices = scores_set,
              selected = "activation.score"
            ),
            selectInput(
              inputId  = "heatmap_y_order",
              label    = "Order y axis by:",
              choices  = y_order_set,
              selected = "cluster"
            ),
            selectInput(
              inputId  = "heatmap_x_order",
              label    = "Order x axis by:",
              choices  = x_order_set,
              selected = "cluster-time"
            ),
            selectInput(
              inputId  = "heatmap_filter_commonDE",
              label    = "Filter genesets to visualise by common DE:",
              choices  = filter_commonDE_set,
              selected = "none"
            ),
            
            # Panel for the "global" or "withinTime" choices
            conditionalPanel(
              condition = 'input.heatmap_filter_commonDE == "global" ||
               input.heatmap_filter_commonDE == "withinTime"',
              numericInput(
                inputId = "heatmap_common_percentage",
                label   = "Filter on which fraction of common DE:",
                value   = 0.25,
                min     = 0,
                max     = 1,
                step    = 0.05
              )
            ),
            
            # Panel for the "score" choice
            conditionalPanel(
              condition = 'input.heatmap_filter_commonDE == "score"',
              wellPanel(
                numericInput(
                  inputId = "heatmap_score_threshold",
                  label   = "Filter Below Sharing Score Threshold:",
                  value   = 5,
                  min     = 0,
                  max     = 13,
                  step    = 1
                )
              )
            ),
            numericInput(
              inputId = "heatmap_quantile_scoreclip",
              label   = "Clip Score At Quantile:",
              value   = 0.995,
              min     = 0,
              max     = 1,
              step    = 0.01
            )
          ),
          # Download section
          wellPanel(
            style = "background-color: #ced4da;",h4("Download Options"),
            # choose format
            selectInput(
              inputId = "heatmap_download_format",
              label   = "Format of download:",
              choices = c("PDF" = "pdf", "PNG" = "png", "JPEG" = "jpeg"),
              selected = "pdf"
            ),
            
            # size controls (in inches)
            numericInput(
              inputId = "heatmap_download_width",
              label   = "Width of download (inches):",
              value   = 35,
              min     = 1
            ),
            numericInput(
              inputId = "heatmap_download_height",
              label   = "Height of download (inches):",
              value   = 15,
              min     = 1
            ),
            
            numericInput(
              inputId = "heatmap_download_dpi",
              label   = "Resolution (DPI):",
              value   = 300,
              min     = 72
            ),
            
            
            
            ## Download plot button
            div(
              style = "width: 100%;",
              downloadButton(
                outputId = "download_heatmap",
                label    = "Download Heatmap",
                class    = "btn-success"
              )
            )
          )
        ),
        
        mainPanel(
          width = 9,
          div(
            style = "text-align: center;",
            plotOutput(
              outputId = "heatmapPlot",
              height = "1000px",
              width  = "160%"
            )
          )
        )
      )
    ),
    
    ## TAB 3 : Expression profiles over time ##
    
    tabPanel(
      title = "Expression Profiles", value = "exprprofiles_tab",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          # Action button
          actionButton(
            inputId = "exprprofiles_Update",
            label   = "Update Expression Profiles",
            class   = "btn-primary",
            width   = '100%'
          ),
          # Data selection group
          wellPanel(
            style = "background-color: #c6d2ed;",
            h4("Data Selection"),
            
            # Results from which method?
            tags$label(
              "DGSA Method:",
              tags$span(
                icon("info-circle"),
                id    = "exprprofiles_method_info",
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            selectInput(
              inputId = "exprprofiles_method",
              label   = NULL,
              choices = methods_set,
              selected = unique(results_df$method)[1],
              multiple = FALSE
            ),
            bsTooltip(
              id      = "exprprofiles_method_info",
              title   = "Choose the method of differential gene set expression analysis whose results are to be viewed.",
              placement = "right",
              trigger   = "hover"
            ),
            
            selectInput(
              inputId = "exprprofiles_conditions",
              label   = "Conditions:",
              choices = conditions_set,
              selected = conditions_set,
              multiple = TRUE
            ),
            selectInput(
              inputId = "exprprofiles_times",
              label   = "Time (days):",
              choices = timepoints_set,
              selected = timepoints_set[1],
              multiple = TRUE
            ),
            selectInput(
              inputId = "exprprofiles_aggregates",
              label   = "Aggregates:",
              choices = aggregates_set,
              selected = aggregates_set,
              multiple = TRUE
            )
          ),
          
          # Visualization settings group
          wellPanel(
            style = "background-color: #ffcad4;",
            h4("Visualization Settings"),
            selectInput(
              inputId  = "exprprofiles_scores",
              label    = "Scores to Display:",
              choices  = scores_set,
              selected = "activation.score"
            ),
            numericInput(
              inputId = "exprprofiles_quantile_scoreclip",
              label   = "Clip Score At Quantile:",
              value   = 0.95,
              min     = 0,
              max     = 1,
              step    = 0.01
            ),
            selectInput(inputId = "exprprofiles_barcolours",
                        label = "Colour bars by...",
                        choices = exprprofiles_barcolours_set,
                        selected = "direction")
          ),
          # Download section
          wellPanel(
            style = "background-color: #ced4da;",h4("Download Options"),
            # choose format
            selectInput(
              inputId = "exprprofiles_download_format",
              label   = "Format of download:",
              choices = c("PDF" = "pdf", "PNG" = "png", "JPEG" = "jpeg"),
              selected = "pdf"
            ),
            
            # size controls (in inches)
            numericInput(
              inputId = "exprprofiles_download_width",
              label   = "Width of download (inches):",
              value   = 35,
              min     = 1
            ),
            numericInput(
              inputId = "exprprofiles_download_height",
              label   = "Height of download (inches):",
              value   = 15,
              min     = 1
            ),
            
            numericInput(
              inputId = "exprprofiles_download_dpi",
              label   = "Resolution (DPI):",
              value   = 300,
              min     = 72
            ),
            
            ## Download plot button
            div(
              style = "width: 100%;",
              downloadButton(
                outputId = "download_exprprofiles",
                label    = "Download exprprofiles",
                class    = "btn-success"
              )
            )
          )
        ),
        
        mainPanel(
          width = 9,
          div(
            style = "text-align: center;",
            plotOutput(
              outputId = "exprprofilesPlot",
              height = "1000px",
              width  = "160%"
            )
          )
        )
      )
    ),
    
    ### TAB 4 : SPIDER PLOTS ###
    tabPanel(
      title = "Spider Plots", value = "spider_tab",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          # Action button
          actionButton(
            inputId = "spider_Update",
            label   = "Update spider Plot",
            class   = "btn-primary",
            width   = '100%'
          ),
          # Data selection group
          wellPanel(
            style = "background-color: #c6d2ed;",
            h4("Data Selection"),
            
            # Results from which method?
            tags$label(
              "DGSA Method:",
              tags$span(
                icon("info-circle"),
                id    = "spider_method_info",
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            selectInput(
              inputId = "spider_method",
              label   = NULL,
              choices = methods_set,
              selected = unique(results_df$method)[1],
              multiple = FALSE
            ),
            bsTooltip(
              id      = "spider_method_info",
              title   = "Choose the method of differential gene set expression analysis whose results are to be viewed.",
              placement = "right",
              trigger   = "hover"
            ),
            
            selectInput(inputId = "spider_grouping",
                        label = "Which plots to show?",
                        choices = spider_grouping_set,
                        selected = "withinTime"),
            
            conditionalPanel(
              condition = "input.spider_grouping == 'withinTime'",
              selectInput(inputId = "spider_conditions_withinTime",
                          label = "Which conditions to display?",
                          choices = conditions_set,
                          selected = conditions_set,
                          multiple = TRUE),
              selectInput(inputId = "spider_times_withinTime",
                          label = "Which time to display?",
                          choices = timepoints_set,
                          selected = timepoints_set[1],
                          multiple = FALSE)
            ),
            
            conditionalPanel(
              condition = "input.spider_grouping == 'withinCondition'",
              selectInput(inputId = "spider_conditions_withinCondition",
                          label = "Which condition to display?",
                          choices = conditions_set,
                          selected = conditions_set[1],
                          multiple = FALSE),
              selectInput(inputId = "spider_times_withinCondition",
                          label = "Which times to display?",
                          choices = timepoints_set,
                          selected = timepoints_set,
                          multiple = TRUE)
            ),
            selectInput(
              inputId = "spider_aggregates",
              label   = "Aggregates:",
              choices = aggregates_set,
              selected = aggregates_set,
              multiple = TRUE
            )
          ),
          # Visualization settings group
          wellPanel(
            style = "background-color: #ffcad4;",
            h4("Visualization Settings"),
            selectInput(
              inputId = "spider_scores",
              label   = "Scores to Display:",
              choices = scores_set,
              selected = "activation.score"
            ),
            numericInput(
              inputId = "spider_quantile_scoreclip",
              label   = "Clip Score At Quantile:",
              value   = 0.95,
              min     = 0,
              max     = 1,
              step    = 0.01
            ),
            selectInput(
              inputId = "spider_grid",
              label = "Fix grid order?",
              choices = spider_grid_set,
              selected = spider_grid_set[2])
          ),
          # Download section
          wellPanel(
            style = "background-color: #ced4da;",h4("Download Options"),
            # choose format
            selectInput(
              inputId = "spider_download_format",
              label   = "Format of download:",
              choices = c("PDF" = "pdf", "PNG" = "png", "JPEG" = "jpeg"),
              selected = "pdf"
            ),
            
            # size controls (in inches)
            numericInput(
              inputId = "spider_download_width",
              label   = "Width of download (inches):",
              value   = 35,
              min     = 1
            ),
            numericInput(
              inputId = "spider_download_height",
              label   = "Height of download (inches):",
              value   = 15,
              min     = 1
            ),
            
            numericInput(
              inputId = "spider_download_dpi",
              label   = "Resolution (DPI):",
              value   = 300,
              min     = 72
            ),
            
            ## Download plot button
            div(
              style = "width: 100%;",
              downloadButton(
                outputId = "download_spider",
                label    = "Download spider",
                class    = "btn-success"
              )
            )
          )
        ),
        
        mainPanel(
          width = 9,
          div(
            style = "text-align: center;",
            plotOutput(
              outputId = "spiderPlot",
              height = "1000px",
              width  = "160%"
            )
          )
        )
      )
    ),
    
    ### TAB 5 : AGGREGATE SENSITIVITY HEATMAP ###
    
    tabPanel(
      title = "Sensitivity Heatmap", value = "sensivity_heatmap_tab",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          # Action button
          actionButton(
            inputId = "sensitivity_heatmap_Update",
            label   = "Update sensitivity heatmap",
            class   = "btn-primary",
            width   = '100%'
          ),
          # Data selection group
          wellPanel(
            style = "background-color: #c6d2ed;",
            h4("Data Selection"),
            
            # Results from which method?
            tags$label(
              "DGSA Method:",
              tags$span(
                icon("info-circle"),
                id    = "sensitivity_heatmap_method_info",
                style = "cursor: help; color: #337ab7; margin-left: 4px;"
              )
            ),
            selectInput(
              inputId = "sensitivity_heatmap_method",
              label   = NULL,
              choices = methods_set,
              selected = unique(results_df$method),
              multiple = TRUE
            ),
            bsTooltip(
              id      = "sensitivity_heatmap_method_info",
              title   = "Choose the method of differential gene set expression analysis whose results are to be viewed.",
              placement = "right",
              trigger   = "hover"
            ),
            
            selectInput(
              inputId = "sensitivity_heatmap_conditions",
              label   = "Conditions:",
              choices = conditions_set,
              selected = conditions_set,
              multiple = TRUE
            ),
            selectInput(
              inputId = "sensitivity_heatmap_time",
              label   = "Time (days):",
              choices = timepoints_set,
              selected = 1,
              multiple = FALSE
            ),
            selectInput(
              inputId = "sensitivity_heatmap_aggregates",
              label   = "Aggregates:",
              choices = aggregates_set,
              selected = aggregates_set,
              multiple = TRUE
            )
          ),
          # Specification parameters
          wellPanel(
            style = "background-color: #b3d89c;",
            h4("Specifications"),
            selectInput(
              inputId = "sensitivity_heatmap_p_approach",
              label   = "Subset of p-values to correct:",
              choices = p_approach_set,
              selected = p_approach_set,
              multiple = TRUE
            ),
            selectInput(
              inputId = "sensitivity_heatmap_p_method",
              label   = "P-value correction method:",
              choices = sensitivity_heatmap_p_correction_set,
              selected = sensitivity_heatmap_p_correction_set,
              multiple = TRUE
            ),
            selectInput(
              inputId = "sensitivity_heatmap_p_threshold",
              label   = "Adjusted p-value significance threshold",
              choices = p_thresholds_set,
              selected = c(0.0001, 0.001, 0.01, 0.05, 0.1),
              multiple = TRUE
            ),
            selectInput(
              inputId = "sensitivity_heatmap_fc_threshold",
              label   = "Absolute Log2(Mean FC) significance threshold",
              choices = fc_thresholds_set,
              selected = seq(0,1,0.2),
              multiple = TRUE
            ),
          ),
          # Visualization settings group
          wellPanel(
            style = "background-color: #ffcad4;",
            h4("Visualization Settings"),
            selectInput(
              inputId = "sensitivity_heatmap_order",
              label   = "Which conditions to display?",
              choices = spider_grid_set,
              selected= "set"
            )
          ),
          # Download section
          wellPanel(
            style = "background-color: #ced4da;",h4("Download Options"),
            # choose format
            selectInput(
              inputId = "sensitivity_heatmap_download_format",
              label   = "Format of download:",
              choices = c("PDF" = "pdf", "PNG" = "png", "JPEG" = "jpeg"),
              selected = "pdf"
            ),
            
            # size controls (in inches)
            numericInput(
              inputId = "sensitivity_heatmap_download_width",
              label   = "Width of download (inches):",
              value   = 24,
              min     = 1
            ),
            numericInput(
              inputId = "sensitivity_heatmap_download_height",
              label   = "Height of download (inches):",
              value   = 20,
              min     = 1
            ),
            
            numericInput(
              inputId = "sensitivity_heatmap_download_dpi",
              label   = "Resolution (DPI):",
              value   = 300,
              min     = 72
            ),
            
            ## Download plot button
            div(
              style = "width: 100%;",
              downloadButton(
                outputId = "download_sensitivity_heatmap",
                label    = "Download sensitivity_heatmap",
                class    = "btn-success"
              )
            )
          )
        ),
        
        mainPanel(
          width = 9,
          div(
            style = "text-align: center;",
            plotOutput(
              outputId = "sensitivityHeatmapPlot",
              height = "1000px",
              width  = "160%"
            )
          )
        )
      )
    ), 
    
    ### END AGGREGATE SENSITIVITY HEATMAP ###
    
    ### TAB 6 : INDIVIDUAL PARAMETER SENSITIVITY ###
    
    tabPanel(
      title = "Individual Parameter Sensitivity", 
      value = "indiv_sensitivity_tab",
      
      sidebarLayout(
        
        # ← sidebarPanel now wraps *both* wellPanels
        sidebarPanel(
          width = 3,
          
          # Action button
          actionButton(
            inputId = "indiv_sensitivity_Update",
            label   = "Update individual parameter sensitivity plot",
            class   = "btn-primary",
            width   = '100%'
          ),
          
          # Data selection group
          wellPanel(
            style = "background-color: #c6d2ed;",
            h4("Data Selection"),
            selectInput(
              inputId = "indiv_sensitivity_conditions",
              label   = "Conditions:",
              choices = conditions_set,
              selected = conditions_set,
              multiple = TRUE
            ),
            selectInput(
              inputId = "indiv_sensitivity_time",
              label   = "Time (days):",
              choices = timepoints_set,
              selected = 1
            ),
            selectInput(
              inputId = "indiv_sensitivity_aggregates",
              label   = "Aggregates:",
              choices = aggregates_set,
              selected = aggregates_set,
              multiple = TRUE
            )
          ),
          
          # Specifications group
          wellPanel(
            style = "background-color: #b3d89c;",
            h4("Specifications"),
            selectInput(
              inputId = "indiv_sensitivity_parameter",
              label   = "Parameter to vary:",
              choices = indiv_sensitivity_parameter_set,
              selected = indiv_sensitivity_parameter_set[1],
              multiple = FALSE
            ),
            conditionalPanel(
              condition = "input.indiv_sensitivity_parameter == 'p_approach'||input.indiv_sensitivity_parameter == 'p_correction'||input.indiv_sensitivity_parameter == 'method'",
              selectInput(
                inputId = "indiv_sensitivity_groupby",
                label   = "Group by which variable for visualisation?",
                choices = indiv_sensitivity_groupby_set,
                selected = "byvaccine",
                multiple = FALSE
              )
            ),
            conditionalPanel(
              condition = "input.indiv_sensitivity_parameter == 'p_threshold'",
              selectInput(
                inputId = "indiv_sensitivity_p_threshold_range",
                label   = "Adjusted p-value significance threshold range",
                choices = c(0.000001, p_thresholds_set, 0.2),
                selected = p_thresholds_set,
                multiple = TRUE
              )
            ),
            conditionalPanel(
              condition = "input.indiv_sensitivity_parameter == 'fc_threshold'",
              selectInput(
                inputId = "indiv_sensitivity_fc_threshold_range",
                label   = "Absolute log2-Fold-Change threshold range",
                choices = fc_thresholds_set,
                selected = fc_thresholds_set,
                multiple = TRUE
              )
            ),
            # Note the "input." prefix in the JS condition
            conditionalPanel(
              condition = "input.indiv_sensitivity_parameter != 'method'",
              selectInput(
                inputId = "indiv_sensitivity_method",
                label   = "Statistical method:",
                choices = methods_set,
                selected = methods_set[1]
              )
            ),
            conditionalPanel(
              condition = "input.indiv_sensitivity_parameter != 'p_correction'",
              selectInput(
                inputId = "indiv_sensitivity_p_correction",
                label   = "P-value correction method:",
                choices = p_correction_set,
                selected = "BH"
              )
            ),
            conditionalPanel(
              condition = "input.indiv_sensitivity_parameter != 'p_approach'",
              selectInput(
                inputId = "indiv_sensitivity_p_approach",
                label   = "P-value correction subset:",
                choices = p_approach_set,
                selected = "withinTime"
              )
            ),
            conditionalPanel(
              condition = "input.indiv_sensitivity_parameter != 'p_threshold'",
              numericInput(
                inputId = "indiv_sensitivity_p_threshold",
                label   = "Adjusted p-value threshold:",
                value   = 0.05,
                min     = 0,
                max     = 1,
                step    = 0.01
              )
            ),
            conditionalPanel(
              condition = "input.indiv_sensitivity_parameter != 'fc_threshold'",
              numericInput(
                inputId = "indiv_sensitivity_fc_threshold",
                label   = "Absolute log2‑fold‑change threshold:",
                value   = 0,
                min     = 0,
                max     = 3,
                step    = 0.1
              )
            )
          ),  # ← end wellPanel
          # Download section
          wellPanel(
            style = "background-color: #ced4da;",h4("Download Options"),
            # choose format
            selectInput(
              inputId = "indiv_sensitivity_download_format",
              label   = "Format of download:",
              choices = c("PDF" = "pdf", "PNG" = "png", "JPEG" = "jpeg"),
              selected = "pdf"
            ),
            
            # size controls (in inches)
            numericInput(
              inputId = "indiv_sensitivity_download_width",
              label   = "Width of download (inches):",
              value   = 20,
              min     = 1
            ),
            numericInput(
              inputId = "indiv_sensitivity_download_height",
              label   = "Height of download (inches):",
              value   = 10,
              min     = 1
            ),
            
            numericInput(
              inputId = "indiv_sensitivity_download_dpi",
              label   = "Resolution (DPI):",
              value   = 300,
              min     = 72
            ),
            
            ## Download plot button
            div(
              style = "width: 100%;",
              downloadButton(
                outputId = "download_indiv_sensitivity",
                label    = "Download Plot",
                class    = "btn-success"
              )
            )
          )
        ), # ← close sidebarPanel
        
        mainPanel(
          width = 9,
          div(
            style = "text-align: center;",
            plotOutput(
              outputId = "indivSensitivityPlot",
              height = "1000px",
              width  = "160%"
            )
          )
        )
      )  # ← close sidebarLayout
    ),  # ← close tabPanel
    
    ### END INDIVIDUAL SENSITIVITY TAB ###
    
    ### TAB 7 : METHOD-SPECIFIC UNIQUE GENESETS ###
    
    tabPanel(
      title = "Method-Specific Unique Genesets", 
      value = "method_unique_tab",
      
      sidebarLayout(
        
        # ← sidebarPanel now wraps *both* wellPanels
        sidebarPanel(
          width = 3,
          
          # Action button
          actionButton(
            inputId = "method_unique_Update",
            label   = "Update method-specific genesets plot",
            class   = "btn-primary",
            width   = '100%'
          ),
          
          # Data selection group
          wellPanel(
            style = "background-color: #c6d2ed;",
            h4("Data Selection"),
            selectInput(
              inputId = "method_unique_conditions",
              label   = "Conditions:",
              choices = conditions_set,
              selected = conditions_set,
              multiple = TRUE
            ),
            selectInput(
              inputId = "method_unique_time",
              label   = "Time (days):",
              choices = timepoints_set,
              selected = 1,
              multiple = TRUE
            ),
            selectInput(
              inputId = "method_unique_aggregates",
              label   = "Aggregates:",
              choices = aggregates_set,
              selected = aggregates_set,
              multiple = TRUE
            ),
            selectInput(
              inputId = "method_unique_reference",
              label   = "Reference DGSA Method:",
              choices = methods_set,
              selected = methods_set[2]
            ),
            selectInput(
              inputId = "method_unique_comparison",
              label   = "DGSA Method to compare with reference:",
              choices = methods_set,
              selected = methods_set[1]
            ),
          ), # End wellPanel
          
          wellPanel(
            style = "background-color: #cce3de;",
            h4("Statistical Significance"),
            selectInput(
              inputId = "method_unique_p_approach",
              label   = "Subset of p-values to correct:",
              choices = p_approach_set,
              selected = "withinTime"
            ),
            selectInput(
              inputId = "method_unique_p_correction",
              label   = "P-value correction method:",
              choices = p_correction_set,
              selected = "BH"
            ),
            numericInput(
              inputId = "method_unique_p_threshold",
              label   = "Adjusted p-value threshold:",
              value   = 0.05,
              min     = 0,
              max     = 1,
              step    = 0.01
            ),
            numericInput(
              inputId = "method_unique_fc_threshold",
              label   = "Absolute log2‑fold‑change threshold:",
              value   = 0,
              min     = 0,
              max     = 3,
              step    = 0.1
            ),
          ), # End wellPanel
          
          wellPanel(
            style = "background-color: #ffcad4;",
            h4("Visualisation Settings"),
            selectInput(
              inputId = "method_unique_filtration",
              label   = "Filter results for visualisation?",
              choices = method_unique_filtration_set,
              selected = "none"
            ),
          ), # End wellPanel
          
          # Download section
          wellPanel(
            style = "background-color: #ced4da;",
            h4("Download Options"),
            # choose format
            selectInput(
              inputId = "method_unique_download_format",
              label   = "Format of download:",
              choices = c("PDF" = "pdf", "PNG" = "png", "JPEG" = "jpeg"),
              selected = "pdf"
            ),
            
            # size controls (in inches)
            numericInput(
              inputId = "method_unique_download_width",
              label   = "Width of download (inches):",
              value   = 20,
              min     = 1
            ),
            numericInput(
              inputId = "method_unique_download_height",
              label   = "Height of download (inches):",
              value   = 10,
              min     = 1
            ),
            
            numericInput(
              inputId = "method_unique_download_dpi",
              label   = "Resolution (DPI):",
              value   = 300,
              min     = 72
            ),
            
            ## Download plot button
            div(
              style = "width: 100%;",
              downloadButton(
                outputId = "download_method_unique",
                label    = "Download Plot",
                class    = "btn-success"
              )
            )
          )
        ), # End sidebarPanel
        
        mainPanel(
          width = 9,
          div(
            style = "text-align: center;",
            plotOutput(
              outputId = "methodUniquePlot",
              height = "1000px",
              width  = "160%"
            )
          )
        )
        
      ) # End sidebarLayout
      
    ) # End tabPanel
    
    ### END METHOD-SPECIFIC UNIQUE GENESETS ###
    
  ) # Close tabsetPanel
) # Close fluidPage

# Example parameters for debugging
input = list()

# Heatmap tab
input$heatmap_method = "dearseq"
input$heatmap_conditions = conditions_set
input$heatmap_times = timepoints_set
input$heatmap_aggregates = aggregates_set
input$heatmap_p_correction = "BH"
input$heatmap_p_approach = "global"
input$heatmap_p_threshold = 0.05
input$heatmap_filter_variable = "none"
input$heatmap_quantile_threshold = 0
input$heatmap_scores = "fc.score"
input$heatmap_y_order = "aggregate"
input$heatmap_x_order = "cluster-time"
input$heatmap_filter_commonDE = "score"
input$heatmap_common_percentage = 0
input$heatmap_score_threshold = 8
input$heatmap_quantile_scoreclip = 0.995

# Spider plot tab
input$spider_method = unique(results_df$method)[1]
input$spider_grouping = "withinTime"
input$spider_conditions_withinTime = conditions_set
input$spider_times_withinTime = timepoints_set[1]
input$spider_conditions_withinCondition = conditions_set[1]
input$spider_conditions_withinCondition = conditions_set[1]
input$spider_aggregates = aggregates_set
input$spider_scores = "fc.score"
input$spider_quantile_scoreclip = 0.95
input$spider_grid = spider_grid_set[2]

# Sensitivity heatmap tab
input$sensitivity_heatmap_method = unique(results_df$method)
input$sensitivity_heatmap_conditions = conditions_set
input$sensitivity_heatmap_time = timepoints_set[4]
input$sensitivity_heatmap_aggregates = aggregates_set
input$sensitivity_heatmap_p_approach =p_approach_set
input$sensitivity_heatmap_p_method = sensitivity_heatmap_p_correction_set
input$sensitivity_heatmap_p_threshold =  c(0.0001, 0.001, 0.01, 0.05, 0.1)
input$sensitivity_heatmap_fc_threshold = seq(0,1,0.2)
input$sensitivity_heatmap_order = "set"

# Individual parameter sensitivity tab
input$indiv_sensitivity_conditions = conditions_set
input$indiv_sensitivity_time = 1
input$indiv_sensitivity_aggregates = aggregates_set
input$indiv_sensitivity_parameter = "method"
input$indiv_sensitivity_method = "dearseq"
input$indiv_sensitivity_p_correction = "BH"
input$indiv_sensitivity_p_approach = "withinTime"
input$indiv_sensitivity_p_threshold = 0.05
input$indiv_sensitivity_fc_threshold = 0
input$indiv_sensitivity_p_threshold_range = c(1e-05, 1e-04, 1e-03, 1e-02, 5e-02, 1e-01)
input$indiv_sensitivity_fc_threshold_range = fc_thresholds_set
input$indiv_sensitivity_groupby = "byvaccine"

# Unique per-method genesets tab
input$method_unique_conditions = conditions_set
input$method_unique_time = c(1,3)
input$method_unique_aggregates = aggregates_set
input$method_unique_p_correction = "BH"
input$method_unique_p_approach = "withinTime"
input$method_unique_p_threshold = 0.05
input$method_unique_fc_threshold = 0
input$method_unique_reference = "qusage"
input$method_unique_comparison = "dearseq"
input$method_unique_filtration = "positive"

server <- function(input, output, session) {
  
  ### BEGIN CIRCOS PLOTS TAB ###
  
  # Store circos plot for downloading
  lastCircos <- reactiveVal(NULL)
  
  output$circosPlot <- renderPlot({
    
    input$circos_Update
    isolate({
      
      # Subset the results by the user selected method, conditions and aggregates
      results_df_circos = results_df %>%
        filter(gs.aggregate %in%  input$circos_aggregates,
               condition %in% input$circos_conditions,
               method == input$circos_method) %>%
        mutate(gs.aggregate = droplevels(gs.aggregate))
      
      # Based on the remaining genesets, define a global ordering
      global_order <- results_df_circos %>%
        select(gs.name, gs.aggregate) %>%
        distinct() %>%
        arrange(gs.aggregate) %>%
        pull(gs.name)
      
      # Apply this global gene set ordering to the dataframe
      results_df_circos <- results_df_circos %>%
        mutate(gs.name = factor(gs.name, levels = global_order))
      
      # Extract a list of conditions to plot
      ## If user selects all to appear (circos_order == set_all), they should all appear regardless of if there are significant results
      if (input$circos_order == "set_available" | input$circos_order == "cluster"){
        condition_names = results_df_circos %>%
          filter(time == input$circos_time) %>%
          pull(condition) %>%
          unique()
      } else if (input$circos_order == "set_all"){
        condition_names <- unique(results_df$condition)
      }
      
      # Drop levels which are not to be included
      condition_order = levels(results_df$condition)
      condition_order = condition_order[condition_order %in% condition_names]
      
      # Create a dataframe containing the metadata for the condition_names
      # This contains the ordering of the conditions in the circos plot, the availability of data, and the availability of the response correlation
      
      if (input$circos_order == "set_all" | input$circos_order == "set_available"){
        # If the user selects a fixed ordering :
        
        condition_circos_metadata <- data.frame(
          order = match(condition_names, condition_order),
          condition = factor(condition_names, levels = condition_order),
          stringsAsFactors = FALSE
        ) %>%
          arrange(condition)
        
      } else if (input$circos_order == "cluster"){
        # Else if the user wants the order to be determined by hierarchical clustering of scores
        sc <- input$circos_scores
        mat_df <- results_df_circos %>%
          filter(time == input$circos_time) %>%
          select(condition, gs.name, score = .data[[sc]]) %>%
          pivot_wider(
            names_from  = gs.name,
            values_from = score
          ) %>%
          # now every missing goes to NA; replace them with 0
          mutate(across(-condition, ~tidyr::replace_na(.x, 0)))
        
        # 2) Turn into a numeric matrix with condition names as row names
        mat <- mat_df %>%
          column_to_rownames("condition") %>%
          as.matrix()
        
        # 3) Compute distances & cluster
        #    Here we use Euclidean; you could also do:
        #      d <- as.dist(1 - cor(t(mat)))    # 1 - Pearson correlation
        d <- dist(mat, method = "euclidean")
        hc <- hclust(d, method = "ward.D2")   # or "complete", "ward.D2", ...
        
        # 4) Extract the condition order
        ordered_conditions <- hc$labels[hc$order]
        
        condition_circos_metadata <- data.frame(
          order = match(condition_names, ordered_conditions),
          condition = factor(condition_names, levels = ordered_conditions),
          stringsAsFactors = FALSE
        ) %>%
          arrange(condition)
      }
      
      # Create columns in the condition metadata to describe the condition availability at each timepoint (gene expression + immune response)
      unique_times <- unique(results_df_circos$time)
      
      condition_colors <- results_df %>%
        select(condition, condition.colour) %>%
        distinct() %>%
        deframe()
      
      for (t in unique_times) {
        # Get condition_names for the expression condition
        conditions_expr_at_time <- results_df_circos %>%
          filter(time == t) %>%
          pull(condition) %>%
          unique()
        
        # Construct dynamic column names for expression availability and color
        expr_col_name <- paste0("condition_expr_available_", t)
        color_col <- paste0("condition_color_", t)
        text_col <- paste0("text_color_", t)
        
        # Get condition_names for the response condition (with non-NA cor_MFC_mean)
        conditions_resp_at_time <- results_df %>%
          filter(time == t, !is.na(corr.mean)) %>%
          pull(condition) %>%
          unique()
        
        # Construct the dynamic column name for response availability
        resp_col_name <- paste0("condition_response_available_", t)
        
        # Update condition_circos_metadata with the new columns
        condition_circos_metadata <- condition_circos_metadata %>%
          mutate(
            !!expr_col_name := condition %in% conditions_expr_at_time,
            !!resp_col_name := condition %in% conditions_resp_at_time,
            !!color_col := if_else(.data[[expr_col_name]],
                                   condition_colors[as.character(condition)],
                                   "#D3D3D3"),
            !!text_col  := if_else(.data[[expr_col_name]],
                                   "black",
                                   "#9999a1")
          )
      }
      
      # Set circos options
      condition_circos_metadata$xmin <- 0
      condition_circos_metadata$xmax <- length(global_order)
      
      # Map gene sets to fixed positions based on global order
      gene_set_positions <- setNames(seq(0.5, length(global_order) - 0.5,
                                         length.out = length(global_order)),
                                     global_order)
      
      # Define colours for each aggregate
      aggregate_colors <- results_df %>%
        select(gs.aggregate, gs.colour) %>%
        distinct() %>%
        { setNames(.$gs.colour, .$gs.aggregate) }
      
      # Compute metadata for the positions of the gene set aggregates on the segments
      aggregate_distinct <- results_df_circos %>%
        select(gs.name, gs.aggregate) %>%
        distinct()
      
      pt <- prop.table(table(aggregate_distinct$gs.aggregate))
      
      aggregate_proportions <- data.frame(
        gs.aggregate = names(pt),
        prop         = as.vector(pt),
        stringsAsFactors = FALSE
      )
      aggregate_proportions$colour = aggregate_colors[aggregate_proportions$gs.aggregate]
      
      # Compute cumulative proportions for each aggregate
      aggregate_proportions$cum_start <- c(0, head(cumsum(aggregate_proportions$prop), -1))
      aggregate_proportions$cum_end <- cumsum(aggregate_proportions$prop)
      
      # FILTRATION : filter significant results based on user inputs
      if (input$circos_p_correction == "none"){
        
        results_df_significant = results_df_circos %>%
          filter(rawPval < input$circos_p_threshold,
                 time == input$circos_time)
      } else {
        
        results_df_significant = results_df_circos %>%
          filter(
            !!sym(paste0(input$circos_p_approach, ".adjPval_", input$circos_p_correction)) < input$circos_p_threshold,
            time == input$circos_time
          )
        
      }
      
      # Filter results by effect size if desired
      if(input$circos_filter_variable != "none" & input$circos_filter_mode == "data"){
        # Global threshold for filtration
        if (input$circos_p_approach == "global"){
          effect_size_threshold <- results_df_circos %>%
            pull(input$circos_filter_variable) %>%
            abs() %>%
            quantile(input$circos_quantile_threshold)
          results_df_significant = results_df_significant %>%
            filter(abs(.data[[input$circos_filter_variable]]) > effect_size_threshold)
          
        } else if (input$circos_p_approach == "withinTime") {
          effect_size_threshold <- results_df_circos %>%
            filter(time == input$circos_time) %>%
            pull(input$circos_filter_variable) %>%
            abs() %>%
            quantile(input$circos_quantile_threshold)
          
          results_df_significant = results_df_significant %>%
            filter(abs(.data[[input$circos_filter_variable]]) > effect_size_threshold)
          
        } else if (input$circos_p_approach == "withinCondition"){
          results_df_significant = results_df_circos %>%
            filter(time == input$circos_time) %>%
            group_by(condition) %>%
            mutate(abs_score = abs(.data[[input$circos_filter_variable]])) %>%
            filter(abs_score >= quantile(abs_score, probs = input$circos_quantile_threshold, na.rm = TRUE)) %>%
            ungroup()
        }
      } else if (input$circos_filter_variable != "none" & input$circos_filter_mode == "user"){
        results_df_significant = results_df_significant %>%
          filter(abs(.data[[input$circos_filter_variable]]) > input$circos_user_threshold)
      }
      
      
      # Identify the shared differentially expressed modules
      if (input$circos_arc == "any"){
        shared_modules <- results_df_significant %>%
          group_by(gs.name) %>%
          filter(n() > 1) %>%
          summarise(condition_names = list(unique(condition)), .groups = "drop")
        
        any_links = (nrow(shared_modules) != 0)
        
        if (any_links) {
          # Prepare the data in the format to plot
          shared_links <- shared_modules %>%
            mutate(pairs = purrr::map(condition_names, ~combn(.x, 2, simplify = FALSE))) %>%
            select(gs.name, pairs) %>%
            unnest(pairs) %>%
            unnest_wider(pairs, names_sep = "_") %>%
            rename(condition1 = pairs_1, condition2 = pairs_2) %>%
            mutate(position = gene_set_positions[gs.name])  # Add position for gene
        }
      } else if (input$circos_arc == "direction"){
        column_name1 <- paste0(input$circos_scores, 1)
        column_name2 <- paste0(input$circos_scores, 2)
        
        # 1. Identify gene‐sets with >1 significant hit and nest their data
        shared_links <- results_df_significant %>%
          # keep only the columns we need
          select(gs.name, condition, .data[[input$circos_scores]]) %>%
          
          # self-join on gs.name to get all condition pairs
          inner_join(
            results_df_significant %>% select(gs.name,
                                              condition,
                                              .data[[input$circos_scores]]),
            by = "gs.name",
            suffix = c("1", "2"),
            relationship = "many-to-many"
          ) %>%
          mutate(
            !!column_name1 := as.numeric(.data[[column_name1]]),
            !!column_name2 := as.numeric(.data[[column_name2]])
          ) %>%
          
          # avoid self‐pairs and duplicate pairs (A–B vs B–A)
          filter(condition1 != condition2) %>%
          
          # require same sign of activation_score (product > 0)
          filter(.data[[column_name1]] * .data[[column_name2]] > 0) %>%
          
          # select / rename for plotting
          transmute(
            gs.name,
            condition1 = condition1,
            condition2 = condition2,
            position = gene_set_positions[gs.name]
          ) %>%
          distinct()
        
        any_links = (nrow(shared_links) != 0)
      } else if (input$circos_arc == "positive"){
        
        shared_modules <- results_df_significant %>%
          filter(.data[[input$circos_scores]] >= 0) %>%
          group_by(gs.name) %>%
          filter(n() > 1) %>%
          summarise(condition_names = list(unique(condition)), .groups = "drop")
        
        any_links = (nrow(shared_modules) != 0)
        
        if (any_links) {
          # Prepare the data in the format to plot
          shared_links <- shared_modules %>%
            mutate(pairs = purrr::map(condition_names, ~combn(.x, 2, simplify = FALSE))) %>%
            select(gs.name, pairs) %>%
            unnest(pairs) %>%
            unnest_wider(pairs, names_sep = "_") %>%
            rename(condition1 = pairs_1, condition2 = pairs_2) %>%
            mutate(position = gene_set_positions[gs.name])  # Add position for gene
        }
        
      } else if (input$circos_arc == "negative"){
        
        shared_modules <- results_df_significant %>%
          filter(.data[[input$circos_scores]] <= 0) %>%
          group_by(gs.name) %>%
          filter(n() > 1) %>%
          summarise(condition_names = list(unique(condition)), .groups = "drop")
        
        any_links = (nrow(shared_modules) != 0)
        
        if (any_links) {
          # Prepare the data in the format to plot
          shared_links <- shared_modules %>%
            mutate(pairs = purrr::map(condition_names, ~combn(.x, 2, simplify = FALSE))) %>%
            select(gs.name, pairs) %>%
            unnest(pairs) %>%
            unnest_wider(pairs, names_sep = "_") %>%
            rename(condition1 = pairs_1, condition2 = pairs_2) %>%
            mutate(position = gene_set_positions[gs.name])  # Add position for gene
        }
        
      }
      
      # Initialize circos plot
      par(mar = rep(0, 4))
      
      # Set layout: 1 row, 2 columns (75% for circos, 25% for legend)
      layout(matrix(1:2, ncol = 2), widths = c(0.75, 0.25))
      
      circos.clear()
      # Set circos dimensions
      # Before you start plotting
      circos.par(cell.padding = c(0, 0, 0, 0),
                 track.margin = c(0, 0.01),
                 start.degree = 81,
                 gap.degree = 2)
      circos.par("canvas.xlim" = c(-1.2, 1.2), "canvas.ylim" = c(-1.2, 1.2))
      # Initialise
      circos.initialize(
        factors = condition_circos_metadata$condition,
        xlim = cbind(condition_circos_metadata$xmin,
                     condition_circos_metadata$xmax)
      )
      
      colour_variable_name = paste0("condition_color_", input$circos_time)
      text_variable_name   = paste0("text_color_", input$circos_time)
      expr_variable_name = paste0("condition_expr_available_", input$circos_time)
      
      
      ### LAYER 1 : SEGMENTS LABELLED BY CONDITION ###
      suppressMessages({
        
        # Plot condition sectors
        circos.trackPlotRegion(
          ylim = c(0, 1),
          factors = condition_circos_metadata$condition,
          track.height = 0.1,
          bg.lwd = 1,
          panel.fun = function(x, y) {
            name <- get.cell.meta.data("sector.index")
            i <- get.cell.meta.data("sector.numeric.index")
            xlim <- get.cell.meta.data("xlim")
            ylim <- get.cell.meta.data("ylim")
            theta <- circlize(mean(xlim), 1.3)[1, 1] %% 360
            dd <- ifelse(
              theta > 200 && theta < 340,
              "outside",
              "inside"
            )
            circos.text(
              x = mean(xlim),
              y = ylim[2] + 0.4,  # Adjusted y position
              labels = name,
              facing = dd,
              cex = 1.5,
              adj = c(0.5, 0.5),
              niceFacing = T,
              col = condition_circos_metadata[[text_variable_name]][i]
            )
            
            # Draw segment with black border
            circos.rect(
              xleft = xlim[1],
              ybottom = ylim[1],
              xright = xlim[2],
              ytop = ylim[2],
              col = condition_circos_metadata[[colour_variable_name]][i],
              border = "black"
            )  # Black border
            
            circos.axis(
              labels = FALSE,
              major.tick = FALSE
            )
          }
        )
      }) # End suppressMessages
      
      # Layer 2 : Correlation with response
      
      if (input$circos_ring == "all") {
        response_variable_name = paste0("condition_response_available_", input$circos_time)
        suppressMessages({
          cor_data <- results_df_circos %>%
            filter(time == input$circos_time,!is.na(.data[[input$circos_correlation]])) %>%
            group_by(condition) %>%
            summarise(
              cor_scores = list(.data[[input$circos_correlation]]),
              gs.names = list(gs.name),
              .groups = "drop"
            ) %>%
            split(.$condition)
          
          cor_data <- cor_data[sapply(cor_data, nrow) > 0]
          
          # Draw the response layer track using the pre-calculated cor_data
          circos.trackPlotRegion(
            ylim = c(-1, 1),
            factors = condition_circos_metadata$condition,
            track.height = 0.1,
            bg.lwd = 1,
            bg.border = condition_circos_metadata[[response_variable_name]] %>% ifelse("black", "white"),
            panel.fun = function(x, y) {
              name <- get.cell.meta.data("sector.index")
              i <- get.cell.meta.data("sector.numeric.index")
              xlim <- get.cell.meta.data("xlim")
              ylim <- get.cell.meta.data("ylim")
              
              # Use pre-calculated data for the current sector
              if (name %in% names(cor_data)) {
                sec_data <- cor_data[[name]]
                # Extract correlation scores and gene set names
                correlation_scores <- sec_data$cor_scores[[1]]
                gs.names <- sec_data$gs.names[[1]]
                
                # Get positions for the gene sets relevant to this condition
                positions <-
                  gene_set_positions[as.character(gs.names)]
                
                # Determine bar colors vectorized over correlation_scores
                bar_colors <-
                  ifelse(correlation_scores > 0, "purple", "orange")
                
                # Normalize score lengths using vectorized computation
                bar_lengths <- correlation_scores
                
                # Draw all bars in one go if circos.rect supports vectorized arguments
                circos.rect(
                  xleft = positions - 0.5,
                  ybottom = 0,
                  xright = positions + 0.5,
                  ytop = bar_lengths,
                  col = bar_colors,
                  border = NA
                )
              }
            }
          )
          
        }) # End suppressMessages
      }
      
      ## LAYER 3 - ACTIVATION SCORES
      if (input$circos_ring %in% c("all","expression")){
        all_scores <- results_df_circos %>%
          filter(time == input$circos_time) %>%
          pull(!!sym(input$circos_scores)) %>%
          as.numeric()
        
        threshold <- quantile(abs(all_scores), input$circos_quantile_scoreclip)
        
        # Clip the activation scores according to the threshold
        max_clipped_score <- max(
          abs(
            pmin(
              pmax(
                all_scores,
                -threshold
              ),
              threshold
            )
          )
        )
        
        # Aggregate the activation scores per condition and gene set.
        # For example, here we use the mean of duplicate entries.
        act_data <- results_df_circos %>%
          filter(time == input$circos_time) %>%
          group_by(condition, gs.name) %>%
          summarise(
            avg_score = mean(!!sym(input$circos_scores), na.rm = TRUE),
            .groups = "drop"
          ) %>%
          group_by(condition) %>%
          summarise(
            gs.names = list(gs.name),
            raw_scores = list(avg_score),
            .groups = "drop"
          ) %>%
          mutate(
            clipped_scores = purrr::map(raw_scores, ~ pmin(pmax(.x, -threshold), threshold))
          ) %>%
          split(.$condition)
        
        act_data <- act_data[sapply(act_data, nrow) > 0]
        
        suppressMessages({
          circos.trackPlotRegion(
            ylim = c(-1, 1),
            factors = condition_circos_metadata$condition,
            track.height = 0.1,
            bg.lwd = 1,
            bg.border = condition_circos_metadata[[ paste0("condition_expr_available_", input$circos_time) ]] %>% ifelse("black", "white"),
            panel.fun = function(x, y) {
              # Get current sector (condition)
              sector_name <- get.cell.meta.data("sector.index")
              xlim <- get.cell.meta.data("xlim")
              ylim <- get.cell.meta.data("ylim")
              
              if (sector_name %in% names(act_data)) {
                sector_act <- act_data[[sector_name]]
                
                # Extract pre-calculated data vectors
                clipped_scores <- sector_act$clipped_scores[[1]]
                gs.names <- sector_act$gs.names[[1]]
                
                # Use match to get positions so that vector lengths are consistent.
                positions <- gene_set_positions[match(gs.names, names(gene_set_positions))]
                
                # Remove any missing positions (if some gene sets are not in the global order)
                keep <- !is.na(positions)
                positions <- positions[keep]
                clipped_scores <- clipped_scores[keep]
                
                # Define colors for each score: red for positive, blue for negative
                bar_colors <- ifelse(clipped_scores > 0, "red", "blue")
                
                # Normalize score lengths to scale within the track
                bar_lengths <- ifelse(
                  clipped_scores < 0,
                  clipped_scores * -(ylim[1] / max_clipped_score),
                  clipped_scores * (ylim[2] / max_clipped_score)
                )
                
                # Draw all bars in a single, vectorized call
                circos.rect(
                  xleft = positions - 0.5,
                  ybottom = 0,
                  xright = positions + 0.5,
                  ytop = bar_lengths,
                  col = bar_colors,
                  border = NA
                )
              }
            }
          )
        }) # End SuppressMessages
      }
      
      # LAYER 4 : GENE SET AGGREGATE POSITIONS
      suppressMessages({
        circos.trackPlotRegion(
          ylim = c(0, 1),
          factors = condition_circos_metadata$condition,
          track.height = 0.05,
          bg.lwd = 1,
          bg.border = condition_circos_metadata[[expr_variable_name]] %>% ifelse("black", "white"),
          panel.fun = function(x, y) {
            # Get sector and drawing parameters once per panel call
            sector_index <- get.cell.meta.data("sector.index")
            xlim <- get.cell.meta.data("xlim")
            
            # Directly extract condition availability using indexing
            condition_available <- condition_circos_metadata[condition_circos_metadata$condition == sector_index, expr_variable_name]
            
            # Vectorize the calculation of start and end positions
            x_start <- xlim[1] + (xlim[2] - xlim[1]) * aggregate_proportions$cum_start
            x_end   <- xlim[1] + (xlim[2] - xlim[1]) * aggregate_proportions$cum_end
            
            # Determine colours based on condition availability
            if (condition_available) {
              col_fill <- aggregate_proportions$colour
              border_col <- rep(NA, length(col_fill))
            } else {
              col_fill <- rep("white", length(x_start))
              border_col <- rep("white", length(x_start))
            }
            
            # Draw all segments in a single, vectorized call
            circos.rect(
              xleft = x_start,
              ybottom = 0,
              xright = x_end,
              ytop = 1,
              col = col_fill,
              border = border_col
            )
          }
        )
      }) # End suppressMessages
      
      # PLOT ARCS
      suppressMessages({
        if (any_links){
          # Pre-calculate the aggregate values for each gene set (if not already present)
          gs.aggregates <- unique(results_df_circos[, c("gs.name", "gs.aggregate")])
          shared_links <- merge(shared_links, gs.aggregates, by = "gs.name")
          
          # make sure gs.aggregate is character
          shared_links$gs.aggregate <- as.character(shared_links$gs.aggregate)
          
          # now name-based lookup “just works”
          shared_links$link_colour <- aggregate_colors[ shared_links$gs.aggregate ]
          
          
          # Use mapply to vectorize the drawing of links
          mapply(function(condition1, condition2, position, link_colour) {
            circos.link(
              sector.index1 = condition1,
              point1 = c(position, position),
              sector.index2 = condition2,
              point2 = c(position, position),
              col = adjustcolor(link_colour, alpha.f = 0.3)
            )
          },
          shared_links$condition1,
          shared_links$condition2,
          shared_links$position,
          shared_links$link_colour)
        }
      }) # End suppressMessages
      
      # PLOT LEGENDS
      if(input$circos_ring == "all"){
        legend_labels <- c(
          expression(bold("Gene Set Aggregate")), input$circos_aggregates,
          expression(bold("Direction of Regulation")), "Upregulated", "Downregulated",
          expression(bold("Correlation with Ab Response")), "Positive", "Negative"
        )
        
        legend_colours <- c(
          NA, aggregate_colors[input$circos_aggregates],  # Gene Set Aggregate
          NA, "red", "blue",     # Gene Set Activity
          NA, "purple", "orange" # Correlation with MFC
        )
        
        # Plot the combined legend with optimized spacing and manually larger titles
        legend("right",
               legend = legend_labels,
               fill = legend_colours,
               border = "white",
               cex = c(1.5, rep(1.2, length(input$circos_aggregates)),  # Gene Set Aggregate
                       1.5, 1.2, 1.2,  # Gene Set Activity
                       1.5, 1.2, 1.2), # Correlation with MFC
               text.width = 1.2 * max(strwidth(legend_labels, cex = 1.5)),  # Ensures alignment
               y.intersp = c(1, rep(0.8, length(input$circos_aggregates)),  # Gene Set Aggregate
                             1.5, 0.8, 0.8,  # Gene Set Activity
                             1.5, 0.8, 0.8), # Correlation with MFC
               inset = c(0, 0.01),
               ncol = 1, bty = "o", xjust = 0,
               text.col = rep("black", length(legend_labels)),
               title = NULL)  # Disable default title
      } else if (input$circos_ring == "expression"){
        
        legend_labels <- c(
          expression(bold("Gene Set Aggregate")), input$circos_aggregates,
          expression(bold("Direction of Regulation")), "Upregulated", "Downregulated"
        )
        
        legend_colours <- c(
          NA, aggregate_colors[input$circos_aggregates],  # Gene Set Aggregate
          NA, "red", "blue"     # Gene Set Activity
        )
        
        # Plot the combined legend with optimized spacing and manually larger titles
        legend("right",
               legend = legend_labels,
               fill = legend_colours,
               border = "white",
               cex = c(1.5, rep(1.2, length(input$circos_aggregates)),  # Gene Set Aggregate
                       1.5, 1.2, 1.2), # Correlation with MFC
               text.width = 1.2 * max(strwidth(legend_labels, cex = 1.5)),  # Ensures alignment
               y.intersp = c(1, rep(0.8, length(input$circos_aggregates)),  # Gene Set Aggregate
                             1.5, 0.8, 0.8), # Correlation with MFC
               inset = c(0, 0.01),
               ncol = 1, bty = "o", xjust = 0,
               text.col = rep("black", length(legend_labels)),
               title = NULL)  # Disable default title
        
        
      } else if (input$circos_ring == "none"){
        
        
        legend_labels <- c(
          expression(bold("Gene Set Aggregate")), input$circos_aggregates
        )
        
        legend_colours <- c(
          NA, aggregate_colors[input$circos_aggregates]
        )
        
        # Plot the combined legend with optimized spacing and manually larger titles
        legend("right",
               legend = legend_labels,
               fill = legend_colours,
               border = "white",
               cex = c(1.5, rep(1.2, length(input$circos_aggregates))), # Correlation with MFC
               text.width = 1.2 * max(strwidth(legend_labels, cex = 1.5)),  # Ensures alignment
               y.intersp = c(1, rep(0.8, length(input$circos_aggregates))), # Correlation with MFC
               inset = c(0, 0.01),
               ncol = 1, bty = "o", xjust = 0,
               text.col = rep("black", length(legend_labels)),
               title = NULL)  # Disable default title
        
      }
      
      # **Record** the finished plot
      lastCircos(recordPlot())
    }) # End isolate
    
  }) # End render plot
  
  # Download the circos plot if clicked
  output$download_circos_plot <- downloadHandler(
    filename = function() {
      ext <- req(input$circos_download_format)
      paste0("circos_plot_day", input$circos_time, "_", Sys.Date(), ".", ext)
    },
    content = function(file) {
      # make sure we have a plot recorded
      req(lastCircos())
      
      fmt <- input$circos_download_format
      w   <- input$circos_download_width
      h   <- input$circos_download_height
      d   <- input$circos_download_dpi
      
      # open the right device
      switch(fmt,
             pdf  = pdf(file, width = w, height = h),
             png  = png(file, width = w, height = h, units = "in", res = d),
             jpeg = jpeg(file, width = w, height = h, units = "in", res = d),
             stop("Unknown format: ", fmt)
      )
      
      # replay and close
      replayPlot(lastCircos())
      dev.off()
    }
  )
  
  
  ### END CIRCOS PLOT TAB ###
  
  ### BEGIN HEATMAP TAB ###
  heatmap_plot <- reactive({
    
    input$heatmap_Update
    isolate({
      
      # Copy data for plotting
      results_df_heatmap <- results_df %>%
        filter(method == input$heatmap_method)
      
      category_order = levels(results_df_heatmap$gs.aggregate)
      
      full_order <- levels(results_df_heatmap$condition)
      
      
      ## STEP 1: Determine significance ---------------------------------------------
      if (input$heatmap_filter_variable == "none") {
        # Only p-value based filtering
        results_df_heatmap <- results_df_heatmap %>%
          mutate(
            significant = !!sym(paste0(input$heatmap_p_approach, ".adjPval_", input$heatmap_p_correction)) < input$heatmap_p_threshold
          )
        
      } else {
        # Effect-size plus p-value filtering
        if (input$heatmap_p_approach == "global") {
          # Global threshold across all data
          effect_size_threshold <- results_df_heatmap %>%
            pull(!!sym(input$heatmap_scores)) %>% abs() %>%
            quantile(input$heatmap_quantile_threshold)
          
          results_df_heatmap <- results_df_heatmap %>%
            mutate(
              significant = (!!sym(paste0(input$heatmap_p_approach, ".adjPval_", input$heatmap_p_correction)) < input$heatmap_p_threshold ) &
                (abs(.data[[input$heatmap_filter_variable]]) > effect_size_threshold)
            )
          
        } else if (input$heatmap_p_approach == "withinTime") {
          # Threshold per time
          thresholds <- results_df_heatmap %>%
            group_by(time) %>%
            summarise(
              threshold_abs = quantile(abs(.data[[input$heatmap_scores]]),
                                       input$heatmap_quantile_threshold, na.rm = TRUE),
              .groups = "drop"
            )
          
          results_df_heatmap <- results_df_heatmap %>%
            left_join(thresholds, by = "time") %>%
            mutate(
              significant = (!!sym(paste0(input$heatmap_p_approach, ".adjPval_", input$heatmap_p_correction)) < input$heatmap_p_threshold ) &
                (abs(.data[[input$heatmap_filter_variable]]) > threshold_abs)
            ) %>%
            select(-threshold_abs)
          
        } else if (input$heatmap_p_approach == "withinComparison") {
          # Threshold per condition-time
          results_df_heatmap <- results_df_heatmap %>%
            group_by(condition, time) %>%
            mutate(
              local_effect = abs(.data[[input$heatmap_scores]]) >
                quantile(abs(.data[[input$heatmap_scores]]),
                         input$heatmap_quantile_threshold, na.rm = TRUE)
            ) %>%
            ungroup() %>%
            mutate(
              significant = (!!sym(paste0(input$heatmap_p_approach, ".adjPval_", input$heatmap_p_correction)) < input$heatmap_p_threshold ) &
                local_effect
            ) %>%
            select(-local_effect)
        }
      }
      
      # Convert to factor for consistent plotting
      results_df_heatmap <- results_df_heatmap %>%
        mutate(significant = factor(significant, levels = c(TRUE, FALSE)))
      
      ## STEP 2: Filter and label ---------------------------------------------------
      results_df_heatmap <- results_df_heatmap %>%
        filter(
          gs.aggregate %in% input$heatmap_aggregates,
          condition %in% input$heatmap_conditions,
          time %in% input$heatmap_times
        ) %>%
        mutate(
          time_label = factor(
            paste0(condition, " - Day ", time),
            levels = unique(paste0(condition, " - Day ", time))
          )
        )
      
      # Clip scores if required (avoids extreme values distorting the colouring)
      clip_threshold  <- results_df_heatmap %>%
        pull(.data[[input$heatmap_scores]]) %>%
        abs() %>%
        quantile(input$heatmap_quantile_scoreclip, na.rm = TRUE)
      
      results_df_heatmap <- results_df_heatmap %>%
        mutate(
          !!input$heatmap_scores := case_when(
            .data[[input$heatmap_scores]] >  clip_threshold  ~  clip_threshold ,
            .data[[input$heatmap_scores]] < -clip_threshold  ~ -clip_threshold ,
            TRUE                                     ~  .data[[input$heatmap_scores]]
          )
        )
      
      ## STEP 3: Filter by commonly DE genesets ------------------------------------
      if (input$heatmap_filter_commonDE == "global"){
        results_df_heatmap <- results_df_heatmap %>%
          group_by(gs.name) %>%
          # compute fraction of non‐significant comparisons per gene set
          mutate(frac_sig = sum(as.logical(significant)) / n()) %>%
          ungroup() %>%
          # keep only those gs.name for which non‐significance is ≤ input$heatmap_common_percentage
          filter(frac_sig >= input$heatmap_common_percentage) %>%
          select(-frac_sig)
      } else if (input$heatmap_filter_commonDE == "withinTime") {
        
        # 1) compute fraction *significant* per gs.name & time
        frac_tbl <- results_df_heatmap %>%
          group_by(gs.name, time) %>%
          summarise(
            frac_sig = sum(significant == "TRUE") / n(),
            .groups = "drop"
          )
        
        # 2) pick up any gene‐set where frac_sig >= threshold in at least one time
        keep_sets <- frac_tbl %>%
          filter(frac_sig >= input$heatmap_common_percentage) %>%
          pull(gs.name) %>% unique()
        
        # 3) filter your full data
        results_df_heatmap <- results_df_heatmap %>%
          filter(gs.name %in% keep_sets)
        
      } else if(input$heatmap_filter_commonDE == "score"){
        # Calculate "sharing scores"
        # Consider up and down-regulation separately
        # Count the number of vaccines for which each gene set is either significantly up- or down-regulated over all times
        # Consider all times, regardless of user input. However, only consider the
        # user inputted vaccines.
        sharing_score_df = results_df %>%
          filter(condition %in% input$heatmap_conditions,
                 method == input$heatmap_method)
        
        calc_sharing_df_multi <- function(x, alpha){
          adjp_cols = paste0(input$heatmap_p_approach, ".adjPval_",input$heatmap_p_correction)
          map_dfr(adjp_cols, function(cur_col) {
            x %>%
              mutate(
                direction      = ifelse(.data[[input$heatmap_scores]] > 0, "up", "down"),
                is_significant = !is.na(.data[[cur_col]]) & .data[[cur_col]] < alpha
              ) %>%
              select(condition, gs.name, direction, is_significant) %>%
              distinct() %>%
              group_by(gs.name, direction) %>%
              summarise(n = sum(is_significant), .groups = "drop") %>%
              group_by(gs.name) %>%
              arrange(-n) %>%
              slice(1) %>%
              ungroup() %>%
              mutate(adjPval_method = cur_col)
          }) %>%
            pivot_wider(
              id_cols    = gs.name,
              names_from  = adjPval_method,
              values_from = n,
              values_fill = 0
            )
        }
        
        sharing_score_df <-
          sharing_score_df %>%
          calc_sharing_df_multi(alpha = input$heatmap_p_threshold)
        
        keep_sets = sharing_score_df %>%
          filter(.data[[paste0(input$heatmap_p_approach, ".adjPval_",input$heatmap_p_correction)]] >= input$heatmap_score_threshold) %>%
          pull(gs.name)
        
        results_df_heatmap <- results_df_heatmap %>%
          filter(gs.name %in% keep_sets)
      }
      
      
      if(nrow(results_df_heatmap) == 0){
        stop("No genesets to plot under current parameters!")
      }
      
      ## STEP 3: Determine axis breaks for separation ------------------------------
      combos <- results_df_heatmap %>%
        select(condition, time) %>%
        distinct() %>%
        arrange(time)
      
      # Vertical lines between input$heatmap_times
      vlines <- which(diff(as.numeric(combos$time)) != 0) + 0.5
      
      # Clustering within-time
      if (input$heatmap_x_order == "cluster-time"){
        # 1. Pivot into a matrix of input$heatmap_conditions × gene‐sets for each time
        mat_list <- results_df_heatmap %>%
          mutate(time = as.numeric(time)) %>%
          select(time_label, time, gs.name.description, !!sym(input$heatmap_scores)) %>%
          pivot_wider(
            names_from  = gs.name.description,
            values_from = !!sym(input$heatmap_scores),
            values_fill = 0
          ) %>%
          group_by(time) %>%
          group_split()
        
        # 2. For each time, cluster input$heatmap_conditions and extract the ordered labels
        condition_order_by_time <- map_dfr(mat_list, function(df) {
          # Extract the matrix and set rownames
          m <- df %>% select(-time_label, -time) %>% as.matrix()
          rownames(m) <- df$time_label
          
          # Decide whether to cluster or to fall back on full_order
          if (nrow(df) > 2) {
            # Hierarchical clustering
            hc <- hclust(dist(m, method = "euclidean"), method = "complete")
            order_labels <- rownames(m)[rev(hc$order)]
          } else {
            # Default to the global full_order of input$heatmap_conditions
            order_labels <- df %>%
              mutate(
                # Extract just the condition name (before " - Day")
                vacc = sub(" - Day.*", "", time_label)
              ) %>%
              arrange(match(vacc, full_order)) %>%
              pull(time_label)
          }
          
          # Return a tibble of time_label order for this time
          tibble(
            time     = unique(df$time),
            condition_order = order_labels
          )
        })
        
        # 3. Combine clustered orders across input$heatmap_times into a single vector
        clustered_levels <- condition_order_by_time %>%
          arrange(time) %>%
          pull(condition_order) %>%
          unlist()
        
        # 4. Re‐factor `time_label` in the plotting data
        results_df_heatmap <- results_df_heatmap %>%
          mutate(
            time_label = factor(
              time_label,
              levels = clustered_levels
            )
          )
        
        
        
      }
      
      ## STEP 4: Optional row ordering ---------------------------------------------
      if (input$heatmap_y_order == "cluster") {
        # Hierarchical clustering of gene sets
        mat <- results_df_heatmap %>%
          select(gs.name.description, time_label, input$heatmap_scores) %>%
          pivot_wider(names_from = time_label,
                      values_from = input$heatmap_scores,
                      values_fill = 0) %>%
          column_to_rownames("gs.name.description") %>%
          as.matrix()
        
        dend  <- hclust(dist(mat), method = "complete")
        order <- rev(rownames(mat)[dend$order])
        
        results_df_heatmap <- results_df_heatmap %>%
          mutate(gs.name.description = factor(gs.name.description, levels = order))
        
      } else if (input$heatmap_y_order == "aggregate") {
        # Order by aggregate category
        results_df_heatmap$gs.aggregate <- factor(
          results_df_heatmap$gs.aggregate,
          levels = category_order
        )
        
        o <- with(results_df_heatmap, order(gs.aggregate, gs.name))
        gene_order <- unique(results_df_heatmap$gs.name.description[o])
        
        results_df_heatmap <- results_df_heatmap %>%
          mutate(
            gs.name.description = factor(gs.name.description, levels = rev(gene_order))
          )
        
        # Horizontal lines between aggregates
        group_tbl   <- results_df_heatmap %>%
          select(gs.name.description, gs.aggregate) %>% distinct()
        sizes       <- rev(table(group_tbl$gs.aggregate))
        hlines      <- cumsum(sizes)[-length(sizes)] + 0.5
      }
      
      ## STEP 5: Build and render heatmap -----------------------------------------
      # make caption for footnote showing parameter values
      if (input$heatmap_filter_commonDE != "none"){
        filtration_string = if (input$heatmap_filter_commonDE == "score"){
          paste0("sharing score above ", input$heatmap_score_threshold) 
        } else if (input$heatmap_filter_commonDE == "global"){
          paste0("across selected times common DE above ", 100*input$heatmap_common_percentage, "%")
        } else if (input$heatmap_filter_commonDE == "withinTime"){
          paste0("within selected times common DE above ", 100*input$heatmap_common_percentage, "%")
        }
        
        if (input$heatmap_filter_variable != "none"){
          
          
          caption = paste0("Parameters:\n",
                           "DGSA method : ", input$heatmap_method, "\n",
                           "- p value correction: ", input$heatmap_p_correction, "     ",
                           "- p value subset: ",              input$heatmap_p_approach,     "     ",
                           "- significance level: ",           input$heatmap_p_threshold,     "\n",
                           "- filtration on ", input$heatmap_filter_variable,      "     ",
                           " at quantile threshold ", input$heatmap_quantile_threshold, "\n",
                           "genesets displayed based on : ", filtration_string) 
        } else if (input$heatmap_filter_variable == "none"){
          
          caption = paste0("Parameters:\n",
                           "DGSA method : ", input$heatmap_method, "\n",
                           "- p value correction: ", input$heatmap_p_correction, "     ",
                           "- p value subset: ",              input$heatmap_p_approach,     "     ",
                           "- significance level: ",           input$heatmap_p_threshold,     "\n",
                           "genesets displayed based on : ", filtration_string)
        }
        
      } else if (input$heatmap_filter_commonDE == "none"){
        
        if (input$heatmap_filter_variable != "none"){
          
          
          caption = paste0("Parameters:\n",
                           "DGSA method : ", input$heatmap_method, "\n",
                           "- p value correction: ", input$heatmap_p_correction, "     ",
                           "- p value subset: ",              input$heatmap_p_approach,     "     ",
                           "- significance level: ",           input$heatmap_p_threshold,     "\n",
                           "- filtration on ", input$heatmap_filter_variable,      "     ",
                           " at quantile threshold ", input$heatmap_quantile_threshold, "\n") 
        } else if (input$heatmap_filter_variable == "none"){
          
          caption = paste0("Parameters:\n",
                           "DGSA method : ", input$heatmap_method,"\n", 
                           "- p value correction: ", input$heatmap_p_correction, "     ",
                           "- p value subset: ",              input$heatmap_p_approach,     "     ",
                           "- significance level: ",           input$heatmap_p_threshold,     "\n")
        }
      }
      
      
      
      
      
      plot <- results_df_heatmap %>%
        ggplot(aes(x = time_label, y = gs.name.description, fill = .data[[input$heatmap_scores]])) +
        geom_tile() +
        scale_fill_gradient2(
          low = "blue", mid = "white", high = "red",
          name = ifelse(input$heatmap_scores == "activation.score",
                        "Activation score",
                        "Mean of Fold Change")
        ) +
        geom_point(
          aes(shape = significant, color = significant), size = 1
        ) +
        scale_shape_manual(values = c("TRUE" = 8, "FALSE" = 0)) +
        scale_color_manual(values = c("TRUE" = "black", "FALSE" = "transparent")) +
        labs(x = "Condition-Time", y = "",
             caption = caption) +
        theme_minimal(base_size = 18) +
        theme(
          axis.text.x    = element_text(angle = 35, hjust = 1, face = "bold", size = 13),
          axis.text.y    = element_text(size = 14),
          axis.title.x   = element_text(size = 20, face = "bold"),
          legend.title   = element_text(size = 18),
          legend.text    = element_text(size = 15),
          panel.grid     = element_blank(),
          axis.ticks     = element_blank(),
          plot.caption          = element_text(
            size       = 14,
            hjust      = 0,
            lineheight = 1.2,
            margin     = margin(t = 15, r = 0, b = 0, l = 0)
          ),
          plot.caption.position = "plot"
        )
      
      # Add horizontal separators if using aggregate ordering
      if (input$heatmap_y_order == "aggregate") {
        # 1) named vector: aggregate → hex
        agg_cols <- results_df_heatmap %>%
          distinct(gs.aggregate, gs.colour) %>%
          arrange(gs.aggregate) %>%
          { setNames(.$gs.colour, .$gs.aggregate) }
        
        # 2) sidebar data
        sidebar_df <- results_df_heatmap %>%
          distinct(gs.name.description, gs.aggregate) %>%
          mutate(
            x = factor("__agg__"),
            gs.name.description = factor(
              gs.name.description,
              levels = levels(results_df_heatmap$gs.name.description)
            )
          )
        
        # 3) number of real heatmap columns
        n_time <- length(levels(results_df_heatmap$time_label))
        
        plot <- plot +
          # 4) draw separators starting at right edge of sidebar (x = 1 + width/2)
          geom_segment(
            data = data.frame(hline = hlines),
            aes(
              x    = 1.25,
              xend = 1.5 + n_time,
              y    = hline,
              yend = hline
            ),
            inherit.aes = FALSE,
            colour    = "black",
            linewidth = 1,
            alpha     = 0.6
          ) +
          # 5) unlock second fill scale for the sidebar
          ggnewscale::new_scale_fill() +
          geom_tile(
            data        = sidebar_df,
            aes(x = x, y = gs.name.description, fill = gs.aggregate),
            width       = 0.25,
            inherit.aes = FALSE
          ) +
          # 6) legend keyed to aggregate names
          scale_fill_manual(
            name   = "Aggregate",
            values = agg_cols,
            guide  = guide_legend(order = 2)
          ) +
          # 7) single discrete x-scale: fake level first, minimal left padding
          scale_x_discrete(
            limits = c("__agg__", levels(results_df_heatmap$time_label)),
            expand = expansion(add = c(0.1, 0)),   # small gap only
            labels = function(x) ifelse(x == "__agg__", "", x)
          ) +
          labs(x = NULL) +
          theme(
            legend.box     = "vertical",
            legend.spacing = unit(0.2, "cm")
          ) +
          geom_vline(xintercept = vlines + 1, linetype = "dashed", linewidth = 1, colour = "black", alpha = 0.5)
      } else {
        plot = plot +
          geom_vline(xintercept = vlines, linetype = "dashed", linewidth = 1, colour = "black", alpha = 0.5)
      }
      print(plot)
      
      
    }) # End Isolate
    
  }) # End Reactive
  
  
  output$heatmapPlot = renderPlot({
    heatmap_plot()
  })
  
  output$download_heatmap <- downloadHandler(
    filename = function() {
      req(input$heatmap_download_format)
      paste0("heatmap", Sys.Date(), ".",
             input$heatmap_download_format)
    },
    content = function(file) {
      # make sure the plot reactive is available
      req(heatmap_plot())
      
      # call ggsave with user inputs
      ggsave(
        filename = file,
        plot     = heatmap_plot(),
        device   = input$heatmap_download_format,
        width    = input$heatmap_download_width,
        height   = input$heatmap_download_height,
        dpi      = input$heatmap_download_dpi,
        units    = "in"
      )
    }
  )
  
  ### END HEATMAP TAB ###
  
  ### BEGIN EXPRESSION PROFILES TAB ###
  
  exprprofiles_plot <- reactive({
    
    input$exprprofiles_Update
    isolate({
      
      # --- 1. Filter & clip data ---
      df_time <- results_df %>%
        filter(
          time %in% input$exprprofiles_times,
          gs.aggregate %in% input$exprprofiles_aggregates,
          condition %in% input$exprprofiles_conditions,
          method == input$exprprofiles_method
        ) %>%
        mutate(
          time = as.numeric(as.character(time)),
          gs.aggregate     = droplevels(gs.aggregate)
        ) %>%
        mutate(
          threshold        = quantile(abs(.data[[input$exprprofiles_scores]]),
                                      input$exprprofiles_quantile_scoreclip,
                                      na.rm = TRUE),
          score = ifelse(
            abs(.data[[input$exprprofiles_scores]]) > threshold,
            sign(.data[[input$exprprofiles_scores]]) * threshold,
            .data[[input$exprprofiles_scores]]
          )
        )
      
      relevant_times = unique(df_time$time)
      relevant_conditions = unique(df_time$condition)
      
      # --- 2. Global gene‐set ordering ---
      global_order <- df_time %>%
        distinct(gs.name, gs.aggregate) %>%
        arrange(gs.aggregate) %>%
        pull(gs.name)
      df_time <- df_time %>%
        mutate(gs.name = factor(gs.name, levels = global_order))
      
      # --- 3. condition ordering (full -> selected) ---
      full_order = levels(df_time$condition)
      plot_order <- full_order[full_order %in% relevant_conditions]
      
      df_time <- df_time %>%
        mutate(condition = factor(condition, levels = plot_order))
      
      # --- 4. Placeholders for missing combos ---
      combo <- expand.grid(
        condition   = plot_order,
        time = relevant_times,
        stringsAsFactors = FALSE
      )
      existing <- df_time %>% distinct(condition, time = time)
      to_add   <- anti_join(combo, existing, by = c("condition","time")) %>%
        mutate(time = time)
      if (nrow(to_add) > 0) {
        ph <- to_add %>%
          mutate(
            gs.name          = global_order[1],
            score = NA,
            gs.aggregate     = NA
          ) %>%
          mutate(
            gs.name  = factor(gs.name, levels = global_order),
            condition  = factor(condition, levels = plot_order)
          )
        df_time <- bind_rows(df_time, ph)
      }
      
      # --- 5. Day‐labels ---
      df_time <- df_time %>%
        mutate(
          day_label = factor(
            paste0("Day ", time),
            levels = paste0("Day ", sort(as.numeric(relevant_times), decreasing = F))
          )
        )
      
      # --- 6. Aggregate stripe data ---
      aggregate_colors <- results_df %>%
        select(gs.aggregate, gs.colour) %>%
        distinct() %>%
        { setNames(.$gs.colour, .$gs.aggregate) }
      
      
      agg_props <- df_time %>%
        distinct(gs.name, gs.aggregate) %>%
        filter(!is.na(gs.aggregate)) %>%     # ⇐ drop the NAs now
        count(gs.aggregate) %>%
        mutate(
          prop      = n / sum(n),            # now sum(n) is only non‐NA counts
          colour    = aggregate_colors[gs.aggregate]
        ) %>%
        arrange(gs.aggregate) %>%
        mutate(
          cum_start = c(0, head(cumsum(prop), -1)),
          cum_end   = cumsum(prop),
          xmin      = cum_start * length(global_order) + 0.5,
          xmax      = cum_end   * length(global_order) + 0.5,
          day_label = factor(
            paste0("Day ", min(relevant_times)),
            levels = paste0("Day ", sort(relevant_times))
          )
        )
      
      
      # --- 7. condition strip colors ---
      condition.colours <- results_df %>%
        distinct(condition, condition.colour) %>%
        filter(condition %in% plot_order) %>%
        deframe()
      condition.colours <- setNames(condition.colours[plot_order], plot_order)
      
      # --- 8. Plot ---
      y_rng       <- range(df_time$score, na.rm = TRUE)
      band_height <- diff(y_rng) * 0.05
      rect_y      <- y_rng[1] - band_height
      
      if (input$exprprofiles_barcolours == "direction"){
        
        p <- ggplot(df_time, aes(gs.name, score)) +
          geom_col(
            data = filter(df_time, !is.na(score)),
            aes(fill = ifelse(score >= 0, "Up-regulated", "Down-regulated")),
            width = 0.8
          ) +
          geom_hline(
            data = filter(df_time, !is.na(score)),
            aes(yintercept = 0), color = "black"
          ) +
          scale_fill_manual(
            values = c("Up-regulated" = "red", "Down-regulated" = "blue"),
            name   = "Direction of regulation"
          ) +
          ggh4x::facet_grid2(
            rows  = vars(fct_rev(day_label)),
            cols  = vars(condition),
            drop  = FALSE,
            strip = strip_themed(
              background_x = elem_list_rect(fill = condition.colours),
              text_x       = elem_list_text(color = "black")
            )
          ) +
          scale_x_discrete(labels = NULL) +
          labs(x = "Gene Sets", y = ifelse(input$exprprofiles_scores == "activation.score",
                                           "Activation Score", "Fold Change")) +
          theme_minimal(base_size = 21) +
          theme(
            strip.text.y.left = element_text(face = "bold", size = 16),
            strip.text.y.right = element_text(face = "bold", size = 25, angle = 0),
            strip.text.x      = element_text(face = "bold", size = 8),
            axis.ticks.x      = element_blank(),
            legend.position   = "right",
            panel.spacing     = unit(0.8, "lines")
          ) +
          new_scale("fill") +
          geom_rect(
            data = agg_props,
            aes(xmin = xmin, xmax = xmax,
                ymin = rect_y, ymax = rect_y + band_height,
                fill = gs.aggregate),
            inherit.aes = FALSE
          ) +
          scale_fill_manual(
            name   = "Gene Set\nAggregate",
            values = aggregate_colors,
            guide  = guide_legend(ncol = 1)
          )
      } else if (input$exprprofiles_barcolours == "Aggregates"){
        p <- ggplot(df_time, aes(gs.name, score, fill = gs.aggregate)) +
          geom_col(width = 0.8) +
          geom_hline(yintercept = 0, color = "black") +
          scale_fill_manual(
            name          = "Gene Set\nAggregate",
            values        = aggregate_colors,
            na.translate  = FALSE            # ← drop the NA entry from the legend
          ) +
          ggh4x::facet_grid2(
            rows  = vars(fct_rev(day_label)),
            cols  = vars(condition),
            drop  = FALSE,
            strip = strip_themed(
              background_x = elem_list_rect(fill = condition.colours),
              text_x       = elem_list_text(color = "black")
            )
          ) +
          scale_x_discrete(labels = NULL) +
          labs(x = "Gene Sets", y = ifelse(input$exprprofiles_scores == "activation.score",
                                           "Activation Score", "Fold Change")) +
          theme_minimal(base_size = 14) +
          theme(
            strip.text.y.left  = element_text(face = "bold", size = 16),
            strip.text.y.right = element_text(face = "bold", size = 25, angle = 0),
            strip.text.x       = element_text(face = "bold", size = 10),
            axis.ticks.x       = element_blank(),
            legend.position    = "right",
            panel.spacing      = unit(1, "lines")
          ) +
          geom_rect(
            data        = agg_props,
            aes(xmin = xmin, xmax = xmax,
                ymin = rect_y, ymax = rect_y + band_height,
                fill = gs.aggregate),
            inherit.aes = FALSE,
            show.legend = FALSE
          )
      }
      
      print(p)
      
    }) # End Isolate
    
  }) # End reactive
  
  
  output$exprprofilesPlot = renderPlot({
    exprprofiles_plot()
  })
  
  output$download_exprprofiles <- downloadHandler(
    filename = function() {
      req(input$exprprofiles_download_format)
      paste0("exprprofiles", Sys.Date(), ".",
             input$exprprofiles_download_format)
    },
    content = function(file) {
      # make sure the plot reactive is available
      req(exprprofiles_plot())
      
      # call ggsave with user inputs
      ggsave(
        filename = file,
        plot     = exprprofiles_plot(),
        device   = input$exprprofiles_download_format,
        width    = input$exprprofiles_download_width,
        height   = input$exprprofiles_download_height,
        dpi      = input$exprprofiles_download_dpi,
        units    = "in"
      )
    }
  )
  
  ### END EXPRESSION PROFILES TAB ###
  
  ### BEGIN SPIDER PLOTS TAB ###
  lastSpider <- reactiveVal(NULL)
  
  output$spiderPlot <- renderPlot({
    input$spider_Update
    isolate({
      
      # Get original levels
      old_levels <- levels(results_df$gs.aggregate)
      
      # Truncate and add "..." only if needed
      new_levels <- ifelse(nchar(old_levels) > 12,
                           paste0(substr(old_levels, 1, 7), "..."),
                           old_levels)
      
      # Create a named vector for mapping old levels to new truncated labels
      names(new_levels) <- old_levels
      
      # Apply transformation and preserve factor level order
      results_df_spider <- results_df %>%
        filter(method == input$spider_method) %>%
        mutate(
          gs.aggregate.short = as.character(gs.aggregate),
          gs.aggregate.short = new_levels[gs.aggregate.short],
          gs.aggregate.short = factor(gs.aggregate.short, levels = unique(new_levels))
        )
      
      # Initialise a list to store spider plots
      radar_plots <- list()
      
      full_order <- levels(results_df_spider$condition)
      
      aggregate_order = levels(results_df_spider$gs.aggregate.short)
      
      all_times = levels(results_df_spider$time)
      
      # If adaptive grid positions, find the conditions which have available data at the given time
      if (input$spider_grid == "available" && input$spider_grouping == "withinTime"){
        spider_conditions_withinTime = results_df_spider %>%
          filter(time == input$spider_times_withinTime,
                 condition %in% input$spider_conditions_withinTime) %>%
          pull(condition) %>%
          unique()
        
        spider_times_withinTime = input$spider_times_withinTime
      } else if (input$spider_grid == "available" && input$spider_grouping == "withinCondition"){
        
        spider_times_withinCondition = results_df_spider %>%
          filter(time %in% input$spider_times_withinCondition,
                 condition == input$spider_conditions_withinCondition) %>%
          pull(time) %>%
          unique()
        
        spider_conditions_withinCondition = input$spider_conditions_withinCondition
        
      } else if (input$spider_grid == "set" && input$spider_grouping == "withinTime"){
        spider_conditions_withinTime = full_order
        spider_times_withinTime = input$spider_times_withinTime
      } else if (input$spider_grid == "set" && input$spider_grouping == "withinCondition"){
        spider_conditions_withinCondition = input$spider_conditions_withinCondition
        spider_times_withinCondition =  all_times
      }
      
      if (input$spider_grouping == "withinTime"){
        
        results_df_spider_time <- results_df_spider %>%
          filter(time == spider_times_withinTime,
                 gs.aggregate %in% input$spider_aggregates,
                 condition %in% spider_conditions_withinTime)
        
        aggregate_scores_all <- results_df_spider_time %>%
          group_by(gs.aggregate.short, condition) %>%
          summarize(avg_score = mean(.data[[input$spider_scores]], na.rm = TRUE), .groups = "drop")
        
        min_score <- min(aggregate_scores_all$avg_score, na.rm = TRUE)
        max_score <- max(aggregate_scores_all$avg_score, na.rm = TRUE)
        
        for (condition_name_temp in spider_conditions_withinTime) {
          
          condition_df <- results_df_spider %>%
            filter(condition == condition_name_temp,
                   time == spider_times_withinTime,
                   gs.aggregate %in% input$spider_aggregates)
          
          if (nrow(condition_df) > 0) {
            aggregate_scores <- condition_df %>%
              group_by(gs.aggregate.short) %>%
              summarize(avg_score = mean(.data[[input$spider_scores]], na.rm = TRUE), .groups = "drop")
            
            radar_data <- aggregate_scores %>%
              pivot_wider(names_from = gs.aggregate.short, values_from = avg_score) %>%
              mutate(group = "Score") %>%
              select(group, everything())
            
            ordered_columns <- c("group", intersect(aggregate_order, colnames(radar_data)))
            radar_data <- radar_data[, ordered_columns]
            
            grid_min <- floor(min_score * 10) / 10
            grid_max <- ceiling(max_score * 10) / 10
            
            colour <- results_df_spider %>%
              filter(condition == condition_name_temp) %>%
              pull(condition.colour) %>%
              unique()
            
            radar_plot <- ggradar(radar_data,
                                  axis.labels = ordered_columns[-1],
                                  grid.min = grid_min,
                                  grid.mid = 0,
                                  grid.max = grid_max,
                                  values.radar = c(grid_min, 0, grid_max),
                                  group.line.width = 0,
                                  group.point.size = 0,
                                  draw.points = FALSE,
                                  fill = TRUE,
                                  group.colours = colour,
                                  axis.label.size = 3,
                                  fill.alpha = 0.8,
                                  grid.label.size = 4
            ) +
              ggtitle(condition_name_temp) +
              theme(plot.title = element_text(hjust = 0.5, size = 15),
                    plot.margin = unit(c(1, 1, 1, 1), "lines"))
            
            radar_plots[[condition_name_temp]] <- radar_plot
          } else {
            # Simulate similar spacing for empty plot
            blank_plot <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
              geom_blank() +
              ggtitle(condition_name_temp) +
              coord_fixed() +
              theme_void() +
              theme(
                plot.title = element_text(hjust = 0.5, size = 15),
                plot.margin = unit(c(1, 1, 1, 1), "lines")
              )
            
            radar_plots[[condition_name_temp]] <- blank_plot
          }
        }
        
        radar_plots_ordered <- radar_plots[intersect(full_order, spider_conditions_withinTime)]
        
        if (length(spider_conditions_withinTime) > 8){
          grid_spider <- grid.arrange(
            grobs = radar_plots_ordered,
            nrow = 3,
            top = textGrob(paste0("Day ", spider_times_withinTime), gp = gpar(fontsize = 30, fontface = "bold"))
          ) } else if(length(spider_conditions_withinTime) < 8 && length(spider_conditions_withinTime) > 4){
            grid_spider <- grid.arrange(
              grobs = radar_plots_ordered,
              nrow = 2,
              top = textGrob(paste0("Day ", spider_times_withinTime), gp = gpar(fontsize = 30, fontface = "bold")))
          } else {
            grid_spider <- grid.arrange(
              grobs = radar_plots_ordered,
              nrow = 1,
              top = textGrob(paste0("Day ", spider_times_withinTime), gp = gpar(fontsize = 30, fontface = "bold")))
          }
      } else if (input$spider_grouping == "withinCondition"){
        
        sorted_spider_times <- spider_times_withinCondition %>% as.character() %>%
          as.numeric() %>%
          sort()
        
        results_df_spider_time <- results_df_spider %>%
          filter(condition == spider_conditions_withinCondition,
                 gs.aggregate %in% input$spider_aggregates,
                 time %in% sorted_spider_times)
        
        aggregate_scores_all <- results_df_spider_time %>%
          group_by(gs.aggregate.short, time) %>%
          summarize(avg_score = mean(.data[[input$spider_scores]], na.rm = TRUE),
                    .groups = "drop")
        
        min_score <- min(aggregate_scores_all$avg_score, na.rm = TRUE)
        max_score <- max(aggregate_scores_all$avg_score, na.rm = TRUE)
        
        for (time_temp in sorted_spider_times) {
          
          condition_df <- results_df_spider %>%
            filter(time == time_temp,
                   condition == spider_conditions_withinCondition,
                   gs.aggregate %in% input$spider_aggregates)
          
          if (nrow(condition_df) > 0) {
            aggregate_scores <- condition_df %>%
              group_by(gs.aggregate.short) %>%
              summarize(avg_score = mean(.data[[input$spider_scores]], na.rm = TRUE), .groups = "drop")
            
            radar_data <- aggregate_scores %>%
              pivot_wider(names_from = gs.aggregate.short, values_from = avg_score) %>%
              mutate(group = "Score") %>%
              select(group, everything())
            
            ordered_columns <- c("group", intersect(aggregate_order, colnames(radar_data)))
            radar_data <- radar_data[, ordered_columns]
            
            grid_min <- floor(min_score * 10) / 10
            grid_max <- ceiling(max_score * 10) / 10
            
            colour <- results_df_spider %>%
              filter(condition == spider_conditions_withinCondition) %>%
              pull(condition.colour) %>%
              unique()
            
            radar_plot <- ggradar(radar_data,
                                  axis.labels = ordered_columns[-1],
                                  grid.min = grid_min,
                                  grid.mid = 0,
                                  grid.max = grid_max,
                                  values.radar = c(grid_min, 0, grid_max),
                                  group.line.width = 0,
                                  group.point.size = 0,
                                  draw.points = FALSE,
                                  fill = TRUE,
                                  group.colours = colour,
                                  axis.label.size = 3,
                                  fill.alpha = 0.8,
                                  grid.label.size = 4
            ) +
              ggtitle(paste0("Day ", time_temp)) +
              theme(plot.title = element_text(hjust = 0.5, size = 15),
                    plot.margin = unit(c(1, 1, 1, 1), "lines"))
            
            radar_plots[[paste0("Day ", time_temp)]] <- radar_plot
          } else {
            # Simulate similar spacing for empty plot
            blank_plot <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
              geom_blank() +
              ggtitle(paste0("Day ", time_temp)) +
              coord_fixed() +
              theme_void() +
              theme(
                plot.title = element_text(hjust = 0.5, size = 15),
                plot.margin = unit(c(1, 1, 1, 1), "lines")
              )
            
            radar_plots[[paste0("Day ", time_temp)]] <- blank_plot
          }
          
        }
        
        radar_plots_ordered <- radar_plots[paste0("Day ", sorted_spider_times)]
        if (length(sorted_spider_times) > 12){
          grid_spider <- grid.arrange(
            grobs = radar_plots_ordered,
            nrow = 4,
            top = textGrob(paste0(spider_conditions_withinCondition), gp = gpar(fontsize = 30, fontface = "bold"))
          )
        } else if (length(sorted_spider_times) < 13 && length(sorted_spider_times) > 8){
          grid_spider <- grid.arrange(
            grobs = radar_plots_ordered,
            nrow = 3,
            top = textGrob(paste0(spider_conditions_withinCondition), gp = gpar(fontsize = 30, fontface = "bold"))
          )
        } else if (length(sorted_spider_times) < 9 && length(sorted_spider_times) > 4){
          grid_spider <- grid.arrange(
            grobs = radar_plots_ordered,
            nrow = 2,
            top = textGrob(paste0(spider_conditions_withinCondition), gp = gpar(fontsize = 30, fontface = "bold"))
          )
        } else {
          grid_spider <- grid.arrange(
            grobs = radar_plots_ordered,
            nrow = 1,
            top = textGrob(paste0(spider_conditions_withinCondition), gp = gpar(fontsize = 30, fontface = "bold"))
          )
        }
        
      }
      
      print(grid_spider)
      
      lastSpider(recordPlot())
      
    }) # End Isolate
    
  }) # End renderPlot
  
  
  # Download the circos plot if clicked
  output$download_spider <- downloadHandler(
    filename = function() {
      paste0("spider_plot_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      # open a PDF device at 12″×8″
      pdf(file, width = 35, height = 15)
      # replay the last recorded Circos plot
      replayPlot(lastSpider())
      dev.off()
    },
    contentType = "application/pdf"
  )
  
  ### BEGIN SENSITIVITY HEATMAP TAB ###
  
  sensitivity_heatmap_plot <- reactive({
    
    input$sensitivity_heatmap_Update
    isolate({
      
      # Convert results to data.table and rename 'method' to 'analysis_method'
      DT <- as.data.table(results_df %>%
                            filter(
                              condition    %in% input$sensitivity_heatmap_conditions,
                              time         == input$sensitivity_heatmap_time,
                              gs.aggregate %in% input$sensitivity_heatmap_aggregates,
                              method       %in% input$sensitivity_heatmap_method
                            ))
      setnames(DT, "method", "analysis_method")
      
      # Extract all desired adjusted p-value columns (e.g., "global.adjPval_BH", etc.)
      combos <- expand.grid(
        approach = input$sensitivity_heatmap_p_approach,
        method   = input$sensitivity_heatmap_p_method,
        stringsAsFactors = FALSE
      )
      adj_cols <- with(combos, paste0(approach, ".adjPval_", method))
      
      # Reshape from wide to long format for each p-value approach × correction method
      dt_long <- data.table::melt(
        DT,
        id.vars      = c(
          "time", "gs.name", "condition", "comparison",
          "fc.score", "gs.colour", "gs.aggregate",
          "analysis_method"
        ),
        measure.vars = adj_cols,
        variable.name= "pval_type",
        value.name   = "adjp"
      )
      dt_long[, c("p_approach_spec", "p_method_spec") :=
                tstrsplit(pval_type, ".adjPval_", fixed = TRUE)]
      
      # Create a grid of reasonable parameter specifications, including analysis_method
      spec_grid <- CJ(
        p_spec           = as.numeric(input$sensitivity_heatmap_p_threshold),
        p_method_spec    = input$sensitivity_heatmap_p_method,
        p_approach_spec  = input$sensitivity_heatmap_p_approach,
        filtration_spec  = as.numeric(input$sensitivity_heatmap_fc_threshold),
        analysis_method  = input$sensitivity_heatmap_method
      )
      n_specifications <- nrow(spec_grid)
      
      # 6) Precompute per‑gene‑set metadata
      gs_meta <- DT[
        , .(
          n.comparisons = uniqueN(comparison),
          gs.colour     = gs.colour[1],
          gs.aggregate  = gs.aggregate[1]
        ),
        by = gs.name
      ]
      
      # 7) Build the full “gs.name × spec” cartesian grid
      full_grid <- CJ(
        gs.name = unique(DT$gs.name),
        spec_id = seq_len(nrow(spec_grid))
      )[
        , spec_id := NULL
      ][
        , cbind(.SD, spec_grid), by = gs.name
      ]
      
      # 8) Count observed DE per spec × gs.name
      obs_counts <- spec_grid[
        dt_long,
        on = .(p_approach_spec, p_method_spec, analysis_method),
        allow.cartesian = TRUE
      ][
        abs(fc.score) >= filtration_spec & adjp <= p_spec,
        .(n.DE = uniqueN(condition)),
        by = .(
          gs.name,
          p_spec, p_method_spec, p_approach_spec,
          filtration_spec, analysis_method
        )
      ]
      
      
      # 9) Stitch everything together and compute percentages
      setkey(gs_meta,   gs.name)
      setkey(full_grid, gs.name)
      setkey(obs_counts, gs.name)
      
      res <- full_grid[gs_meta]
      res[, n.DE := 0L]
      
      # join in observed counts
      key_cols <- c(
        "gs.name", "p_spec", "p_method_spec", "p_approach_spec",
        "filtration_spec", "analysis_method"
      )
      setkeyv(res,        key_cols)
      setkeyv(obs_counts, key_cols)
      res[obs_counts, n.DE := i.n.DE]
      
      # finalize sensitivity_gs_results
      sensitivity_gs_results <- res[
        , percent.DE := fifelse(n.DE > 0, 100L * n.DE / n.comparisons, 0L)
      ][
        , .(
          gs.name, analysis_method,
          p_spec, p_method_spec, p_approach_spec, filtration_spec,
          n.DE, n.comparisons, percent.DE,
          gs.colour, gs.aggregate
        )
      ]
      
      # Summarise and rank specifications
      summary_df <- sensitivity_gs_results %>%
        group_by(
          analysis_method,
          p_method_spec, p_approach_spec,
          filtration_spec, p_spec, gs.aggregate
        ) %>%
        summarise(percent.DE = mean(percent.DE, na.rm = TRUE), .groups = "drop")
      
      ranked_specs <- summary_df %>%
        group_by(
          analysis_method, p_spec,
          p_method_spec, p_approach_spec, filtration_spec
        ) %>%
        summarise(avg_percent_DE = mean(percent.DE, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(avg_percent_DE), desc(filtration_spec)) %>%
        mutate(rank = row_number()) %>%
        mutate(
          spec_id = paste0(
            "p=", p_spec,
            "|m=", p_method_spec,
            "|a=", p_approach_spec,
            "|f=", filtration_spec,
            "|d=", analysis_method
          )
        )
      
      # 1) Unique conditions & gene‐sets
      conds    <- DT[, unique(condition)]
      gs_names <- DT[, unique(gs.name)]
      
      # 2) Build full grid of (condition × gs.name × specs)
      full_grid2 <- CJ(
        condition      = conds,
        gs.name        = gs_names,
        spec_id        = seq_len(nrow(spec_grid)),
        unique         = TRUE
      )[
        , spec_id := NULL
      ][
        , cbind(.SD, spec_grid), by = .(condition, gs.name)
      ]
      
      # 3) Flag DE in the long table per spec
      de_flags <- spec_grid[
        dt_long,
        on = .(p_approach_spec, p_method_spec, analysis_method),
        allow.cartesian = TRUE
      ][
        abs(fc.score) >= filtration_spec & adjp <= p_spec,
        .(DE = 1L),
        by = .(
          condition, gs.name,
          p_spec, p_method_spec, p_approach_spec,
          filtration_spec, analysis_method
        )
      ]
      
      # 4) Left‑join flags into full_grid2, fill missing → DE = 0
      setkeyv(full_grid2, c(
        "condition", "gs.name",
        "p_spec", "p_method_spec", "p_approach_spec",
        "filtration_spec", "analysis_method"
      ))
      setkeyv(de_flags, key(full_grid2))
      
      full_de <- de_flags[full_grid2]
      full_de[is.na(DE), DE := 0L]
      
      # 5) Bring in aggregate labels
      agg_map <- gs_meta[, .(gs.name, gs.aggregate)]
      setkeyv(agg_map, "gs.name")
      setkeyv(full_de, "gs.name")
      full_de <- agg_map[full_de]
      
      # 6) Compute proportion DE per spec × condition × aggregate
      agg_de_summary <- full_de[
        , .(prop.DE = sum(DE) / .N),
        by = .(
          condition,
          p_spec, p_method_spec, p_approach_spec,
          filtration_spec, analysis_method,
          gs.aggregate
        )
      ] %>%
        mutate(
          spec_id = paste0(
            "p=", p_spec,
            "|m=", p_method_spec,
            "|a=", p_approach_spec,
            "|f=", filtration_spec,
            "|d=", analysis_method
          )
        ) %>%
        left_join(
          ranked_specs %>% select(spec_id, rank),
          by = "spec_id"
        ) %>%
        mutate(
          spec_id = factor(spec_id,
                           levels = ranked_specs$spec_id[order(ranked_specs$rank)])
        ) %>%
        left_join(
          DT %>% distinct(gs.aggregate, gs.colour),
          by = "gs.aggregate"
        )
      
      # Determine condition ordering
      if (input$sensitivity_heatmap_order == "available") {
        condition_order <- agg_de_summary %>%
          group_by(condition) %>%
          summarise(mean_prop_DE = mean(prop.DE, na.rm = TRUE), .groups = "drop") %>%
          arrange(desc(mean_prop_DE)) %>%
          pull(condition)
      } else {
        condition_order <- input$sensitivity_heatmap_conditions
      }
      
      agg_de_summary <- agg_de_summary %>%
        mutate(condition = factor(condition, levels = rev(condition_order)))
      
      # Prepare facet strip colours
      strip_colours <- agg_de_summary %>%
        distinct(gs.aggregate, gs.colour) %>%
        arrange(gs.aggregate) %>%
        deframe()
      facet_levels <- input$sensitivity_heatmap_aggregates
      strip_background_elements <- lapply(facet_levels, function(f) element_rect(fill = strip_colours[[f]]))
      strip_text_elements       <- lapply(facet_levels, function(f) element_text(face = "bold", color = "white"))
      
      # Build plotting grid and merge in prop.DE
      full_plot_grid <- CJ(
        condition    = condition_order,
        spec_id      = levels(agg_de_summary$spec_id),
        gs.aggregate = input$sensitivity_heatmap_aggregates,
        unique       = TRUE
      )
      heatmap_df <- full_plot_grid %>%
        left_join(
          agg_de_summary %>% select(condition, spec_id, gs.aggregate, prop.DE),
          by = c("condition", "spec_id", "gs.aggregate")
        ) %>%
        mutate(
          condition    = factor(condition,    levels = rev(condition_order)),
          spec_id      = factor(spec_id,      levels = rev(levels(agg_de_summary$spec_id))),
          gs.aggregate = factor(gs.aggregate, levels = input$sensitivity_heatmap_aggregates)
        )
      
      # Identify conditions with all NA
      na_conditions <- heatmap_df %>%
        group_by(condition) %>%
        summarize(all_na = all(is.na(prop.DE)), .groups = "drop") %>%
        filter(all_na) %>%
        pull(condition)
      label_map <- setNames(vapply(levels(heatmap_df$condition), function(lvl) {
        if (lvl %in% na_conditions) {
          # sprintf('<span style="color:#ced4da">%s</span>', lvl)
          lvl
        } else {
          lvl
        }
      }, FUN.VALUE = character(1)),
      levels(heatmap_df$condition))
      
      
      # Plot
      ggplot(heatmap_df, aes(x = spec_id, y = condition, fill = prop.DE)) +
        geom_tile() +
        scale_y_discrete(labels = label_map) +
        scale_fill_gradient(
          low      = "white",
          high     = "#70e000",
          na.value = "white",
          limits   = c(0, 1),
          name     = str_wrap("Proportion of DE genesets within aggregate", 20),
          guide    = guide_colorbar(
            title.position = "top",
            title.hjust    = 0.5,
            title.theme    = element_text(size = 20, margin = margin(b = 30)),
            label.theme    = element_text(size = 16),
            barwidth       = unit(4.5, "lines"),
            barheight      = unit(10, "lines")
          )
        ) +
        labs(
          x     = "Specifications (ranked by mean common differential expression)",
          y     = "conditions",
          title = paste0(
            "Differential Expression Across Analysis Specifications\n",
            "(Day ", input$sensitivity_heatmap_time, ")"
          )
        ) +
        facet_wrap2(
          ~ gs.aggregate,
          scales = "free_y",
          axes   = "y",
          strip  = strip_themed(
            background_x = strip_background_elements,
            text_x       = strip_text_elements
          )
        ) +
        theme_minimal(base_size = 16) +
        theme(
          plot.title   = element_text(hjust = 0.5, size = 28, face = "bold"),
          axis.title   = element_text(size = 30),
          strip.text.x = element_text(size = 12, face = "bold", color = "black"),
          axis.text.x  = element_blank(),
          axis.text.y = ggtext::element_markdown(size = 10),
        )
    }) # End Isolate
    
  }) # End reactive
  
  output$sensitivityHeatmapPlot <- renderPlot({
    sensitivity_heatmap_plot()
  })
  
  # Download handler for high-resolution plot
  output$download_sensitivity_heatmap <- downloadHandler(
    filename = function() {
      req(input$sensitivity_heatmap_download_format)
      paste0("sensitivity_heatmap", Sys.Date(), ".",
             input$sensitivity_heatmap_download_format)
    },
    content = function(file) {
      # make sure the plot reactive is available
      req(sensitivity_heatmap_plot())
      
      # call ggsave with user inputs
      ggsave(
        filename = file,
        plot     = sensitivity_heatmap_plot(),
        device   = input$sensitivity_heatmap_download_format,
        width    = input$sensitivity_heatmap_download_width,
        height   = input$sensitivity_heatmap_download_height,
        dpi      = input$sensitivity_heatmap_download_dpi,
        units    = "in"
      )
    }
  )
  
  ### BEGIN INDIVIDUAL PARAMETER SENSITIVITY PLOT TAB ### 
  indiv_sensitivity_plot <- reactive({
    
    input$indiv_sensitivity_Update
    isolate({
      
      # First create a new data table with the data of interest
      DT <- as.data.table(results_df) %>% 
        filter(time == input$indiv_sensitivity_time,
               condition %in% input$indiv_sensitivity_conditions,
               gs.aggregate %in% input$indiv_sensitivity_aggregates,
               !is.na(rawPval))
      
      # Now convert this to long format to have one row per specification
      dt_long <- DT %>%
        pivot_longer(
          cols = matches("\\.adjPval_"),             # all columns like “approach.adjPval_method”
          names_to  = c("approach", "p_method"),     # two new cols
          names_pattern = "(.*)\\.adjPval_(.*)",     # regex: split at “.adjPval_”
          values_to = "adj_pval"                     # the cell value
        )
      
      # Now, according to the fixed parameter values, we can calculate the percentage of DE genesets under each of the varying specifications. 
      if (input$indiv_sensitivity_parameter == "method"){
        
        dt_long_filtered = dt_long %>% 
          filter(approach == input$indiv_sensitivity_p_approach,
                 p_method == input$indiv_sensitivity_p_correction) %>% 
          mutate(significant = (adj_pval <= input$indiv_sensitivity_p_threshold &
                                  abs(fc.score) >= input$indiv_sensitivity_fc_threshold)) 
        
        df_prop_de <- dt_long_filtered %>%
          group_by(condition, method, condition.colour) %>% 
          summarise(
            prop_signif = mean(significant),
            .groups = "drop"
          ) %>% 
          # reorder 'condition' by its mean prop_signif (across methods)
          mutate(
            condition = fct_reorder(condition, prop_signif, .fun = mean)
          )
        
        if (input$indiv_sensitivity_groupby == "byparameter"){
          plot <- df_prop_de %>% 
            ggplot(aes(x = method, y = prop_signif, fill = condition)) +
            
            scale_fill_manual(
              name   = "Vaccine",
              values = set_names(df_prop_de$condition.colour, df_prop_de$condition),
              breaks = levels(df_prop_de$condition),
              labels = levels(df_prop_de$condition)
            ) +
            
            geom_col(position = position_dodge(width = 0.7), width = 0.6) +
            
            labs(
              title = glue::glue(
                "Proportion of DE Genesets by DGSA Method\n",
                "(Day {input$indiv_sensitivity_time})"
              ),
              x       = "DGSA Method",
              y       = "Proportion of DE genesets",
              caption = paste0(
                "Fixed parameters:\n",
                "- p value correction: ", input$indiv_sensitivity_p_correction, "     ",
                "- p value subset: ",              input$indiv_sensitivity_p_approach,     "\n",
                "- significance level: ",           input$indiv_sensitivity_p_threshold,     "     ",
                "- abs log2 FC threshold: ",         input$indiv_sensitivity_fc_threshold
              )
            ) +
            
            theme_minimal(base_size = 14) +
            theme(
              plot.title            = element_text(hjust = 0.5, face = "bold", size = 40),
              axis.title            = element_text(face = "bold", size = 30),
              axis.text             = element_text(size = 25),
              legend.title          = element_text(face = "bold", size = 25),
              legend.text           = element_text(size = 20),
              
              # caption styling
              plot.caption          = element_text(
                size       = 14,
                hjust      = 0,
                lineheight = 1.2,
                margin     = margin(t = 15, r = 0, b = 0, l = 0)
              ),
              plot.caption.position = "plot"
            ) +
            
            ylim(0, 1)
        } else if (input$indiv_sensitivity_groupby == "byvaccine"){
          plot <- df_prop_de %>%
            ggplot(aes(
              x      = condition,
              y      = prop_signif,
              fill   = method
            )) +
            
            # use your own named vector of method‐colours here
            scale_fill_manual(
              name   = "DGSA Method",
              values = set_names(
                c("#1b9e77", "#d95f02"),            # e.g. one hex per method
                unique(df_prop_de$method)                       # matching the factor levels
              ),
              breaks = unique(df_prop_de$method),
              labels = unique(df_prop_de$method)
            ) +
            
            geom_col(
              position = position_dodge(width = 0.7),
              width    = 0.6
            ) +
            
            labs(
              title = glue::glue(
                "Proportion of DE Genesets by DGSA Method\n",
                "(Day {input$indiv_sensitivity_time})"
              ),
              x       = "Vaccine",
              y       = "Proportion of DE genesets",
              caption = paste0(
                "Fixed parameters:\n",
                "- p value correction: ", input$indiv_sensitivity_p_correction, "    ",
                "- p value subset: ",    input$indiv_sensitivity_p_approach,     "\n",
                "- significance level: ", input$indiv_sensitivity_p_threshold,    "    ",
                "- abs log2 FC threshold: ",           input$indiv_sensitivity_fc_threshold
              )
            ) +
            
            theme_minimal(base_size = 14) +
            theme(
              plot.title    = element_text(hjust = 0.5, face = "bold", size = 40),
              axis.title    = element_text(face = "bold", size = 30),
              axis.text     = element_text(size = 25),
              axis.text.x   = element_text(angle = 45, hjust = 1),    # ← rotated labels
              legend.title  = element_text(face = "bold", size = 25),
              legend.text   = element_text(size = 20),
              plot.caption  = element_text(size = 14, hjust = 0, lineheight = 1.2),
              plot.caption.position = "plot"
            ) +
            
            ylim(0, 1)
        }
        
      } else if (input$indiv_sensitivity_parameter == "p_correction"){
        
        dt_long_filtered = dt_long %>% 
          filter(approach == input$indiv_sensitivity_p_approach,
                 method == input$indiv_sensitivity_method) %>% 
          mutate(significant = (adj_pval <= input$indiv_sensitivity_p_threshold &
                                  abs(fc.score) >= input$indiv_sensitivity_fc_threshold)) 
        
        df_prop_de <- dt_long_filtered %>%
          group_by(condition, p_method, condition.colour) %>% 
          summarise(
            prop_signif = mean(significant),
            .groups = "drop"
          ) %>% 
          # reorder 'condition' by its mean prop_signif (across methods)
          mutate(
            condition = fct_reorder(condition, prop_signif, .fun = mean)
          )
        
        if (input$indiv_sensitivity_groupby == "byparameter"){
          plot <- df_prop_de %>% 
            ggplot(aes(x = p_method, y = prop_signif, fill = condition)) +
            
            # use the hex‐codes in condition.colour, in the same order as levels(condition)
            scale_fill_manual(
              name   = "Vaccine",
              values = set_names(df_prop_de$condition.colour, df_prop_de$condition),
              breaks = levels(df_prop_de$condition),
              labels = levels(df_prop_de$condition)
            ) +
            
            geom_col(position = position_dodge(width = 0.7), width = 0.6) +
            
            labs(
              title = paste0(
                "Proportion of DE Genesets by p-value Correction Method \n",
                "(Day ", input$indiv_sensitivity_time, ")"
              ),
              x     = "p-value Correction Method",
              y     = "Proportion of DE genesets",
              caption = paste0(
                "Fixed parameters:\n",
                "- DGSA Method: ", input$indiv_sensitivity_method, "     ",
                "- p value subset: ",              input$indiv_sensitivity_p_approach,     "\n",
                "- significance level: ",           input$indiv_sensitivity_p_threshold,     "     ",
                "- abs log2 FC threshold: ",         input$indiv_sensitivity_fc_threshold
              )
            ) +
            
            theme_minimal(base_size = 14) +
            theme(
              plot.title   = element_text(hjust = 0.5, face = "bold", size = 40),
              axis.title   = element_text(face = "bold", size = 30),
              axis.text    = element_text(size = 25),
              legend.title = element_text(face = "bold", size = 25),
              legend.text  = element_text(size = 20),
              plot.caption          = element_text(
                size       = 14,
                hjust      = 0,
                lineheight = 1.2,
                margin     = margin(t = 15, r = 0, b = 0, l = 0)
              ),  
              plot.caption.position = "plot"
            ) +
            ylim(0, 1)
        } else if (input$indiv_sensitivity_groupby == "byvaccine"){
          plot <- df_prop_de %>%
            ggplot(aes(
              x    = condition,
              y    = prop_signif,
              fill = p_method
            )) +
            
            # map each correction method to a colour
            scale_fill_manual(
              name   = "p‑value Correction Method",
              values = set_names(
                c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"),
                unique(df_prop_de$p_method)          # match factor levels
              ),
              breaks = unique(df_prop_de$p_method),
              labels = unique(df_prop_de$p_method)
            ) +
            
            geom_col(
              position = position_dodge(width = 0.7),
              width    = 0.6
            ) +
            
            labs(
              title = glue::glue(
                "Proportion of DE Genesets by Correction Method \n",
                "(Day {input$indiv_sensitivity_time})"
              ),
              x       = "Vaccine",
              y       = "Proportion of DE genesets",
              caption = paste0(
                "Fixed parameters:\n",
                "- DGSA Method: ", input$indiv_sensitivity_method, "    ",
                "- p‑value subset: ", input$indiv_sensitivity_p_approach, "\n",
                "- significance level: ", input$indiv_sensitivity_p_threshold,    "    ",
                "- abs log2‑FC threshold: ", input$indiv_sensitivity_fc_threshold
              )
            ) +
            
            theme_minimal(base_size = 14) +
            theme(
              plot.title    = element_text(hjust = 0.5, face = "bold", size = 40),
              axis.title    = element_text(face = "bold", size = 30),
              axis.text     = element_text(size = 25),
              axis.text.x   = element_text(angle = 45, hjust = 1),    # rotated labels
              legend.title  = element_text(face = "bold", size = 25),
              legend.text   = element_text(size = 20),
              plot.caption  = element_text(size = 14, hjust = 0, lineheight = 1.2),
              plot.caption.position = "plot"
            ) +
            
            ylim(0, 1)
        }
        
      } else if (input$indiv_sensitivity_parameter == "p_approach"){
        
        dt_long_filtered = dt_long %>% 
          filter(p_method == input$indiv_sensitivity_p_correction,
                 method == input$indiv_sensitivity_method) %>% 
          mutate(significant = (adj_pval <= input$indiv_sensitivity_p_threshold &
                                  abs(fc.score) >= input$indiv_sensitivity_fc_threshold)) 
        
        df_prop_de <- dt_long_filtered %>%
          group_by(condition, approach, condition.colour) %>% 
          summarise(
            prop_signif = mean(significant),
            .groups = "drop"
          ) %>% 
          # reorder 'condition' by its mean prop_signif (across methods)
          mutate(
            condition = fct_reorder(condition, prop_signif, .fun = mean)
          )
        
        if (input$indiv_sensitivity_groupby == "byparameter"){
          plot <- df_prop_de %>% 
            ggplot(aes(x = approach, y = prop_signif, fill = condition)) +
            
            # use the hex‐codes in condition.colour, in the same order as levels(condition)
            scale_fill_manual(
              name   = "Vaccine",
              values = set_names(df_prop_de$condition.colour, df_prop_de$condition),
              breaks = levels(df_prop_de$condition),
              labels = levels(df_prop_de$condition)
            ) +
            
            geom_col(position = position_dodge(width = 0.7), width = 0.6) +
            
            labs(
              title = paste0(
                "Proportion of DE Genesets by p-value adjustment subset \n",
                "(Day ", input$indiv_sensitivity_time, ")"
              ),
              x     = "p-value adjustment subset",
              y     = "Proportion of DE genesets",
              caption = paste0(
                "Fixed parameters:\n",
                "- DGSA Method: ", input$indiv_sensitivity_method, "     ",
                "- p value correction: ",              input$indiv_sensitivity_p_correction,     "\n",
                "- significance level: ",           input$indiv_sensitivity_p_threshold,     "     ",
                "- abs log2 FC threshold: ",         input$indiv_sensitivity_fc_threshold
              )
            ) +
            
            theme_minimal(base_size = 14) +
            theme(
              plot.title   = element_text(hjust = 0.5, face = "bold", size = 40),
              axis.title   = element_text(face = "bold", size = 30),
              axis.text    = element_text(size = 25),
              legend.title = element_text(face = "bold", size = 25),
              legend.text  = element_text(size = 20),
              plot.caption          = element_text(
                size       = 14,
                hjust      = 0,
                lineheight = 1.2,
                margin     = margin(t = 15, r = 0, b = 0, l = 0)
              ),  
              plot.caption.position = "plot"
            ) +
            ylim(0, 1)
        } else if (input$indiv_sensitivity_groupby == "byvaccine"){
          plot <- df_prop_de %>%
            ggplot(aes(
              x    = condition,
              y    = prop_signif,
              fill = approach
            )) +
            
            # map each p‑value adjustment subset to a colour
            scale_fill_manual(
              name   = "p‑value Subset",
              values = set_names(
                c("#1b9e77", "#d95f02", "#7570b3"),   # one hex per approach
                unique(df_prop_de$approach)           # match factor levels
              ),
              breaks = unique(df_prop_de$approach),
              labels = unique(df_prop_de$approach)
            ) +
            
            geom_col(
              position = position_dodge(width = 0.7),
              width    = 0.6
            ) +
            
            labs(
              title = glue::glue(
                "Proportion of DE Genesets by p-value adjustment subset\n",
                "(Day {input$indiv_sensitivity_time})"
              ),
              x       = "Vaccine",
              y       = "Proportion of DE genesets",
              caption = paste0(
                "Fixed parameters:\n",
                "- DGSA Method: ",           input$indiv_sensitivity_method,      "    ",
                "- p value correction: ",    input$indiv_sensitivity_p_correction, "\n",
                "- significance level: ",    input$indiv_sensitivity_p_threshold,  "    ",
                "- abs log2‑FC threshold: ", input$indiv_sensitivity_fc_threshold
              )
            ) +
            
            theme_minimal(base_size = 14) +
            theme(
              plot.title    = element_text(hjust = 0.5, face = "bold", size = 40),
              axis.title    = element_text(face = "bold", size = 30),
              axis.text     = element_text(size = 25),
              axis.text.x   = element_text(angle = 45, hjust = 1),
              legend.title  = element_text(face = "bold", size = 25),
              legend.text   = element_text(size = 20),
              plot.caption  = element_text(size = 14, hjust = 0, lineheight = 1.2),
              plot.caption.position = "plot"
            ) +
            
            ylim(0, 1)
        }
        
      } else if (input$indiv_sensitivity_parameter == "p_threshold"){
        
        dt_long_filtered = dt_long %>% 
          filter(p_method == input$indiv_sensitivity_p_correction,
                 method == input$indiv_sensitivity_method,
                 approach == input$indiv_sensitivity_p_approach) 
        
        dt_long_filtered <- bind_cols(
          dt_long_filtered,
          map_dfc(input$indiv_sensitivity_p_threshold_range, function(thresh_char) {
            thresh <- as.numeric(thresh_char)
            indicator <- with(dt_long_filtered,
                              abs(fc.score) >= input$indiv_sensitivity_fc_threshold &
                                adj_pval <= thresh
            )
            colname <- paste0("significant_threshold_", gsub("\\.", "_", thresh_char))
            tibble(!!colname := indicator)
          })
        )
        
        df_prop_de <- dt_long_filtered %>%
          # 1) pivot the threshold‐columns into long format
          pivot_longer(
            cols      = starts_with("significant_threshold_"),
            names_to  = "threshold",
            names_prefix = "significant_threshold_",
            values_to = "significant"
          ) %>%
          # 2) parse the threshold name into a numeric
          mutate(
            threshold = as.numeric(str_replace(threshold, "_", "."))
          ) %>%
          # 3) compute the mean(significant) per group
          group_by(condition, condition.colour, threshold) %>%
          summarise(
            prop_signif = mean(significant),
            .groups     = "drop"
          ) %>% 
          mutate(condition = fct_reorder(condition, prop_signif, .fun = mean))
        
        custom_pct <- function(x) {
          pct <- x * 100
          vapply(pct, function(p) {
            if (p < 1) {
              # format with exactly two decimals, then drop trailing zeros
              txt <- formatC(p, format = "f", digits = 5)
              txt <- sub("\\.?0+$", "", txt)
            } else {
              # no decimals for 1% and above
              txt <- as.character(round(p))
            }
            paste0(txt, "%")
          }, FUN.VALUE = character(1))
        }
        
        plot <- df_prop_de %>% 
          ggplot(aes(x = threshold, 
                     y = prop_signif, 
                     colour = condition, 
                     group = condition)) +
          scale_color_manual(
            name   = "Vaccine",
            values = set_names(df_prop_de$condition.colour, df_prop_de$condition),
            breaks = levels(df_prop_de$condition),
            labels = levels(df_prop_de$condition)
          ) +
          geom_point(size = 2) +
          geom_line(linewidth = 2, alpha = 0.7) +
          
          # <— Fix: coerce breaks to numeric
          scale_x_log10(
            breaks = as.numeric(input$indiv_sensitivity_p_threshold_range),
            labels = custom_pct,
          ) +
          
          labs(
            title = paste0(
              "Proportion of DE Genesets by significance level\n",
              "(Day ", input$indiv_sensitivity_time, ")"
            ),
            x = "Significance level",
            y = "Proportion of DE genesets",
            caption = paste0(
              "Fixed parameters:\n",
              "- DGSA Method: ", input$indiv_sensitivity_method, "     ",
              "- p value correction: ",              input$indiv_sensitivity_p_correction,     "\n",
              "- p-value subset: ",           input$indiv_sensitivity_p_approach,     "     ",
              "- abs log2 FC threshold: ",         input$indiv_sensitivity_fc_threshold
            )
          ) +
          theme_minimal(base_size = 14) +
          theme(
            plot.title   = element_text(hjust = 0.5, face = "bold", size = 40),
            axis.title   = element_text(face = "bold", size = 30),
            axis.text    = element_text(size = 25),
            legend.title = element_text(face = "bold", size = 25),
            legend.text  = element_text(size = 20),
            plot.caption          = element_text(
              size       = 14,
              hjust      = 0,
              lineheight = 1.2,
              margin     = margin(t = 15, r = 0, b = 0, l = 0)
            ),  
            plot.caption.position = "plot"
          ) +
          ylim(0, 1)
        
      } else if (input$indiv_sensitivity_parameter == "fc_threshold"){
        
        dt_long_filtered = dt_long %>% 
          filter(p_method == input$indiv_sensitivity_p_correction,
                 method == input$indiv_sensitivity_method,
                 approach == input$indiv_sensitivity_p_approach) 
        
        dt_long_filtered <- bind_cols(
          dt_long_filtered,
          map_dfc(input$indiv_sensitivity_fc_threshold_range, function(thresh) {
            # Create a logical vector
            indicator <- with(dt_long_filtered,
                              abs(fc.score) >= thresh &
                                adj_pval <= input$indiv_sensitivity_p_threshold)
            # Name the column like significant_threshold_0.01 etc.
            colname <- paste0("significant_threshold_", gsub("\\.", "_", format(thresh)))
            tibble(!!colname := indicator)
          })
        )
        
        df_prop_de <- dt_long_filtered %>%
          # 1) pivot the threshold‐columns into long format
          pivot_longer(
            cols      = starts_with("significant_threshold_"),
            names_to  = "threshold",
            names_prefix = "significant_threshold_",
            values_to = "significant"
          ) %>%
          # 2) parse the threshold name into a numeric
          mutate(
            threshold = as.numeric(str_replace(threshold, "_", "."))
          ) %>%
          # 3) compute the mean(significant) per group
          group_by(condition, condition.colour, threshold) %>%
          summarise(
            prop_signif = mean(significant),
            .groups     = "drop"
          )
        
        plot <- df_prop_de %>% 
          ggplot(aes(x = threshold, y = prop_signif, colour = condition, group = condition)) +
          
          scale_color_manual(
            name   = "Vaccine",
            values = set_names(df_prop_de$condition.colour, df_prop_de$condition),
            breaks = levels(df_prop_de$condition),
            labels = levels(df_prop_de$condition)
          ) +
          
          geom_point(size = 2) +
          geom_line(linewidth = 2,
                    alpha = 0.7) +
          
          labs(
            title = paste0(
              "Proportion of DE Genesets by Absolute log2-fold-change threshold\n",
              "(Day ", input$indiv_sensitivity_time, ")"
            ),
            x = "Absolute log2-FC threshold",
            y = "Proportion of DE genesets",
            caption = paste0(
              "Fixed parameters:\n",
              "- DGSA Method: ", input$indiv_sensitivity_method, "     ",
              "- p value correction: ",              input$indiv_sensitivity_p_correction,     "\n",
              "- p-value subset: ",           input$indiv_sensitivity_p_approach,     "     ",
              "- significance level: ",         input$indiv_sensitivity_p_threshold
            )
          ) +
          
          theme_minimal(base_size = 14) +
          theme(
            plot.title   = element_text(hjust = 0.5, face = "bold", size = 40),
            axis.title   = element_text(face = "bold", size = 30),
            axis.text    = element_text(size = 25),
            legend.title = element_text(face = "bold", size = 25),
            legend.text  = element_text(size = 20),
            plot.caption          = element_text(
              size       = 14,
              hjust      = 0,
              lineheight = 1.2,
              margin     = margin(t = 15, r = 0, b = 0, l = 0)
            ),  
            plot.caption.position = "plot"
          ) +
          
          ylim(0, 1)
      }
      
      plot
      
    }) # End isolate
    
  }) # End reactive
  
  output$indivSensitivityPlot <- renderPlot({
    indiv_sensitivity_plot()
  })
  
  # Download handler for high-resolution plot
  output$download_indiv_sensitivity <- downloadHandler(
    filename = function() {
      req(input$indiv_sensitivity_download_format)
      paste0("indiv_sensitivity_", 
             input$indiv_sensitivity_parameter, 
             "_day",
             input$indiv_sensitivity_time,
             ".",
             input$indiv_sensitivity_download_format)
    },
    content = function(file) {
      # make sure the plot reactive is available
      req(indiv_sensitivity_plot())
      
      # call ggsave with user inputs
      ggsave(
        filename = file,
        plot     = indiv_sensitivity_plot(),
        device   = input$indiv_sensitivity_download_format,
        width    = input$indiv_sensitivity_download_width,
        height   = input$indiv_sensitivity_download_height,
        dpi      = input$indiv_sensitivity_download_dpi,
        units    = "in"
      )
    }
  )
  
  ### END INDIVIDUAL PARAMETER SENSITIVITY PLOT TAB ###
  
  
  ### BEGIN METHOD-UNIQUE GENESETS TAB ###
  
  method_unique_plot <- reactive({
    
    input$method_unique_Update
    isolate({
      
      results_df <- results_df %>%
        mutate(time = factor(time, levels = sort(unique(time))))
      
      # 1. build your “significant” table across all selected times
      #    (note we now keep time in the data, and allow length(input$...time) ≥ 1)
      pval_col <- paste0(
        input$method_unique_p_approach, 
        ".adjPval_", 
        input$method_unique_p_correction
      )
      
      sig_df <- results_df %>%
        filter(
          condition     %in% input$method_unique_conditions,
          time          %in% input$method_unique_time,
          gs.aggregate  %in% input$method_unique_aggregates
        ) %>%
        mutate(
          significant = (.data[[pval_col]] <= input$method_unique_p_threshold) &
            (abs(fc.score)    >= input$method_unique_fc_threshold)
        ) %>%
        filter(significant) %>%
        select(time, condition, method, gs.name, gs.aggregate, gs.colour, fc.score)
      
      # 1b. optional positive/negative filtration
      if (input$method_unique_filtration == "positive") {
        sig_df <- filter(sig_df, fc.score > 0)
      } else if (input$method_unique_filtration == "negative") {
        sig_df <- filter(sig_df, fc.score < 0)
      }
      
      # 2. build ref- & cmp‐lists *by* time & condition
      ref_sets <- sig_df %>%
        filter(method == input$method_unique_reference) %>%
        group_by(time, condition) %>%
        summarise(ref_list = list(gs.name), .groups = "drop")
      
      cmp_sets <- sig_df %>%
        filter(method == input$method_unique_comparison) %>%
        group_by(time, condition) %>%
        summarise(cmp_list = list(gs.name), .groups = "drop")
      
      # 3. compute unique sets (cmp − ref), keeping time
      unique_df <- inner_join(cmp_sets, ref_sets, by = c("time","condition")) %>%
        mutate(unique = map2(cmp_list, ref_list, setdiff)) %>%
        select(time, condition, unique) %>%
        unnest(unique) %>%
        rename(gs.name = unique)
      
      # 4. attach aggregate & colour; count & proportion *by* time & condition
      plot_df <- unique_df %>%
        left_join(
          sig_df %>% distinct(gs.name, gs.aggregate, gs.colour),
          by = "gs.name"
        ) %>%
        group_by(time, condition, gs.aggregate, gs.colour) %>%
        summarise(count = n(), .groups = "drop") %>%
        group_by(time, condition) %>%
        mutate(prop = count / sum(count)) %>%
        ungroup()
      
      # 5. plot: horizontal bars, faceted by time. 
      #    free_y so each facet only shows its own vaccines; x‐scale fixed across facets
      ggplot(plot_df, aes(
        x    = count,
        y    = fct_reorder(condition, count, .fun = sum), 
        fill = gs.aggregate
      )) +
        geom_col(position = "stack") +
        facet_wrap(
          ~ time, 
          ncol         = 1, 
          scales       = "free_y", 
          strip.position = "top",
          labeller     = labeller(time = function(x) paste0("Day ", x))
        ) +
        scale_fill_manual(
          name   = "Gene Set Aggregate",
          values = aggregate_colors     # assumes that named vector is in your environment
        ) +
        labs(
          title    = glue(
            "Unique Genesets: {input$method_unique_comparison} vs. {input$method_unique_reference}"
          ),
          subtitle = glue("Day(s): {paste0(input$method_unique_time, collapse = ', ')}"),
          x        = "Number of uniquely identified genesets",
          y        = NULL,
          caption  = paste0(
            "Fixed specs:\n",
            "- p value subset: ",    input$method_unique_p_approach,    "    ",
            "- p value correction: ", input$method_unique_p_correction, "\n",
            "- significance level: ", input$method_unique_p_threshold,   "    ",
            "- abs log2 FC threshold: ", input$method_unique_fc_threshold
          )
        ) +
        theme_minimal(base_size = 20) +
        theme(
          strip.text      = element_text(face = "bold", size = 24),
          axis.text.y     = element_text(size = 18),
          axis.title.x    = element_text(face = "bold", size = 20),
          legend.title    = element_text(face = "bold", size = 20),
          legend.text     = element_text(size = 16),
          plot.title      = element_text(face = "bold", size = 30, hjust = 0.5),
          plot.subtitle   = element_text(size = 22, hjust = 0.5),
          plot.caption    = element_text(size = 12, hjust = 0)
        )
      
    }) # End isolate
    
  }) # End reactive
  
  output$methodUniquePlot <- renderPlot({
    method_unique_plot()
  })
  
  # Download handler for high-resolution plot
  output$download_method_unique <- downloadHandler(
    filename = function() {
      req(input$method_unique_download_format)
      paste0("method_unique_genesets_day", 
             input$method_unique_time,
             ".",
             input$method_unique_download_format)
    },
    content = function(file) {
      # make sure the plot reactive is available
      req(method_unique_plot())
      
      # call ggsave with user inputs
      ggsave(
        filename = file,
        plot     = method_unique_plot(),
        device   = input$method_unique_download_format,
        width    = input$method_unique_download_width,
        height   = input$method_unique_download_height,
        dpi      = input$method_unique_download_dpi,
        units    = "in"
      )
    }
  )
  
  ### END METHOD-UNIQUE GENESETS TAB ###
  
} # End server

shinyApp(ui, server)
