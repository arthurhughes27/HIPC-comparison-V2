suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(fs)
  library(shiny)
  library(shinyBS)
  library(patchwork)
  library(grid)
  library(gridExtra)
})

# Load in all the possible results that shiny may need access to

# Cumulative approach, with clinical variables, with TBA
p_prediction_results_all_cumulative_withClinical_withTBA = fs::path(
  "output",
  "results",
  "prediction",
  "prediction_results_all_cumulative_withClinical_withTBA.rds"
)

prediction_results_all_cumulative_withClinical_withTBA = readRDS(p_prediction_results_all_cumulative_withClinical_withTBA)

# Cumulative approach, without clinical variables, with TBA
p_prediction_results_all_cumulative_withoutClinical_withTBA = fs::path(
  "output",
  "results",
  "prediction",
  "prediction_results_all_cumulative_withoutClinical_withTBA.rds"
)

prediction_results_all_cumulative_withoutClinical_withTBA = readRDS(p_prediction_results_all_cumulative_withoutClinical_withTBA)

# Sequential approach, with clinical variables, with TBA
p_prediction_results_all_sequential_withClinical_withTBA = fs::path(
  "output",
  "results",
  "prediction",
  "prediction_results_all_sequential_withClinical_withTBA.rds"
)

prediction_results_all_sequential_withClinical_withTBA = readRDS(p_prediction_results_all_sequential_withClinical_withTBA)

# Sequential approach, without clinical variables, with TBA
p_prediction_results_all_sequential_withoutClinical_withTBA = fs::path(
  "output",
  "results",
  "prediction",
  "prediction_results_all_sequential_withoutClinical_withTBA.rds"
)

prediction_results_all_sequential_withoutClinical_withTBA = readRDS(p_prediction_results_all_sequential_withoutClinical_withTBA)

# Cumulative approach, with clinical variables, without TBA
p_prediction_results_all_cumulative_withClinical_withoutTBA = fs::path(
  "output",
  "results",
  "prediction",
  "prediction_results_all_cumulative_withClinical_withoutTBA.rds"
)

prediction_results_all_cumulative_withClinical_withoutTBA = readRDS(p_prediction_results_all_cumulative_withClinical_withoutTBA)

# Cumulative approach, without clinical variables, without TBA
p_prediction_results_all_cumulative_withoutClinical_withoutTBA = fs::path(
  "output",
  "results",
  "prediction",
  "prediction_results_all_cumulative_withoutClinical_withoutTBA.rds"
)

prediction_results_all_cumulative_withoutClinical_withoutTBA = readRDS(p_prediction_results_all_cumulative_withoutClinical_withoutTBA)

# Sequential approach, with clinical variables, without TBA
p_prediction_results_all_sequential_withClinical_withoutTBA = fs::path(
  "output",
  "results",
  "prediction",
  "prediction_results_all_sequential_withClinical_withoutTBA.rds"
)

prediction_results_all_sequential_withClinical_withoutTBA = readRDS(p_prediction_results_all_sequential_withClinical_withoutTBA)

# Sequential approach, without clinical variables, without TBA
p_prediction_results_all_sequential_withoutClinical_withoutTBA = fs::path(
  "output",
  "results",
  "prediction",
  "prediction_results_all_sequential_withoutClinical_withoutTBA.rds"
)

prediction_results_all_sequential_withoutClinical_withoutTBA = readRDS(p_prediction_results_all_sequential_withoutClinical_withoutTBA)

### DEFINE USER PARAMETERS ###

### RESULTS COMPARISON TAB PARAMETERS ###

approaches_set = c("Sequential" = "sequential", "Cumulative" = "cumulative")

include_clinical_set = c(
  "Include clinical variables" = "withClinical",
  "Do not include clinical variables" = "withoutClinical"
)

include_tba_set = c("Include TBA" = "withTBA", "Exclude TBA" = "withoutTBA")

highlight_best_set = c("No highlighting" = "none",
                       "Spearman Correlation" = "spearman",
                       "Standardised RMSE" = "sRMSE",
                       "sRMSE + (1-Spearman)" = "combined")

### CV PREDICTIONS TAB ###

vaccine_set = names(prediction_results_all_sequential_withoutClinical_withTBA)

all_sets <- unique(unlist(
  lapply(
    prediction_results_all_sequential_withClinical_withTBA,
    names
  )
))

# Extract numeric part from strings like "day 7"
day_sets <- all_sets[grepl("^Day\\s*[0-9]+$", all_sets)]
day_nums <- as.numeric(sub("Day\\s*", "", day_sets))

feature_set_set <- c("clinical", paste0("Day ", sort(day_nums)))

### VARIABLE IMPORTANCE TAB ###


### DEFINE R SHINY PAGE LAYOUT ###
ui <- fluidPage(
  # Aesthetic style options
  tags$head(tags$style(
    HTML(
      "
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
      "
    )
  )),
  
  # Title of application
  titlePanel("Vaccine-Transcriptomics Predictions Application"),
  
  # Tab layouts
  tabsetPanel(
    id = "main_tabs",
    
    # 0) Home/about page
    tabPanel(title = "About", value = "about_tab", fluidRow(column(
      width = 12,
      h2("Application Background"),
      p(
        "This RShiny app allows you to explore results of an analysis of the prediction of vaccine response with transcriptomic data."
      ),
      p(
        "The 'Results Comparison' tab allows for the comparison of cross-validated model prediction metrics across vaccines and features."
      ),
      p(
        "The 'Cross-Validation Predictions' tab allows for the visualisation of the cross-validation predictions for a particular vaccine and feature set."
      ),
      p(
        "The 'Variable Importance' tab allows for the visualisation of the mean Random Forest variable importance across cross-validation folds for a particular vaccine and feature set."
      )
    ))),
    # Close tabPanel for Home page
    
    ### TAB 1 : RESULTS COMPARISON ###
    tabPanel(
      title = "Results Comparison",
      value = "results_comparison_tab",
      sidebarLayout(
        sidebarPanel(
          # Specify width of panel
          width = 3,
          
          # Update plot button
          actionButton(
            inputId = "results_comparison_update",
            label   = "Update Figure",
            class   = "btn-primary",
            width   = '100%'
          ),
          # close actionButton
          
          # Data selection group
          wellPanel(
            style = "background-color: #c6d2ed;",
            h4("Results selection"),
            
            # Results from which approach?
            tags$label(
              "Temporal feature selection approach:",
              tags$span(icon("info-circle"), id    = "results_comparison_approach_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
            ),
            selectInput(
              inputId = "results_comparison_approach",
              label   = NULL,
              choices = approaches_set,
              selected = approaches_set[1],
              multiple = FALSE
            ),
            bsTooltip(
              id      = "results_comparison_approach_info",
              title   = "Choose the feature selection approach.",
              placement = "right",
              trigger   = "hover"
            ),
            
            # Include clinical variables?
            tags$label(
              "Include baseline clinical variables?",
              tags$span(icon("info-circle"), id    = "results_comparison_include_clinical_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
            ),
            selectInput(
              inputId = "results_comparison_include_clinical",
              label   = NULL,
              choices = include_clinical_set,
              selected = include_clinical_set[1],
              multiple = FALSE
            ),
            bsTooltip(
              id      = "results_comparison_include_clinical_info",
              title   = "Include baseline clinical variables as predictors?",
              placement = "right",
              trigger   = "hover"
            ),
            
            # Include TBA modules?
            tags$label(
              "Include TBA modules?",
              tags$span(icon("info-circle"), id    = "results_comparison_include_tba_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
            ),
            selectInput(
              inputId = "results_comparison_include_tba",
              label   = NULL,
              choices = include_tba_set,
              selected = include_tba_set[2],
              multiple = FALSE
            ),
            bsTooltip(
              id      = "results_comparison_include_tba_info",
              title   = "Include modules without labels as predictors?",
              placement = "right",
              trigger   = "hover"
            ),
            
            # Include TBA modules?
            tags$label(
              "Include which feature sets?",
              tags$span(icon("info-circle"), id    = "results_comparison_feature_set_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
            ),
            selectInput(
              inputId = "results_comparison_feature_set",
              label   = NULL,
              choices = feature_set_set,
              selected = feature_set_set,
              multiple = TRUE
            ),
            bsTooltip(
              id      = "results_comparison_feature_set_info",
              title   = "Include which feature sets for the evaluation?",
              placement = "right",
              trigger   = "hover"
            )
          ),
          # Close first wellPanel (data selection)
          
          # ---------- Download section ----------
          wellPanel(
            style = "background-color: #ced4da;",
            h4("Download Options"),
            # choose format
            selectInput(
              inputId = "results_comparison_download_format",
              label   = "Format of download:",
              choices = c(
                "PDF" = "pdf",
                "PNG" = "png",
                "JPEG" = "jpeg"
              ),
              selected = "pdf"
            ),
            
            # size controls (in inches)
            numericInput(
              inputId = "results_comparison_download_width",
              label   = "Width of download (inches):",
              value   = 15,
              min     = 1
            ),
            numericInput(
              inputId = "results_comparison_download_height",
              label   = "Height of download (inches):",
              value   = 10,
              min     = 1
            ),
            
            numericInput(
              inputId = "results_comparison_download_dpi",
              label   = "Resolution (DPI):",
              value   = 300,
              min     = 72
            ),
            
            ## Download plot button
            div(
              style = "width: 100%;",
              downloadButton(
                outputId = "download_results_comparison_plot",
                label    = "Download Results Comparison Plot",
                class    = "btn-success"
              )
            )
          ) # Close download wellPanel
          
        ),
        # Close sidebarPanel
        
        # Plot settings (main panel)
        mainPanel(
          width = 9,
          plotOutput(outputId = "results_comparison_plot", height = "600px")
        ) # Close mainPanel
        
      ) # Close sidebarLayout
    ),
    # Close tabPanel for results comparison tab
    
    ### TAB 2 : CROSS-VALIDATION PREDICTIONS ###
    tabPanel(
      title = "CV Predictions",
      value = "cv_predictions_tab",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          actionButton(
            inputId = "cv_predictions_update",
            label   = "Update Figure",
            class   = "btn-primary",
            width   = '100%'
          ),
          
          wellPanel(
            style = "background-color: #c6d2ed;",
            h4("Results selection"),
            
            tags$label(
              "Vaccine",
              tags$span(icon("info-circle"), id = "cv_predictions_vaccine_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
            ),
            selectInput(
              inputId = "cv_predictions_vaccine",
              label   = NULL,
              choices = vaccine_set,
              selected = vaccine_set[1],
              multiple = FALSE
            ),
            bsTooltip(
              id      = "cv_predictions_vaccine_info",
              title   = "Choose the vaccine for which to display variable importance",
              placement = "right",
              trigger   = "hover"
            ),
            
            tags$label(
              "Feature set",
              tags$span(icon("info-circle"), id = "cv_predictions_feature_set_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
            ),
            ## <-- DYNAMIC UI here: server will render a selectInput for available feature sets
            uiOutput("cv_predictions_feature_set_ui"),
            bsTooltip(
              id      = "cv_predictions_feature_set_info",
              title   = "Choose the feature set for which to display CV predictions",
              placement = "right",
              trigger   = "hover"
            ),
            
            tags$label(
              "Temporal feature selection approach:",
              tags$span(icon("info-circle"), id = "cv_predictions_approach_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
            ),
            selectInput(
              inputId = "cv_predictions_approach",
              label   = NULL,
              choices = approaches_set,
              selected = approaches_set[1],
              multiple = FALSE
            ),
            bsTooltip(
              id      = "cv_predictions_approach_info",
              title   = "Choose the feature selection approach.",
              placement = "right",
              trigger   = "hover"
            ),
            
            tags$label(
              "Include baseline clinical variables?",
              tags$span(icon("info-circle"), id = "cv_predictions_include_clinical_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
            ),
            selectInput(
              inputId = "cv_predictions_include_clinical",
              label   = NULL,
              choices = include_clinical_set,
              selected = include_clinical_set[1],
              multiple = FALSE
            ),
            bsTooltip(
              id      = "cv_predictions_include_clinical_info",
              title   = "Include baseline clinical variables as predictors?",
              placement = "right",
              trigger   = "hover"
            ),
            
            # Include TBA modules?
            tags$label(
              "Include TBA modules?",
              tags$span(icon("info-circle"), id    = "cv_predictions_include_tba_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
            ),
            selectInput(
              inputId = "cv_predictions_include_tba",
              label   = NULL,
              choices = include_tba_set,
              selected = include_tba_set[2],
              multiple = FALSE
            ),
            bsTooltip(
              id      = "cv_predictions_include_tba_info",
              title   = "Include modules without labels as predictors?",
              placement = "right",
              trigger   = "hover"
            ),
          
          # Highlight best predictor set?
          tags$label(
            "Highlight best predictor set?",
            tags$span(icon("info-circle"), 
                      id  = "cv_predictions_highlight_best_info",
                      style = "cursor: help; color: #337ab7; margin-left: 4px;")
          ),
          selectInput(
            inputId = "cv_predictions_highlight_best",
            label   = NULL,
            choices = highlight_best_set,
            selected = highlight_best_set[1],
            multiple = FALSE
          ),
          bsTooltip(
            id      = "cv_predictions_highlight_best_info",
            title   = "Highlight the best set of predictors?",
            placement = "right",
            trigger   = "hover"
          )
        ),
          # close wellPanel
          
          # ---------- Download section ----------
          wellPanel(
            style = "background-color: #ced4da;",
            h4("Download Options"),
            selectInput(
              inputId = "cv_predictions_download_format",
              label   = "Format of download:",
              choices = c(
                "PDF" = "pdf",
                "PNG" = "png",
                "JPEG" = "jpeg"
              ),
              selected = "pdf"
            ),
            numericInput(
              inputId = "cv_predictions_download_width",
              label   = "Width of download (inches):",
              value   = 15,
              min     = 1
            ),
            numericInput(
              inputId = "cv_predictions_download_height",
              label   = "Height of download (inches):",
              value   = 4,
              min     = 1
            ),
            numericInput(
              inputId = "cv_predictions_download_dpi",
              label   = "Resolution (DPI):",
              value   = 300,
              min     = 72
            ),
            div(
              style = "width: 100%;",
              downloadButton(
                outputId = "download_cv_predictions_plot",
                label    = "Download CV Predictions Plot",
                class    = "btn-success"
              )
            )
          ) # close download wellPanel
          
        ),
        # close sidebarPanel
        
        mainPanel(
          width = 9,
          plotOutput(outputId = "cv_predictions_plot", height = "600px")
        ) # close mainPanel
        
      ) # close sidebarLayout
    ),
    # close tabPanel
    # Close tabPanel for cv predictions tab
    
    ### TAB 3 : VARIABLE IMPORTANCE ###
    tabPanel(
      title = "Variable Importance",
      value = "variable_importance_tab",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          actionButton(
            inputId = "variable_importance_update",
            label   = "Update Figure",
            class   = "btn-primary",
            width   = '100%'
          ),
          
          wellPanel(
            style = "background-color: #c6d2ed;",
            h4("Results selection"),
            
            tags$label(
              "Vaccine",
              tags$span(icon("info-circle"), id = "variable_importance_vaccine_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
            ),
            selectInput(
              inputId = "variable_importance_vaccine",
              label   = NULL,
              choices = vaccine_set,
              selected = vaccine_set[1],
              multiple = FALSE
            ),
            bsTooltip(
              id      = "variable_importance_vaccine_info",
              title   = "Choose the vaccine for which to display variable importance",
              placement = "right",
              trigger   = "hover"
            ),
            
            tags$label(
              "Feature set",
              tags$span(icon("info-circle"), id = "variable_importance_feature_set_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
            ),
            ## <-- DYNAMIC UI here: server will render a selectInput for available feature sets
            uiOutput("variable_importance_feature_set_ui"),
            bsTooltip(
              id      = "variable_importance_feature_set_info",
              title   = "Choose the feature set for which to display CV predictions",
              placement = "right",
              trigger   = "hover"
            ),
            
            tags$label(
              "Temporal feature selection approach:",
              tags$span(icon("info-circle"), id = "variable_importance_approach_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
            ),
            selectInput(
              inputId = "variable_importance_approach",
              label   = NULL,
              choices = approaches_set,
              selected = approaches_set[1],
              multiple = FALSE
            ),
            bsTooltip(
              id      = "variable_importance_approach_info",
              title   = "Choose the feature selection approach.",
              placement = "right",
              trigger   = "hover"
            ),
            
            tags$label(
              "Include baseline clinical variables?",
              tags$span(icon("info-circle"), id = "variable_importance_include_clinical_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
            ),
            selectInput(
              inputId = "variable_importance_include_clinical",
              label   = NULL,
              choices = include_clinical_set,
              selected = include_clinical_set[1],
              multiple = FALSE
            ),
            bsTooltip(
              id      = "variable_importance_include_clinical_info",
              title   = "Include baseline clinical variables as predictors?",
              placement = "right",
              trigger   = "hover"
            ),
            
            # Include TBA modules?
            tags$label(
              "Include TBA modules?",
              tags$span(icon("info-circle"), id    = "variable_importance_include_tba_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
            ),
            selectInput(
              inputId = "variable_importance_include_tba",
              label   = NULL,
              choices = include_tba_set,
              selected = include_tba_set[2],
              multiple = FALSE
            ),
            bsTooltip(
              id      = "variable_importance_include_tba_info",
              title   = "Include modules without labels as predictors?",
              placement = "right",
              trigger   = "hover"
            ),
            # Number of variables to include
            tags$label(
              "Top N variables to include",
              tags$span(icon("info-circle"), id = "variable_importance_topN_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
            ),
            numericInput(
              inputId = "variable_importance_topN",
              label   = NULL,
              value   = 20,
              min     = 10,
              max = 100,
              step = 5
            ),
            bsTooltip(
              id      = "variable_importance_topN_info",
              title   = "Top N variables to include on plot by mean importance",
              placement = "right",
              trigger   = "hover"
            ),
            
            
          ),
          # close wellPanel
          
          # ---------- Download section ----------
          wellPanel(
            style = "background-color: #ced4da;",
            h4("Download Options"),
            selectInput(
              inputId = "variable_importance_download_format",
              label   = "Format of download:",
              choices = c(
                "PDF" = "pdf",
                "PNG" = "png",
                "JPEG" = "jpeg"
              ),
              selected = "pdf"
            ),
            numericInput(
              inputId = "variable_importance_download_width",
              label   = "Width of download (inches):",
              value   = 15,
              min     = 1
            ),
            numericInput(
              inputId = "variable_importance_download_height",
              label   = "Height of download (inches):",
              value   = 8,
              min     = 1
            ),
            numericInput(
              inputId = "variable_importance_download_dpi",
              label   = "Resolution (DPI):",
              value   = 300,
              min     = 72
            ),
            div(
              style = "width: 100%;",
              downloadButton(
                outputId = "download_variable_importance_plot",
                label    = "Download CV Predictions Plot",
                class    = "btn-success"
              )
            )
          ) # close download wellPanel
          
        ),
        # close sidebarPanel
        
        mainPanel(
          width = 9,
          plotOutput(outputId = "variable_importance_plot", height = "600px")
        ) # close mainPanel
        
      ) # close sidebarLayout
    ),
    # Close tabPanel for variable importance tab
    
  ) # Close tabsetPanel
  
) # Close fluidPage

### EXAMPLE INPUTS FOR DEBUGGING ###

input = list()
input$results_comparison_approach = "sequential"
input$results_comparison_include_clinical = "withoutClinical"
input$results_comparison_include_tba = "withoutTBA"
input$results_comparison_feature_set = feature_set_set

input$cv_predictions_approach = "sequential"
input$cv_predictions_include_clinical = "withClinical"
input$cv_predictions_vaccine = "Yellow Fever (LV)"
input$cv_predictions_feature_set = c("clinical", "Day 0", "Day 3", "Day 7", "Day 14")
input$cv_predictions_include_tba = "withoutTBA"
input$cv_predictions_highlight_best = "combined"

input$variable_importance_approach = "sequential"
input$variable_importance_include_clinical = "withClinical"
input$variable_importance_include_tba = "withTBA"
input$variable_importance_vaccine = "Yellow Fever (LV)"
input$variable_importance_feature_set = c("Day 0", "Day 3")
input$variable_importance_topN = 30

### SERVER FUNCTION ###

server <- function(input, output, session) {
  #### BEGIN HEATMAP TAB ###
  results_comparison_plot <- reactive({
    input$results_comparison_update
    isolate({
      # Select the correct dataframe
      results_name <- paste0(
        "prediction_results_all_",
        input$results_comparison_approach,
        "_",
        input$results_comparison_include_clinical,
        "_",
        input$results_comparison_include_tba
      )
      
      # Fetch the object from the environment
      selected_results <- get(results_name, envir = .GlobalEnv)
      
      # Now plot the results
      
      # --- 1. Build vaccine x set grid and extract Rspearman + sRMSE in one pass ----
      vaccines <- names(selected_results)
      
      # gather all predictor-set names that appear anywhere
      all_sets <- input$results_comparison_feature_set
      
      # expand grid of all combinations
      grid <- expand.grid(vaccine = vaccines,
                          set     = all_sets,
                          stringsAsFactors = FALSE)
      
      # pre-allocate numeric vectors for efficiency
      n <- nrow(grid)
      Rs_vec    <- rep(NA_real_, n)
      sRMSE_vec <- rep(NA_real_, n)
      
      # single-pass extraction (avoids rowwise / repeated list traversal overhead)
      for (i in seq_len(n)) {
        vac  <- grid$vaccine[i]
        setn <- grid$set[i]
        # guard against missing vaccine or set
        if (is.null(selected_results[[vac]]))
          next
        res_set <- selected_results[[vac]][[setn]]
        if (is.null(res_set) || !is.list(res_set))
          next
        if (!("metrics" %in% names(res_set)))
          next
        metrics <- res_set[["metrics"]]
        # extract if present (coerce to numeric, keep NA if missing)
        if ("Rspearman" %in% names(metrics)) {
          Rs_vec[i] <- as.numeric(metrics[["Rspearman"]])
        }
        if ("sRMSE" %in% names(metrics)) {
          sRMSE_vec[i] <- as.numeric(metrics[["sRMSE"]])
        }
      }
      
      # combine into data frame
      plot_df <- bind_cols(as_tibble(grid), tibble(Rs = Rs_vec, sRMSE = sRMSE_vec))
      
      # formatted labels: empty string for NA (so no annotation on missing tiles)
      plot_df <- plot_df %>%
        mutate(
          Rs_label    = ifelse(is.na(Rs), "", sprintf("%.2f", Rs)),
          sRMSE_label = ifelse(is.na(sRMSE), "", sprintf("%.2f", sRMSE))
        )
      
      # factors for plotting order (adjust ordering if desired)
      plot_df <- plot_df %>%
        mutate(
          set     = factor(set, levels = all_sets),
          vaccine = factor(vaccine, levels = vaccines)
        )
      
      # --- 2. Spearman R heatmap (white = 0 -> blue = 1; NA -> grey) ----
      heatmap_plot_R <- ggplot(plot_df, aes(x = set, y = vaccine, fill = Rs)) +
        geom_tile(colour = "white") +
        # annotate only non-NA values (labels are empty for NA)
        geom_text(aes(label = Rs_label),
                  colour = "black",
                  size = 5) +
        scale_fill_gradient(
          name = expression(rho),
          low = "white",
          high = "#0072B2",
          # blue
          na.value = "grey80",
          limits = c(0, 1),
          # anchor scale to 0..1
          oob = scales::squish,
          # values outside 0..1 will be squished to limits
          breaks = seq(0, 1, by = 0.25)
        ) +
        coord_fixed(ratio = 1) +   # square tiles
        labs(x = "Predictor set", y = "Vaccine") +
        theme_minimal(base_size = 20) +
        ggtitle(expression(paste("Spearman ", rho))) +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 18, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 25),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.margin = margin(
            t = 10,
            r = 2,
            b = 5,
            l = 10
          ),
          legend.position = "none"   # remove individual legend
        )
      
      heatmap_plot_sRMSE <- ggplot(plot_df, aes(x = set, y = vaccine, fill = sRMSE)) +
        geom_tile(colour = "white") +
        # annotate only non-NA values
        geom_text(aes(label = sRMSE_label),
                  colour = "black",
                  size = 5) +
        scale_fill_gradient(
          name = "sRMSE",
          low = "#2ca02c",
          # green at 0
          high = "white",
          # white at 1 and above (squished)
          na.value = "grey80",
          limits = c(0, 1),
          oob = scales::squish,
          breaks = seq(0, 1, by = 0.25)
        ) +
        coord_fixed(ratio = 1) +   # square tiles
        labs(x = "Predictor set", y = "Vaccine") +
        theme_minimal(base_size = 20) +
        ggtitle("Standardised RMSE") +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 18, hjust = 0.5),
          axis.title.y = element_text(size = 25),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.margin = margin(
            t = 5,
            r = 2,
            b = 10,
            l = 10
          ),
          legend.position = "none"   # remove individual legend
        )
      
      
      plot_subtitle = paste0(
        ifelse(
          input$results_comparison_approach == "sequential",
          "Sequential",
          "Cumulative"
        ),
        " prediction approach, clinical variables ",
        ifelse(
          input$results_comparison_include_clinical == "withClinical",
          "included",
          "not included"
        )
      )
      
      
      # --- Combine with shared legend and common title ---
      combined <- heatmap_plot_R / heatmap_plot_sRMSE +
        plot_layout(
          ncol = 1,
          heights = c(1, 1),
          guides = "collect"   # collect legends into one
        ) +
        plot_annotation(
          title = "Evaluation metrics of CV predictions",
          subtitle = plot_subtitle,
          theme = theme(
            plot.title = element_text(
              size = 26,
              face = "bold",
              hjust = 0.5
            ),
            plot.subtitle = element_text(
              size = 20,
              face = "bold",
              hjust = 0.5
            )
          )
        ) & theme(
          legend.position = "right",
          legend.box.margin = margin(0, 0, 0, 0),
          legend.key.width = unit(1, "cm")  # optional: reduce legend width
        )
      
      print(combined)
      
    }) # End Isolate
    
  }) # End Reactive
  
  output$results_comparison_plot = renderPlot({
    results_comparison_plot()
  }) # End renderPlot
  
  output$download_results_comparison_plot <- downloadHandler(
    filename = function() {
      req(input$results_comparison_download_format)
      paste0(
        "results_comparison_",
        input$results_comparison_approach,
        "_",
        input$results_comparison_include_clinical,
        "_",
        input$results_comparison_include_tba,
        ".",
        input$results_comparison_download_format
      )
    },
    content = function(file) {
      # make sure the plot reactive is available
      req(results_comparison_plot())
      
      # call ggsave with user inputs
      ggsave(
        filename = file,
        plot     = results_comparison_plot(),
        device   = input$results_comparison_download_format,
        width    = input$results_comparison_download_width,
        height   = input$results_comparison_download_height,
        dpi      = input$results_comparison_download_dpi,
        units    = "in"
      )
    }
  )
  ### END RESULTS COMPARISON TAB ###
  
  ### BEGIN CV PREDICTIONS TAB ###
  
  # 1) Dynamic UI for feature sets (timepoints) depending on vaccine / approach / include_clinical
  output$cv_predictions_feature_set_ui <- renderUI({
    req(
      input$cv_predictions_vaccine,
      input$cv_predictions_approach,
      input$cv_predictions_include_clinical,
      input$cv_predictions_include_tba
    )
    
    approach_key <- tolower(gsub("\\s+", "", input$cv_predictions_approach))
    
    include_label <- tolower(input$cv_predictions_include_clinical)
    include_suffix <- if (grepl("withClinical", include_label, ignore.case = TRUE)) {
      "withClinical"
    } else {
      "withoutClinical"
    }
    
    include_label_tba <- tolower(input$cv_predictions_include_tba)
    include_suffix_tba <- if (grepl("withTBA", include_label_tba, ignore.case = TRUE)) {
      "withTBA"
    } else {
      "withoutTBA"
    }
    
    results_name <- paste0(
      "prediction_results_all_",
      approach_key,
      "_",
      include_suffix,
      "_",
      include_suffix_tba
    )
    
    results_obj <- get0(results_name, envir = .GlobalEnv)
    
    choices <- character(0)
    if (!is.null(results_obj) &&
        input$cv_predictions_vaccine %in% names(results_obj)) {
      vaccine_entry <- results_obj[[input$cv_predictions_vaccine]]
      if (is.list(vaccine_entry) && length(vaccine_entry) > 0) {
        choices <- names(vaccine_entry)
      }
    }
    
    if (length(choices) == 0) {
      selectInput(
        inputId = "cv_predictions_feature_set",
        label = NULL,
        choices = c("No feature sets available" = ""),
        selected = "",
        multiple = TRUE
      )
    } else {
      previous <- isolate(input$cv_predictions_feature_set)
      
      # Keep only previous selections that are still present in choices,
      # preserving the user's selection order.
      valid_prev <- intersect(as.character(previous), choices)
      
      # If nothing valid remains, fall back to the first choice (or NULL to select none)
      selected_choice <- if (length(valid_prev) > 0)
        valid_prev
      else
        choices[1]
      
      selectInput(
        inputId = "cv_predictions_feature_set",
        label = NULL,
        choices = choices,
        selected = selected_choice,
        multiple = TRUE
      )
    }
  })
  
  # 2) Reactive that builds (and returns) the ggplot object for the selected options
  cv_predictions_plot <- reactive({
    # This line makes the reactive re-run when the user clicks the Update Figure button
    input$cv_predictions_update
    
    isolate({
      # build the object name the same way we do in renderUI
      approach_key <- tolower(gsub("\\s+", "", input$cv_predictions_approach))
      include_label <- tolower(input$cv_predictions_include_clinical)
      include_suffix <- if (grepl("withClinical", include_label, ignore.case = TRUE)) {
        "withClinical"
      } else {
        "withoutClinical"
      }
      include_label_tba <- tolower(input$cv_predictions_include_tba)
      include_suffix_tba <- if (grepl("withTBA", include_label_tba, ignore.case = TRUE)) {
        "withTBA"
      } else {
        "withoutTBA"
      }
      results_name <- paste0(
        "prediction_results_all_",
        approach_key,
        "_",
        include_suffix,
        "_",
        include_suffix_tba
      )
      
      # Safely fetch the results object
      selected_results <- get0(results_name, envir = .GlobalEnv)
      
      if (length(input$cv_predictions_feature_set) == 1) {
        # SINGLE FEATURE SET SELECTION CASE
        selected_results_vaccine <- selected_results[[input$cv_predictions_vaccine]]
        
        selected_results_vaccine_feature_set <- selected_results_vaccine[[input$cv_predictions_feature_set]]
        
        # Try to extract the plot object; expect it under [["plots"]][["cv_results"]]
        plot_obj <- NULL
        if (is.list(selected_results_vaccine_feature_set) &&
            !is.null(selected_results_vaccine_feature_set[["plots"]]) &&
            !is.null(selected_results_vaccine_feature_set[["plots"]][["cv_results"]])) {
          plot_obj <- selected_results_vaccine_feature_set[["plots"]][["cv_results"]]
        }
        
        validate(
          need(
            !is.null(plot_obj),
            "Plot not found for the selected vaccine / feature set / approach."
          )
        )
        
        # replace the letter "R" with the greek letter rho
        for (i in seq_along(plot_obj$layers)) {
          lyr <- plot_obj$layers[[i]]
          
          # Static annotation, e.g. annotate("text", label = "R = 0.85")
          if (!is.null(lyr$aes_params$label) &&
              is.character(lyr$aes_params$label)) {
            lbl <- lyr$aes_params$label
            if (grepl("\\bR\\b", lbl)) {
              if (grepl("=", lbl)) {
                val <- trimws(sub(".*=", "", lbl))
                plot_obj$layers[[i]]$aes_params$label <- bquote(rho == .(val))
              } else {
                plot_obj$layers[[i]]$aes_params$label <- expression(rho)
              }
            }
          }
        }
        
        print(plot_obj)
        
        # END SINGLE FEATURE SET SELECTION CASE
      } else if (length(input$cv_predictions_feature_set) > 1) {
        # MULTIPLE FEATURE SETS CASE
        
        selected_results_vaccine <- selected_results[[input$cv_predictions_vaccine]]
        # Extract metrics for selected timepoints
        selected_results_vaccine_feature_set_filtered = selected_results_vaccine[input$cv_predictions_feature_set]
        # Compute the combined metric for each set
        
        if (input$cv_predictions_highlight_best == "none"){
          best_set = "none"
        } else if (input$cv_predictions_highlight_best == "spearman"){
          set_metrics <- sapply(input$cv_predictions_feature_set, function(set) {
            res <- selected_results_vaccine_feature_set_filtered[[set]]
            1 - res[["metrics"]][["Rspearman"]]
          })
          
          # Identify the set with the minimum value
          best_set <- names(which.min(set_metrics))
        } else if (input$cv_predictions_highlight_best == "sRMSE"){
          set_metrics <- sapply(input$cv_predictions_feature_set, function(set) {
            res <- selected_results_vaccine_feature_set_filtered[[set]]
            res[["metrics"]][["sRMSE"]]
          })
          
          # Identify the set with the minimum value
          best_set <- names(which.min(set_metrics))
        } else if (input$cv_predictions_highlight_best == "combined"){
          set_metrics <- sapply(input$cv_predictions_feature_set, function(set) {
            res <- selected_results_vaccine_feature_set_filtered[[set]]
            res[["metrics"]][["sRMSE"]] + (1 - res[["metrics"]][["Rspearman"]])
          })
          
          # Identify the set with the minimum value
          best_set <- names(which.min(set_metrics))
        }
        
        #jump
        day_sets <- input$cv_predictions_feature_set[grepl("^Day\\s*[0-9]+$", input$cv_predictions_feature_set)]
        day_nums <- as.numeric(sub("Day\\s*", "", day_sets))
        
        if ("clinical" %in% input$cv_predictions_feature_set) {
          sorted_feature_sets <- c("clinical", paste0("Day ", sort(day_nums)))
        } else {
          sorted_feature_sets <- c(paste0("Day ", sort(day_nums)))
        }
        
        # --- annotation tuning (change these if you want smaller/larger annotations) ---
        ann_size  <- 4   # geom_text/geom_label size (ggplot2 units) - smaller number => smaller text
        ann_alpha <- 0.9   # transparency for annotation text
        
        # --- blank placeholder (pure white, no border) --------------------------------
        blank_placeholder <- function() {
          ggplot() + geom_blank() + theme_void() +
            theme(
              plot.background = element_rect(fill = "white", colour = NA),
              plot.margin = margin(6, 6, 6, 6)
            )
        }
        
        # --- build grobs: top row headers, bottom row plots ----------------------------
        n_col <- length(sorted_feature_sets)
        
        # We'll create 2 rows (headers + plots) and n_col columns
        total_rows <- 2
        total_cols <- n_col
        
        grob_list <- vector("list", total_rows * total_cols)
        
        # top row: feature set headers
        for (j in seq_len(n_col)) {
          h <- textGrob(
            label = paste0(sorted_feature_sets[j], " (n = ",
                           nrow(selected_results_vaccine_feature_set_filtered[[j]][["cv_results"]]) 
                           ,")"),
            gp = gpar(fontsize = 16, fontface = "bold"),
            just = "center"
          )
          grob_list[[j]] <- h
        }
        
        # bottom row: the plots (with panel border inside each plot); blanks if missing
        for (j in seq_len(n_col)) {
          feat <- sorted_feature_sets[j]
          
          is_best <- identical(feat, best_set)
          
          got_plot <- NULL
          
          if (!is.null(selected_results[[input$cv_predictions_vaccine]]) &&
              !is.null(selected_results[[input$cv_predictions_vaccine]][[feat]])) {
            inner <- selected_results[[input$cv_predictions_vaccine]][[feat]]
            if (!is.null(inner[["plots"]]) &&
                !is.null(inner[["plots"]][["cv_results"]])) {
              got_plot <- inner[["plots"]][["cv_results"]]
            }
          }
          
          if (!is.null(got_plot) && inherits(got_plot, "ggplot")) {
            # Keep metric annotations but reduce their size/alpha
            p_mod <- got_plot + labs(title = NULL, subtitle = NULL)
            
            # shrink any geom_text / geom_label layers
            if (length(p_mod$layers) > 0) {
              for (k in seq_along(p_mod$layers)) {
                if (inherits(p_mod$layers[[k]]$geom, "GeomText") ||
                    inherits(p_mod$layers[[k]]$geom, "GeomLabel")) {
                  # set fixed text size & alpha
                  p_mod$layers[[k]]$aes_params$size  <- ann_size
                  p_mod$layers[[k]]$geom_params$size <- ann_size
                  p_mod$layers[[k]]$aes_params$alpha  <- ann_alpha
                  p_mod$layers[[k]]$geom_params$alpha <- ann_alpha
                  
                  # for labels, reduce padding/rounding so boxes are smaller
                  if (inherits(p_mod$layers[[k]]$geom, "GeomLabel")) {
                    p_mod$layers[[k]]$geom_params$label.padding <- grid::unit(0.05, "lines")
                    p_mod$layers[[k]]$geom_params$label.r <- grid::unit(0.05, "lines")
                  }
                }
              }
            }
            
            # choose border appearance depending on whether this is the best set
            panel_border_col <-  "grey70"
            panel_border_size <- 0.8
            plot_bg_col <- if (is_best) "#5DC53B" else NA
            plot_bg_size <- if (is_best) 1.5 else 0
            
            # apply compact theme and panel border
            p_mod <- p_mod +
              theme(
                plot.title = element_blank(),
                plot.subtitle = element_blank(),
                plot.margin = margin(6, 6, 6, 6),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(size = 7),
                axis.text.y = element_text(size = 7),
                panel.border = element_rect(
                  colour = panel_border_col,
                  fill = NA,
                  size = panel_border_size
                ),
                panel.background = element_rect(fill = "white", colour = NA),
                plot.background = element_rect(
                  fill = "white",
                  colour = plot_bg_col,
                  size = plot_bg_size
                )
              )
            
            grob_list[[n_col + j]] <- ggplotGrob(p_mod) # second row index = n_col + j
          } else {
            # missing -> blank white placeholder (no border)
            grob_list[[n_col + j]] <- ggplotGrob(blank_placeholder())
          }
        }
        
        # --- layout matrix ------------------------------------------------------------
        layout_mat <- matrix(
          seq_len(total_rows * total_cols),
          nrow = total_rows,
          ncol = total_cols,
          byrow = TRUE
        )
        
        # sizing: tweak these if the result looks cramped
        col_widths <- rep(4, n_col)
        row_heights <- c(1, 5)  # 1st row small for headers, second row big for plots
        
        # Create the arranged grob
        final_grob <- arrangeGrob(
          grobs = grob_list,
          layout_matrix = layout_mat,
          widths = col_widths,
          heights = row_heights,
          top = textGrob(
            paste0(
              "Cross validation predictions - ",
              input$cv_predictions_vaccine
            ),
            gp = gpar(fontsize = 18, fontface = "bold")
          )
        )
        
        final_grob_with_axes <- gridExtra::grid.arrange(
          final_grob,
          left = textGrob(
            "Predicted",
            rot = 90,
            gp = gpar(fontsize = 20)
          ),
          bottom = textGrob("Observed", gp = gpar(fontsize = 20))
        )
        
        return(final_grob_with_axes)
        
      }
      
    }) # end isolate
  }) # end reactive
  
  # 3) Render the plot in the UI
  output$cv_predictions_plot <- renderPlot({
    p <- cv_predictions_plot()
    req(p)
    
    if (inherits(p, "ggplot")) {
      print(p)
    } else if (inherits(p, c("gtable", "grob"))) {
      grid::grid.newpage()
      grid::grid.draw(p)
    } else {
      tryCatch({
        print(p)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, "Unable to render the selected plot object.")
      })
    }
  })
  
  # 4) Download handler (unchanged except using get0-safe reactive above)
  output$download_cv_predictions_plot <- downloadHandler(
    filename = function() {
      fmt <- req(input$cv_predictions_download_format)
      vaccine_clean <- gsub("[()]", "", gsub(" ", "", input$cv_predictions_vaccine))
      
      day_sets <- input$cv_predictions_feature_set[grepl("^Day\\s*[0-9]+$", input$cv_predictions_feature_set)]
      day_nums <- as.numeric(sub("Day\\s*", "", day_sets))
      
      if ("clinical" %in% input$cv_predictions_feature_set) {
        sorted_feature_sets <- c("clinical", paste0("Day ", sort(day_nums)))
      } else {
        sorted_feature_sets <- c(paste0("Day ", sort(day_nums)))
      }
      
      feature_sets_collapsed <- paste(sorted_feature_sets, collapse = "-")
      feature_clean <- gsub("[()]", "", gsub(" ", "", feature_sets_collapsed))
      
      paste0(
        "cv_predictions_",
        input$cv_predictions_approach,
        "_",
        input$cv_predictions_include_clinical,
        "_",
        input$results_comparison_include_tba,
        "_",
        vaccine_clean,
        "_",
        feature_clean,
        ".",
        fmt
      )
    },
    content = function(file) {
      p <- req(cv_predictions_plot())
      fmt <- tolower(req(input$cv_predictions_download_format))
      width_in  <- req(input$cv_predictions_download_width)
      height_in <- req(input$cv_predictions_download_height)
      dpi <- req(input$cv_predictions_download_dpi)
      
      # If ggplot, use ggsave
      if (inherits(p, "ggplot")) {
        ggplot2::ggsave(
          filename = file,
          plot     = p,
          device   = fmt,
          width    = width_in,
          height   = height_in,
          dpi      = dpi,
          units    = "in"
        )
        return()
      }
      
      # If grob/gtable: draw to device manually
      if (fmt %in% c("pdf", "svg")) {
        if (fmt == "pdf")
          grDevices::pdf(file = file,
                         width = width_in,
                         height = height_in)
        if (fmt == "svg")
          grDevices::svg(file = file,
                         width = width_in,
                         height = height_in)
        grid::grid.draw(p)
        grDevices::dev.off()
        return()
      }
      
      # raster formats: convert to pixels
      px_width  <- ceiling(width_in * dpi)
      px_height <- ceiling(height_in * dpi)
      
      if (fmt %in% c("png")) {
        grDevices::png(
          filename = file,
          width = px_width,
          height = px_height,
          res = dpi,
          units = "px"
        )
        grid::grid.draw(p)
        grDevices::dev.off()
        return()
      } else if (fmt %in% c("jpg", "jpeg")) {
        grDevices::jpeg(
          filename = file,
          width = px_width,
          height = px_height,
          res = dpi,
          units = "px"
        )
        grid::grid.draw(p)
        grDevices::dev.off()
        return()
      } else if (fmt %in% c("tiff")) {
        grDevices::tiff(
          filename = file,
          width = px_width,
          height = px_height,
          res = dpi,
          units = "px"
        )
        grid::grid.draw(p)
        grDevices::dev.off()
        return()
      } else {
        stop("Unsupported format for saving grob: ", fmt)
      }
    }
  )
  
  ### END CV PREDICTIONS TAB ###
  
  ### BEGIN VARIABLE IMPORTANCE TAB ###
  
  # 1) Dynamic UI for feature sets (timepoints) depending on vaccine / approach / include_clinical
  output$variable_importance_feature_set_ui <- renderUI({
    req(
      input$variable_importance_vaccine,
      input$variable_importance_approach,
      input$variable_importance_include_clinical,
      input$variable_importance_include_tba
    )
    
    approach_key <- tolower(gsub("\\s+", "", input$variable_importance_approach))
    
    include_label <- tolower(input$variable_importance_include_clinical)
    include_suffix <- if (grepl("withClinical", include_label, ignore.case = TRUE)) {
      "withClinical"
    } else {
      "withoutClinical"
    }
    include_label_tba <- tolower(input$variable_importance_include_tba)
    include_suffix_tba <- if (grepl("withTBA", include_label_tba, ignore.case = TRUE)) {
      "withTBA"
    } else {
      "withoutTBA"
    }
    
    results_name <- paste0(
      "prediction_results_all_",
      approach_key,
      "_",
      include_suffix,
      "_",
      include_suffix_tba
    )
    
    results_obj <- get0(results_name, envir = .GlobalEnv)
    
    choices <- character(0)
    if (!is.null(results_obj) &&
        input$variable_importance_vaccine %in% names(results_obj)) {
      vaccine_entry <- results_obj[[input$variable_importance_vaccine]]
      if (is.list(vaccine_entry) && length(vaccine_entry) > 0) {
        choices <- names(vaccine_entry)
      }
    }
    
    if (length(choices) == 0) {
      selectInput(
        inputId = "variable_importance_feature_set",
        label = NULL,
        choices = c("No feature sets available" = ""),
        selected = "",
        multiple = TRUE
      )
    } else {
      previous <- isolate(input$variable_importance_feature_set)
      
      # keep any previous selections that remain valid; otherwise pick the first choice
      valid_prev <- intersect(as.character(previous), choices)
      selected_choice <- if (length(valid_prev) > 0)
        valid_prev
      else
        choices[1]
      
      selectInput(
        inputId = "variable_importance_feature_set",
        label = NULL,
        choices = choices,
        selected = selected_choice,
        multiple = TRUE
      )
    }
  })
  
  # 2) Reactive that builds (and returns) the ggplot object for the selected options
  variable_importance_plot <- reactive({
    # This line makes the reactive re-run when the user clicks the Update Figure button
    input$variable_importance_update
    
    isolate({
      # build the object name the same way we do in renderUI
      approach_key <- tolower(gsub("\\s+", "", input$variable_importance_approach))
      include_label <- tolower(input$variable_importance_include_clinical)
      include_suffix <- if (grepl("withClinical", include_label, ignore.case = TRUE)) {
        "withClinical"
      } else {
        "withoutClinical"
      }
      include_label_tba <- tolower(input$variable_importance_include_tba)
      include_suffix_tba <- if (grepl("withTBA", include_label_tba, ignore.case = TRUE)) {
        "withTBA"
      } else {
        "withoutTBA"
      }
      results_name <- paste0(
        "prediction_results_all_",
        approach_key,
        "_",
        include_suffix,
        "_",
        include_suffix_tba
      )
      
      # Safely fetch the results object
      selected_results <- get0(results_name, envir = .GlobalEnv)
      
      selected_results_vaccine <- selected_results[[input$variable_importance_vaccine]]
      
      # IF ONLY ONE TIMEPOINT SELECTED #
      if (length(input$variable_importance_feature_set) == 1) {
        selected_results_vaccine_feature_set <- selected_results_vaccine[[input$variable_importance_feature_set]]
        
        # Try to extract the plot object; expect it under [["var_imp"]]
        
        metrics = selected_results_vaccine_feature_set[["metrics"]]
        
        plot_obj = selected_results_vaccine_feature_set[["var_imp"]]
        
        validate(
          need(
            !is.null(plot_obj),
            "Plot not found for the selected vaccine / feature set / approach."
          )
        )
        
        # plot topN mean standardised feature importances with coloured errorbars and points
        
        # Prepare plotting data
        plot_df <- plot_obj %>%
          slice_min(order_by = mean_rank,
                    n = input$variable_importance_topN) %>%
          mutate(
            mean_imp = as.numeric(mean_imp),
            sd_imp   = as.numeric(sd_imp),
            feature  = forcats::fct_reorder(feature, mean_imp),
            xmin = pmax(0, mean_imp - sd_imp),
            xmax = pmin(100, mean_imp + sd_imp)
          )
        
        colour_df = plot_obj %>%
          dplyr::select(feature_group, feature_colour) %>%
          distinct()
        
        pred_cols = colour_df$feature_colour
        
        names(pred_cols) = colour_df$feature_group
        
        plot_title = paste0(
          "Variable Importance : " ,
          input$variable_importance_approach,
          " prediction approach, clinical variables ",
          ifelse(
            input$variable_importance_include_clinical == "withClinical",
            "included",
            "not included"
          )
        )
        
        # Plot average standardised variable importance
        vi_plot <- plot_df %>%
          ggplot(aes(x = mean_imp, y = feature)) +
          geom_errorbar(
            aes(
              xmin = xmin,
              xmax = xmax,
              colour = feature_group
            ),
            width = 0.25,
            size = 1,
            orientation = "y"
          ) +
          geom_point(aes(colour = feature_group), size = 3) +
          scale_colour_manual(
            name = "Feature set",
            values = pred_cols,
            na.value = "#777777",
            guide = guide_legend(
              override.aes = list(
                linetype = 1,
                shape = 16,
                size = 5
              ),
              keywidth  = unit(2.2, "cm"),
              keyheight = unit(0.9, "cm"),
              ncol = 1
            )
          ) +
          scale_x_continuous(limits = c(0, 100),
                             expand = expansion(mult = c(0.02, 0.12))) +
          labs(
            x = "Mean standardised training-set variable importance",
            y = NULL,
            title = plot_title,
            subtitle = bquote("sRMSE =" ~ .(sprintf(
              "%.2f", metrics$sRMSE
            )) * ", " ~ rho == .(
              sprintf("%.2f", metrics$Rspearman)
            ))
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(
              size = 20,
              face = "bold",
              hjust = 0.5
            ),
            plot.subtitle = element_text(size = 20, hjust = 0.5),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.y = element_blank(),
            legend.title = element_text(size = 16),
            legend.text  = element_text(size = 14),
            legend.key.size = unit(1.1, "cm")
          )
        
        print(vi_plot) # END IF CLAUSE FOR SINGULAR FEATURE SET SELECTION
      } else if (length(input$variable_importance_feature_set) > 1) {
        # BEGIN IF CLAUSE FOR MULTIPLE FEATURE SET SELECTIONS
        
        
        # --- blank placeholder (pure white, no border) --------------------------------
        blank_placeholder <- function() {
          ggplot() + geom_blank() + theme_void() +
            theme(
              plot.background = element_rect(fill = "white", colour = NA),
              plot.margin = margin(6, 6, 6, 6)
            )
        }
        
        # --- helper function to build a single VI plot ---------------------------------
        make_vi_plot <- function(vac,
                                 feature_set,
                                 topN,
                                 wrap_width = 20) {
          selected_results_vaccine_feature_set <- selected_results_vaccine[[feature_set]]
          
          if (is.null(selected_results_vaccine_feature_set))
            return(blank_placeholder())
          
          metrics <- selected_results_vaccine_feature_set[["metrics"]]
          plot_obj <- selected_results_vaccine_feature_set[["var_imp"]]
          
          # Prepare topN data and wrap feature names
          plot_df <- plot_obj %>%
            slice_min(order_by = mean_rank, n = topN) %>%
            mutate(
              mean_imp = as.numeric(mean_imp),
              sd_imp   = as.numeric(sd_imp),
              feature  = forcats::fct_reorder(
                stringr::str_wrap(feature, width = wrap_width),
                mean_imp
              ),
              xmin = pmax(0, mean_imp - sd_imp),
              xmax = pmin(100, mean_imp + sd_imp)
            )
          
          # colours
          colour_df <- plot_obj %>%
            dplyr::select(feature_group, feature_colour) %>%
            distinct()
          pred_cols <- colour_df$feature_colour
          names(pred_cols) <- colour_df$feature_group
          
          # Build plot
          p <- ggplot(plot_df, aes(x = feature, y = mean_imp)) +
            geom_errorbar(
              aes(
                ymin = xmin,
                ymax = xmax,
                colour = feature_group
              ),
              width = 0.25,
              size = 1,
              orientation = "x"
            ) +
            geom_point(aes(colour = feature_group), size = 2) +
            scale_colour_manual(values = pred_cols, na.value = "#777777") +
            scale_y_continuous(limits = c(0, 100),
                               expand = expansion(mult = c(0.02, 0.12))) +
            labs(
              x = NULL,
              y = "Mean standardized VI",
              title = feature_set,
              subtitle = bquote(sRMSE == .(
                sprintf("%.2f", metrics$sRMSE)
              ) ~ ", " ~ rho == .(
                sprintf("%.2f", metrics$Rspearman)
              ))
            ) +
            coord_flip() +  # flip coordinates
            theme_minimal() +
            theme(
              plot.title = element_text(
                size = 12,
                face = "bold",
                hjust = 0.5
              ),
              plot.subtitle = element_text(size = 10, hjust = 0.5),
              axis.text = element_text(size = 7),
              axis.title = element_text(size = 10),
              panel.grid.major.y = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none",
              panel.border = element_rect(
                colour = "grey70",
                fill = NA,
                size = 0.8
              ),
              plot.margin = margin(4, 4, 4, 4)
            )
          
          return(p)
        }
        
        
        # --- build grobs for all feature sets ----------------------------------------
        # First sort the feature sets
        day_sets <- input$variable_importance_feature_set[grepl("^Day\\s*[0-9]+$",
                                                                input$variable_importance_feature_set)]
        day_nums <- as.numeric(sub("Day\\s*", "", day_sets))
        
        if ("clinical" %in% input$variable_importance_feature_set) {
          sorted_feature_sets <- c("clinical", paste0("Day ", sort(day_nums)))
        } else {
          sorted_feature_sets <- c(paste0("Day ", sort(day_nums)))
        }
        
        grob_list <- lapply(sorted_feature_sets, function(fs)
          make_vi_plot(
            input$variable_importance_vaccine,
            fs,
            input$variable_importance_topN
          ))
        
        # Assign the arranged grob to a variable
        final_vi_grob <- arrangeGrob(
          grobs = grob_list,
          nrow = 1,
          top = textGrob(
            paste0(
              "Variable Importance - ",
              input$variable_importance_vaccine
            ),
            gp = gpar(fontsize = 18, fontface = "bold")
          )
        )
        
        return(final_vi_grob)
        
      }
      
    }) # end isolate
  }) # end reactive
  
  output$variable_importance_plot <- renderPlot({
    p <- variable_importance_plot()
    req(p) # ensures we have something to draw
    
    if (inherits(p, "ggplot")) {
      print(p)
    } else if (inherits(p, c("gtable", "grob"))) {
      grid::grid.newpage()
      grid::grid.draw(p)
    } else {
      # best-effort fallback
      tryCatch({
        print(p)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, "Unable to render the selected plot object.")
      })
    }
  })
  
  # 4) Download handler (unchanged except using get0-safe reactive above)
  output$download_variable_importance_plot <- downloadHandler(
    filename = function() {
      fmt <- req(input$variable_importance_download_format)
      vaccine_clean <- gsub("[()]",
                            "",
                            gsub(" ", "", input$variable_importance_vaccine))
      feature_sets_collapsed <- paste(input$variable_importance_feature_set, collapse = "-")
      feature_clean <- gsub("[()]", "", gsub(" ", "", feature_sets_collapsed))
      
      paste0(
        "variable_importance_",
        input$variable_importance_approach,
        "_",
        input$variable_importance_include_clinical,
        "_",
        input$results_comparison_include_tba,
        "_",
        vaccine_clean,
        "_",
        feature_clean,
        "_",
        input$variable_importance_topN,
        "variables.",
        fmt
      )
    },
    content = function(file) {
      p <- req(variable_importance_plot())  # get the returned plot/grob
      fmt <- tolower(req(input$variable_importance_download_format))
      width_in  <- req(input$variable_importance_download_width)
      height_in <- req(input$variable_importance_download_height)
      dpi <- req(input$variable_importance_download_dpi)
      
      # If ggplot, use ggsave (handles many formats)
      if (inherits(p, "ggplot")) {
        ggplot2::ggsave(
          filename = file,
          plot     = p,
          device   = fmt,
          width    = width_in,
          height   = height_in,
          dpi      = dpi,
          units    = "in"
        )
        return()
      }
      
      # If grob/gtable, draw to a graphics device manually
      if (fmt %in% c("pdf", "svg")) {
        if (fmt == "pdf")
          grDevices::pdf(file = file,
                         width = width_in,
                         height = height_in)
        if (fmt == "svg")
          grDevices::svg(file = file,
                         width = width_in,
                         height = height_in)
        grid::grid.draw(p)
        grDevices::dev.off()
        return()
      }
      
      # raster formats - width/height must be in pixels
      px_width  <- ceiling(width_in * dpi)
      px_height <- ceiling(height_in * dpi)
      
      if (fmt %in% c("png")) {
        grDevices::png(
          filename = file,
          width = px_width,
          height = px_height,
          res = dpi,
          units = "px"
        )
        grid::grid.draw(p)
        grDevices::dev.off()
        return()
      } else if (fmt %in% c("jpg", "jpeg")) {
        grDevices::jpeg(
          filename = file,
          width = px_width,
          height = px_height,
          res = dpi,
          units = "px"
        )
        grid::grid.draw(p)
        grDevices::dev.off()
        return()
      } else if (fmt %in% c("tiff")) {
        grDevices::tiff(
          filename = file,
          width = px_width,
          height = px_height,
          res = dpi,
          units = "px"
        )
        grid::grid.draw(p)
        grDevices::dev.off()
        return()
      } else {
        stop("Unsupported format for saving grob: ", fmt)
      }
    }
  )
  
  ### END VARIABLE IMPORTANCE TAB ###
  
  
} # close server function

# Run the shiny application
shinyApp(ui, server)
