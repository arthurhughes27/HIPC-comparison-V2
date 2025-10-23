suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(fs)
  library(shiny)
  library(shinyBS)
  library(patchwork)
})


# Load in all the possible results that shiny may need access to

# Cumulative approach, with clinical variables
p_prediction_results_all_cumulative_withClinical = fs::path("output",
                                                            "results",
                                                            "prediction_results_all_cumulative_withClinical.rds")

prediction_results_all_cumulative_withClinical = readRDS(p_prediction_results_all_cumulative_withClinical)

# Cumulative approach, without clinical variables
p_prediction_results_all_cumulative_withoutClinical = fs::path("output",
                                                               "results",
                                                               "prediction_results_all_cumulative_withoutClinical.rds")

prediction_results_all_cumulative_withoutClinical = readRDS(p_prediction_results_all_cumulative_withoutClinical)

# Sequential approach, with clinical variables
p_prediction_results_all_sequential_withClinical = fs::path("output",
                                                            "results",
                                                            "prediction_results_all_sequential_withClinical.rds")

prediction_results_all_sequential_withClinical = readRDS(p_prediction_results_all_sequential_withClinical)

# Sequential approach, without clinical variables
p_prediction_results_all_sequential_withoutClinical = fs::path("output",
                                                               "results",
                                                               "prediction_results_all_sequential_withoutClinical.rds")

prediction_results_all_sequential_withoutClinical = readRDS(p_prediction_results_all_sequential_withoutClinical)

### DEFINE USER PARAMETERS ###

### RESULTS COMPARISON TAB PARAMETERS ###

approaches_set = c("Sequential" = "sequential", "Cumulative" = "cumulative")

include_clinical_set = c(
  "Include clinical variables" = "withClinical",
  "Do not include clinical variables" = "withoutClinical"
)

### CV PREDICTIONS TAB ###

vaccine_set = names(prediction_results_all_sequential_withoutClinical)

feature_set_set = unique(unlist(
  lapply(prediction_results_all_sequential_withClinical, names)
))

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
              value   = 35,
              min     = 1
            ),
            numericInput(
              inputId = "results_comparison_download_height",
              label   = "Height of download (inches):",
              value   = 15,
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
        mainPanel(width = 9, div(
          style = "text-align: center;",
          plotOutput(
            outputId = "results_comparison_plot",
            height = "1000px",
            width  = "160%"
          )
        )) # Close mainPanel
        
      ) # Close sidebarLayout
    ),
    # Close tabPanel for results comparison tab
    
    ### TAB 2 : CROSS-VALIDATION PREDICTIONS ###
    tabPanel(
      title = "CV Predictions",
      value = "cv_predictions_tab",
      sidebarLayout(
        sidebarPanel(
          # Specify width of panel
          width = 3,
          
          # Update plot button
          actionButton(
            inputId = "cv_predictions_update",
            label   = "Update Figure",
            class   = "btn-primary",
            width   = '100%'
          ),
          # close actionButton
          
          # Data selection group
          wellPanel(
            style = "background-color: #c6d2ed;",
            h4("Results selection"),
            
            # Which vaccine?
            tags$label(
              "Vaccine",
              tags$span(icon("info-circle"), id    = "cv_predictions_vaccine_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
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
            
            # Which feature_set?
            tags$label(
              "Feature set",
              tags$span(icon("info-circle"), id    = "cv_predictions_feature_set_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
            ),
            selectInput(
              inputId = "cv_predictions_feature_set",
              label   = NULL,
              choices = feature_set_set,
              selected = feature_set_set[2],
              multiple = FALSE
            ),
            bsTooltip(
              id      = "cv_predictions_feature_set_info",
              title   = "Choose the feature set for which to display CV predictions",
              placement = "right",
              trigger   = "hover"
            ),
            
            # Results from which approach?
            tags$label(
              "Temporal feature selection approach:",
              tags$span(icon("info-circle"), id    = "cv_predictions_approach_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
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
            
            # Include clinical variables?
            tags$label(
              "Include baseline clinical variables?",
              tags$span(icon("info-circle"), id    = "cv_predictions_include_clinical_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
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
            )
          ),
          # Close first wellPanel (data selection)
          
          # ---------- Download section ----------
          wellPanel(
            style = "background-color: #ced4da;",
            h4("Download Options"),
            # choose format
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
            
            # size controls (in inches)
            numericInput(
              inputId = "cv_predictions_download_width",
              label   = "Width of download (inches):",
              value   = 35,
              min     = 1
            ),
            numericInput(
              inputId = "cv_predictions_download_height",
              label   = "Height of download (inches):",
              value   = 15,
              min     = 1
            ),
            
            numericInput(
              inputId = "cv_predictions_download_dpi",
              label   = "Resolution (DPI):",
              value   = 300,
              min     = 72
            ),
            
            ## Download plot button
            div(
              style = "width: 100%;",
              downloadButton(
                outputId = "download_cv_predictions_plot",
                label    = "Download CV Predictions Plot",
                class    = "btn-success"
              )
            )
          ) # Close download wellPanel
          
        ),
        # Close sidebarPanel
        
        # Plot settings (main panel)
        mainPanel(width = 7, div(
          style = "text-align: center;",
          plotOutput(
            outputId = "cv_predictions_plot",
            height = "1000px",
            width  = "160%"
          )
        )) # Close mainPanel
        
      ) # Close sidebarLayout
    ),
    # Close tabPanel for cv predictions tab
    
    ### TAB 3 : VARIABLE IMPORTANCE ###
    tabPanel(
      title = "Variable importance",
      value = "variable_importance_tab",
      sidebarLayout(
        sidebarPanel(
          # Specify width of panel
          width = 3,
          
          # Update plot button
          actionButton(
            inputId = "variable_importance_update",
            label   = "Update Figure",
            class   = "btn-primary",
            width   = '100%'
          ),
          # close actionButton
          
          # Data selection group
          wellPanel(
            style = "background-color: #c6d2ed;",
            h4("Results selection"),
            
            # Which vaccine?
            tags$label(
              "Vaccine",
              tags$span(icon("info-circle"), id    = "variable_importance_vaccine_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
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
              title   = "Choose the vaccine for which to display CV predictions",
              placement = "right",
              trigger   = "hover"
            ),
            
            # Which feature_set?
            tags$label(
              "Feature set",
              tags$span(icon("info-circle"), id    = "variable_importance_feature_set_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
            ),
            selectInput(
              inputId = "variable_importance_feature_set",
              label   = NULL,
              choices = feature_set_set,
              selected = feature_set_set[2],
              multiple = FALSE
            ),
            bsTooltip(
              id      = "variable_importance_feature_set_info",
              title   = "Choose the feature set for which to display variable importance",
              placement = "right",
              trigger   = "hover"
            ),
            
            # Results from which approach?
            tags$label(
              "Temporal feature selection approach:",
              tags$span(icon("info-circle"), id    = "variable_importance_approach_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
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
            
            # Include clinical variables?
            tags$label(
              "Include baseline clinical variables?",
              tags$span(icon("info-circle"), id    = "variable_importance_include_clinical_info", style = "cursor: help; color: #337ab7; margin-left: 4px;")
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
            )
          ),
          # Close first wellPanel (data selection)
          
          # ---------- Download section ----------
          wellPanel(
            style = "background-color: #ced4da;",
            h4("Download Options"),
            # choose format
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
            
            # size controls (in inches)
            numericInput(
              inputId = "variable_importance_download_width",
              label   = "Width of download (inches):",
              value   = 35,
              min     = 1
            ),
            numericInput(
              inputId = "variable_importance_download_height",
              label   = "Height of download (inches):",
              value   = 15,
              min     = 1
            ),
            
            numericInput(
              inputId = "variable_importance_download_dpi",
              label   = "Resolution (DPI):",
              value   = 300,
              min     = 72
            ),
            
            ## Download plot button
            div(
              style = "width: 100%;",
              downloadButton(
                outputId = "download_variable_importance_plot",
                label    = "Download Variable Importance Plot",
                class    = "btn-success"
              )
            )
          ) # Close download wellPanel
          
        ),
        # Close sidebarPanel
        
        # Plot settings (main panel)
        mainPanel(width = 9, div(
          style = "text-align: center;",
          plotOutput(
            outputId = "variable_importance_plot",
            height = "1000px",
            width  = "160%"
          )
        )) # Close mainPanel
        
      ) # Close sidebarLayout
    ) # Close tabPanel for variable importance tab
    
  ) # Close tabsetPanel
  
) # Close fluidPage


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
        input$results_comparison_include_clinical
      )
      
      # Fetch the object from the environment
      selected_results <- get(results_name, envir = .GlobalEnv)
      
      # Now plot the results
      
      # --- 1. Build vaccine x set grid and extract Rspearman + sRMSE in one pass ----
      vaccines <- names(selected_results)
      
      # gather all predictor-set names that appear anywhere
      all_sets <- unique(unlist(lapply(selected_results, names)))
      
      # Extract numeric part from strings like "day 7"
      day_sets <- all_sets[grepl("^Day\\s*[0-9]+$", all_sets)]
      day_nums <- as.numeric(sub("Day\\s*", "", day_sets))
      
      # Combine into the desired order
      all_sets <- c("clinical", paste0("Day ", sort(day_nums)))
      
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
          name = "R",
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
        ggtitle("Spearman R") +
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
      
      
      # --- Combine with shared legend and common title ---
      combined <- heatmap_plot_R / heatmap_plot_sRMSE +
        plot_layout(
          ncol = 1,
          heights = c(1, 1),
          guides = "collect"   # collect legends into one
        ) +
        plot_annotation(
          title = "Evaluation metrics of CV predictions",
          subtitle = "Sequential prediction set approach, clinical variables included",
          theme = theme(
            plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 20, face = "bold", hjust = 0.5)
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
        "results_comparison",
        Sys.Date(),
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
  

} # close server function

# Run the shiny application
shinyApp(ui, server)
