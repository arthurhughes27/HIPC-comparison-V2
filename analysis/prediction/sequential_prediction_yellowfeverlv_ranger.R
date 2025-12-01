# A script to perform within-vaccine predictions based on sequential transcriptomics data over time
# The vaccines were selected based on number of available samples and heterogeneity in immune response
# The timepoints were selected per-vaccine based on the number of available samples
# The justification for these decisions can be found in the files analysis/transcriptomic_sample_descriptions.R
# and analysis/immResp_sample_descriptions.R
# In this file, we use random forest models to make predictions in a nested cross-validation framework within each vaccine's data

# Packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(caret)
library(doParallel)
library(stringr)
library(ggtext)
library(forcats)
library(scales)
library(patchwork)

# Directory to store engineered data
processed_data_folder = "data"

# Directory to store figures
prediction_figures_folder = fs::path("output", "figures", "prediction")

# Paths of transformed data
## With TBA modules
# p_load_sequential_list_yellowfeverlv_withTBA = fs::path(processed_data_folder,
#                                                            "sequential_list_yellowfeverlv_withTBA.rds")

## Without TBA modules
p_load_sequential_list_yellowfeverlv_withoutTBA = fs::path(processed_data_folder,
                                                        "sequential_list_yellowfeverlv_withoutTBA.rds")



# Load sequential prediction sets lists
# sequential_list_yellowfeverlv_withTBA = readRDS(p_load_sequential_list_yellowfeverlv_withTBA)
sequential_prediction_sets = readRDS(p_load_sequential_list_yellowfeverlv_withoutTBA)

# Assign each prediction set a colour code for later visualisation
prediction_set_colours = list(
  "clinical" = "#001219",
  "Day 0" = "#0a9396",
  "Day 1" = "#e9d8a6",
  "Day 3" = "#ee9b00",
  "Day 7" = "#bb3e03",
  "Day 10" = "#ae2012",
  "Day 14" = "#780000"
)

# prediction_set_list = sequential_prediction_sets
# transformation = "mean-dearseq"
# set_name = "Day 0"
# n_outer = 10
# n_inner = 5
# response_name = "immResp_MFC_anyAssay_log2_MFC"
# clinical_cols = c("age_imputed", "gender")
# include_clinical = TRUE
# seed = 22072025
# n_cores = 8
# verbose = FALSE


# Write a function to perform prediction for a given vaccine with information sequential at a given timepoint
sequential_prediction_function = function(prediction_set_list,
                                          set_name = "clinical",
                                          transformation = c("mean", "median", "max", "mean-z", "median-z", "max-z", "mean-rank", "median-rank", "max-rank", "ssgsea"),
                                          n_outer,
                                          n_inner,
                                          response_name,
                                          clinical_cols,
                                          include_clinical = TRUE,
                                          seed = 22072025,
                                          n_cores = 1,
                                          verbose = TRUE) {
  
  
  # Set seeds
  set.seed(seed)
  RNGkind("L'Ecuyer-CMRG")
  
  # Extract the number corresponding to the prediction set e.g. X from "Day X"
  if (!set_name == "clinical") {
    # Make an exception if the set is the clinical variables only
    timepoint_number = as.numeric(sub("Day ", "", set_name))
  } else {
    timepoint_number = "clinical"
  }
  
  # Extract the predictor set from the predictor set list
  df_predict = prediction_set_list[[set_name]][[transformation]] 
  
  if (is.null(df_predict)) {
    return(NULL)
  } else {
    df_predict = df_predict %>%
      dplyr::select(-participant_id)
  }
  
  # Extract the response from the prediction set
  y_vec = df_predict[[response_name]]
  
  # Include or exclude clinical variables
  if (include_clinical) {
    X_mat = df_predict %>%
      dplyr::select(-all_of(response_name), -vaccine_name, -vaccine_colour) %>%
      as.data.frame()
  } else {
    X_mat = df_predict %>%
      dplyr::select(
        -all_of(response_name),
        -all_of(clinical_cols),
        -vaccine_name,
        -vaccine_colour
      ) %>%
      as.data.frame()
  }
  
  # Randomly assign outer folds
  folds <- createFolds(y_vec, k = n_outer, list = TRUE)
  
  # initialise parallel computation procedure if desired
  if (n_cores > 1) {
    cl <- makePSOCKcluster(n_cores)
    parallel::clusterSetRNGStream(cl, iseed = seed)
    registerDoParallel(cl)
  }
  
  # Initialise vector to store out-of-fold predictions
  oof_pred <- numeric(nrow(df_predict))
  
  # initialise list for variable importance per-fold
  vi_list <- list()
  vi_list_names <- character(length(folds))
  
  for (fold_idx in seq_along(folds)) {
    # outer cross-validation loop
    test_idx  <- folds[[fold_idx]] # Test folds id's
    train_idx <-
      setdiff(seq_len(nrow(df_predict)), test_idx) # Train fold id's
    
    # Training set
    x_train <- as.data.frame(X_mat[train_idx, ]) # predictors
    y_train <- y_vec[train_idx] # response
    
    # Testing set
    x_test  <-
      as.data.frame(X_mat[test_idx, , drop = FALSE]) # predictors
    
    # inner cross-validation loop setup
    ctrl <- trainControl(
      method = "cv",
      number = n_inner,
      allowParallel = ifelse(n_cores > 1, TRUE, FALSE)
    )
    
    # Register time before fitting model
    time_before = Sys.time()
    # Fit model
    fit <- train(
      x = x_train,
      y = y_train,
      method = "ranger",
      trControl = ctrl,
      tuneLength = 10,
      metric = "RMSE",
      num.trees = 250,
      num.threads = 1,
      importance = "permutation",
    )
    
    # Time after fitting
    time_after = Sys.time()
    # Fit time
    time_elapsed = round(as.numeric(time_after - time_before, units = "secs"), 0)
    
    # Training variable importance
    vi_train_fold <- varImp(fit, scale = T)$importance
    names(vi_train_fold) <- "importance"
    
    # Store VI information in a list
    vi_list[[fold_idx]] <- vi_train_fold %>%
      tibble::rownames_to_column(var = "feature") %>%
      mutate(fold = fold_idx) %>%
      select(fold, feature, importance)
    
    vi_list_names[fold_idx] <- paste0("fold", fold_idx)
    
    # If verbose mode desired, print a message for the outer fold completion
    if (verbose) {
      print(
        paste0(
          "Fold ",
          fold_idx,
          " of ",
          n_outer,
          " completed in ",
          time_elapsed,
          " seconds."
        )
      )
    }
    
    # Predict on the held‑out fold
    oof_pred[test_idx] <- predict(fit, newdata = x_test)
    
    # Clear the unused memory in R
    gc()
  }
  
  # Stop parallel backend
  stopCluster(cl)
  registerDoSEQ()
  
  # Now derive the results to be saved as a list
  
  # Store the CV predictions in a dataframe
  cv_results <- data.frame("observed"  = y_vec, "predicted" = oof_pred)
  
  # Store derived CV evaluation metrics in a list
  metrics = list(
    "Rsquared" = cor(cv_results$observed, cv_results$predicted, method = "pearson")^2,
    "Rspearman" = cor(cv_results$observed, cv_results$predicted, method = "spearman"),
    "RMSE" = sqrt(mean((cv_results$observed - cv_results$predicted)^2, na.rm = TRUE
    )),
    "sRMSE" = (sqrt(mean((cv_results$observed - cv_results$predicted)^2, na.rm = TRUE
    ))) / sd(cv_results$observed)
  )
  
  # Save the vaccine colour for later plotting
  vaccine_colour = df_predict %>%
    pull(vaccine_colour) %>%
    unique()
  
  # Store cv versus predicted values plot
  
  # compute limits
  lims <- range(c(cv_results$observed, cv_results$predicted), na.rm = TRUE)
  # add some padding to the limits
  pad  <- diff(lims) * 0.02
  # store limits with padding
  lims <- c(lims[1] - pad, lims[2] + pad)
  
  # compute metrics annotation and place it near top-left
  annot_df <- data.frame("label" = sprintf("R = %.2f\nsRMSE = %.2f", metrics$Rspearman, metrics$sRMSE)) %>%
    mutate(x = lims[1] + 0.02 * diff(lims),
           y = lims[2] - 0.02 * diff(lims))
  
  if (transformation == "ssgsea"){
    transformation_string = "ssGSEA scores"
  } else if (transformation == "pc1"){
    transformation_string = "PC1 scores"
  } else if(transformation == "mean-dearseq"){
    transformation_string = "Mean dearseq scores"
  } else if(transformation == "median-dearseq"){
    transformation_string = "Median dearseq scores"
  } else if(transformation == "max-dearseq"){
    transformation_string = "Max dearseq scores"
  } else if (transformation %in% c("mean", "median", "max", "mean-z", "median-z", "max-z", "mean-rank", "median-rank", "max-rank")){
    transformation_string = paste0(paste(toupper(substr(transformation, 1, 1)), substr(transformation, 2, nchar(transformation)), sep=""), " of gene-set expression")
  }
  
  # Create a subtitle describing the prediction parameters
  if (timepoint_number == "clinical") {
    subtitle = paste0("Baseline demographic information only")
  } else if (include_clinical && timepoint_number != "clinical") {
    subtitle = paste0(transformation_string, 
      " at day ",
      timepoint_number,
      " post-vaccination (demographic variables included)"
    )
  } else {
    subtitle = paste0(transformation_string, 
                      " at day ",
                      timepoint_number,
                      " post-vaccination (demographic variables not included)"
    )
  }
  
  # cv predicted versus observed values plot
  cv_prediction_plot = cv_results %>%
    ggplot(aes(x = observed, y = predicted)) +
    geom_point(alpha = 0.8,
               size = 2,
               colour = vaccine_colour) +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = "red",
      linewidth = 1,
      linetype = "dashed"
    ) +
    # same limits on both axes
    scale_x_continuous(limits = lims, expand = c(0, 0)) +
    scale_y_continuous(limits = lims, expand = c(0, 0)) +
    coord_fixed() +
    labs(
      x = "Observed",
      y = "Predicted",
      title = "Yellow Fever (LV): cross-validation predicted versus observed values",
      subtitle =  subtitle
    ) +
    # add annotation in top-left
    geom_label(
      data = annot_df,
      aes(x = x, y = y, label = label),
      hjust = 0,
      vjust = 1,
      size = 10,
      fill = "white",
      alpha = 0.8,
      linewidth = 0     # no border line
    ) +
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 28),
      plot.subtitle = element_text(hjust = 0.5, size = 19),
      axis.title = element_text(size = 21),
      axis.text  = element_text(size = 14)
    )
  
  # Store average normalised variable importance across outer folds
  vi_all <- bind_rows(vi_list)
  
  pred_cols <- unlist(prediction_set_colours)
  
  # prepare summary of variable importance across folds
  vi_summary <- vi_all %>%
    group_by(feature) %>%
    summarise(
      n_folds = n(),
      mean_imp = mean(importance, na.rm = TRUE),
      sd_imp = sd(importance),
      .groups = "drop"
    ) %>%
    mutate(mean_rank = dplyr::dense_rank(dplyr::desc(mean_imp))) %>%
    arrange(mean_rank) %>%
    mutate(
      feature_group = case_when(
        feature %in% clinical_cols ~ "clinical",
        # <- titlecase to match pred_cols
        str_detect(feature, "_Day\\d+") ~ str_extract(feature, "Day\\d+"),
        TRUE ~ "other"
      ),
      feature_group = str_replace(feature_group, "Day", "Day "),
      feature_colour = pred_cols[feature_group]
    ) %>%
    mutate(feature_colour = ifelse(is.na(feature_colour), "#777777", feature_colour))
  
  # plot top 20 mean standardised feature importances with coloured errorbars and points
  topN <- 20
  
  # Prepare plotting data
  plot_df <- vi_summary %>%
    slice_min(order_by = mean_rank, n = topN) %>%
    mutate(
      mean_imp = as.numeric(mean_imp),
      sd_imp   = as.numeric(sd_imp),
      feature  = forcats::fct_reorder(feature, mean_imp),
      xmin = pmax(0, mean_imp - sd_imp),
      xmax = pmin(100, mean_imp + sd_imp)
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
    scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0.02, 0.12))) +
    labs(
      x = "Mean standardised training-set variable importance",
      y = NULL,
      title = paste0(
        "Mean standardised variable importance (top ",
        min(topN, length(unique(
          vi_summary$feature
        ))),
        " features)"
      ),
      subtitle = bquote("sRMSE =" ~ .(sprintf(
        "%.2f", metrics$sRMSE
      )) * ", " ~ rho == .(sprintf(
        "%.2f", metrics$Rspearman
      )))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(
        size = 25,
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
  
  # Final output of function as a list
  res = list(
    "cv_results" = cv_results,
    "var_imp" = vi_summary,
    "vaccine_colour" = vaccine_colour,
    "metrics" = metrics,
    "plots" = list("cv_results" = cv_prediction_plot, "var_imp" = vi_plot)
  )
  
}

# Now we need to apply this function to all of our vaccines and sequential predictor sets

# prepare a results list with the same structure as the predictor set list
prediction_results_sequential_list_yellowfeverlv_withoutTBA = vector("list", length(sequential_prediction_sets))
names(prediction_results_sequential_list_yellowfeverlv_withoutTBA) <- names(sequential_prediction_sets)


# For each predictor set
for (set in names(sequential_prediction_sets)){
  transformation_names <- names(sequential_prediction_sets[[set]])
  prediction_results_sequential_list_yellowfeverlv_withoutTBA[[set]] <- vector("list", length(transformation_names))
  names(prediction_results_sequential_list_yellowfeverlv_withoutTBA[[set]]) <- transformation_names
  # For each transformation
  for (trans in transformation_names) {
    
    # Progress message
    message(sprintf("Running: set = '%s'   transformation = '%s'", set, trans))
    # Run prediction function
    res <-
      sequential_prediction_function(
        prediction_set_list = sequential_prediction_sets,
        transformation = trans,
        set_name = set,
        n_outer = 10,
        n_inner = 5,
        response_name = "immResp_MFC_anyAssay_log2_MFC",
        clinical_cols = c("age_imputed", "gender"),
        include_clinical = TRUE,
        seed = 22072025,
        n_cores = 10,
        verbose = FALSE
      )
    
    # store the result
    prediction_results_sequential_list_yellowfeverlv_withoutTBA[[set]][[trans]] <- res
    
    gc()
    
  }
}


# Save these results
p_prediction_results_sequential_list_yellowfeverlv_withoutTBA = fs::path(
  "output",
  "results",
  "prediction",
  "prediction_results_sequential_list_yellowfeverlv_withoutTBA.rds"
)


saveRDS(prediction_results_sequential_list_yellowfeverlv_withoutTBA,
        file = p_prediction_results_sequential_list_yellowfeverlv_withoutTBA)

# Plot - bar chart of sRMSE across predictor sets, stratified by transformation

prediction_results_sequential_list_yellowfeverlv_withoutTBA = readRDS(p_prediction_results_sequential_list_yellowfeverlv_withoutTBA)

# your list (use the real object from your environment)
res_list <- prediction_results_sequential_list_yellowfeverlv_withoutTBA

# 1) baseline sRMSE (clinical)
baseline_srmse <- res_list[["clinical"]][["clinical"]][["metrics"]][["sRMSE"]]

# 2) get predictor-set names in original order, but exclude "clinical" for the bars
all_sets <- names(res_list)
sets_for_bars <- all_sets[all_sets != "clinical"]  # preserves original order

# 3) build a data.frame of (set, transformation, sRMSE)
df <- do.call(rbind, lapply(sets_for_bars, function(set_name) {
  trans_names <- names(res_list[[set_name]])
  # if there are no transformations, return NULL to skip
  if (length(trans_names) == 0) return(NULL)
  data.frame(
    set = set_name,
    trans = trans_names,
    sRMSE = vapply(trans_names, function(tr) {
      # safe extraction; returns NA if something unexpected happens
      v <- tryCatch(res_list[[set_name]][[tr]][["metrics"]][["sRMSE"]],
                    error = function(e) NA_real_)
      as.numeric(v)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
}))

# ensure factor order of sets follows the list order
df$set <- factor(df$set, levels = sets_for_bars)

# optional: order transformations within legend (alphabetical here; change if you want)
df$trans <- factor(df$trans, levels = unique(df$trans))

# 4) plot: grouped bars per predictor set, one bar per transformation
# Create a dummy dataframe for the baseline line
# baseline line as a dummy data frame
my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(14)

baseline_df <- data.frame(
  label = "Baseline Demographic Model",
  y = baseline_srmse
)

p_srmse <- ggplot(df, aes(x = set, y = sRMSE, fill = trans)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  
  geom_hline(
    data = baseline_df,
    aes(yintercept = y, linetype = label),
    color = "black",
    size = 0.8,
    show.legend = TRUE
  ) +
  
  scale_fill_manual(values = my_colors) +
  scale_linetype_manual(
    name = "", 
    values = c("Baseline Demographic Model" = "dashed")
  ) +
  
  guides(
    linetype = guide_legend(
      order = 1,                       # appears first (top)
      override.aes = list(fill = NA, colour = "black")
    ),
    fill = guide_legend(
      order = 2,                       # appears below (bottom)
      override.aes = list(linetype = 0, colour = NA)
    )
  ) +
  
  labs(
    x = "Predictor set",
    y = "Standardised RMSE (sRMSE)",
    fill = "Transformation",
    title = "Yellow Fever : Random Forest C.V. prediction standardised RMSE by predictor set and transformation",
    subtitle = "Baseline demographic model performance shown as dashed line (smaller = better)"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    legend.box = "vertical",           # stack legends vertically
    legend.spacing.y = unit(0.5, "cm"), # vertical space between boxes
    legend.key.width = unit(3, "lines"),
    plot.margin = margin(10, 20, 20, 10)
  ) +
  coord_cartesian(clip = "off") +
  ylim(0, 1)

print(p_srmse)


# Same for R

# your list (use the real object from your environment)
res_list <- prediction_results_sequential_list_yellowfeverlv_withoutTBA

# 1) baseline sRMSE (clinical)
baseline_Rspearman <- res_list[["clinical"]][["clinical"]][["metrics"]][["Rspearman"]]

# 2) get predictor-set names in original order, but exclude "clinical" for the bars
all_sets <- names(res_list)
sets_for_bars <- all_sets[all_sets != "clinical"]  # preserves original order

# 3) build a data.frame of (set, transformation, sRMSE)
df <- do.call(rbind, lapply(sets_for_bars, function(set_name) {
  trans_names <- names(res_list[[set_name]])
  # if there are no transformations, return NULL to skip
  if (length(trans_names) == 0) return(NULL)
  data.frame(
    set = set_name,
    trans = trans_names,
    Rspearman = vapply(trans_names, function(tr) {
      # safe extraction; returns NA if something unexpected happens
      v <- tryCatch(res_list[[set_name]][[tr]][["metrics"]][["Rspearman"]],
                    error = function(e) NA_real_)
      as.numeric(v)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
}))

# ensure factor order of sets follows the list order
df$set <- factor(df$set, levels = sets_for_bars)

# optional: order transformations within legend (alphabetical here; change if you want)
df$trans <- factor(df$trans, levels = unique(df$trans))

# 4) plot: grouped bars per predictor set, one bar per transformation
# Create a dummy dataframe for the baseline line
# baseline line as a dummy data frame
baseline_df <- data.frame(
  label = "Baseline Demographic Model",
  y = baseline_Rspearman
)

p_Rspearman <- ggplot(df, aes(x = set, y = Rspearman, fill = trans)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  
  geom_hline(
    data = baseline_df,
    aes(yintercept = y, linetype = label),
    color = "black",
    size = 0.8,
    show.legend = TRUE
  ) +
  
  scale_fill_manual(values = my_colors) +
  scale_linetype_manual(
    name = "", 
    values = c("Baseline Demographic Model" = "dashed")
  ) +
  
  guides(
    linetype = guide_legend(
      order = 1,                       # appears first (top)
      override.aes = list(fill = NA, colour = "black")
    ),
    fill = guide_legend(
      order = 2,                       # appears below (bottom)
      override.aes = list(linetype = 0, colour = NA)
    )
  ) +
  
  labs(
    x = "Predictor set",
    y = "Spearman R",
    fill = "Transformation",
    title = "Yellow Fever : Random Forest C.V. prediction Spearman R by predictor set and transformation",
    subtitle = "Baseline demographic model performance shown as dashed line (higher = better)"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    legend.box = "vertical",           # stack legends vertically
    legend.spacing.y = unit(0.5, "cm"), # vertical space between boxes
    legend.key.width = unit(3, "lines"),
    plot.margin = margin(10, 20, 20, 10)
  ) +
  coord_cartesian(clip = "off") +
  ylim(0, 1)

print(p_Rspearman)


# Save evaluation plots

# Directory to store figures
prediction_figures_folder = fs::path("output", "figures", "prediction", "transformation")

ggsave(
  filename = "evaluation_transformation_srmse_yellowfeverlv_ranger.pdf",
  path = prediction_figures_folder,
  plot = p_srmse,
  width = 45,
  height = 18,
  units = "cm"
)

ggsave(
  filename = "evaluation_transformation_rspearman_yellowfeverlv_ranger.pdf",
  path = prediction_figures_folder,
  plot = p_Rspearman,
  width = 45,
  height = 18,
  units = "cm"
)
