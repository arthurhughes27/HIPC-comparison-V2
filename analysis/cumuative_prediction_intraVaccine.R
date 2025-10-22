# A script to perform within-vaccine predictions based on accumulated transcriptomics data over time
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

# Directory to store engineered data
processed_data_folder = "data"

# Directory to store figures
prediction_figures_folder = fs::path("output", "figures", "prediction")

# Path to cumulative prediction sets list
p_load_cumulative_prediction_sets <- fs::path(processed_data_folder,
                                              "cumulative_prediction_sets_list.rds")

# Load cumulative prediction sets list
cumulative_prediction_sets = readRDS(p_load_cumulative_prediction_sets)

# Assign each prediction set a colour code for later visualisation
prediction_set_colours = list(
  "clinical" = "#001219",
  "Day 0" = "#005f73",
  "Day 1" = "#0a9396",
  "Day 3" = "#94d2bd",
  "Day 7" = "#e9d8a6",
  "Day 10" = "#ee9b00",
  "Day 14" = "#bb3e03",
  "Day 28" = "#ae2012"
)


# Write a function to perform prediction for a given vaccine with information accumulated up to a given timepoint
cumulative_prediction_function = function(prediction_set_list,
                                          vaccine,
                                          set_name = "clinical",
                                          n_outer,
                                          n_inner,
                                          response_name,
                                          clinical_cols,
                                          include_clinical = TRUE,
                                          seed = 22072025,
                                          n_cores = 1,
                                          verbose = TRUE) {
  if (!set_name == "clinical") {
    timepoint_number = as.numeric(sub("Day ", "", set_name))
  } else {
    timepoint_number = "clinical"
  }
  
  df_predict = prediction_set_list[[vaccine]][[set_name]] %>%
    dplyr::select(-participant_id)
  
  y_vec = df_predict[[response_name]]
  
  if (include_clinical) {
    X_mat = df_predict %>%
      dplyr::select(-all_of(response_name), -vaccine_name, -vaccine_colour) %>%
      as.data.frame()
  } else {
    X_mat = df_predict %>%
      dplyr::select(
        -all_of(response_name),-all_of(clinical_cols),-vaccine_name,-vaccine_colour
      ) %>%
      as.data.frame()
  }
  
  folds <- createFolds(y_vec, k = n_outer, list = TRUE)
  
  # initialise parallel computation procedure
  if (n_cores > 1) {
    cores <- parallel::detectCores()
    cl <- makePSOCKcluster(1 + cores / 2)
    registerDoParallel(cl)
  }
  
  # Initialise vector for out-of-fold predictions
  oof_pred <- numeric(nrow(df_predict))
  
  # initialise list for variable importance
  vi_list <- list()
  vi_list_names <- character(length(folds))
  
  for (fold_idx in seq_along(folds)) {
    # 10-fold Cross-validation
    test_idx  <- folds[[fold_idx]] # Test folds id's
    train_idx <-
      setdiff(seq_len(nrow(df_predict)), test_idx) # Train fold id's
    
    # Training set
    x_train <- as.data.frame(X_mat[train_idx, ])
    y_train <- y_vec[train_idx]
    
    # Testing set
    x_test  <-
      as.data.frame(X_mat[test_idx, , drop = FALSE])
    
    # Repeated CV: 5 folds
    ctrl <- trainControl(
      method = "cv",
      number = n_inner,
      allowParallel = ifelse(n_cores > 1, TRUE, FALSE)
    )
    
    time_before = Sys.time()
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
    
    time_after = Sys.time()
    time_elapsed = round(as.numeric(time_after - time_before, units = "secs"), 0)
    
    # Training variable importance
    vi_train_fold <- varImp(fit, scale = T)$importance
    names(vi_train_fold) <- "importance"
    
    vi_list[[fold_idx]] <- vi_train_fold %>%
      tibble::rownames_to_column(var = "feature") %>%
      mutate(fold = fold_idx) %>%
      select(fold, feature, importance)
    
    vi_list_names[fold_idx] <- paste0("fold", fold_idx)
    if (verbose) {
      print(paste0(
        "Fold ",
        fold_idx,
        " of ",
        n_outer,
        " completed in ",
        time_elapsed,
        " seconds."
      ))
    }
    
    # Predict on the held‑out fold
    oof_pred[test_idx] <- predict(fit, newdata = x_test)
    
    gc()
  }
  
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
  
  
  # Store cv versus predicted values plot
  
  # Save the vaccine colour for plotting
  vaccine_colour = df_predict %>%
    pull(vaccine_colour) %>%
    unique()
  
  # compute limits
  lims <- range(c(cv_results$observed, cv_results$predicted), na.rm = TRUE)
  # add some padding to the limits
  pad  <- diff(lims) * 0.02
  # store limits with padding
  lims <- c(lims[1] - pad, lims[2] + pad)
  
  # compute annotation (Rspearman R and standardized RMSE) and place it near top-left
  annot_df <- data.frame("label" = sprintf("R = %.2f\nsRMSE = %.2f", metrics$Rspearman, metrics$sRMSE)) %>%
    mutate(x = lims[1] + 0.02 * diff(lims),
           y = lims[2] - 0.02 * diff(lims))
  
  # Create a subtitle
  if (timepoint_number == "clinical") {
    subtitle = paste0("Baseline demographic information only")
  } else if (include_clinical && timepoint_number != "clinical") {
    subtitle = paste0(
      "Gene-set information up to day ",
      timepoint_number,
      " post-vaccination (demographic variables included)"
    )
  } else {
    subtitle = paste0(
      "Gene-set information up to day ",
      timepoint_number,
      " post-vaccination (demographic variables NOT included)"
    )
  }
  
  # plot
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
      title = paste0(vaccine, " : cross-validation predicted versus observed values"),
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
      label.size = 0     # no border line
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
  
  
  pred_cols <- unlist(prediction_set_colours) # named character vector
  
  # prepare summary (your existing code, but ensure clinical capitalisation)
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
  
  # plot topN with coloured errorbars and points
  topN <- 20
  
  plot_df <- vi_summary %>%
    slice_min(order_by = mean_rank, n = topN) %>%
    mutate(
      mean_imp = as.numeric(mean_imp),
      sd_imp   = as.numeric(sd_imp),
      feature  = forcats::fct_reorder(feature, mean_imp),
      xmin = pmax(0, mean_imp - sd_imp),
      xmax = pmin(100, mean_imp + sd_imp)
    )
  
  vi_plot <- plot_df %>%
    ggplot(aes(x = mean_imp, y = feature)) +
    # colour the error bars by feature_group
    geom_errorbarh(
      aes(
        xmin = xmin,
        xmax = xmax,
        colour = feature_group
      ),
      width = 0.25,
      size = 1
    ) +
    # colour the points by the same group
    geom_point(aes(colour = feature_group), size = 3) +
    scale_colour_manual(
      name = "Feature set",
      values = pred_cols,
      na.value = "#777777",
      guide = guide_legend(
        override.aes = list(
          linetype = 1,
          shape = 16,
          size = 5    # make the legend points/lines larger
        ),
        keywidth  = unit(2.2, "cm"),
        # wider legend keys
        keyheight = unit(0.9, "cm"),
        # taller legend keys
        ncol = 1
      )
    ) +
    scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0.02, 0.12))) +
    labs(
      x = "Mean standardised training-set variable importance",
      y = NULL,
      title = "Mean standardised training variable importance (training)",
      subtitle = paste0("Top ", min(topN, length(
        unique(vi_summary$feature)
      )), " features across folds")
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
      # larger legend title
      legend.text  = element_text(size = 14),
      # larger legend labels
      legend.key.size = unit(1.1, "cm")         # spacing around each legend key
    )
  
  
  res = list(
    "cv_results" = cv_results,
    "var_imp" = vi_summary,
    "vaccine_colour" = vaccine_colour,
    "metrics" = metrics,
    "plots" = list("cv_results" = cv_prediction_plot, "var_imp" = vi_plot)
  )
  
}

# Now we need to apply this function to all of our vaccines and cumulative predictor sets

# prepare a same-structured results list
prediction_results_all <- vector("list", length(cumulative_prediction_sets))
names(prediction_results_all) <- names(cumulative_prediction_sets)

for (vac in names(cumulative_prediction_sets)) {
  # initialize inner list with same names as the predictor-sets for this vaccine
  set_names <- names(cumulative_prediction_sets[[vac]])
  prediction_results_all[[vac]] <- vector("list", length(set_names))
  names(prediction_results_all[[vac]]) <- set_names
  
  for (set_name in set_names) {
    message(sprintf("Running: vaccine = '%s'   set = '%s'", vac, set_name))
    # call your function (sequential)
    res <-
      cumulative_prediction_function(
        prediction_set_list = cumulative_prediction_sets,
        vaccine = vac,
        set_name = set_name,
        n_outer = 10,
        n_inner = 5,
        response_name = "immResp_MFC_anyAssay_log2_MFC",
        clinical_cols = c("age_imputed", "gender", "race"),
        include_clinical = TRUE,
        seed = 22072025,
        n_cores = 10,
        verbose = FALSE
      )
    
    # store the result (replacing the predictor-set df with the result list/object)
    prediction_results_all[[vac]][[set_name]] <- res
  }
}
