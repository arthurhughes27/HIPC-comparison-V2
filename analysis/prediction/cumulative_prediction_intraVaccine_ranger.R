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
library(scales)
library(patchwork)

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
  "Day 0" = "#0a9396",
  "Day 1" = "#e9d8a6",
  "Day 3" = "#ee9b00",
  "Day 7" = "#bb3e03",
  "Day 10" = "#ae2012",
  "Day 14" = "#780000"
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
  df_predict = prediction_set_list[[vaccine]][[set_name]] %>%
    dplyr::select(-participant_id)
  
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
        -all_of(response_name),-all_of(clinical_cols),-vaccine_name,-vaccine_colour
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
  
  # Create a subtitle describing the prediction parameters
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

# Now we need to apply this function to all of our vaccines and cumulative predictor sets

# prepare a results list with the same structure as the predictor set list
prediction_results_all_cumulative_withClinical_withTBA_ranger <- vector("list", length(cumulative_prediction_sets))
names(prediction_results_all_cumulative_withClinical_withTBA_ranger) <- names(cumulative_prediction_sets)

# For each vaccine
for (vac in names(cumulative_prediction_sets)) {
  # initialize inner list with same names as the predictor-sets for this vaccine
  set_names <- names(cumulative_prediction_sets[[vac]])
  prediction_results_all_cumulative_withClinical_withTBA_ranger[[vac]] <- vector("list", length(set_names))
  names(prediction_results_all_cumulative_withClinical_withTBA_ranger[[vac]]) <- set_names
  
  # For each predictor set
  for (set_name in set_names) {
    # Progress message
    message(sprintf("Running: vaccine = '%s'   set = '%s'", vac, set_name))
    # Run prediction function
    res <-
      cumulative_prediction_function(
        prediction_set_list = cumulative_prediction_sets,
        vaccine = vac,
        set_name = set_name,
        n_outer = 10,
        n_inner = 5,
        response_name = "immResp_MFC_anyAssay_log2_MFC",
        clinical_cols = c("age_imputed", "gender"),
        include_clinical = TRUE,
        seed = 22072025,
        n_cores = 8,
        verbose = FALSE
      )
    
    # store the result
    prediction_results_all_cumulative_withClinical_withTBA_ranger[[vac]][[set_name]] <- res
    
    gc()
  }
}

# Save these results
p_prediction_results_all_cumulative_withClinical_withTBA_ranger = fs::path(
  "output",
  "results",
  "prediction",
  "prediction_results_all_cumulative_withClinical_withTBA_ranger.rds"
)

saveRDS(prediction_results_all_cumulative_withClinical_withTBA_ranger,
        file = p_prediction_results_all_cumulative_withClinical_withTBA_ranger)

# Do the same prediction task but remove clinical variables from the prediction

prediction_results_all_cumulative_withoutClinical_withTBA_ranger <- vector("list", length(cumulative_prediction_sets))
names(prediction_results_all_cumulative_withoutClinical_withTBA_ranger) <- names(cumulative_prediction_sets)
# For each vaccine
for (vac in names(cumulative_prediction_sets)) {
  # initialize inner list with same names as the predictor-sets for this vaccine
  set_names <- names(cumulative_prediction_sets[[vac]])[-1]
  prediction_results_all_cumulative_withoutClinical_withTBA_ranger[[vac]] <- vector("list", length(set_names))
  names(prediction_results_all_cumulative_withoutClinical_withTBA_ranger[[vac]]) <- set_names
  
  # For each predictor set
  for (set_name in set_names) {
    # Progress message
    message(sprintf("Running: vaccine = '%s'   set = '%s'", vac, set_name))
    # Run prediction function
    res <-
      cumulative_prediction_function(
        prediction_set_list = cumulative_prediction_sets,
        vaccine = vac,
        set_name = set_name,
        n_outer = 10,
        n_inner = 5,
        response_name = "immResp_MFC_anyAssay_log2_MFC",
        clinical_cols = c("age_imputed", "gender"),
        include_clinical = FALSE,
        seed = 22072025,
        n_cores = 10,
        verbose = FALSE
      )
    
    # store the result
    prediction_results_all_cumulative_withoutClinical_withTBA_ranger[[vac]][[set_name]] <- res
    
    gc()
  }
}

# Save these results
p_prediction_results_all_cumulative_withoutClinical_withTBA_ranger = fs::path(
  "output",
  "results",
  "prediction",
  "prediction_results_all_cumulative_withoutClinical_withTBA_ranger.rds"
)

saveRDS(prediction_results_all_cumulative_withoutClinical_withTBA_ranger,
        file = p_prediction_results_all_cumulative_withoutClinical_withTBA_ranger)

#-------------------------------
# Plot results
#-------------------------------

# Load the 'with clinical' results

p_prediction_results_all_cumulative_withClinical_withTBA_ranger = fs::path(
  "output",
  "results",
  "prediction",
  "prediction_results_all_cumulative_withClinical_withTBA_ranger.rds"
)

prediction_results_all_cumulative_withClinical_withTBA_ranger = readRDS(p_prediction_results_all_cumulative_withClinical_withTBA_ranger)

# --- 1. Build vaccine x set grid and extract Rspearman + sRMSE in one pass ----
vaccines <- names(prediction_results_all_cumulative_withClinical_withTBA_ranger)

# gather all predictor-set names that appear anywhere
all_sets <- c("clinical", "Day 0", "Day 1", "Day 3", "Day 7", "Day 10", "Day 14")

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
  if (is.null(prediction_results_all_cumulative_withClinical_withTBA_ranger[[vac]]))
    next
  res_set <- prediction_results_all_cumulative_withClinical_withTBA_ranger[[vac]][[setn]]
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
  labs(x = "Predictor set",
       y = "Vaccine",
       title = expression(paste(
         "Spearman ", rho, " by vaccine and predictor set - Random Forest"
       ))) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 25, hjust = 0.5),
  )

# Print Spearman R heatmap
print(heatmap_plot_R)

# --- 3. sRMSE heatmap (0 = green -> 1 = white; values > 1 shown as white; NA -> grey) ----
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
  labs(x = "Predictor set", y = "Vaccine", title = "Standardised RMSE by vaccine and predictor set - Random Forest") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 25, hjust = 0.5),
  )

# Print sRMSE heatmap
print(heatmap_plot_sRMSE)

title = expression(paste("Spearman ", rho, " by vaccine and predictor set - Random Forest"))

# --- Modify heatmap plots ---
heatmap_plot_R_mod <- heatmap_plot_R +
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

heatmap_plot_sRMSE_mod <- heatmap_plot_sRMSE +
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
combined <- heatmap_plot_R_mod / heatmap_plot_sRMSE_mod +
  plot_layout(
    ncol = 1,
    heights = c(1, 1),
    guides = "collect"   # collect legends into one
  ) +
  plot_annotation(
    title = "Evaluation metrics of CV predictions - Random Forest",
    subtitle = "Cumulative prediction set approach, clinical variables included",
    theme = theme(
      plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  ) & theme(
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.key.width = unit(1, "cm")  # optional: reduce legend width
  )

# --- Display combined plot ---
print(combined)

# --- Save figure ---
ggsave(
  filename = "evaluation_combined_cumulative_withClinical_withTBA_ranger.pdf",
  path = prediction_figures_folder,
  plot = combined,
  width = 40,
  height = 25,
  units = "cm"
)

# Now repeat for the 'without Clinical' results

# Load the 'without clinical' results

p_prediction_results_all_cumulative_withoutClinical_withTBA_ranger = fs::path(
  "output",
  "results",
  "prediction",
  "prediction_results_all_cumulative_withoutClinical_withTBA_ranger.rds"
)

prediction_results_all_cumulative_withoutClinical_withTBA_ranger = readRDS(p_prediction_results_all_cumulative_withoutClinical_withTBA_ranger)

# --- 1. Build vaccine x set grid and extract Rspearman + sRMSE in one pass ----
vaccines <- names(prediction_results_all_cumulative_withoutClinical_withTBA_ranger)

# gather all predictor-set names that appear anywhere
all_sets <- c("Day 0", "Day 1", "Day 3", "Day 7", "Day 10", "Day 14")

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
  if (is.null(prediction_results_all_cumulative_withoutClinical_withTBA_ranger[[vac]]))
    next
  res_set <- prediction_results_all_cumulative_withoutClinical_withTBA_ranger[[vac]][[setn]]
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
  labs(x = "Predictor set",
       y = "Vaccine",
       title = expression(paste(
         "Spearman ", rho, " by vaccine and predictor set - Random Forest"
       ))) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 25, hjust = 0.5),
  )

# Print Spearman R heatmap
print(heatmap_plot_R)


# --- 3. sRMSE heatmap (0 = green -> 1 = white; values > 1 shown as white; NA -> grey) ----
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
  labs(x = "Predictor set", y = "Vaccine", title = "Standardised RMSE by vaccine and predictor set - Random Forest") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 25, hjust = 0.5),
  )

# Print sRMSE heatmap
print(heatmap_plot_sRMSE)


# --- Modify heatmap plots ---
heatmap_plot_R_mod <- heatmap_plot_R +
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

heatmap_plot_sRMSE_mod <- heatmap_plot_sRMSE +
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
combined <- heatmap_plot_R_mod / heatmap_plot_sRMSE_mod +
  plot_layout(
    ncol = 1,
    heights = c(1, 1),
    guides = "collect"   # collect legends into one
  ) +
  plot_annotation(
    title = "Evaluation metrics of CV predictions - Random Forest",
    subtitle = "Cumulative prediction set approach, clinical variables not included",
    theme = theme(
      plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  ) & theme(
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.key.width = unit(1, "cm")  # optional: reduce legend width
  )

# --- Display combined plot ---
print(combined)

# --- Save figure ---
ggsave(
  filename = "evaluation_combined_cumulative_withoutClinical_withTBA_ranger.pdf",
  path = prediction_figures_folder,
  plot = combined,
  width = 40,
  height = 25,
  units = "cm"
)
