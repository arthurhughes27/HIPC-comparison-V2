# R script to engineer preprocessed expression data into gene-set level fold-change features.
# We adopt a "cumulative information" approach
# That is, we aim to see how accumulating information from multiple timepoints changes predictive quality
# Here, we start with a baseline set of predictors, then add on predictors derived at each one of some selected timepoints


# Directory to store engineered data
processed_data_folder = "data"

# Path to processed gene-level data and processed genesets
p_load_expr_all_norm <- fs::path(processed_data_folder, "hipc_merged_all_norm.rds")
p_load_btm <- fs::path(processed_data_folder, "BTM_processed.rds")

# Load merged gene-level data
hipc_merged_all_norm = readRDS(p_load_expr_all_norm)
# Load genesets
BTM = readRDS(p_load_btm)

# Filter out observations which do not have an immune response value
hipc_merged_all_norm_response = hipc_merged_all_norm %>%
  filter(!is.na(immResp_MFC_anyAssay_log2_MFC))

# Now code a helper functions to calculate gene-set features at baseline (i.e. mean gene-set expression at baseline)

compute_geneset_baseline <- function(df,
                                     # expression dataframe
                                     genesets,
                                     # list of genesets (i.e. each element is a vector of gene names)
                                     geneset_names,
                                     # vector of geneset names
                                     id_col,
                                     # column name of the participant identifier in df
                                     time_col) {
  ## Basic checks
  if (!is.data.frame(df)) {
    stop("df must be a data.frame") # check df is a dataframe
  }
  if (!is.list(genesets)) {
    stop("genesets must be a list") # check genesets are given as a list
  }
  if (length(genesets) != length(geneset_names)) {
    stop("genesets and geneset_names must have same length") # check each geneset corresponds to a name
  }
  if (!(id_col %in% names(df))) {
    stop("id_col not found in df") # check the participant identifier is present
  }
  if (!(time_col %in% names(df))) {
    stop("time_col not found in df") # check the timepoint identifier is present
  }
  
  # Ensure timepoint is numeric
  df[[time_col]] <- as.numeric(df[[time_col]])
  
  # Extract gene names from df
  gene_cols <- df %>%
    dplyr::select(a1cf:zzz3) %>%
    colnames()
  
  # Since some participants have multiple pre-vaccination timepoints, we need to select the most recent one
  # Extract the data corresponding to all pre-vaccination measurements
  baseline_candidates <- df[!is.na(df[[time_col]]) &
                              df[[time_col]] <= 0, , drop = FALSE]
  
  # order by id ascending, time descending so the first row per id is the closest-to-0 baseline
  ord_bl <- order(baseline_candidates[[id_col]], -baseline_candidates[[time_col]])
  bl_sorted <- baseline_candidates[ord_bl, , drop = FALSE]
  baseline_df <- bl_sorted[!duplicated(bl_sorted[[id_col]]), , drop = FALSE]
  baseline_df <- baseline_df[, c(id_col, gene_cols), drop = FALSE]
  
  # ---- compute geneset baseline means ----
  # convert baseline gene values to numeric matrix
  bl_mat <- as.matrix(baseline_df[, gene_cols, drop = FALSE])
  storage.mode(bl_mat) <- "numeric"
  
  # preallocate result dataframes
  result <-
    data.frame(matrix(nrow = length(rownames(baseline_df)), ncol = (length(geneset_names) + 1)))
  
  # set the rownames
  rownames(result) = baseline_df[[id_col]]
  
  # set the column names
  colnames(result) = c(id_col, paste0(geneset_names, "_Day0"))
  
  # set the participant ids
  result[[id_col]] <- baseline_df[[id_col]]
  
  # fill geneset columns with row-wise means across present genes
  for (i in seq_along(genesets)) {
    # For each geneset
    # Which are the genes in this set
    set_genes <- genesets[[i]]
    # Which of these genes are present in our data?
    present <- intersect(set_genes, gene_cols)
    # What is the mean value of these genes' expression
    set_means <- rowMeans(bl_mat[, present, drop = FALSE], na.rm = TRUE)
    # Set this as the result for the baseline matrix
    result[[paste0(geneset_names[i], "_Day0")]] <- as.numeric(set_means)
  }
  
  rownames(result) <- NULL
  return(result)
}

# The following is a helper function which computes gene-set fold change at a given timepoint.

compute_geneset_fc <- function(df,
                               # expression dataframe
                               genesets,
                               # list of genesets (i.e. each element is a vector of gene names)
                               geneset_names,
                               # vector of geneset names
                               id_col,
                               # column name of the participant identifier in df
                               time_col,
                               # column name of the timepoint identifier in df
                               timepoint) {
  ## Basic checks
  if (!is.data.frame(df)) {
    stop("df must be a data.frame") # check df is a dataframe
  }
  if (!is.list(genesets)) {
    stop("genesets must be a list") # check genesets are given as a list
  }
  if (length(genesets) != length(geneset_names)) {
    stop("genesets and geneset_names must have same length") # check each geneset corresponds to a name
  }
  if (!(id_col %in% names(df))) {
    stop("id_col not found in df") # check the participant identifier is present
  }
  if (!(time_col %in% names(df))) {
    stop("time_col not found in df") # check the timepoint identifier is present
  }
  
  # Ensure timepoint is numeric
  df[[time_col]] <- as.numeric(df[[time_col]])
  
  # Extract gene names from df
  gene_cols <- df %>%
    dplyr::select(a1cf:zzz3) %>%
    colnames()
  
  # Since some participants have multiple pre-vaccination timepoints, we need to select the most recent one
  # Extract the data corresponding to all pre-vaccination measurements
  baseline_candidates <- df[!is.na(df[[time_col]]) &
                              df[[time_col]] <= 0, , drop = FALSE]
  
  # order by id ascending, time descending so the first row per id is the closest-to-0 baseline
  ord_bl <- order(baseline_candidates[[id_col]], -baseline_candidates[[time_col]])
  bl_sorted <- baseline_candidates[ord_bl, , drop = FALSE]
  baseline_df <- bl_sorted[!duplicated(bl_sorted[[id_col]]), , drop = FALSE]
  baseline_df <- baseline_df[, c(id_col, gene_cols), drop = FALSE]
  
  # convert baseline gene values to numeric matrix
  bl_mat <- as.matrix(baseline_df[, gene_cols, drop = FALSE])
  storage.mode(bl_mat) <- "numeric"
  
  # Extract post-vaccination measurements at arugment timepoint
  tp_candidates <-
    df[!is.na(df[[time_col]]) &
         df[[time_col]] == timepoint, , drop = FALSE]
  
  # assume individuals do not have multiple measurements per timepoint
  tp_df <- tp_candidates[, c(id_col, gene_cols), drop = FALSE]
  
  # ---- keep only participants with one baseline AND one post measurement ----
  bl_ids <- as.character(baseline_df[[id_col]])
  tp_ids <- as.character(tp_df[[id_col]])
  common_ids <-
    intersect(bl_ids, tp_ids)  # order follows baseline_df
  
  # Make two dataframes, one containing the baseline expression and
  # the other containing the post-vaccination expression
  # we should make sure these contain only the relevant participants, and they are in the same order
  baseline_sub <-
    baseline_df[match(common_ids, as.character(baseline_df[[id_col]])), , drop = FALSE]
  tp_sub       <-
    tp_df[match(common_ids, as.character(tp_df[[id_col]])), , drop = FALSE]
  
  # ---- convert to numeric matrices and compute difference (post - baseline) ----
  bl_mat <- as.matrix(baseline_sub[, gene_cols, drop = FALSE])
  tp_mat <- as.matrix(tp_sub[, gene_cols, drop = FALSE])
  storage.mode(bl_mat) <- "numeric"
  storage.mode(tp_mat) <- "numeric"
  
  # Compute fold-change
  diff_mat <- tp_mat - bl_mat
  rownames(diff_mat) <- as.character(baseline_sub[[id_col]])
  diff_mat = diff_mat %>%
    as.data.frame()
  
  # ---- compute geneset fold-change features
  # preallocate result dataframes
  result <-
    data.frame(matrix(nrow = length(rownames(diff_mat)), ncol = (length(geneset_names) + 1)))
  
  # set the rownames
  rownames(result) = baseline_sub[[id_col]]
  
  # set the column names
  colnames(result) = c(id_col, paste0(geneset_names, "_Day", timepoint))
  
  # set the participant ids
  result[[id_col]] <- baseline_sub[[id_col]]
  
  # fill geneset columns with row-wise means across present genes
  for (i in seq_along(genesets)) {
    # For each geneset
    # Which are the genes in this set
    set_genes <- genesets[[i]]
    # Which of these genes are present in our data?
    present <- intersect(set_genes, gene_cols)
    # What is the mean value of these genes' expression
    set_means <- rowMeans(diff_mat[, present, drop = FALSE], na.rm = TRUE)
    # Set this as the result for the baseline matrix
    result[[paste0(geneset_names[i], "_Day", timepoint)]] <- as.numeric(set_means)
  }
  
  rownames(result) <- NULL
  return(result)
}

# Now we can calculate geneset-level features for each timepoint of interest
# The timepoints we consider are vaccine-dependent
# These timepoints were selected based on the (cumulative) number of samples available for each vaccine
# The timepoints of interest are c(0, 1, 3, 7, 10, 14)

# Select a clinical baseline dataframe containing vaccine, age, gender, and race
d0_clinical = hipc_merged_all_norm_response %>%
  select(
    participant_id,
    immResp_MFC_anyAssay_log2_MFC,
    vaccine_name,
    vaccine_colour,
    age_imputed,
    gender
  ) %>%
  distinct()

d0 = compute_geneset_baseline(
  df = hipc_merged_all_norm_response,
  genesets = BTM[["genesets"]],
  geneset_names = BTM[["geneset.names.descriptions"]],
  id_col = "participant_id",
  time_col = "study_time_collected"
)

d0.17 = compute_geneset_fc(
  df = hipc_merged_all_norm_response,
  genesets = BTM[["genesets"]],
  geneset_names = BTM[["geneset.names.descriptions"]],
  id_col = "participant_id",
  time_col = "study_time_collected",
  timepoint = 0.17
)

d1 = compute_geneset_fc(
  df = hipc_merged_all_norm_response,
  genesets = BTM[["genesets"]],
  geneset_names = BTM[["geneset.names.descriptions"]],
  id_col = "participant_id",
  time_col = "study_time_collected",
  timepoint = 1
)

d2 = compute_geneset_fc(
  df = hipc_merged_all_norm_response,
  genesets = BTM[["genesets"]],
  geneset_names = BTM[["geneset.names.descriptions"]],
  id_col = "participant_id",
  time_col = "study_time_collected",
  timepoint = 2
)


d3 = compute_geneset_fc(
  df = hipc_merged_all_norm_response,
  genesets = BTM[["genesets"]],
  geneset_names = BTM[["geneset.names.descriptions"]],
  id_col = "participant_id",
  time_col = "study_time_collected",
  timepoint = 3
)

d7 = compute_geneset_fc(
  df = hipc_merged_all_norm_response,
  genesets = BTM[["genesets"]],
  geneset_names = BTM[["geneset.names.descriptions"]],
  id_col = "participant_id",
  time_col = "study_time_collected",
  timepoint = 7
)

d10 = compute_geneset_fc(
  df = hipc_merged_all_norm_response,
  genesets = BTM[["genesets"]],
  geneset_names = BTM[["geneset.names.descriptions"]],
  id_col = "participant_id",
  time_col = "study_time_collected",
  timepoint = 10
)

d14 = compute_geneset_fc(
  df = hipc_merged_all_norm_response,
  genesets = BTM[["genesets"]],
  geneset_names = BTM[["geneset.names.descriptions"]],
  id_col = "participant_id",
  time_col = "study_time_collected",
  timepoint = 14
)

# For the sequential approach, we save each dataframe per vaccine per timepoint
# containing only the information present at that timepoint

sequential_merge_per_vaccine <- function(timepoints_of_interest, vaccine) {
  # List to store sequential dataframes
  sequential_list <- list()
  
  # Start with baseline clinical data filtered for the selected vaccine
  clinical_df <- d0_clinical %>%
    dplyr::filter(vaccine_name == vaccine)
  
  sequential_list[[1]] <- clinical_df  # first element is baseline
  
  # Keep track of names
  list_names <- c("clinical")
  
  # Loop over each timepoint
  for (tp in timepoints_of_interest) {
    # Construct the dataframe name dynamically, e.g., "d3", "d7"
    df_name <- paste0("d", tp)
    
    # Check if the dataframe exists
    if (!exists(df_name)) {
      warning(paste0(
        "Dataframe ",
        df_name,
        " does not exist. Skipping timepoint ",
        tp,
        "."
      ))
      next
    }
    
    # Get the dataframe object
    tp_df <- get(df_name)
    
    # Keep participant_id first
    tp_df <- tp_df %>%
      dplyr::select(participant_id, dplyr::everything())
    
    # Cumulative merge: keep only participants present in both
    current_df <- dplyr::inner_join(clinical_df, tp_df, by = "participant_id")
    
    # Append to the list
    sequential_list[[length(sequential_list) + 1]] <- current_df
    
    # Append name
    list_names <- c(list_names, paste0("Day ", tp))
  }
  
  # Set names for the list
  names(sequential_list) <- list_names
  
  return(sequential_list)
}


# For the cumulative approach, we have to merge these data together in a cumulative manner
# That is, we should have one dataframe per vaccine per timepoint.
# This dataframe should contain all of the features up to a given timepoint.
# We have to do this vaccine-by-vaccine since not all vaccines share the same desired timepoints.
# Here is a helper function to derive cumulatively merged dataframes for a given vaccine
# It outputs a list for a given vaccine containing dataframes at each timepoint

cumulative_merge_per_vaccine <- function(timepoints_of_interest, vaccine) {
  # List to store cumulative dataframes
  cumulative_list <- list()
  
  # Start with baseline clinical data filtered for the selected vaccine
  current_df <- d0_clinical %>%
    dplyr::filter(vaccine_name == vaccine)
  
  cumulative_list[[1]] <- current_df  # first element is baseline
  
  # Keep track of names
  list_names <- c("clinical")
  
  # Loop over each timepoint
  for (tp in timepoints_of_interest) {
    # Construct the dataframe name dynamically, e.g., "d3", "d7"
    df_name <- paste0("d", tp)
    
    # Check if the dataframe exists
    if (!exists(df_name)) {
      warning(paste0(
        "Dataframe ",
        df_name,
        " does not exist. Skipping timepoint ",
        tp,
        "."
      ))
      next
    }
    
    # Get the dataframe object
    tp_df <- get(df_name)
    
    # Keep participant_id first
    tp_df <- tp_df %>%
      dplyr::select(participant_id, dplyr::everything())
    
    # Cumulative merge: keep only participants present in both
    current_df <- dplyr::inner_join(current_df, tp_df, by = "participant_id")
    
    # Append to the list
    cumulative_list[[length(cumulative_list) + 1]] <- current_df
    
    # Append name
    list_names <- c(list_names, paste0("Day ", tp))
  }
  
  # Set names for the list
  names(cumulative_list) <- list_names
  
  return(cumulative_list)
}

# Specify the vaccines of interst and the corresponding timepoints for each in the
# sequential and cumulative approaches.
# These vaccines were chosen based on sample sizes and heterogeneity in immune response values
specified_timepoints_list_sequential <- list(
  "Meningococcus (PS)" = c(0, 3, 7),
  "Meningococcus (CJ)" = c(0, 3, 7),
  "Influenza (IN)" = c(0, 1, 3, 7, 14),
  "Hepatitis A/B (IN/RP)" = c(0, 7),
  "Yellow Fever (LV)" = c(0, 0.17, 3, 7, 10, 14)
)

specified_timepoints_list_cumulative <- list(
  "Meningococcus (PS)" = c(0, 3, 7),
  "Meningococcus (CJ)" = c(0, 3, 7),
  "Influenza (IN)" = c(0, 1, 3, 7, 14),
  "Hepatitis A/B (IN/RP)" = c(0, 7),
  "Yellow Fever (LV)" = c(0, 3, 7, 10, 14)
)

# Now apply the sequential merge function to each 
sequential_prediction_sets_list <- lapply(names(specified_timepoints_list_sequential), function(vac) {
  timepoints <- specified_timepoints_list_sequential[[vac]]
  sequential_merge_per_vaccine(timepoints_of_interest = timepoints, vaccine = vac)
})

# Name the top-level list elements with the vaccine names
names(sequential_prediction_sets_list) <- names(specified_timepoints_list_cumulative)

# Specify the path to save the list
p_save_prediction_sets_list <- fs::path(processed_data_folder, "sequential_prediction_sets_list.rds")

# Save the data
saveRDS(sequential_prediction_sets_list, p_save_prediction_sets_list)

# Now apply the cumulative merge function to each
cumulative_prediction_sets_list <- lapply(names(specified_timepoints_list_cumulative), function(vac) {
  timepoints <- specified_timepoints_list_cumulative[[vac]]
  cumulative_merge_per_vaccine(timepoints_of_interest = timepoints, vaccine = vac)
})

# Name the top-level list elements with the vaccine names
names(cumulative_prediction_sets_list) <- names(specified_timepoints_list_cumulative)

# Specify the path to save the list
p_save_prediction_sets_list <- fs::path(processed_data_folder, "cumulative_prediction_sets_list.rds")

# Save the data
saveRDS(cumulative_prediction_sets_list, p_save_prediction_sets_list)


# Now we do the same but for predictor sets without NA modules 
# First find the genesets to remove
# Load genesets
BTM_withoutTBA = readRDS(p_load_btm)

idx <- which(BTM_withoutTBA[["geneset.descriptions"]] == "TBA")

BTM_withoutTBA[["genesets"]] <- BTM_withoutTBA[["genesets"]][-idx]
BTM_withoutTBA[["geneset.descriptions"]] <- BTM_withoutTBA[["geneset.descriptions"]][-idx]
BTM_withoutTBA[["geneset.names"]] <- BTM_withoutTBA[["geneset.names"]][-idx]
BTM_withoutTBA[["geneset.aggregates"]] = BTM_withoutTBA[["geneset.aggregates"]][-idx]
BTM_withoutTBA[["geneset.names.descriptions"]] = BTM_withoutTBA[["geneset.names.descriptions"]][-idx]

# Select a clinical baseline dataframe containing vaccine, age, gender, and race
d0_clinical = hipc_merged_all_norm_response %>%
  select(
    participant_id,
    immResp_MFC_anyAssay_log2_MFC,
    vaccine_name,
    vaccine_colour,
    age_imputed,
    gender
  ) %>%
  distinct()

d0 = compute_geneset_baseline(
  df = hipc_merged_all_norm_response,
  genesets = BTM_withoutTBA[["genesets"]],
  geneset_names = BTM_withoutTBA[["geneset.names.descriptions"]],
  id_col = "participant_id",
  time_col = "study_time_collected"
)

d0.17 = compute_geneset_fc(
  df = hipc_merged_all_norm_response,
  genesets = BTM_withoutTBA[["genesets"]],
  geneset_names = BTM_withoutTBA[["geneset.names.descriptions"]],
  id_col = "participant_id",
  time_col = "study_time_collected",
  timepoint = 0.17
)

d1 = compute_geneset_fc(
  df = hipc_merged_all_norm_response,
  genesets = BTM_withoutTBA[["genesets"]],
  geneset_names = BTM_withoutTBA[["geneset.names.descriptions"]],
  id_col = "participant_id",
  time_col = "study_time_collected",
  timepoint = 1
)

d2 = compute_geneset_fc(
  df = hipc_merged_all_norm_response,
  genesets = BTM_withoutTBA[["genesets"]],
  geneset_names = BTM_withoutTBA[["geneset.names.descriptions"]],
  id_col = "participant_id",
  time_col = "study_time_collected",
  timepoint = 2
)


d3 = compute_geneset_fc(
  df = hipc_merged_all_norm_response,
  genesets = BTM_withoutTBA[["genesets"]],
  geneset_names = BTM_withoutTBA[["geneset.names.descriptions"]],
  id_col = "participant_id",
  time_col = "study_time_collected",
  timepoint = 3
)

d7 = compute_geneset_fc(
  df = hipc_merged_all_norm_response,
  genesets = BTM_withoutTBA[["genesets"]],
  geneset_names = BTM_withoutTBA[["geneset.names.descriptions"]],
  id_col = "participant_id",
  time_col = "study_time_collected",
  timepoint = 7
)

d10 = compute_geneset_fc(
  df = hipc_merged_all_norm_response,
  genesets = BTM_withoutTBA[["genesets"]],
  geneset_names = BTM_withoutTBA[["geneset.names.descriptions"]],
  id_col = "participant_id",
  time_col = "study_time_collected",
  timepoint = 10
)

d14 = compute_geneset_fc(
  df = hipc_merged_all_norm_response,
  genesets = BTM_withoutTBA[["genesets"]],
  geneset_names = BTM_withoutTBA[["geneset.names.descriptions"]],
  id_col = "participant_id",
  time_col = "study_time_collected",
  timepoint = 14
)

# Now apply the sequential merge function to each 
sequential_prediction_sets_list_withoutTBA <- lapply(names(specified_timepoints_list_sequential), function(vac) {
  timepoints <- specified_timepoints_list_sequential[[vac]]
  sequential_merge_per_vaccine(timepoints_of_interest = timepoints, vaccine = vac)
})

# Name the top-level list elements with the vaccine names
names(sequential_prediction_sets_list_withoutTBA) <- names(specified_timepoints_list_cumulative)

# Specify the path to save the list
p_save_prediction_sets_list_withoutTBA <- fs::path(processed_data_folder, "sequential_prediction_sets_list_withoutTBA.rds")

# Save the data
saveRDS(sequential_prediction_sets_list_withoutTBA, p_save_prediction_sets_list_withoutTBA)

# Now apply the cumulative merge function to each
cumulative_prediction_sets_list_withoutTBA <- lapply(names(specified_timepoints_list_cumulative), function(vac) {
  timepoints <- specified_timepoints_list_cumulative[[vac]]
  cumulative_merge_per_vaccine(timepoints_of_interest = timepoints, vaccine = vac)
})

# Name the top-level list elements with the vaccine names
names(cumulative_prediction_sets_list_withoutTBA) <- names(specified_timepoints_list_cumulative)

# Specify the path to save the list
p_save_prediction_sets_list_withoutTBA <- fs::path(processed_data_folder, "cumulative_prediction_sets_list_withoutTBA.rds")

# Save the data
saveRDS(cumulative_prediction_sets_list_withoutTBA, p_save_prediction_sets_list_withoutTBA)
