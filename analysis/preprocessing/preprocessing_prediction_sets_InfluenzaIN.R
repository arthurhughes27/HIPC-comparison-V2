# R script to engineer Influenza expression data into gene-set level features
# - read inputs, define timepoints, and filter to relevant samples

library(dplyr)
library(GSVA)
library(fs)   # used for path construction (fs::path used below)

# Directory containing engineered / processed data files
processed_data_folder <- "data"

# Paths to processed gene-level data and gene-set objects
p_load_expr_all_norm <- fs::path(processed_data_folder, "hipc_merged_all_norm.rds")
p_load_btm           <- fs::path(processed_data_folder, "BTM_processed.rds")

# Load data objects
hipc_merged_all_norm <- readRDS(p_load_expr_all_norm)
BTM                  <- readRDS(p_load_btm)

# Timepoints of interest (numeric)
timepoints_of_interest <- c(0, 3, 7, 14)

# Filter to samples with non-missing immune response, Influenza vaccine,
# and collected at one of the specified timepoints.
hipc_merged_all_norm_filtered <- hipc_merged_all_norm %>%
  filter(
    !is.na(immResp_MFC_anyAssay_log2_MFC),
    vaccine_name == "Influenza (IN)",
    study_time_collected %in% timepoints_of_interest
  )

gene_names = hipc_merged_all_norm_filtered %>% 
  dplyr::select(a1cf:zzz3) %>% 
  colnames()

# df = hipc_merged_all_norm_filtered
# genesets = BTM[["genesets"]]
# geneset_names = BTM[["geneset.names.descriptions"]]
# id_col = "participant_id"
# time_col = "study_time_collected"
# transformation = "ssgsea"
# timepoint = 7

compute_geneset_baseline_transformation <- function(df,
                                                    genesets,
                                                    geneset_names,
                                                    id_col,
                                                    time_col,
                                                    transformation = c("mean",
                                                                       "median",
                                                                       "max",
                                                                       "mean-rank",
                                                                       "median-rank",
                                                                       "max-rank",
                                                                       "pc1",
                                                                       "ssgsea")) {
  # ---- basic checks --------------------------------------------------------
  if (!is.data.frame(df)) stop("df must be a data.frame")
  if (!is.list(genesets)) stop("genesets must be a list (each element: vector of gene names)")
  if (length(genesets) != length(geneset_names)) stop("genesets and geneset_names must have same length")
  if (!(id_col %in% names(df))) stop("id_col not found in df")
  if (!(time_col %in% names(df))) stop("time_col not found in df")
  if (!exists("gene_names")) stop("gene_names vector not found in environment (required).")
  
  # pick exactly one transformation (validates the input)
  transformation <- match.arg(transformation)
  split_strings <- strsplit(transformation, "-")[[1]]
  
  # detect column-wise transformation
  if ("rank" %in% split_strings) {
    col_transformation <- "rank"
  } else if ("pc1" %in% split_strings) {
    col_transformation <- "pc1"
  } else if ("ssgsea" %in% split_strings) {
    col_transformation <- "ssgsea"
  } else {
    col_transformation <- NA
  }
  
  # row-wise aggregator
  if (split_strings[1] == "mean") {
    row_transformation <- function(x) mean(x, na.rm = TRUE)
  } else if (split_strings[1] == "median") {
    row_transformation <- function(x) median(x, na.rm = TRUE)
  } else if (split_strings[1] == "max") {
    row_transformation <- function(x) max(x, na.rm = TRUE)
  }
  
  # ---- coerce time column -------------------------------------------------
  df[[time_col]] <- as.numeric(df[[time_col]])
  if (any(is.na(df[[time_col]]))) warning("Some time_col values could not be coerced to numeric (NAs present).")
  
  # ---- determine gene columns from provided gene_names ---------------------
  gene_cols <- intersect(gene_names, colnames(df))
  if (length(gene_cols) == 0) stop("No gene columns found in df matching gene_names.")
  
  # ---- select most recent baseline (<= 0) per participant ------------------
  baseline_candidates <- df[!is.na(df[[time_col]]) & df[[time_col]] <= 0, , drop = FALSE]
  if (nrow(baseline_candidates) == 0) stop("No baseline (<= 0) rows found in df.")
  ord_bl <- order(baseline_candidates[[id_col]], -baseline_candidates[[time_col]])
  bl_sorted <- baseline_candidates[ord_bl, , drop = FALSE]
  baseline_df <- bl_sorted[!duplicated(bl_sorted[[id_col]]), , drop = FALSE]
  baseline_df <- baseline_df[, c(id_col, gene_cols), drop = FALSE]
  
  # numeric matrix: rows = samples, cols = genes
  bl_mat <- as.matrix(baseline_df[, gene_cols, drop = FALSE])
  storage.mode(bl_mat) <- "numeric"
  
  # prepare result skeleton
  out_colnames <- paste0(geneset_names, "_Day0")
  result <- data.frame(matrix(NA_real_, nrow = nrow(baseline_df), ncol = length(out_colnames)),
                       stringsAsFactors = FALSE)
  colnames(result) <- out_colnames
  result[[id_col]] <- baseline_df[[id_col]]
  result <- result[, c(id_col, out_colnames), drop = FALSE]
  
  # ---- ssGSEA branch (use your ssgseaParam + gsva(params) pattern) ---------
  if (!is.na(col_transformation) && col_transformation == "ssgsea") {
    if (!requireNamespace("GSVA", quietly = TRUE)) {
      stop("Package 'GSVA' is required for ssgsea transformation. Please install it.")
    }
    # prepare gene set list and names from inputs
    BTM_list <- genesets
    names(BTM_list) <- geneset_names
    # use ssgseaParam + gsva(params) pattern (matches your current GSVA API)
    params <- ssgseaParam(
      exprData = t(bl_mat),    # genes x samples
      geneSets = BTM_list,
      minSize = 1,
      maxSize = Inf,
      normalize = TRUE
    )
    ssgsea_results <- t(gsva(params))    # samples x genesets
    colnames(ssgsea_results) <- paste0(colnames(ssgsea_results), "_Day0")
    ssgsea_df <- as.data.frame(ssgsea_results, stringsAsFactors = FALSE)
    ssgsea_df <- cbind(id = baseline_df[[id_col]], ssgsea_df)
    names(ssgsea_df)[1] <- id_col
    # reorder columns to match result and return
    result <- ssgsea_df[, colnames(result), drop = FALSE]
    rownames(result) <- NULL
    return(result)
  }
  
  # ---- compute geneset-level features for each geneset ---------------------
  for (i in seq_along(genesets)) {
    set_genes <- genesets[[i]]
    present <- intersect(set_genes, gene_cols)
    if (length(present) == 0) {
      result[[out_colnames[i]]] <- NA_real_
      next
    }
    bl_mat_filtered <- bl_mat[, present, drop = FALSE] # samples x present_genes
    
    # apply transformation
    if (is.na(col_transformation)) {
      aggregated <- apply(bl_mat_filtered, MARGIN = 1, FUN = row_transformation)
    } else if (col_transformation == "rank") {
      col_ranked <- apply(bl_mat_filtered, 2, function(col) rank(-col, ties.method = "average"))
      aggregated <- apply(col_ranked, MARGIN = 1, FUN = row_transformation)
    } else if (col_transformation == "pc1") {
      if (ncol(bl_mat_filtered) == 1) {
        aggregated <- as.numeric(scale(bl_mat_filtered)[, 1])
      } else {
        pca_res <- stats::prcomp(bl_mat_filtered, center = TRUE, scale. = TRUE)
        aggregated <- as.numeric(pca_res$x[, 1])
      }
    }
    
    result[[out_colnames[i]]] <- as.numeric(aggregated)
  }
  
  rownames(result) <- NULL
  return(result)
}



# The following is a helper function which computes gene-set post-vaccination fold-change transformation at a given timepoint.
compute_geneset_postvax_fc_transformation <- function(df,
                                                      genesets,
                                                      geneset_names,
                                                      id_col,
                                                      time_col,
                                                      timepoint,
                                                      transformation = c("mean",
                                                                         "median",
                                                                         "max",
                                                                         "mean-rank",
                                                                         "median-rank",
                                                                         "max-rank",
                                                                         "pc1",
                                                                         "ssgsea")) {
  # ---- input checks --------------------------------------------------------
  if (!is.data.frame(df)) stop("df must be a data.frame")
  if (!is.list(genesets)) stop("genesets must be a list (each element: vector of gene names)")
  if (length(genesets) != length(geneset_names)) stop("genesets and geneset_names must have same length")
  if (!(id_col %in% names(df))) stop("id_col not found in df")
  if (!(time_col %in% names(df))) stop("time_col not found in df")
  if (!exists("gene_names")) stop("gene_names vector not found in environment (required).")
  
  # pick exactly one transformation value
  transformation <- match.arg(transformation)
  split_strings <- strsplit(transformation, "-")[[1]]
  
  # detect column-wise transformation
  if ("rank" %in% split_strings) {
    col_transformation <- "rank"
  } else if ("pc1" %in% split_strings) {
    col_transformation <- "pc1"
  } else if ("ssgsea" %in% split_strings) {
    col_transformation <- "ssgsea"
  } else {
    col_transformation <- NA
  }
  
  # row-wise aggregator
  if (split_strings[1] == "mean") {
    row_transformation <- function(x) mean(x, na.rm = TRUE)
  } else if (split_strings[1] == "median") {
    row_transformation <- function(x) median(x, na.rm = TRUE)
  } else if (split_strings[1] == "max") {
    row_transformation <- function(x) max(x, na.rm = TRUE)
  }
  
  # ---- time coercion ------------------------------------------------------
  df[[time_col]] <- as.numeric(df[[time_col]])
  if (any(is.na(df[[time_col]]))) warning("Some time_col values could not be coerced to numeric (NAs present).")
  
  # ---- gene columns via provided gene_names --------------------------------
  # use gene_names (defined earlier) to determine gene columns in df
  gene_cols <- intersect(gene_names, colnames(df))
  if (length(gene_cols) == 0) stop("No gene columns found in df matching gene_names.")
  
  # ---- choose most recent baseline (<= 0) per participant ------------------
  baseline_candidates <- df[!is.na(df[[time_col]]) & df[[time_col]] <= 0, , drop = FALSE]
  if (nrow(baseline_candidates) == 0) stop("No baseline (<= 0) rows found in df.")
  ord_bl <- order(baseline_candidates[[id_col]], -baseline_candidates[[time_col]])
  bl_sorted <- baseline_candidates[ord_bl, , drop = FALSE]
  baseline_df <- bl_sorted[!duplicated(bl_sorted[[id_col]]), , drop = FALSE]
  baseline_df <- baseline_df[, c(id_col, gene_cols), drop = FALSE]
  
  # ---- select post-vax rows at specified timepoint ------------------------
  tp_candidates <- df[!is.na(df[[time_col]]) & df[[time_col]] == timepoint, , drop = FALSE]
  if (nrow(tp_candidates) == 0) stop("No rows found for requested timepoint: ", timepoint)
  tp_df <- tp_candidates[, c(id_col, gene_cols), drop = FALSE]
  
  # ---- keep participants with both baseline AND post-vax -------------------
  bl_ids <- as.character(baseline_df[[id_col]])
  tp_ids <- as.character(tp_df[[id_col]])
  common_ids <- intersect(bl_ids, tp_ids)   # order follows baseline_df
  if (length(common_ids) == 0) stop("No participants have both baseline and post-vax at timepoint ", timepoint)
  
  baseline_sub <- baseline_df[match(common_ids, as.character(baseline_df[[id_col]])), , drop = FALSE]
  tp_sub       <- tp_df[match(common_ids, as.character(tp_df[[id_col]])), , drop = FALSE]
  
  # ---- compute difference matrix (post - baseline) -------------------------
  bl_mat <- as.matrix(baseline_sub[, gene_cols, drop = FALSE])
  tp_mat <- as.matrix(tp_sub[, gene_cols, drop = FALSE])
  storage.mode(bl_mat) <- "numeric"
  storage.mode(tp_mat) <- "numeric"
  diff_mat <- tp_mat - bl_mat
  rownames(diff_mat) <- as.character(baseline_sub[[id_col]])
  # keep diff_mat as numeric matrix for transformations
  
  # ---- prepare result skeleton ---------------------------------------------
  out_colnames <- paste0(geneset_names, "_Day", timepoint)
  result <- data.frame(matrix(NA_real_, nrow = nrow(diff_mat), ncol = length(out_colnames)),
                       stringsAsFactors = FALSE)
  colnames(result) <- out_colnames
  result[[id_col]] <- baseline_sub[[id_col]]
  result <- result[, c(id_col, out_colnames), drop = FALSE]
  
  # ---- ssGSEA branch (use your ssgseaParam + gsva(params) pattern) ---------
  if (!is.na(col_transformation) && col_transformation == "ssgsea") {
    if (!requireNamespace("GSVA", quietly = TRUE)) {
      stop("Package 'GSVA' is required for ssgsea transformation. Please install it.")
    }
    # prepare gene set list (named)
    BTM_list <- genesets
    names(BTM_list) <- geneset_names
    # note: using the ssgseaParam + gsva(params) style you prefer
    params <- ssgseaParam(
      exprData = t(diff_mat),   # genes x samples
      geneSets = BTM_list,
      minSize = 1,
      maxSize = Inf,
      normalize = TRUE
    )
    ssgsea_results <- t(gsva(params))   # keep same structure as your original code
    colnames(ssgsea_results) <- paste0(colnames(ssgsea_results), "_Day", timepoint)
    # bind ids (ensure naming consistent with id_col)
    ssgsea_df <- as.data.frame(ssgsea_results, stringsAsFactors = FALSE)
    ssgsea_df <- cbind(id = rownames(diff_mat), ssgsea_df)
    names(ssgsea_df)[1] <- id_col
    # reorder columns to match result
    result <- ssgsea_df[, colnames(result), drop = FALSE]
    rownames(result) <- NULL
    return(result)
  }
  
  # ---- compute geneset-level features for each geneset ---------------------
  for (i in seq_along(genesets)) {
    set_genes <- genesets[[i]]
    present <- intersect(set_genes, colnames(diff_mat))
    if (length(present) == 0) {
      result[[out_colnames[i]]] <- NA_real_
      next
    }
    diff_mat_filtered <- diff_mat[, present, drop = FALSE]  # samples x genes
    
    # apply requested transformations
    if (is.na(col_transformation)) {
      aggregated <- apply(diff_mat_filtered, MARGIN = 1, FUN = row_transformation)
    } else if (col_transformation == "rank") {
      col_ranked <- apply(diff_mat_filtered, 2, function(col) rank(-col, ties.method = "average"))
      aggregated <- apply(col_ranked, MARGIN = 1, FUN = row_transformation)
    } else if (col_transformation == "pc1") {
      if (ncol(diff_mat_filtered) == 1) {
        aggregated <- as.numeric(scale(diff_mat_filtered)[, 1])
      } else {
        pca_res <- stats::prcomp(diff_mat_filtered, center = TRUE, scale. = TRUE)
        aggregated <- as.numeric(pca_res$x[, 1])
      }
    }
    
    result[[out_colnames[i]]] <- as.numeric(aggregated)
  }
  
  rownames(result) <- NULL
  return(result)
}


# Now we can calculate geneset-level features for each timepoint of interest

transformations_of_interest = c("mean",
                                "median" ,
                                "max",
                                "mean-rank",
                                "median-rank",
                                "max-rank",
                                "pc1",
                                "ssgsea")

timepoints_of_interest <- c(0, 3, 7, 14)


# First, select a clinical baseline dataframe containing vaccine, age, gender
d0_clinical <- hipc_merged_all_norm_filtered %>%
  dplyr::select(
    participant_id,
    immResp_MFC_anyAssay_log2_MFC,
    vaccine_name,
    vaccine_colour,
    age_imputed,
    gender
  ) %>%
  dplyr::distinct()

d0 <- lapply(transformations_of_interest, FUN = function(trans) {
  compute_geneset_baseline_transformation(
    df = hipc_merged_all_norm_filtered,
    genesets = BTM[["genesets"]],
    geneset_names = BTM[["geneset.names.descriptions"]],
    id_col = "participant_id",
    time_col = "study_time_collected",
    transformation = trans
  )
})

# name the list elements by the transformation used
names(d0) <- transformations_of_interest

d3 = lapply(
  transformations_of_interest,
  FUN = function(trans) {
    compute_geneset_postvax_fc_transformation(
      df = hipc_merged_all_norm_filtered,
      timepoint = 3,
      genesets = BTM[["genesets"]],
      geneset_names = BTM[["geneset.names.descriptions"]],
      id_col = "participant_id",
      time_col = "study_time_collected",
      transformation = trans
    )
  }
)

names(d3) = transformations_of_interest


d7 = lapply(
  transformations_of_interest,
  FUN = function(trans) {
    compute_geneset_postvax_fc_transformation(
      df = hipc_merged_all_norm_filtered,
      timepoint = 7,
      genesets = BTM[["genesets"]],
      geneset_names = BTM[["geneset.names.descriptions"]],
      id_col = "participant_id",
      time_col = "study_time_collected",
      transformation = trans
    )
  }
)

names(d7) = transformations_of_interest

d14 = lapply(
  transformations_of_interest,
  FUN = function(trans) {
    compute_geneset_postvax_fc_transformation(
      df = hipc_merged_all_norm_filtered,
      timepoint = 14,
      genesets = BTM[["genesets"]],
      geneset_names = BTM[["geneset.names.descriptions"]],
      id_col = "participant_id",
      time_col = "study_time_collected",
      transformation = trans
    )
  }
)

names(d14) = transformations_of_interest


# --- organize clinical and geneset data into a nested list -------------------
# Top-level: timepoints; second-level: transformations
# Each dataframe combines baseline clinical info with geneset features

# baseline clinical dataframe
clinical_df <- d0_clinical

# initialize main list with clinical info as first element
sequential_list <- list(clinical = list(clinical = clinical_df))

# loop over timepoints of interest
for (tp in timepoints_of_interest) {
  time_name <- paste0("Day ", tp)
  time_numeric = paste0("d", tp)
  
  # retrieve the list of transformations for this timepoint
  # NOTE: ensure that objects d0, d3, d7, d14 exist in environment
  tp_df_list <- get(time_numeric)
  
  # create a sublist for this timepoint
  sequential_list[[time_name]] <- list()
  
  # loop over transformations
  for (trans in transformations_of_interest) {
    tp_df <- tp_df_list[[trans]]
    
    # reorder columns: participant_id first
    tp_df <- tp_df %>% dplyr::select(participant_id, dplyr::everything())
    
    # merge baseline clinical info with geneset features
    current_df <- dplyr::inner_join(clinical_df, tp_df, by = "participant_id")
    
    # save into nested list
    sequential_list[[time_name]][[trans]] <- current_df
  }
}


# --- cumulative (information accrual) merging -------------------------------
# Aim: produce one dataframe per timepoint per transformation that contains
#      all features available up to that timepoint.

# ensure timepoints are sorted ascending
timepoints_of_interest <- sort(timepoints_of_interest)

# check presence of clinical dataframe and join key
if (!exists("clinical_df")) stop("clinical_df not found in environment.")
if (!("participant_id" %in% colnames(clinical_df))) stop("'participant_id' column missing from clinical_df.")

# initialize result and cumulative storage
cumulative_list <- list(clinical = list(clinical = clinical_df))
# create a named list of cumulative dataframes (one per transformation),
# starting from the clinical baseline
cumulative_by_trans <- setNames(
  lapply(transformations_of_interest, function(x) clinical_df),
  transformations_of_interest
)

# iterate through timepoints and accumulate features
for (tp in timepoints_of_interest) {
  time_name <- paste0("Day ", tp)
  time_numeric = paste0("d", tp)
  
  # attempt to retrieve the object named "d{tp}" (may be NULL if not created)
  tp_df_list <- tryCatch(get(time_numeric), error = function(e) NULL)
  
  # create container for this timepoint in the cumulative list
  cumulative_list[[time_name]] <- list()
  
  for (trans in transformations_of_interest) {
    # if no data available for this timepoint/transformation, keep the running cumulative df
    if (is.null(tp_df_list) || is.null(tp_df_list[[trans]])) {
      cumulative_list[[time_name]][[trans]] <- cumulative_by_trans[[trans]]
      next
    }
    
    # ensure the tp dataframe has participant_id and place that column first
    tp_df <- tp_df_list[[trans]] %>%
      dplyr::select(participant_id, dplyr::everything())
    if (!("participant_id" %in% colnames(tp_df))) {
      stop("participant_id not found in data for ", time_name, " / ", trans)
    }
    
    # perform cumulative join using inner join
    cumulative_by_trans[[trans]] <- dplyr::inner_join(cumulative_by_trans[[trans]], tp_df, by = "participant_id")
    
    # store the snapshot for this timepoint/transformation
    cumulative_list[[time_name]][[trans]] <- cumulative_by_trans[[trans]]
  }
}

# Specify the path to save the lists
p_save_sequential_list_influenzain_withTBA <- fs::path(processed_data_folder, "sequential_list_influenzain_withTBA.rds")
p_save_cumulative_list_influenzain_withTBA <- fs::path(processed_data_folder, "cumulative_list_influenzain_withTBA.rds")

# Save the data
saveRDS(sequential_list, p_save_sequential_list_influenzain_withTBA)
saveRDS(cumulative_list, p_save_cumulative_list_influenzain_withTBA)






# Now we do the same but for predictor sets without TBA modules 
# Reload genesets
BTM_withoutTBA = readRDS(p_load_btm)

idx <- which(BTM_withoutTBA[["geneset.descriptions"]] == "TBA")

BTM_withoutTBA[["genesets"]] <- BTM_withoutTBA[["genesets"]][-idx]
BTM_withoutTBA[["geneset.descriptions"]] <- BTM_withoutTBA[["geneset.descriptions"]][-idx]
BTM_withoutTBA[["geneset.names"]] <- BTM_withoutTBA[["geneset.names"]][-idx]
BTM_withoutTBA[["geneset.aggregates"]] = BTM_withoutTBA[["geneset.aggregates"]][-idx]
BTM_withoutTBA[["geneset.names.descriptions"]] = BTM_withoutTBA[["geneset.names.descriptions"]][-idx]

transformations_of_interest = c("mean",
                                "median" ,
                                "max",
                                "mean-rank",
                                "median-rank",
                                "max-rank",
                                "pc1",
                                "ssgsea")

timepoints_of_interest <- c(0, 3, 7, 14)


# First, select a clinical baseline dataframe containing vaccine, age, gender
d0_clinical <- hipc_merged_all_norm_filtered %>%
  dplyr::select(
    participant_id,
    immResp_MFC_anyAssay_log2_MFC,
    vaccine_name,
    vaccine_colour,
    age_imputed,
    gender
  ) %>%
  dplyr::distinct()

d0 <- lapply(transformations_of_interest, FUN = function(trans) {
  compute_geneset_baseline_transformation(
    df = hipc_merged_all_norm_filtered,
    genesets = BTM_withoutTBA[["genesets"]],
    geneset_names = BTM_withoutTBA[["geneset.names.descriptions"]],
    id_col = "participant_id",
    time_col = "study_time_collected",
    transformation = trans
  )
})

# name the list elements by the transformation used
names(d0) <- transformations_of_interest

d3 = lapply(
  transformations_of_interest,
  FUN = function(trans) {
    compute_geneset_postvax_fc_transformation(
      df = hipc_merged_all_norm_filtered,
      timepoint = 3,
      genesets = BTM_withoutTBA[["genesets"]],
      geneset_names = BTM_withoutTBA[["geneset.names.descriptions"]],
      id_col = "participant_id",
      time_col = "study_time_collected",
      transformation = trans
    )
  }
)

names(d3) = transformations_of_interest


d7 = lapply(
  transformations_of_interest,
  FUN = function(trans) {
    compute_geneset_postvax_fc_transformation(
      df = hipc_merged_all_norm_filtered,
      timepoint = 7,
      genesets = BTM_withoutTBA[["genesets"]],
      geneset_names = BTM_withoutTBA[["geneset.names.descriptions"]],
      id_col = "participant_id",
      time_col = "study_time_collected",
      transformation = trans
    )
  }
)

names(d7) = transformations_of_interest

d14 = lapply(
  transformations_of_interest,
  FUN = function(trans) {
    compute_geneset_postvax_fc_transformation(
      df = hipc_merged_all_norm_filtered,
      timepoint = 14,
      genesets = BTM_withoutTBA[["genesets"]],
      geneset_names = BTM_withoutTBA[["geneset.names.descriptions"]],
      id_col = "participant_id",
      time_col = "study_time_collected",
      transformation = trans
    )
  }
)

names(d14) = transformations_of_interest


# --- organize clinical and geneset data into a nested list -------------------
# Top-level: timepoints; second-level: transformations
# Each dataframe combines baseline clinical info with geneset features

# baseline clinical dataframe
clinical_df <- d0_clinical

# initialize main list with clinical info as first element
sequential_list <- list(clinical = list(clinical = clinical_df))

# loop over timepoints of interest
for (tp in timepoints_of_interest) {
  time_name <- paste0("Day ", tp)
  time_numeric = paste0("d", tp)
  
  # retrieve the list of transformations for this timepoint
  # NOTE: ensure that objects d0, d3, d7, d14 exist in environment
  tp_df_list <- get(time_numeric)
  
  # create a sublist for this timepoint
  sequential_list[[time_name]] <- list()
  
  # loop over transformations
  for (trans in transformations_of_interest) {
    tp_df <- tp_df_list[[trans]]
    
    # reorder columns: participant_id first
    tp_df <- tp_df %>% dplyr::select(participant_id, dplyr::everything())
    
    # merge baseline clinical info with geneset features
    current_df <- dplyr::inner_join(clinical_df, tp_df, by = "participant_id")
    
    # save into nested list
    sequential_list[[time_name]][[trans]] <- current_df
  }
}


# --- cumulative (information accrual) merging -------------------------------
# Aim: produce one dataframe per timepoint per transformation that contains
#      all features available up to that timepoint.

# ensure timepoints are sorted ascending
timepoints_of_interest <- sort(timepoints_of_interest)

# check presence of clinical dataframe and join key
if (!exists("clinical_df")) stop("clinical_df not found in environment.")
if (!("participant_id" %in% colnames(clinical_df))) stop("'participant_id' column missing from clinical_df.")

# initialize result and cumulative storage
cumulative_list <- list(clinical = list(clinical = clinical_df))
# create a named list of cumulative dataframes (one per transformation),
# starting from the clinical baseline
cumulative_by_trans <- setNames(
  lapply(transformations_of_interest, function(x) clinical_df),
  transformations_of_interest
)

# iterate through timepoints and accumulate features
for (tp in timepoints_of_interest) {
  time_name <- paste0("Day ", tp)
  time_numeric = paste0("d", tp)
  
  # attempt to retrieve the object named "d{tp}" (may be NULL if not created)
  tp_df_list <- tryCatch(get(time_numeric), error = function(e) NULL)
  
  # create container for this timepoint in the cumulative list
  cumulative_list[[time_name]] <- list()
  
  for (trans in transformations_of_interest) {
    # if no data available for this timepoint/transformation, keep the running cumulative df
    if (is.null(tp_df_list) || is.null(tp_df_list[[trans]])) {
      cumulative_list[[time_name]][[trans]] <- cumulative_by_trans[[trans]]
      next
    }
    
    # ensure the tp dataframe has participant_id and place that column first
    tp_df <- tp_df_list[[trans]] %>%
      dplyr::select(participant_id, dplyr::everything())
    if (!("participant_id" %in% colnames(tp_df))) {
      stop("participant_id not found in data for ", time_name, " / ", trans)
    }
    
    # perform cumulative join using inner join
    cumulative_by_trans[[trans]] <- dplyr::inner_join(cumulative_by_trans[[trans]], tp_df, by = "participant_id")
    
    # store the snapshot for this timepoint/transformation
    cumulative_list[[time_name]][[trans]] <- cumulative_by_trans[[trans]]
  }
}

# Specify the path to save the lists
p_save_sequential_list_influenzain_withoutTBA <- fs::path(processed_data_folder, "sequential_list_influenzain_withoutTBA.rds")
p_save_cumulative_list_influenzain_withoutTBA <- fs::path(processed_data_folder, "cumulative_list_influenzain_withoutTBA.rds")

# Save the data
saveRDS(sequential_list, p_save_sequential_list_influenzain_withoutTBA)
saveRDS(cumulative_list, p_save_cumulative_list_influenzain_withoutTBA)