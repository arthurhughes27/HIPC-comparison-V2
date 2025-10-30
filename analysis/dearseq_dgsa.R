# File to perform differential gene-set analysis with dearseq across multiple vaccines and preprocess the results
## We perform one differential gene-set analysis per vaccine (pathogen + vaccine type)
## We adjust for the effects of age, sex, and study in the linear mixed model framework of dearseq
## To make the comparison with Hagan et al.'s results as fair as possible, we use the young, non-cross-study-normalised dataset

# Load the necessary packages
library(fs)
library(tidyverse)
library(dearseq)

# Specify folder within folder root where the raw data lives
processed_data_folder = "data"

# Use fs::path() to specify the data paths robustly
p_load_hipc_merged_young_noNorm <- fs::path(processed_data_folder, "hipc_merged_young_noNorm.rds")
p_load_BTM <- fs::path(processed_data_folder, "BTM_processed.rds")

# Read in the files
hipc_merged_young_noNorm <- readRDS(p_load_hipc_merged_young_noNorm)
BTM <- readRDS(p_load_BTM)

# Extract the gene names
gene_names_hipc = hipc_merged_young_noNorm %>%
  select(a1cf:zzz3) %>%
  colnames()

# Results list to store outputs
results_list <- list()

# Make a code to indicate active treatment (post-vaccination) and control (pre-vaccination) samples
hipc_merged_young_noNorm <- hipc_merged_young_noNorm %>%
  mutate(vaccine_code = as.factor(ifelse(time_post_last_vax > 0, 2, 1)),
         study_accession = study_accession %>% as.factor()) %>%
  relocate(vaccine_code, .after = vaccine_name)

# Identify all post-vaccination timepoints (greater than zero) and sort
timepoints <- hipc_merged_young_noNorm %>%
  filter(time_post_last_vax > 0) %>%
  pull(time_post_last_vax) %>%
  unique() %>%
  sort()

# Number of distinct timepoints and initialize index vector
n_timepoints <- length(timepoints)
index <- integer(n_timepoints + 1)  # integer vector of zeros


# Iterate over each timepoint and perform comparisons against control
for (i in seq_along(timepoints)) {
  day <- timepoints[i]
  
  # Identify vaccines
  valid_vaccines <- hipc_merged_young_noNorm %>%
    group_by(vaccine_name, participant_id) %>%
    summarise(
      has_post = any(time_post_last_vax == day),
      has_pre = any(time_post_last_vax <= 0),
      .groups = "drop"
    ) %>%
    filter(has_post & has_pre) %>%
    pull(vaccine_name) %>%
    unique()
  
  # Skip if no valid vaccines found
  if (length(valid_vaccines) == 0)
    next
  
  # Record how many comparisons this day
  index[i + 1] <- length(valid_vaccines)
  
  # Loop through each vaccine for this day
  for (vax in valid_vaccines) {
    # Progress message
    message(paste0("Analysing ", vax, " at day ", day, " against self-controls."))
    
    # Keep only participants with BOTH pre and post-vaccine measurements
    df <- hipc_merged_young_noNorm %>%
      filter(vaccine_name == vax,
             time_post_last_vax == day |
               time_post_last_vax <= 0) %>%
      group_by(participant_id) %>%
      filter(any(time_post_last_vax == day),
             any(time_post_last_vax <= 0)) %>%
      mutate(max_pre_time = if (any(time_post_last_vax <= 0)) {
        max(study_time_collected[time_post_last_vax <= 0], na.rm = TRUE)
      } else {
        NA_real_
      }) %>%
      filter(time_post_last_vax == day |
               study_time_collected == max_pre_time) %>%
      select(-max_pre_time) %>%
      ungroup()
    
    # Convert participant IDs to numeric factor
    df <- df %>%
      mutate(participant_id = as.factor(participant_id) %>% as.numeric() %>% as.factor())
    subject_ids <- df$participant_id
    
    # Test variable: vaccine code (numeric)
    phi <- as.numeric(df$vaccine_code)
    
    # Build the design matrix for covariates, removing if collinearity is found
    X_full <-
      model.matrix(~ vaccine_code + age_imputed + gender + study_accession, data = df)
    
    # Remove linearly dependent columns
    qrX    <- qr(X_full)
    keep   <-
      sort(qrX$pivot[seq_len(qrX$rank)])    # indices of independent columns
    X_indep <- X_full[, keep, drop = FALSE]
    
    # 3. Drop all the vaccine_code columns
    is_vax_col <- grepl("^vaccine_code", colnames(X_indep))
    x <- X_indep[, !is_vax_col, drop = FALSE]
    
    message(
      paste0(
        "After resolving collinearity issues, adjusting for the following covariates: ",
        paste(colnames(x), collapse = ", ")
      )
    )
    
    # Prepare expression matrix (genes × samples)
    exprmat <- df %>%
      select(all_of(gene_names_hipc)) %>%
      t()
    
    exprmat = na.omit(exprmat)
    
    exprmat_postvax = df %>%
      filter(time_post_last_vax == day) %>%
      select(all_of(gene_names_hipc)) %>%
      t()
    
    exprmat_postvax = na.omit(exprmat_postvax)
    
    response = df %>%
      filter(time_post_last_vax == day) %>%
      pull(immResp_MFC_anyAssay_log2_MFC)
    
    # Do permutation testing if sample size < 200
    if (ncol(exprmat) < 200) {
      # P-values
      test_res <- dgsa_seq(
        exprmat = exprmat,
        covariates = x,
        variables2test = phi %>% as.matrix(),
        weights_var2test_condi = FALSE,
        genesets = BTM[["genesets"]],
        sample_group = subject_ids,
        which_weights = 'loclin',
        which_test = 'permutation',
        progressbar = TRUE,
        parallel_comp = TRUE,
        preprocessed = TRUE,
        gene_based_weights = FALSE,
        transform = TRUE,
        padjust_methods = "BH",
        bw = "nrd",
        kernel = "gaussian",
        homogen_traj = FALSE,
        na.rm_gsaseq = TRUE,
        verbose = FALSE,
        nb_cores = 1,
        n_perm = 1000,
        adaptive = TRUE
      )
    } else {
      test_res <- dgsa_seq(
        exprmat = exprmat,
        covariates = x,
        variables2test = phi %>% as.matrix(),
        weights_var2test_condi = TRUE,
        genesets = BTM[["genesets"]],
        sample_group = subject_ids,
        which_weights = 'loclin',
        which_test = 'asymptotic',
        progressbar = TRUE,
        parallel_comp = TRUE,
        preprocessed = TRUE,
        gene_based_weights = FALSE,
        transform = TRUE,
        padjust_methods = "BH",
        bw = "nrd",
        kernel = "gaussian",
        homogen_traj = FALSE,
        na.rm_gsaseq = TRUE,
        verbose = FALSE,
        nb_cores = 7
      )
    }
    
    #  Scores
    score_res = calculate_scores(
      y = exprmat,
      x,
      BTM = BTM,
      phi = phi %>% as.matrix(),
      use_phi = TRUE,
      preprocessed = TRUE,
      gene_based = FALSE,
      bw = "nrd",
      active_level = 2,
      sample_group = subject_ids
    )
    
    # Correlation
    cor_res = calculate_gs_correlation(y = exprmat_postvax,
                                       response = response,
                                       BTM = BTM)
    
    # Store and name results
    idx <- sum(index[1:i]) + which(valid_vaccines == vax)
    
    results_list[[idx]] <- list(pvals = test_res$pvals,
                                score = score_res,
                                cor = cor_res)
    
    names(results_list)[idx] <- sprintf("%s vs %s - Day %s", vax, "Control", day)
    
    # Log percentage of significant gene sets
    sig_pct <- 100 * mean(test_res$pvals$adjPval < 0.05)
    
    # Output progress message
    message(
      sprintf(
        "%s vs %s: %.1f%% significant gene sets at day %s",
        "Control",
        vax,
        sig_pct,
        day
      )
    )
    
    # Clean up
    rm(exprmat,
       test_res,
       exprmat_postvax,
       cor_res,
       score_res,
       df,
       x,
       qrX,
       X_full,
       X_indep)
    gc()
  }
}

# Specify path for saving list of results
p_results = fs::path("output", "results", "dearseq_dgsa_list.rds")

# Save results
saveRDS(results_list, file = p_results)

# If necessary, load back in results
dearseq_dgsa_list = readRDS(p_results)

# Now, we do some preprocessing of these results to collect them into a dataframe

# We write a helper function to extract the information from the list and process them

#'Pre-processing of many gene-set differential analyses for comparison
#'
#'A function which takes lists of results from multiple gene-set analyses with
#'the \code{dgsa_seq()} function and preprocesses it for visualisation
#'
#'@param raw_pvals_list a list of raw-values from a \code{dgsa_seq()} analysis.
#' Each element of the list is the output of one differential analysis i.e. a vector of raw p-values
#' where each entry of the vector corresponds to one geneset.
#'
#'@param scores_list a list of outputs from the \code{compute_scores()} function.
#' Each list entry should be a list with at least two elements \code{activation.scores} and \code{fc.scores}.
#' The ordering and names of the elements should match exactly that of the \code{dgsa_seq_list} argument.
#'
#'@param correlations_list optional list giving correlations between each geneset and a given downstream
#' response to the condition under investigation, given by the \code{calculate_gs_correlation()} function.
#' Each list entry should be a list with at least two elements \code{mean.corr} and \code{corr.mean}.
#'
#'@param conditions a vector of strings giving the active condition being investigated in each element.
#' If not given, assume all elements are from distinct conditions and will give generic names to each.
#'
#'@param condition_colors an optional vector of strings containing hex codes corresponding to each
#'of the \code{conditions}. If not given, will auto-generate colours for each condition with
#'\code{Rcolorbrewer}.
#'
#'@param times a vector of numeric times (in days post-perturbation) corresponding to each element of
#'the results list. If not given, assume all elements are from the same timepoint.
#'
#'@param genesets a \code{GSA.genesets} object containing the following named elements \itemize{
#' \item \code{genesets} : a list where each element gives the names of the genes within a set
#' \item \code{geneset.names} : a vector of strings containing the names of the genesets
#' \item \code{geneset.descriptions} : a vector of strings containing descriptions for the genesets
#' \item \code{geneset.aggregates} : a vector of strings containing aggregate groups for the genesets.
#' If this is given as a factor variable, the ordering will be preserved for visualisation. Else,
#' the default ordering will be alphabetical.
#'}
#'
#'@param aggregate_colors an optional vector of strings containing hex codes corresponding to each of
#'the \code{geneset.aggregates}. If given, the ordering that these colours will be assigned is given
#'by the factor ordering of \code{geneset.aggregates}. If the factor ordering is not given, this
#'defaults to alphabetical. If not given, colours will be assigned automatically for each aggregate
#'with \code{Rcolorbrewer}.
#'
#'@returns comparison_dataframe : a dataframe of results ready to be passed for visualisation.
#'
#'@import dplyr, Rcolorbrewer, purrr
#'@author Arthur Hughes
#'

dgsa_comparison_preprocessing = function(raw_pvals_list,
                                         scores_list,
                                         correlations_list = NULL,
                                         conditions = NULL,
                                         condition_colors = NULL,
                                         times = NULL,
                                         genesets,
                                         aggregate_colors = NULL) {
  n_gene_sets = length(raw_pvals_list[[1]])
  n_comparisons = length(raw_pvals_list)
  
  ### INPUT TRANSFORMATION ###
  # If no conditions are specified, assume each element is a disinct condition and assign generic names
  if (is.null(conditions)) {
    conditions = paste0('Condition ', 1:n_comparisons)
  }
  
  # If no times are specified, assume each element is from the same time
  if (is.null(times)) {
    times = rep("NA", n_comparisons)
  }
  
  condition_times = paste0(conditions, " - Day ", times)
  
  ### VALIDITY CHECKS ###
  # Make sure everything has the same number of elements
  
  if (length(scores_list[[1]][["activation.scores"]]) != n_gene_sets |
      length(correlations_list[[1]][["mean.corr"]]) != n_gene_sets |
      length(genesets[["geneset.names"]]) != n_gene_sets) {
    stop("Number of genesets does not match across inputs!")
  }
  
  if (length(conditions) != n_comparisons |
      length(times) != n_comparisons |
      length(scores_list) != n_comparisons |
      length(correlations_list) != n_comparisons) {
    stop("Number of comparisons does not match across inputs!")
  }
  
  # Transform genesets object into geneset metadata
  gs_metadata <- list(
    names       = genesets$geneset.names,
    descriptions = genesets$geneset.descriptions,
    aggregates  = genesets$geneset.aggregates,
    sets        = genesets$genesets
  )
  
  # Define a colour for each distinct condition
  # If conditions is given as a factor, conserve the ordering
  # Else, order alphabetically
  if (is.factor(conditions)) {
    unique_conditions <- levels(conditions)
  } else {
    unique_conditions <- sort(unique(conditions))
    conditions = factor(conditions, levels = unique_conditions)
  }
  
  n_conditions = length(unique_conditions)
  
  # Assign colours to conditions
  ## If colours given, use these
  if (!is.null(condition_colors)) {
    condition_colours = setNames(condition_colors, unique_conditions)
  } else {
    # Else generate with Rcolorbrewer
    colours <-
      brewer.pal(min(n_conditions, brewer.pal.info["Paired", "maxcolors"]), "Paired")
    
    # If more unique conditions than available colors, you can interpolate
    if (n_conditions > length(colours)) {
      colours <- colorRampPalette(colours)(n_conditions)
    }
    
    # Assign colours to condition names
    condition_colours <- setNames(colours, unique_conditions)
  }
  
  # Now define colours for geneset aggregates
  # If genesets is given as a factor, conserve the ordering
  # Else, order alphabetically
  if (is.factor(gs_metadata$aggregates)) {
    unique_aggregates <- levels(gs_metadata$aggregates)
  } else {
    unique_aggregates <- sort(unique(gs_metadata$aggregates))
    
    gs_metadata$aggregates = factor(gs_metadata$aggregates, levels = unique_aggregates)
  }
  n_aggregates = length(unique_aggregates)
  
  # If colours are provided, use them
  if (!is.null(aggregate_colors)) {
    aggregate_colours <- setNames(aggregate_colors, unique_aggregates)
  } else {
    # Else generate with Rcolorbrewer
    colours <-
      brewer.pal(min(n_aggregates, brewer.pal.info["Spectral", "maxcolors"]), "Spectral")
    
    # If more unique aggregates than available colors, you can interpolate
    if (n_aggregates > length(colours)) {
      colours <- colorRampPalette(colours)(n_aggregates)
    }
    # Assign colours to condition names
    aggregate_colours <- setNames(colours, unique_aggregates)
  }
  
  
  # Assemble results into a single data frame
  if (is.null(correlations_list)) {
    results_df <-
      purrr::map_dfr(seq_along(raw_pvals_list), function(idx) {
        # Build base tibble
        tibble(
          comparison          = condition_times[idx],
          condition           = conditions[idx],
          time                = times[idx],
          condition.colour    = condition_colours[conditions][idx],
          gs.name             = gs_metadata$names,
          gs.description      = gs_metadata$descriptions,
          gs.name.description = paste0(gs_metadata$names, " - ", gs_metadata$descriptions),
          gs.aggregate        = gs_metadata$aggregates,
          gs.colour           = aggregate_colours[gs_metadata$aggregates],
          activation.score    = scores_list[[idx]][["activation.scores"]] %>% unlist(),
          fc.score            = scores_list[[idx]][["fc.scores"]] %>% unlist(),
          corr.mean           = NA_real_,
          mean.corr           = NA_real_,
          rawPval             = raw_pvals_list[[idx]]
        )
      })
  } else {
    results_df <-
      purrr::map_dfr(seq_along(raw_pvals_list), function(idx) {
        # Build base tibble
        tibble(
          comparison          = condition_times[idx],
          condition           = conditions[idx],
          time                = times[idx],
          condition.colour    = condition_colours[conditions][idx],
          gs.name             = gs_metadata$names,
          gs.description      = gs_metadata$descriptions,
          gs.name.description = paste0(gs_metadata$names, " - ", gs_metadata$descriptions),
          gs.aggregate        = gs_metadata$aggregates,
          gs.colour           = aggregate_colours[gs_metadata$aggregates],
          activation.score    = scores_list[[idx]][["activation.scores"]] %>% unlist(),
          fc.score            = scores_list[[idx]][["fc.scores"]] %>% unlist(),
          corr.mean           = correlations_list[[idx]][["corr.mean"]],
          mean.corr           = correlations_list[[idx]][["mean.corr"]],
          rawPval             = raw_pvals_list[[idx]]
        )
      })
  }
  
  # Now perform different p-value correction methods/approaches
  correction_methods <-
    c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY")
  for (method in correction_methods) {
    # Global correction across all gene sets
    results_df <- results_df %>%
      mutate(!!paste0("global.adjPval_", method) := p.adjust(rawPval, method = method))
    
    # Timepoint-wise correction
    results_df <- results_df %>%
      group_by(time) %>%
      mutate(!!paste0("withinTime.adjPval_", method) := p.adjust(rawPval, method = method)) %>%
      ungroup()
    
    # Comparison-wise correction
    results_df <- results_df %>%
      group_by(comparison) %>%
      mutate(
        !!paste0("withinComparison.adjPval_", method) := p.adjust(rawPval, method = method)
      ) %>%
      ungroup()
  }
  
  
  return(results_df)
}

# Now we extract the arguments from the list to pass to the function
# The raw p-values
raw_pvals_list = lapply(dearseq_dgsa_list, function(x)
  x[["pvals"]][["rawPval"]])

# The list of scores
scores_list <- lapply(dearseq_dgsa_list, function(x) {
  activation <- x[["score"]][["activation.scores"]]
  fc <- x[["score"]][["fc.scores"]]
  
  list(activation.scores = activation, fc.scores = fc)
})

# The list of correlations
correlations_list = lapply(dearseq_dgsa_list, function(x) {
  # Extract the two score vectors
  mean.corr <- x[["cor"]][["mean.corr"]]
  corr.mean <- x[["cor"]][["corr.mean"]]
  
  # Combine them into a matrix
  list(mean.corr = mean.corr, corr.mean = corr.mean)
})

# The names of the experimental conditions (vaccines)
conditions <- sub("^(.*?)\\s+vs.*", "\\1", names(dearseq_dgsa_list))

# The factor ordering for the conditions
conditions_order <- c(
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

conditions = factor(conditions, levels = conditions_order)

# The ordered timepoints
times <-
  as.numeric(sub(".*Day\\s+([0-9.]+)$", "\\1", names(dearseq_dgsa_list)))

times = times %>%
  factor(levels = unique(times))

# A colour to associate with each vaccine
condition_colors = c(
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

# A colour to associate with each gene-set aggregate
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

# Now apply the function to process the list
results_df = dgsa_comparison_preprocessing(
  raw_pvals_list = raw_pvals_list,
  scores_list = scores_list,
  correlations_list = correlations_list,
  conditions = conditions,
  condition_colors = condition_colors,
  times = times,
  genesets = BTM,
  aggregate_colors = aggregate_colors
)

# Do some post-processing modifications
results_df = results_df %>%
  mutate(
    activation.score = as.numeric(activation.score),
    fc.score = as.numeric(fc.score),
    method = "dearseq"
  )

# Save the processed dgsa results dataframe
p_results_df = fs::path("output", "results", "dearseq_dgsa_results_processed.rds")
saveRDS(results_df, file = p_results_df)
