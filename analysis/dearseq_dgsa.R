# File to perform differential gene-set analysis with dearseq across multiple vaccines and preprocess the results
## We perform one differential gene-set analysis per vaccine (pathogen + vaccine type)
## We adjust for the effects of age, sex, and study in the linear mixed model framework of dearseq
## To make the comparison with Hagan et al.'s results as fair as possible, we use the young, non-cross-study-normalised dataset

# Load the necessary packages
library(fs)
library(dplyr)
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
    
    message(
      sprintf(
        "%s vs %s: %.1f%% significant gene sets at day %s",
        "Control",
        vax,
        sig_pct,
        day
      )
    )
    
    # # Clean up
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

dgsa_results_directory = "./output/results"

p_results = fs::path("output", "results", "dearseq_reanalysis_list.rds")

saveRDS(results_list, file = p_results)
