# R script to perform meta-analysis with rise

library(tidyverse)
library(SurrogateRank)
library(grid)    # for unit()
library(scales)  # pretty formatting

tp = 1

# Directory containing engineered / processed data files
processed_data_folder <- "data"

# Paths to processed gene-level data and gene-set objects
p_load_expr_all_norm <- fs::path(processed_data_folder, "hipc_merged_all_norm.rds")

# Load data objects
hipc_merged_all_norm <- readRDS(p_load_expr_all_norm)

# Extract the gene names
gene_names = hipc_merged_all_norm %>%
  dplyr::select(a1cf:zzz3) %>%
  colnames()

# Timepoints of interest (numeric)
timepoints_to_keep <- c(0, tp)

# Filter to samples with non-missing immune response, Influenza vaccine,
# and collected at one of the specified timepoints.
hipc_merged_all_norm_filtered <- hipc_merged_all_norm %>%
  filter(
    !is.na(immResp_MFC_anyAssay_log2_MFC),
    vaccine_name == "Influenza (IN)",
    study_time_collected %in% timepoints_to_keep
  ) %>%
  select(
    participant_id,
    study_accession,
    study_time_collected,
    immResp_MFC_anyAssay_pre_value,
    immResp_MFC_anyAssay_post_value,
    all_of(gene_names)
  )

# Extract the study names
study_names = hipc_merged_all_norm_filtered %>%
  filter(study_time_collected > 0) %>%
  pull(study_accession) %>%
  unique()

timepoints_of_interest = c(tp)

epsilon = 0.2

for (sdy in study_names) {
  left_out_study = sdy
  
  study_names_filtered = study_names[study_names != left_out_study]
  
  for (sdy_alt in study_names_filtered) {
    hipc_merged_all_norm_filtered_sdy_alt = hipc_merged_all_norm_filtered %>%
      filter(study_time_collected %in% c(0, tp),
             study_accession == sdy_alt) %>%
      group_by(participant_id) %>%
      filter(all(c(0, tp) %in% study_time_collected)) %>%
      ungroup()
    
    # Pre-vaccination response
    yzero_sdy_alt = hipc_merged_all_norm_filtered_sdy_alt %>%
      filter(study_time_collected == 0) %>%
      pull(immResp_MFC_anyAssay_pre_value)
    
    # Post-vaccination response
    yone_sdy_alt = hipc_merged_all_norm_filtered_sdy_alt %>%
      filter(study_time_collected == tp) %>%
      pull(immResp_MFC_anyAssay_post_value)
    
    # Pre-vaccination gene expression
    szero_sdy_alt = hipc_merged_all_norm_filtered_sdy_alt %>%
      filter(study_time_collected == 0) %>%
      select(all_of(gene_names))
    
    # Post-vaccination gene expression
    sone_sdy_alt = hipc_merged_all_norm_filtered_sdy_alt %>%
      filter(study_time_collected == tp) %>%
      select(all_of(gene_names))
    
    # Extract the estimated effect on the response on the relevant samples
    yresult_sdy_alt_main = SurrogateRank::test.surrogate.extension(
      yone = yone_sdy_alt,
      yzero = yzero_sdy_alt,
      sone = yone_sdy_alt,
      szero = yzero_sdy_alt,
      power.want.s = 0.8,
      paired = T,
      alternative = "two.sided"
    )[["u.y"]]
    
    test_result_sdy_alt = SurrogateRank::rise.screen(
      yone = yone_sdy_alt,
      yzero = yzero_sdy_alt,
      sone = sone_sdy_alt,
      szero = szero_sdy_alt,
      epsilon = epsilon,
      alternative = "two.sided",
      paired = T,
      alpha = 0.05,
      p.correction = "none",
      weight.mode = "diff.epsilon"
    )
    
    
    test_result_df = test_result_sdy_alt[["screening.metrics"]]
    test_result_df$study_accession = sdy_alt
    test_result_df$n = 2 * length(yone_sdy_alt)
    test_result_df$u.y = yresult_sdy_alt_main
    test_result_df$u.s = yresult_sdy_alt_main - test_result_df$delta
    
    test_result_df = test_result_df %>%
      select(
        marker,
        study_accession,
        n,
        epsilon,
        u.y,
        u.s,
        delta,
        ci_upper,
        ci_lower,
        sd,
        p_unadjusted
      )
    
    if (sdy_alt == study_names_filtered[1]) {
      all_results = test_result_df
    } else {
      all_results = rbind(all_results, test_result_df)
    }
  }
  
  # ---- compute per-marker fixed-effect pooled results -----------------------
  summary_df <- all_results %>%
    # keep only required columns (tolerant to extra cols)
    transmute(marker,
              study_accession,
              n,
              delta,
              ci_upper,
              ci_lower,
              sd,
              p_unadjusted) %>%
    
    # group by marker and compute pooled FE; exclude rows with sd == 0 or NA for pooling
    group_by(marker) %>%
    group_modify(~{
      dat <- .x
      
      # exclude studies with sd == 0 or NA for the pooling step
      dat_use <- dat %>% filter(!is.na(sd) & sd > 0)
      
      if(nrow(dat_use) == 0) {
        # no usable studies -> return NA row
        return(tibble(
          pooled_delta = NA_real_,
          se_pooled = NA_real_,
          var_pooled = NA_real_,
          ci_pooled_lower = NA_real_,
          ci_pooled_upper = NA_real_,
          p_two_sided = NA_real_,
          p1_pooled = NA_real_,
          p2_pooled = NA_real_,
          p_tost_pooled = NA_real_,
          k_used = 0L,
          total_N = sum(dat$n, na.rm = TRUE),
          W_sum = NA_real_
        ))
      }
      
      # study-level variances and weights
      dat_use <- dat_use %>%
        mutate(var = sd^2,
               weight = 1 / var)
      
      W_sum <- sum(dat_use$weight, na.rm = TRUE)
      
      # pooled fixed-effect estimate (inverse-variance)
      pooled_delta <- sum(dat_use$weight * dat_use$delta, na.rm = TRUE) / W_sum
      var_pooled <- 1 / W_sum
      se_pooled <- sqrt(var_pooled)
      
      # 95% CI (normal approximation)
      ci_pooled_lower <- pooled_delta - qnorm(0.975) * se_pooled
      ci_pooled_upper <- pooled_delta + qnorm(0.975) * se_pooled
      
      # two-sided test vs 0 (for reporting)
      p_two_sided <- 2 * (1 - pnorm(abs(pooled_delta / se_pooled)))
      
      # pooled one-sided p-values for TOST (normal approx)
      Z1 <- (pooled_delta - epsilon) / se_pooled   # test delta < epsilon
      Z2 <- (pooled_delta + epsilon) / se_pooled   # test delta > -epsilon
      
      p1_pooled <- pnorm(Z1)           # one-sided p for delta < epsilon
      p2_pooled <- 1 - pnorm(Z2)       # one-sided p for delta > -epsilon
      
      p_tost_pooled <- max(p1_pooled, p2_pooled)  # TOST pooled p-value
      
      # return a single-row tibble for this marker
      tibble(
        pooled_delta = pooled_delta,
        se_pooled = se_pooled,
        var_pooled = var_pooled,
        ci_pooled_lower = ci_pooled_lower,
        ci_pooled_upper = ci_pooled_upper,
        p_two_sided = p_two_sided,
        p1_pooled = p1_pooled,
        p2_pooled = p2_pooled,
        p_tost_pooled = p_tost_pooled,
        k_used = nrow(dat_use),
        total_N = sum(dat$n, na.rm = TRUE),
        W_sum = W_sum
      )
    }, .keep = TRUE) %>%
    ungroup()
  
  # ---- optional: reorder columns and view -----------------------------------
  summary_df <- summary_df %>%
    mutate(p_adjusted = p.adjust(p_tost_pooled, method = "BH")) %>% 
    select(marker,
           k_used,
           total_N,
           W_sum,
           pooled_delta, se_pooled, var_pooled,
           ci_pooled_lower, ci_pooled_upper,
           p_two_sided,
           p1_pooled, p2_pooled, p_tost_pooled, p_adjusted)
  
  signature_genes = summary_df %>% 
    filter(p_adjusted < 0.05) %>% 
    pull(marker)
  
  weight_df = summary_df %>% 
    filter(marker %in% signature_genes) %>% 
    mutate(weight = 1/abs(pooled_delta)) %>% 
    select(marker, weight)
  
  hipc_merged_all_norm_filtered_sdy_eval = hipc_merged_all_norm_filtered %>%
    filter(study_time_collected %in% c(0, tp), study_accession == sdy) %>%
    group_by(participant_id) %>%
    filter(all(c(0, tp) %in% study_time_collected)) %>%
    ungroup()
  
  # Pre-vaccination response
  yzero_sdy = hipc_merged_all_norm_filtered_sdy_eval %>%
    filter(study_time_collected == 0) %>%
    pull(immResp_MFC_anyAssay_pre_value)
  
  # Post-vaccination response
  yone_sdy = hipc_merged_all_norm_filtered_sdy_eval %>%
    filter(study_time_collected == tp) %>%
    pull(immResp_MFC_anyAssay_post_value)
  
  # Pre-vaccination gene expression
  szero_sdy = hipc_merged_all_norm_filtered_sdy_eval %>%
    filter(study_time_collected == 0) %>%
    select(all_of(signature_genes))
  
  # Post-vaccination gene expression
  sone_sdy = hipc_merged_all_norm_filtered_sdy_eval %>%
    filter(study_time_collected == tp) %>%
    select(all_of(signature_genes))
  
  # Extract the estimated effect on the response on the relevant samples
  yresult_sdy_main = SurrogateRank::test.surrogate.extension(
    yone = yone_sdy,
    yzero = yzero_sdy,
    sone = yone_sdy,
    szero = yzero_sdy,
    power.want.s = 0.8,
    paired = T,
    alternative = "two.sided"
  )[["u.y"]]
  
  evaluate_result_sdy = SurrogateRank::rise.evaluate(
    yone = yone_sdy,
    yzero = yzero_sdy,
    sone = sone_sdy,
    szero = szero_sdy,
    epsilon = epsilon,
    alternative = "two.sided",
    paired = T,
    alpha = 0.05,
    p.correction = "none",
    evaluate.weights = T,
    screening.weights = weight_df,
    markers = signature_genes
  )
  
  evaluate_result_df = as.data.frame(as.list(evaluate_result_sdy[["gamma.s.evaluate"]]))
  evaluate_result_df$study_accession = sdy
  evaluate_result_df$n = 2 * length(yone_sdy)
  evaluate_result_df$u.y = yresult_sdy_main
  evaluate_result_df$u.s = yresult_sdy_main - evaluate_result_df$delta
  evaluate_result_df$n_genes = length(signature_genes)
  
  evaluate_result_df = evaluate_result_df %>%
    select(study_accession,
           n,
           n_genes,
           epsilon,
           u.y,
           u.s,
           delta,
           ci_upper,
           ci_lower,
           sd,
           p_unadjusted)
  
  if (sdy == study_names[1]) {
    all_results_eval = evaluate_result_df
  } else {
    all_results_eval = rbind(all_results_eval, evaluate_result_df)
  }
}

min_n = min(all_results_eval$n)
max_n = max(all_results_eval$n)

br <- pretty(c(min_n, max_n), n = 3)
br <- unique(c(min_n, br[br >= min_n & br <= max_n], max_n))
# if pretty returned only one value (rare), create sensible breaks
if (length(br) == 1) br <- unique(round(c(min_n, (min_n + max_n)/2, max_n)))
# ensure increasing numeric vector
br <- sort(as.numeric(br))

p <- ggplot(all_results_eval, aes(x = u.y, y = u.s)) +
  geom_point(aes(size = n), alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 0.7) +
  coord_fixed(xlim = c(0.7,1.05), ylim = c(0.7,1.05), ratio = 1, expand = FALSE) +
  labs(
    title = "Leave-One-Study-Out Cross Validation",
    x = expression(U[Y]),
    y = expression(U[S])
  ) +
  # scale point area so area ∝ n_plot; use limits and breaks based on data
  scale_size_area(
    name = "n",
    max_size = 13,
    limits = c(min_n, max_n),
    breaks = br,
    labels = br
  ) +
  theme_minimal() +
  theme(
    # big, centered title/subtitle
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    # larger axis titles & text
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    # legend styling
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    # comfortable margins
    plot.margin = margin(10, 10, 10, 10)
  )

p
