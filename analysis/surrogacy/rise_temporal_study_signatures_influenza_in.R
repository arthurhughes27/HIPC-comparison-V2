# File to find trial-level gene expression surrogates of influenza vaccination
# At each timepoint, we derive a signature from the study with the largest sample size
# We then apply the signature to each study one by one and perform trial-level meta-analysis

library(tidyverse)
library(SurrogateRank)

# Directory containing engineered / processed data files
processed_data_folder <- "data"

# Paths to processed gene-level data and gene-set objects
p_load_expr_all_norm <- fs::path(processed_data_folder, "hipc_merged_all_norm.rds")

# Load data objects
hipc_merged_all_norm <- readRDS(p_load_expr_all_norm)

# Timepoints of interest (numeric)
timepoints_of_interest <- c(0, 1, 3, 7, 10, 14)

# Filter to samples with non-missing immune response, Influenza vaccine,
# and collected at one of the specified timepoints.
hipc_merged_all_norm_filtered <- hipc_merged_all_norm %>%
  filter(
    !is.na(immResp_MFC_anyAssay_log2_MFC),
    vaccine_name == "Influenza (IN)",
    study_time_collected %in% timepoints_of_interest
  )

# Extract the gene names
gene_names = hipc_merged_all_norm_filtered %>%
  dplyr::select(a1cf:zzz3) %>%
  colnames()

# Extract the study names
study_names = hipc_merged_all_norm_filtered %>%
  filter(study_time_collected > 0) %>%
  pull(study_accession) %>%
  unique()

timepoints_of_interest = c(1)

# Initialise results vector
rise_studysignatures_results_list = vector("list", length(timepoints_of_interest))
names(rise_studysignatures_results_list) = as.character(timepoints_of_interest)

n_markers = 200

for (tp in timepoints_of_interest) {
  message(paste0("Analysing timepoint ", tp, "..."))
  
  # Filter the data by timepoint
  
  hipc_merged_all_norm_filtered_tp = hipc_merged_all_norm_filtered %>%
    filter(study_time_collected %in% c(0, tp)) %>%
    group_by(participant_id) %>%
    filter(all(c(0, tp) %in% study_time_collected)) %>%
    ungroup()
  
  study_names_tp = hipc_merged_all_norm_filtered_tp %>%
    pull(study_accession) %>%
    unique()
  
  n_studies_tp = study_names_tp %>%
    length()
  
  # Find the largest study at this timepoint
  
  largest_study = hipc_merged_all_norm_filtered_tp$study_accession %>%
    table() %>%
    which.max() %>%
    names()
  # largest_study = "SDY80"
  
  largest_study_n = hipc_merged_all_norm_filtered_tp %>% 
    filter(study_accession == largest_study) %>% 
    pull(study_accession) %>%
    table() %>%
    max()
  
  message(
    paste0(
      "Largest study is ",
      largest_study,
      " with ",
      largest_study_n,
      " observations."
    )
  )
  
  # Derive a signature based on the largest study
  hipc_merged_all_norm_filtered_tp_largest = hipc_merged_all_norm_filtered_tp %>%
    filter(study_accession == largest_study)
  
  # Arrange data
  rise_df = hipc_merged_all_norm_filtered_tp_largest %>%
    dplyr::select(
      participant_id,
      study_time_collected,
      immResp_MFC_anyAssay_pre_value,
      immResp_MFC_anyAssay_post_value,
      all_of(gene_names)
    ) %>%
    arrange(participant_id)
  
  # Pre-vaccination response
  yzero = rise_df %>%
    filter(study_time_collected == 0) %>%
    pull(immResp_MFC_anyAssay_pre_value)
  
  # Post-vaccination response
  yone = rise_df %>%
    filter(study_time_collected == tp) %>%
    pull(immResp_MFC_anyAssay_post_value)
  
  # Pre-vaccination gene expression
  szero = rise_df %>%
    filter(study_time_collected == 0) %>%
    dplyr::select(all_of(gene_names))
  
  # Post-vaccination gene expression
  sone = rise_df %>%
    filter(study_time_collected == tp) %>%
    dplyr::select(all_of(gene_names))
  
  # Extract the estimated effect on the response on the relevant samples
  yresult = SurrogateRank::test.surrogate.extension(
    yone = yone,
    yzero = yzero,
    sone = yone,
    szero = yzero,
    power.want.s = 0.8,
    paired = T,
    alternative = "two.sided"
  )[["u.y"]]
  
  # Screen all markers
  rise_screen_result = rise.screen(yone = yone,
                                   yzero = yzero,
                                   sone = sone,
                                   szero = szero,
                                   alpha = 0.05,
                                   power.want.s = 0.8,
                                   # epsilon = 0.3,
                                   p.correction = "BH",
                                   n.cores = 1,
                                   alternative = "two.sided",
                                   paired = TRUE)
  
  rise_screen_result_df = rise_screen_result[["screening.metrics"]] %>% 
    arrange(abs(delta))
  
  rise_screen_genes = rise_screen_result_df %>% 
    slice_head(n = n_markers) %>% 
    pull(marker)
  
  # Pre-vaccination gene expression
  szero_screened = szero %>%
    dplyr::select(all_of(rise_screen_genes))
  
  # Post-vaccination gene expression
  sone_screened = sone %>%
    dplyr::select(all_of(rise_screen_genes))
  
  weight_df = rise_screen_result_df %>% 
    slice_head(n = n_markers) %>% 
    mutate(weight = abs(1/delta)) %>% 
    select(marker, weight)
  
  rise_result = rise.evaluate(
    yone,
    yzero,
    sone = sone_screened,
    szero = szero_screened,
    alpha = 0.05,
    power.want.s = 0.8,
    # epsilon = 0.3,
    p.correction = "BH",
    n.cores = 1,
    alternative = "two.sided",
    paired = TRUE,
    evaluate.weights = T,
    screening.weights = weight_df
  )
  
  message(paste0(length(rise_result[["screening.results"]][["significant.markers"]]), " significant markers identified in signature."))
  
  rise_study_df = data.frame(as.list(rise_result[["gamma.s.evaluate"]]))
  
  # Store the treatment effects on each marker
  rise_study_df$u_y = yresult
  rise_study_df$u_s = yresult - rise_study_df$delta
  
  # Store the number of observations
  rise_study_df$n = 2 * length(yone)
  
  rise_study_df$study_accession = largest_study
  
  
  char_tp = as.character(tp)
  # Store the result
  rise_studysignatures_results_list[[char_tp]] = vector("list", 1)
  names(rise_studysignatures_results_list[[char_tp]]) = largest_study
  
  rise_studysignatures_results_list[[char_tp]][[largest_study]] = vector("list", n_studies_tp)
  names(rise_studysignatures_results_list[[char_tp]][[largest_study]]) = study_names_tp
  
  rise_studysignatures_results_list[[char_tp]][[largest_study]][[largest_study]] = rise_study_df
  
  remaining_studies = hipc_merged_all_norm_filtered_tp %>%
    filter(study_accession != largest_study) %>%
    pull(study_accession) %>%
    unique()
  
  for (study in remaining_studies) {
    
    # Derive a signature based on the largest study
    hipc_merged_all_norm_filtered_tp_study = hipc_merged_all_norm_filtered_tp %>%
      filter(study_accession == study)
    
    # Arrange data
    rise_df = hipc_merged_all_norm_filtered_tp_study %>%
      dplyr::select(
        participant_id,
        study_time_collected,
        immResp_MFC_anyAssay_pre_value,
        immResp_MFC_anyAssay_post_value,
        all_of(gene_names)
      ) %>%
      arrange(participant_id)
    
    # Pre-vaccination response
    yzero = rise_df %>%
      filter(study_time_collected == 0) %>%
      pull(immResp_MFC_anyAssay_pre_value)
    
    # Post-vaccination response
    yone = rise_df %>%
      filter(study_time_collected == tp) %>%
      pull(immResp_MFC_anyAssay_post_value)
    
    # Pre-vaccination gene expression
    szero = rise_df %>%
      filter(study_time_collected == 0) %>%
      dplyr::select(all_of(gene_names))
    
    # Post-vaccination gene expression
    sone = rise_df %>%
      filter(study_time_collected == tp) %>%
      dplyr::select(all_of(gene_names))
    
    # Extract the estimated effect on the response on the relevant samples
    yresult = SurrogateRank::test.surrogate.extension(
      yone = yone,
      yzero = yzero,
      sone = yone,
      szero = yzero,
      power.want.s = 0.8,
      paired = T,
      alternative = "two.sided"
    )[["u.y"]]
    
    study_evaluate = rise.evaluate(
      yone,
      yzero,
      sone,
      szero,
      alpha = 0.05,
      power.want.s = 0.8,
      alternative = "two.sided",
      paired = T,
      evaluate.weights = T,
      screening.weights = weight_df,
      markers = weight_df$marker
    )
    
    rise_study_df_evaluate = data.frame(as.list(study_evaluate[["gamma.s.evaluate"]]))
    
    # Store the treatment effects on each marker
    rise_study_df_evaluate$u_y = yresult
    rise_study_df_evaluate$u_s = yresult - rise_study_df_evaluate$delta
    
    # Store the number of observations
    rise_study_df_evaluate$n = 2 * length(yone)
    
    rise_study_df_evaluate$study_accession = study
    
    rise_studysignatures_results_list[[char_tp]][[largest_study]][[study]] = rise_study_df_evaluate
    
  }
}

df_list <- rise_studysignatures_results_list[[as.character(timepoints_of_interest)]][[largest_study]]
combined_df <- do.call(rbind, df_list)

plot_min = min(combined_df$u_y,
               combined_df$u_s) - 0.5
plot_min = 0

# 1) Get min and max of n
n_min <- min(combined_df$n, na.rm = TRUE)
n_max <- max(combined_df$n, na.rm = TRUE)

# 2) Round min up and max down to nearest 10
legend_min <- ceiling(n_min / 10) * 10
legend_max <- floor(n_max / 10) * 10

# 3) Split interval into 4 equally spaced points: min, value1, value2, max
legend_ns <- seq(legend_min, legend_max, length.out = 4)

# 4) Round intermediate values to nearest 10
legend_ns[c(2, 3)] <- round(legend_ns[c(2, 3)] / 50) * 50


ggplot(combined_df, aes(x = u_y, y = u_s)) +
  # point size mapped to n so radius ~ n (area ~ n)
  geom_point(aes(size = n), alpha = 0.8) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "red",
    linetype = "dashed",
    size = 0.8
  ) +
  scale_x_continuous(limits = c(plot_min, 1.01), expand = c(0, 0)) +
  scale_y_continuous(limits = c(plot_min, 1.01), expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  labs(
    x = expression(U[Y]),
    y = expression(U[S]),
    title = paste0(
      "Influenza (IN) : Treatment effects on response vs marker across trials (Day ",
      tp,
      ": N trials = ",
      length(unique(
        combined_df$study_accession
      )),
      ")"
    ),
    size = "Trial N"  # legend title
  ) +
  # size scale: adjust range to taste; breaks/labels show natural n
  scale_size_continuous(
    range = c(1.5, 7),
    # min/max point sizes (tweak if needed)
    breaks = legend_ns,
    # positions in mapped space
    labels = as.character(legend_ns)              # show natural-scale labels in legend
  ) +
  theme_minimal(base_size = 25) +
  theme(
    plot.title = element_text(
      size = 30,
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 20)
  )
