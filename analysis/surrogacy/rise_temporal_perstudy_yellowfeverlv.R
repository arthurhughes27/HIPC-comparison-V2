# File to find trial-level gene expression surrogates of Yellow Fever vaccination
# We do a paired-data analysis per-study to account for study effects

library(tidyverse)
library(SurrogateRank)
library(purrr)
library(stringr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Directory containing engineered / processed data files
processed_data_folder <- "data"

# Paths to processed gene-level data and gene-set objects
p_load_expr_all_norm <- fs::path(processed_data_folder, "hipc_merged_all_norm.rds")

# Load data objects
hipc_merged_all_norm <- readRDS(p_load_expr_all_norm)

# Timepoints of interest (numeric)
timepoints_of_interest <- c(0, 1, 3, 7, 10, 14)

# Filter to samples with non-missing immune response, Yellow Fever vaccine,
# and collected at one of the specified timepoints.
hipc_merged_all_norm_filtered <- hipc_merged_all_norm %>%
  filter(
    !is.na(immResp_MFC_anyAssay_log2_MFC),
    vaccine_name == "Yellow Fever (LV)",
    study_time_collected %in% timepoints_of_interest
  )

gene_names = hipc_merged_all_norm_filtered %>% 
  dplyr::select(a1cf:zzz3) %>% 
  colnames()


study_names = hipc_merged_all_norm_filtered %>% 
  filter(study_time_collected > 0) %>% 
  pull(study_accession) %>% 
  unique()

rise_results_list = vector("list", length(study_names))
names(rise_results_list) = study_names

for(study in study_names) {
  message(paste0("Analysing study ", study, "..."))
  
  hipc_merged_all_norm_filtered_study = hipc_merged_all_norm_filtered %>%
    filter(study_accession == study)
  
  timepoints <- hipc_merged_all_norm_filtered_study %>%
    filter(study_time_collected > 0) %>%
    count(study_time_collected) %>%     
    filter(n >= 4) %>%                   
    pull(study_time_collected) %>%
    sort()
  
  rise_results_list[[study]] = vector("list", length(timepoints))
  names(rise_results_list[[study]]) = paste0("Day ", timepoints)
  
  for (tp in timepoints) {
    message(paste0("Day ", tp, " in progress."))
    
    
    hipc_merged_all_norm_filtered_study_tp = hipc_merged_all_norm_filtered_study %>%
      filter(study_time_collected %in% c(0, tp)) %>%
      group_by(participant_id) %>%
      filter(all(c(0, tp) %in% study_time_collected)) %>%
      ungroup()
    
    rise_df = hipc_merged_all_norm_filtered_study_tp %>%
      dplyr::select(
        participant_id,
        study_time_collected,
        immResp_MFC_anyAssay_pre_value,
        immResp_MFC_anyAssay_post_value,
        all_of(gene_names)
      ) %>%
      arrange(participant_id)
    
    yzero = rise_df %>%
      filter(study_time_collected == 0) %>%
      pull(immResp_MFC_anyAssay_pre_value)
    
    yone = rise_df %>%
      filter(study_time_collected == tp) %>%
      pull(immResp_MFC_anyAssay_post_value)
    
    szero = rise_df %>%
      filter(study_time_collected == 0) %>%
      dplyr::select(all_of(gene_names))
    
    sone = rise_df %>%
      filter(study_time_collected == tp) %>%
      dplyr::select(all_of(gene_names))
    
    yresult = SurrogateRank::test.surrogate.extension(yone = yone, yzero = yzero,
                                                      sone = yone, szero = yzero,
                                                      power.want.s = 0.8,
                                                      paired = T,
                                                      alternative = "two.sided")[["u.y"]]
    
    rise_result = rise.screen(
      yone,
      yzero,
      sone,
      szero,
      alpha = 0.05,
      power.want.s = 0.8,
      p.correction = "BH",
      n.cores = 1,
      alternative = "two.sided",
      paired = TRUE
    )
    
    rise_result[["screening.metrics"]]$u_y = yresult
    rise_result[["screening.metrics"]]$u_s = yresult - rise_result[["screening.metrics"]]$delta
    rise_result[["screening.metrics"]]$n = 2*length(yone)
    
    rise_results_list[[study]][[paste0("Day ", tp)]] = rise_result[["screening.metrics"]]
    
  }
  
}

p_rise_results_list_yellowfeverlv = fs::path(
  "output",
  "results",
  "surrogacy",
  "rise_results_list_yellowfeverlv.rds"
)


saveRDS(rise_results_list,
        file = p_rise_results_list_yellowfeverlv)

rise_results_list = readRDS(p_rise_results_list_yellowfeverlv)

# Bind the results together including study and timepoint identifiers

rise_results_df <- purrr::imap_dfr(rise_results_list, ~ {
  study_list <- .x
  study_name <- .y
  
  # inner: for each timepoint data.frame in this study
  purrr::imap_dfr(study_list, ~ {
    df <- .x
    tp_name <- .y
    
    # skip non-data.frame or empty elements
    if (!is.data.frame(df) ||
        nrow(df) == 0)
      return(tibble::tibble())
    
    # extract numeric from "Day x" (handles "Day 7", " Day 7 ", etc.)
    tp_num <- as.numeric(stringr::str_extract(tp_name, "\\d+"))
    
    df %>%
      mutate(study_accession = study_name,
             study_time_collected = tp_num)
  })
}) %>%
  # put the new columns first
  relocate(study_accession, study_time_collected)

tp = 3

rise_results_df_timepoint = rise_results_df_filtered <- rise_results_df %>%
  filter(study_time_collected == tp)


# summary table: number of non-missing pairs and spearman correlation per marker
res_by_marker <- rise_results_df_timepoint %>%
  group_by(marker) %>%
  summarize(
    n_pairs = sum(!is.na(u_y) & !is.na(u_s)),
    mae = ifelse(
      n_pairs > 0,
      mean(abs(u_y - u_s), na.rm = TRUE),
      NA_real_
    ),
    spearman = ifelse(
      n_pairs > 2,
      cor(u_y, u_s, method = "spearman", use = "pairwise.complete.obs"),
      NA_real_
    ),
    .groups = "drop"
  )



gene_name = "ifi44"


rise_results_df_filtered <- rise_results_df %>%
  filter(marker == gene_name,
         study_time_collected == tp)



p <- ggplot(rise_results_df_filtered, aes(x = u_y, y = u_s)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +
  scale_x_continuous(limits = c(-0.1, 1.1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.1, 1.1), expand = c(0, 0)) +
  coord_fixed(ratio = 1) +                    # keep the same scale on both axes
  labs(x = "UY", y = "US", title = paste0("U_Y vs U_S for gene ", gene_name)) +
  theme_minimal(base_size = 18) +
  theme(
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 16),
    plot.title = element_text(size = 20, hjust = 0.5)
  )

print(p)
