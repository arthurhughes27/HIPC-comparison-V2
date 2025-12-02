# File to find trial-level gene expression surrogates of influenza vaccination
# We do a paired-data analysis per-study to account for study effects

library(tidyverse)
library(SurrogateRank)
library(purrr)
library(stringr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(forcats)
library(tidytext)

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

# Initialise results vector 
rise_results_list = vector("list", length(study_names))
names(rise_results_list) = study_names

for(study in study_names) { # For each study
  message(paste0("Analysing study ", study, "...")) 
  
  # Filter data by study
  hipc_merged_all_norm_filtered_study = hipc_merged_all_norm_filtered %>%
    filter(study_accession == study)
  
  # Extract possible timepoints (must have at least 4 observations i.e. 2 pairs)
  timepoints <- hipc_merged_all_norm_filtered_study %>%
    filter(study_time_collected > 0) %>%
    count(study_time_collected) %>%     
    filter(n >= 4) %>%                   
    pull(study_time_collected) %>%
    sort()
  
  # Initialise within-study results 
  rise_results_list[[study]] = vector("list", length(timepoints))
  names(rise_results_list[[study]]) = paste0("Day ", timepoints)
  
  for (tp in timepoints) { # For each available timepoint
    message(paste0("Day ", tp, " in progress."))
    
    # Filter by timepoint 
    hipc_merged_all_norm_filtered_study_tp = hipc_merged_all_norm_filtered_study %>%
      filter(study_time_collected %in% c(0, tp)) %>%
      group_by(participant_id) %>%
      filter(all(c(0, tp) %in% study_time_collected)) %>%
      ungroup()
    
    # Arrange data
    rise_df = hipc_merged_all_norm_filtered_study_tp %>%
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
    yresult = SurrogateRank::test.surrogate.extension(yone = yone, yzero = yzero,
                                            sone = yone, szero = yzero,
                                            power.want.s = 0.8,
                                            paired = T,
                                            alternative = "two.sided")[["u.y"]]
    
    # Screen all markers
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
    
    # Store the treatment effects on each marker
    rise_result[["screening.metrics"]]$u_y = yresult
    rise_result[["screening.metrics"]]$u_s = yresult - rise_result[["screening.metrics"]]$delta
    
    # Store the number of observations
    rise_result[["screening.metrics"]]$n = 2*length(yone)
    
    # Store the results in the list
    rise_results_list[[study]][[paste0("Day ", tp)]] = rise_result[["screening.metrics"]]
    
  }
  
}

# Specify the path to save the results
p_rise_results_list_influenzain = fs::path(
  "output",
  "results",
  "surrogacy",
  "rise_results_list_influenzain.rds"
)

# Save the results
saveRDS(rise_results_list,
        file = p_rise_results_list_influenzain)

# Read back in the results if necessary
rise_results_list = readRDS(p_rise_results_list_influenzain)

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



### FIRST EXPLORATION OF RESULTS ###
# Calculate a weighted average of delta values across studies for each timepoint
# The weights will be the square root of the number of observations
# such that larger samples provide more evidence

# First, derive weighted averages
rise_summary <- rise_results_df %>%
  group_by(marker, study_time_collected) %>%
  summarise(
    w_sum = sum(sqrt(n), na.rm = TRUE),
    weighted_delta = if (w_sum > 0) sum(delta * sqrt(n), na.rm = TRUE) / w_sum else NA_real_,
    n_studies = n_distinct(study_accession),
    .groups = "drop"
  )

# Parameters
timepoints_of_interest_heatmap <- c(1, 3, 7, 14)
n_markers <- 150
cutoff <- 0.9  # define the cutoff

# Prepare plotting data
df_plot <- rise_summary %>%
  filter(study_time_collected %in% timepoints_of_interest_heatmap) %>%
  mutate(plot_value = 1 - weighted_delta) %>%
  group_by(study_time_collected) %>%
  arrange(desc(plot_value)) %>%
  slice(1:n_markers) %>%
  ungroup() %>%
  mutate(
    timepoint = factor(study_time_collected, levels = timepoints_of_interest_heatmap),
    marker_reordered = reorder_within(marker, plot_value, timepoint, .desc = TRUE),
    color_flag = ifelse(plot_value > cutoff, "above_cutoff", "below_cutoff")
  )

plot_min <- min(df_plot$plot_value, na.rm = TRUE) - 0.05

ggplot(df_plot, aes(x = fct_rev(marker_reordered), y = plot_value, color = color_flag)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_hline(yintercept = cutoff, linetype = "dashed", color = "red", size = 0.7) +
  scale_color_manual(values = c("above_cutoff" = "red", "below_cutoff" = "#2C7BB6"), guide = "none") +
  labs(
    title = "Markers ranked by average trial-level surrogacy",
    x = "Marker",
    y = "1 - weighted average of delta"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  ylim(plot_min, 1) +
  facet_wrap(~ timepoint, ncol = 1, scales = "free_x",
             labeller = labeller(timepoint = function(x) paste0("Day ", x))) +
  scale_x_reordered()

# params
tp <- 7

# Define a vector of genes which have high average trial-level 
# Surrogacy at any timepoint
markers_to_keep <- df_plot %>%
  filter(plot_value > cutoff,
         study_time_collected == tp) %>%
  pull(marker) %>%
  unique()

markers_to_keep

# Filter data for the chosen timepoint and markers
rise_results_df_filtered <- rise_results_df %>%
  filter(marker %in% markers_to_keep,
         study_time_collected == tp) %>%
  # ensure marker factor ordering matches markers_to_keep
  mutate(marker = factor(marker, levels = markers_to_keep))

# decide number of columns for facet grid
n_markers <- length(markers_to_keep)
ncol_facets <- min(5, n_markers)  # adjust 5 -> any default you prefer

plot_min = min(rise_results_df_filtered$u_y,rise_results_df_filtered$u_s) - 0.01

# plotting
p_facets <- ggplot(rise_results_df_filtered, aes(x = u_y, y = u_s)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 0.8) +
  scale_x_continuous(limits = c(plot_min, 1.01), expand = c(0, 0)) +
  scale_y_continuous(limits = c(plot_min, 1.01), expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  labs(
    x = "U_Y",
    y = "U_S",
    title = paste0("U_Y vs U_S — timepoint = ", tp)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12),
    plot.title = element_text(size = 16, hjust = 0.5),
    strip.text = element_text(size = 12)   # facet label text size
  ) +
  facet_wrap(~ marker, ncol = ncol_facets, scales = "fixed")

print(p_facets)
