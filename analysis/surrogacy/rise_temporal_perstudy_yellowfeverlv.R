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

# Filter to samples with non-missing immune response, Yellow Fever vaccine,
# and collected at one of the specified timepoints.
hipc_merged_all_norm_filtered <- hipc_merged_all_norm %>%
  filter(
    !is.na(immResp_MFC_anyAssay_log2_MFC),
    vaccine_name == "Yellow Fever (LV)",
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
p_rise_results_list_yellowfeverlv = fs::path(
  "output",
  "results",
  "surrogacy",
  "rise_results_list_yellowfeverlv.rds"
)

# Save the results
readRDS(p_rise_results_list_yellowfeverlv)

# Save the results
saveRDS(rise_results_list,
        file = p_rise_results_list_yellowfeverlv)

# Read back in the results if necessary
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



### FIRST EXPLORATION OF RESULTS ###
# Calculate a weighted average of absolute delta values across studies for each timepoint
# The weights will be the square root of the number of observations
# such that larger samples provide more evidence

# First, derive weighted averages
rise_summary <- rise_results_df %>%
  group_by(marker, study_time_collected) %>%
  summarise(
    w_sum = sum(n, na.rm = TRUE),
    weighted_delta = if (w_sum > 0) sum(abs(delta) * n, na.rm = TRUE) / w_sum else NA_real_,
    n_studies = n_distinct(study_accession),
    .groups = "drop"
  )

# Parameters
timepoints_to_retain = c(3,7,14)
n_markers <- 150
cutoff <- 0.9  # define the cutoff

# Prepare plotting data
df_plot <- rise_summary %>%
  filter(study_time_collected %in% timepoints_to_retain) %>%
  mutate(plot_value = 1 - weighted_delta) %>%
  group_by(study_time_collected) %>%
  arrange(desc(plot_value)) %>%
  slice(1:n_markers) %>%
  ungroup() %>%
  mutate(
    timepoint = factor(study_time_collected, levels = timepoints_to_retain),
    marker_reordered = reorder_within(marker, plot_value, timepoint, .desc = TRUE),
    color_flag = ifelse(plot_value > cutoff, "above_cutoff", "below_cutoff")
  )

# Build labels for facets: "Day x - n studies"
# This assumes rise_summary has a column `n_studies` giving the number of studies per timepoint.
# If not available, replace the summarise(...) line below with an alternative (e.g. count distinct study IDs).
labels_vec <- rise_summary %>%
  filter(study_time_collected %in% timepoints_to_retain) %>%
  group_by(study_time_collected) %>%
  summarise(n_studies = first(n_studies)) %>%    # replace with an appropriate summary if needed
  ungroup() %>%
  mutate(label = paste0("Day ", study_time_collected, " (", n_studies, " studies)")) %>%
  { setNames(.$label, as.character(.$study_time_collected)) }

plot_min <- min(df_plot$plot_value, na.rm = TRUE) - 0.05

p1 <- ggplot(df_plot, aes(x = fct_rev(marker_reordered), y = plot_value, color = color_flag)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_hline(yintercept = cutoff, linetype = "dashed", color = "red", size = 0.7) +
  scale_color_manual(values = c("above_cutoff" = "red", "below_cutoff" = "#2C7BB6"), guide = "none") +
  labs(
    title = "Markers ranked by weighted average trial-level surrogacy",
    x = "Marker",
    ## y label with Greek Delta and subscript s:
    y = expression(1 - "weighted average of " * delta[s])
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
  coord_cartesian(ylim = c(plot_min, 1)) +
  facet_wrap(~ timepoint, ncol = 1, scales = "free_x",
             labeller = labeller(timepoint = labels_vec)) +
  scale_x_reordered()

# Print plot
print(p1)

surrogacy_figures_folder = fs::path("output", "figures", "surrogacy")

ggsave(
  filename = "yellowfeverlv_surrogacy_weightedranks.pdf",
  path = surrogacy_figures_folder,
  plot = p1,
  width = 40,
  height = 30,
  units = "cm"
)

# EXPLORING TRIAL-LEVEL SURROGACY WITH CROSS-TRIAL METRICS#

# params
for (tp in timepoints_to_retain) {
  top_markers = 15
  
  # Define a vector of genes which have high average trial-level
  # Surrogacy at any timepoint
  markers_to_keep <- df_plot %>%
    filter(study_time_collected == tp) %>%
    arrange(weighted_delta) %>%
    slice_head(n = top_markers) %>%
    pull(marker) %>%
    unique()
  
  markers_to_keep
  
  # Filter data for the chosen timepoint and markers
  rise_results_df_filtered <- rise_results_df %>%
    filter(marker %in% markers_to_keep, study_time_collected == tp) %>%
    # ensure marker factor ordering matches markers_to_keep
    mutate(marker = factor(marker, levels = markers_to_keep))
  
  # decide number of columns for facet grid
  n_markers <- length(markers_to_keep)
  ncol_facets <- min(5, n_markers)  # adjust 5 -> any default you prefer
  
  plot_min = min(rise_results_df_filtered$u_y,
                 rise_results_df_filtered$u_s) - 0.1
  # plot_min = 0.5
  
  # 1) Get min and max of n
  n_min <- min(rise_results_df_filtered$n, na.rm = TRUE)
  n_max <- max(rise_results_df_filtered$n, na.rm = TRUE)
  
  # 2) Round min up and max down to nearest 10
  legend_min <- ceiling(n_min / 10) * 10
  legend_max <- floor(n_max / 10) * 10
  
  # 3) Split interval into 4 equally spaced points: min, value1, value2, max
  legend_ns <- seq(legend_min, legend_max, length.out = 4)
  
  # 4) Round intermediate values to nearest 10
  legend_ns[c(2, 3)] <- round(legend_ns[c(2, 3)] / 50) * 50
  
  p2 <- ggplot(rise_results_df_filtered, aes(x = u_y, y = u_s)) +
    # point size mapped to n so radius ~ n (area ~ n)
    geom_point(aes(size = n), alpha = 0.5) +
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
        "Treatment effects on response vs marker across trials (Day ",
        tp,
        ": N trials = ",
        length(unique(
          rise_results_df_filtered$study_accession
        )),
        ")"
      ),
      subtitle = paste0(
        "Top ",
        n_markers,
        " markers by weighted average of trial-level delta"
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
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(
        size = 20,
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 20),
      axis.text  = element_text(size = 12),
      
      # increase spacing between facets so each panel has more room
      panel.spacing = unit(1.5, "lines"),
      
      # draw a border around each facet panel and ensure a white background
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 0.5
      ),
      panel.background = element_rect(fill = "white", colour = NA),
      
      # make facet strip noticeable
      strip.background = element_rect(
        fill = "grey95",
        colour = "black",
        size = 0.3
      ),
      strip.text = element_text(size = 13, face = "bold"),
      
      # legend styling
      legend.key = element_rect(fill = "white", colour = NA),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    facet_wrap( ~ marker, ncol = ncol_facets, scales = "fixed")
  
  print(p2)
  
  
  ggsave(
    paste0("yellowfeverlv_gene_grid_day", tp, ".pdf"),
    path = surrogacy_figures_folder,
    plot = p2,
    device = cairo_pdf,
    width = 35,
    height = 25,
    units = "cm"
  )
}
