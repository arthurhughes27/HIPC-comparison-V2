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

# Paths to processed gene-level data
p_load_expr_all_norm <- fs::path(processed_data_folder, "hipc_merged_all_norm.rds")

# Load data object
hipc_merged_all_norm <- readRDS(p_load_expr_all_norm)

# Timepoints of interest (numeric)
timepoints_of_interest <- c(0, 1, 3, 7)

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
# The top level elements will be the studies
# Inside each study will be results from 1 timepoint
rise_results_list = vector("list", length(study_names))
names(rise_results_list) = study_names

for(study in study_names) {
  # For each study
  # Print progress message
  message(paste0("Analysing study ", study, "..."))
  
  # Filter data by study
  hipc_merged_all_norm_filtered_study = hipc_merged_all_norm_filtered %>%
    filter(study_accession == study)
  
  # Extract timepoints for the study (must have at least 4 observations i.e. 2 pairs)
  timepoints <- hipc_merged_all_norm_filtered_study %>%
    filter(study_time_collected > 0) %>%
    count(study_time_collected) %>%
    filter(n >= 4) %>%
    pull(study_time_collected) %>%
    sort()
  
  # Initialise within-study results
  rise_results_list[[study]] = vector("list", length(timepoints))
  names(rise_results_list[[study]]) = paste0("Day ", timepoints)
  
  for (tp in timepoints) {
    # For each available timepoint
    # Progress message
    message(paste0("Day ", tp, " in progress."))
    
    # Filter by timepoint (post-vaccination and baseline)
    hipc_merged_all_norm_filtered_study_tp = hipc_merged_all_norm_filtered_study %>%
      filter(study_time_collected %in% c(0, tp)) %>%
      group_by(participant_id) %>%
      filter(all(c(0, tp) %in% study_time_collected)) %>%
      ungroup()
    
    # Arrange data by participant identifier and select relevant columns
    rise_df = hipc_merged_all_norm_filtered_study_tp %>%
      dplyr::select(
        participant_id,
        study_time_collected,
        immResp_MFC_anyAssay_pre_value,
        immResp_MFC_anyAssay_post_value,
        all_of(gene_names)
      ) %>%
      arrange(participant_id)
    
    # Extract pre-vaccination primary response
    yzero = rise_df %>%
      filter(study_time_collected == 0) %>%
      pull(immResp_MFC_anyAssay_pre_value)
    
    # Extract post-vaccination primary response
    yone = rise_df %>%
      filter(study_time_collected == tp) %>%
      pull(immResp_MFC_anyAssay_post_value)
    
    # Extract pre-vaccination surrogate candidates
    szero = rise_df %>%
      filter(study_time_collected == 0) %>%
      dplyr::select(all_of(gene_names))
    
    # Extract post-vaccination surrogate candidates
    sone = rise_df %>%
      filter(study_time_collected == tp) %>%
      dplyr::select(all_of(gene_names))
    
    # Extract the estimated effect on the primary response
    yresult = SurrogateRank::test.surrogate.extension(
      yone = yone,
      yzero = yzero,
      sone = yone,
      szero = yzero,
      power.want.s = 0.8,
      paired = T,
      alternative = "two.sided"
    )[["u.y"]]
    
    # Screen all markers with the paired, two one-sided test
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
    
    # Store the estimated treatment effects in the results
    rise_result[["screening.metrics"]]$u_y = yresult
    rise_result[["screening.metrics"]]$u_s = yresult - rise_result[["screening.metrics"]]$delta
    
    # Store the number of observations
    rise_result[["screening.metrics"]]$n = 2 * length(yone)
    
    # Store the results in the list element
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
# Calculate a weighted average of absolute delta values across studies for each timepoint
# The weights will be the number of observations
# We call this metric weighted mean absolute error (wMAE)

# First, derive weighted averages
rise_summary <- rise_results_df %>%
  group_by(marker, study_time_collected) %>%
  summarise(
    w_sum = sum(n, na.rm = TRUE),
    wmae = sum(abs(delta) * n, na.rm = TRUE) / w_sum,
    mae = mean(abs(delta)),
    n_studies = n_distinct(study_accession),
    .groups = "drop"
  )

# Parameters
timepoints_to_retain <- c(1, 3, 7)
n_markers <- 150
cutoff <- 0.9  # define the cutoff

# Prepare plotting data
df_plot <- rise_summary %>%
  filter(study_time_collected %in% timepoints_to_retain) %>%
  mutate(plot_value = 1 - wmae) %>%
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
    title = paste0("Top ", n_markers, " markers ranked by weighted average trial-level surrogacy"),
    x = "Marker",
    ## y label with Greek Delta and subscript s:
    y = expression(1 - "wMAE")
  ) +
  theme_bw(base_size = 25) +
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
  filename = "influenzain_surrogacy_weightedranks.pdf",
  path = surrogacy_figures_folder,
  plot = p1,
  width = 40,
  height = 30,
  units = "cm"
)

# EXPLORING TRIAL-LEVEL SURROGACY WITH CROSS-TRIAL METRICS#

# params
for (tp in timepoints_to_retain) {
  top_markers <- 20
  
  # Define a vector of genes ordered by wmae (ascending). use distinct() to keep order.
  markers_to_keep <- df_plot %>%
    filter(study_time_collected == tp) %>%
    arrange(wmae) %>%
    slice_head(n = top_markers) %>%
    distinct(marker, .keep_all = TRUE) %>%   # ensure unique but preserve order
    pull(marker)
  
  # Compute wMAE per marker for annotation
  label_df <- df_plot %>%
    filter(study_time_collected == tp, marker %in% markers_to_keep) %>%
    group_by(marker) %>%
    summarise(wmae = mean(wmae, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(wmae_fmt = paste0("wMAE = ", formatC(wmae, format = "f", digits = 3)))
  
  # Filter data and explicitly set factor levels to enforce ordering in facets
  rise_results_df_filtered <- rise_results_df %>%
    filter(marker %in% markers_to_keep, study_time_collected == tp) %>%
    mutate(marker = factor(marker, levels = markers_to_keep, ordered = TRUE))
  
  # decide number of columns for facet grid
  n_markers <- length(markers_to_keep)
  ncol_facets <- min(5, n_markers)
  
  plot_min <- min(rise_results_df_filtered$u_y, rise_results_df_filtered$u_s) - 0.1
  
  # legend n-breaks
  n_min <- min(rise_results_df_filtered$n, na.rm = TRUE)
  n_max <- max(rise_results_df_filtered$n, na.rm = TRUE)
  legend_min <- ceiling(n_min / 10) * 10
  legend_max <- floor(n_max / 10) * 10
  legend_ns <- seq(legend_min, legend_max, length.out = 4)
  legend_ns[c(2, 3)] <- round(legend_ns[c(2, 3)] / 50) * 50
  
  p2 <- ggplot(rise_results_df_filtered, aes(x = u_y, y = u_s)) +
    geom_point(aes(size = n), alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 0.8) +
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
        length(unique(rise_results_df_filtered$study_accession)),
        ")"
      ),
      subtitle = paste0("Top ", n_markers, " markers by weighted average of trial-level delta"),
      size = "Trial N"
    ) +
    scale_size_continuous(range = c(1.5, 7), breaks = legend_ns, labels = as.character(legend_ns)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 12),
      panel.spacing = unit(1.5, "lines"),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      panel.background = element_rect(fill = "white", colour = NA),
      strip.background = element_rect(fill = "grey95", colour = "black", size = 0.3),
      strip.text = element_text(size = 13, face = "bold"),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    facet_wrap(~ marker, ncol = ncol_facets, scales = "fixed")
  
  # Add red wMAE text inside each facet (top-left).
  # Make sure the marker column in text_df is a factor with the same levels.
  text_df <- label_df %>%
    filter(marker %in% levels(rise_results_df_filtered$marker)) %>%
    mutate(
      marker = factor(marker, levels = markers_to_keep, ordered = TRUE),
      x = plot_min + 0.03 * (1.01 - plot_min),
      y = 1.01 - 0.03 * (1.01 - plot_min)
    )
  
  p2 <- p2 +
    geom_text(
      data = text_df,
      aes(x = x, y = y, label = wmae_fmt, group = marker),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 4,
      color = "red"
    )
  
  print(p2)
  
  ggsave(
    paste0("influenzain_gene_grid_day", tp, ".pdf"),
    path = surrogacy_figures_folder,
    plot = p2,
    device = cairo_pdf,
    width = 35,
    height = 25,
    units = "cm"
  )
}


# ---------------------------
# Parameters
# ---------------------------
top_markers <- 30
timepoints_to_retain <- c(1,3,7)
highlight_threshold <- 0.2
scale_power <- 1.8
sep <- "___"

# ---------------------------
# Select top markers
# ---------------------------
markers_to_keep_heatmap <- df_plot %>%
  filter(study_time_collected %in% timepoints_to_retain) %>% 
  arrange(wmae) %>%
  slice_head(n = top_markers) %>%
  pull(marker) %>%
  unique()

rise_results_df_filtered <- rise_results_df %>% 
  filter(study_time_collected %in% timepoints_to_retain,
         marker %in% markers_to_keep_heatmap)

# ---------------------------
# Pivot to wide and numeric matrix
# ---------------------------
df_wide <- rise_results_df_filtered %>%
  filter(!is.na(delta)) %>%
  mutate(col_id = paste0(study_accession, sep, study_time_collected)) %>%
  select(marker, col_id, delta, p_adjusted) %>%
  pivot_wider(names_from = col_id, values_from = c(delta, p_adjusted))

# Extract matrices
mat <- as.matrix(df_wide %>% select(starts_with("delta_")))
rownames(mat) <- df_wide$marker

p_mat <- as.matrix(df_wide %>% select(starts_with("p_adjusted_")))
rownames(p_mat) <- df_wide$marker

# Keep only non-empty columns
keep_cols <- colSums(!is.na(mat)) > 0
mat <- mat[, keep_cols, drop = FALSE]
p_mat <- p_mat[, keep_cols, drop = FALSE]
if (ncol(mat) == 0) stop("No columns with data remain after filtering.")

# Column metadata
col_full <- colnames(mat)
col_study <- sub(paste0(sep, ".*$"), "", col_full)  # remove suffix after ___
col_study <- sub("^delta_", "", col_study)           # remove delta_ prefix
col_time_raw <- sub(paste0("^.*", sep), "", col_full)
col_time_num <- as.integer(col_time_raw)
col_time_factor <- factor(paste0("Day ", col_time_num),
                          levels = paste0("Day ", timepoints_to_retain))

# ----- NEW: lookup sample sizes (n) per study × time and build labeled study names -----
# use rise_results_df as source of n; choose the maximum n if multiple rows exist
n_lookup <- rise_results_df %>%
  dplyr::distinct(study_accession, study_time_collected, n) %>%
  dplyr::group_by(study_accession, study_time_collected) %>%
  dplyr::summarise(n_avail = max(n, na.rm = TRUE), .groups = "drop")

# create a small data.frame for columns and join to get n
cols_df <- tibble::tibble(
  col_full = col_full,
  study = col_study,
  time = col_time_num
) %>%
  dplyr::left_join(n_lookup, by = c("study" = "study_accession", "time" = "study_time_collected")) %>%
  dplyr::mutate(
    # display "N = NA" if missing; you can change this to e.g. omit (N = ...) if preferred
    study_label = ifelse(is.na(n_avail),
                         paste0(study, " (N = NA)"),
                         paste0(study, " (N = ", n_avail, ")"))
  )
# -------------------------------------------------------------------------------

# Reorder columns
ord_cols <- order(match(col_time_num, timepoints_to_retain), col_study)
mat <- mat[, ord_cols, drop = FALSE]
p_mat <- p_mat[, ord_cols, drop = FALSE]
col_study <- col_study[ord_cols]
col_time_factor <- factor(paste0("Day ", col_time_num[ord_cols]), levels = paste0("Day ", timepoints_to_retain))

# ----- Replace displayed column names with labeled versions that include N -----
# `cols_df$study_label` correspond to col_full order, so reindex and apply ordering
col_labels_ordered <- cols_df$study_label[ord_cols]
colnames(mat) <- col_labels_ordered
colnames(p_mat) <- col_labels_ordered
# -----------------------------------------------------------------------------

# ---------------------------
# Color mapping
# ---------------------------
color_transform <- function(x, threshold = highlight_threshold, power = scale_power) {
  y <- exp(- (abs(x)/threshold)^power)
  return(y)
}
col_fun <- function(x) {
  x_trans <- color_transform(x)
  colorRamp2(c(0,1), c("white", "#C40014"))(x_trans)
}

# ---------------------------
# Create significance symbols
# ---------------------------
sig_mat <- ifelse(p_mat < 0.05, "*", "")

ht_rise <- Heatmap(
  mat,
  name = "Delta",
  col = col_fun,
  na_col = "grey95",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_split = col_time_factor,
  cluster_column_slices = FALSE,
  column_dend_side = "top",
  column_dend_height = unit(3, "cm"),
  column_gap = unit(5, "mm"),
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 8),
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  row_dend_width = unit(3, "cm"),
  row_title = "Gene",
  show_heatmap_legend = F,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(sig_mat[i,j] == "*") {
      grid.text("*", x = x, y = y, gp = gpar(col = "black", fontsize = 14))
    }
  }
)

# ---------------------------
# 2) Create the color legend manually
# ---------------------------
max_abs <- 1
legend_at <- c(-max_abs, -highlight_threshold, 0, highlight_threshold, max_abs)
legend_labels <- formatC(legend_at, digits = 2, format = "f")

color_legend <- Legend(
  col_fun = col_fun,
  title = expression(delta[s]),
  at = legend_at,
  labels = legend_labels,
  legend_height = unit(3, "cm"),
  legend_width = unit(4, "cm"),
  border = TRUE,
  title_gp = gpar(fontsize = 16, fontface = "bold"),
  labels_gp = gpar(fontsize = 12)
)

# ---------------------------
# 3) Create the star legend
# ---------------------------
star_legend <- Legend(
  labels = expression(p[adj] < 0.05),,
  legend_gp = gpar(col = "black", fontsize = 14),
  title = "Significance",
  type = "points",
  pch = "*",
  size = unit(5, "mm")
)

# ---------------------------
# 4) Combine legends into a single object
# ---------------------------
combined_legend <- packLegend(star_legend, color_legend)

# ---------------------------
# 5) Draw heatmap with combined legend and a separate main title row
# ---------------------------

main_title <- paste0("Surrogacy heatmap - Influenza (IN): top ",  n_markers ," markers by wMAE across time")

pdf_file <- surrogacy_figures_folder / "heatmap_surrogacy_influenzain.pdf"
pdf(pdf_file, width = 14, height = 10)

# create a 2-row layout: top row for title, bottom row for the heatmap
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2,
                                           ncol = 1,
                                           heights = unit.c(unit(0.9, "in"), unit(1, "npc") - unit(0.9, "in")))))

# draw title in top row
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.text(main_title, gp = gpar(fontsize = 16, fontface = "bold"))
popViewport()

# draw heatmap in the second row without starting a new page (preserves split titles)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(ht_rise,
     annotation_legend_list = list(combined_legend),
     heatmap_legend_side = "right",
     newpage = FALSE)
popViewport()

dev.off()
