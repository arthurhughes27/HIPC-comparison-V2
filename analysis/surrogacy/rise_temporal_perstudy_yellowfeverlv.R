# File to find trial-level gene expression surrogates of Yellow Fever vaccination
# We do a paired-data analysis per-study to account for study effects

library(tidyverse)
library(SurrogateRank)
library(purrr)
library(stringr)

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
  
  timepoints = hipc_merged_all_norm_filtered_study %>%
    filter(study_time_collected > 0) %>%
    pull(study_time_collected) %>%
    unique() %>%
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


library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

## -----------------------
## User params
## -----------------------
p_limit <- 0.05            # significance threshold for p_adjusted
sharing_thresh <- 2        # keep genes with sharing_score >= this
sigma <- 0.21               # how tightly red is concentrated near zero; smaller => tighter
power <- 2                # exponent; >2 sharper near zero

## -----------------------
## Basic checks
## -----------------------
required_cols <- c("marker", "study_accession", "study_time_collected", "delta", "p_adjusted")
if (!all(required_cols %in% names(rise_results_df))) {
  stop("rise_results_df must contain columns: ", paste(required_cols, collapse = ", "))
}

## -----------------------
## Compute sharing score:
## - For each marker × study_accession: is there any timepoint with p_adjusted < p_limit?
## - For each marker: count the number of studies with any_sig == TRUE
## - Keep markers with sharing_score >= sharing_thresh
## -----------------------
sharing_df <- rise_results_df %>%
  filter(!is.na(p_adjusted)) %>%
  group_by(marker, study_accession) %>%
  summarise(any_sig_in_study = any(p_adjusted < p_limit, na.rm = TRUE), .groups = "drop") %>%
  group_by(marker) %>%
  summarise(sharing_score = sum(as.integer(any_sig_in_study), na.rm = TRUE), .groups = "drop")

kept_markers <- sharing_df %>%
  filter(sharing_score >= sharing_thresh) %>%
  pull(marker)

if (length(kept_markers) == 0) {
  stop("No markers pass the sharing_score filter (sharing_thresh = ", sharing_thresh, ").")
}

## -----------------------
## Filter the main results to the kept markers
## -----------------------
rise_results_df_filt <- rise_results_df %>%
  filter(marker %in% kept_markers)

## -----------------------
## 1) Build delta and p matrices aligned by (marker, col_id)
## -----------------------
processed <- rise_results_df_filt %>%
  mutate(
    study_time_collected = as.numeric(as.character(study_time_collected)),
    col_id = paste0("Day", study_time_collected, "_", study_accession)
  )

mat_df <- processed %>%
  group_by(marker, col_id) %>%
  summarise(delta = mean(delta, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = col_id, values_from = delta, values_fill = NA)

p_df <- processed %>%
  group_by(marker, col_id) %>%
  summarise(p_adj = min(p_adjusted, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = col_id, values_from = p_adj, values_fill = NA)

if (ncol(mat_df) <= 1) {
  stop("After filtering/aggregation no study×timepoint columns remain to plot.")
}

mat <- mat_df %>% column_to_rownames("marker") %>% as.matrix()
p_mat <- p_df %>% column_to_rownames("marker") %>% as.matrix()

## -----------------------
## 2) Column metadata from col_id
## -----------------------
col_meta <- tibble(col_id = colnames(mat)) %>%
  tidyr::separate(
    col_id,
    into = c("day_label", "study"),
    sep = "_",
    extra = "merge",
    remove = FALSE
  ) %>%
  mutate(study_time_collected = suppressWarnings(as.numeric(sub("Day", "", day_label))))

if (any(is.na(col_meta$study_time_collected))) {
  warning("Some column names did not parse to a numeric day. Check that col_id = 'DayX_study' is correct.")
}

## -----------------------
## 3) Order columns by numeric day ascending then study name
## -----------------------
col_meta <- col_meta %>% arrange(study_time_collected, study)
col_ord <- col_meta$col_id
mat <- mat[, col_ord, drop = FALSE]
p_mat <- p_mat[, col_ord, drop = FALSE]
col_meta <- col_meta %>% arrange(study_time_collected, study)

## -----------------------
## 4) Column names = study labels (shown under the heatmap)
## -----------------------
colnames(mat) <- col_meta$study
colnames(p_mat) <- col_meta$study

## -----------------------
## 5) Non-linear continuous colour mapping (map |delta| -> strength -> white/red)
## -----------------------
strength_to_col <- colorRamp2(c(0, 1), c("white", "red"))

col_fun <- function(x) {
  na_idx <- is.na(x)
  if (all(na_idx)) return(rep(NA_character_, length(x)))
  s <- exp(- (abs(x) / sigma) ^ power)  # 1 at x=0, decays with |x|
  cols <- strength_to_col(s)
  cols[na_idx] <- NA_character_
  cols
}

legend_breaks <- seq(-1, 1, length.out = 5)
legend_cols <- col_fun(legend_breaks)

## -----------------------
## 6) Significance mask for stars
## -----------------------
sig_mat <- (p_mat < p_limit)
sig_mat[is.na(sig_mat)] <- FALSE

## -----------------------
## 7) column_split with labels "Day x" (slice titles will be shown above column dendrogram)
## -----------------------
timepoint_levels <- sort(unique(col_meta$study_time_collected))
timepoint_levels <- timepoint_levels[!is.na(timepoint_levels)]
column_split_factor <- factor(col_meta$study_time_collected,
                              levels = timepoint_levels,
                              labels = paste0("Day ", timepoint_levels))

## -----------------------
## 8) Build heatmap (show_column_split = TRUE to display Day slice titles above dendrogram)
## -----------------------
ht <- Heatmap(
  mat,
  name = "delta",
  col = col_fun,
  na_col = "grey90",
  cluster_rows = F,
  show_row_dend = FALSE,
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  show_column_dend = TRUE,
  column_dend_side = "top",           # dendrogram above heatmap
  column_dend_height = unit(2, "cm"),
  column_split = column_split_factor,
  column_gap = unit(2, "mm"),
  column_names_side = "bottom",       # study labels under the heatmap
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  heatmap_legend_param = list(
    title = "delta",
    at = legend_breaks,
    labels = as.character(round(legend_breaks, 2))
  ),
  column_title = "Study Name",
  column_title_side = "bottom",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  row_title = "Gene",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (isTRUE(sig_mat[i, j])) {
      # compute a fontsize in points proportional to the cell height (in mm -> approximate pts)
      # convertHeight(height, "mm", valueOnly=TRUE) returns mm, 1 mm ≈ 2.83465 points
      cell_h_mm <- grid::convertHeight(height, "mm", valueOnly = TRUE)
      fontsize_pts <- round(pmin(10, pmax(4, cell_h_mm * 5)))    # tune multiplier if needed
      
      # small vertical offset (fraction of cell height) to compensate for font baseline
      v_off <- height * 0.06  # 6% of cell height; reduce/increase if needed
      
      grid::grid.text(
        "*",
        x = x,
        y = y - v_off,            # nudge downward slightly so it visually centers
        gp = grid::gpar(fontsize = fontsize_pts, fontface = "bold"),
        just = "centre"           # ensure centred justification
      )
    }
  }
)

# Draw the heatmap
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
