# Libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)


p_dearseq_dgsa_results_processed = fs::path("output",
                                            "results",
                                            "dearseq",
                                            "dearseq_dgsa_results_processed.rds")
p_qusage_dgsa_results_processed = fs::path("output",
                                           "results",
                                           "qusage",
                                           "qusage_dgsa_results_processed.rds")

results_dearseq = readRDS(p_dearseq_dgsa_results_processed)
results_qusage = readRDS(p_qusage_dgsa_results_processed)

results_df = bind_rows(results_dearseq, results_qusage)

rm(results_dearseq, results_qusage)

all_row_names <- unique(results_df$gs.name.description)
fixed_row_names_width <- max_text_width(all_row_names) + unit(2, "mm")

plot_circos = function(method_name = NULL,
                       conditions,
                       times,
                       aggregates,
                       fixed_row_names_width = 150,
                       p_correction = c("BH", "bonferroni", "holm", "hommel", "hochberg", "BY"),
                       p_approach = c("global", "withinTime", "withinComparison"),
                       p_threshold = 0.05,
                       filter_mode = c("none", "user", "data"),
                       user_threshold = 0.5,
                       quantile_threshold = 0.5,
                       scores = c("fc.score", "activation.score"),
                       y_order = c("cluster", "aggregate"),
                       x_order = c("cluster", "set"),
                       filter_commonDE = c(
                         "No filtration" = "none",
                         "Filter on across-time common DE" = "global",
                         "Filter on within-time common DE" = "withinTime",
                         "Filter by sharing score" = "score"
                       ),
                       common_proportion = 0,
                       score_threshold = 0,
                       quantile_scoreclip = 0.995,
                       legend_max = NULL) {
  # Copy results data into a new object
  results_df_heatmap <- results_df
  
  # Filter by chosen DGSA method
  if (!is.null(method_name)) {
    results_df_heatmap = results_df_heatmap %>%
      filter(method == method_name)
  }
  
  
  
  ## STEP 1 : DETERMINATION OF STATISTICAL SIGNIFICANCE ##
  ## We use two criteria for determining significance - p-value and effect size ##
  ## The thresholds and strategies for their calculation are determined by the user inputs ##
  ## This needs to be done on the entire data (i.e. regardless of which data are to be visualised later) ##
  
  # First, get names for the adjusted p-value and effect size columns of interest
  adj_Pval_column_name = paste0(p_approach, ".adjPval_", p_correction)
  
  filter_variable_column_name = scores
  
  if (filter_mode == "none") {
    # NO EFFECT SIZE FILTRATION CASE
    
    # Significance is uniquely determined by adjusted p-values under the threshold
    results_df_heatmap <- results_df_heatmap %>%
      mutate(significant = .data[[adj_Pval_column_name]] < p_threshold)
    
  } else if (filter_mode != "none") {
    # EFFECT SIZE FILTRATION CASES
    
    if (filter_mode == "data") {
      # DATA-DRIVEN EFFECT SIZE FILTRATION
      
      if (p_approach == "global") {
        # GLOBAL EFFECT SIZE FILTRATION CASE
        
        # Find the global effect size filtration threshold
        # In this case, this is the user-specified quantile of the absolute scores
        effect_size_threshold <- results_df_heatmap %>%
          pull(.data[[filter_variable_column_name]]) %>%
          abs() %>%
          quantile(quantile_threshold)
        
        # Filter the results by this effect size threshold
        results_df_heatmap <- results_df_heatmap %>%
          mutate(significant = (.data[[adj_Pval_column_name]] < p_threshold) &
                   (abs(.data[[filter_variable_column_name]]) > effect_size_threshold))
        
      } else if (p_approach == "withinTime") {
        # WITHIN-TIME EFFECT SIZE FILTRATION CASE
        
        # Find the within-time effect size filtration thresholds
        # In this case, this is the user-specified quantile of the absolute scores within each time
        thresholds <- results_df_heatmap %>%
          group_by(time) %>%
          summarise(threshold_abs = quantile(abs(.data[[filter_variable_column_name]]), quantile_threshold),
                    .groups = "drop")
        
        # Filter by within-time thresholds
        results_df_heatmap <- results_df_heatmap %>%
          left_join(thresholds, by = "time") %>%
          mutate(significant = (.data[[adj_Pval_column_name]] < p_threshold) &
                   (abs(.data[[filter_variable_column_name]]) > threshold_abs)) %>%
          select(-threshold_abs)
        
      } else if (p_approach == "withinComparison") {
        # WITHIN-COMPARISON EFFECT SIZE FILTRATION CASE
        
        # Find the within-comparison effect size filtration thresholds
        # In this case, this is the user-specified quantile of the absolute scores within each condition-time
        results_df_heatmap <- results_df_heatmap %>%
          group_by(condition, time) %>%
          mutate(local_effect = abs(.data[[filter_variable_column_name]]) >
                   quantile(abs(.data[[filter_variable_column_name]]), quantile_threshold, na.rm = TRUE)) %>%
          ungroup() %>%
          mutate(significant = (.data[[adj_Pval_column_name]] < p_threshold) &
                   local_effect) %>%
          select(-local_effect)
      }
      # END DATA-DRIVEN FILTRATION
    } else if (filter_mode == "user") {
      # USER-DEFINED THRESHOLD CASE
      
      # Here, we simply filter on the user-defined score thresholds
      results_df_heatmap <- results_df_heatmap %>%
        mutate(significant = (.data[[adj_Pval_column_name]] < p_threshold) &
                 (abs(.data[[filter_variable_column_name]]) > user_threshold))
    }
  } # END SIGNIFICANCE DETERMINATION
  
  # Convert significance variable to factor
  results_df_heatmap <- results_df_heatmap %>%
    mutate(significant = factor(significant, levels = c(TRUE, FALSE)))
  
  ## STEP 2: DATA SELECTION ##
  ## Here, we select the user-specified data (i.e. timepoints, conditions, aggregates)
  results_df_heatmap <- results_df_heatmap %>%
    filter(gs.aggregate %in% aggregates,
           condition %in% conditions,
           time %in% times) %>%
    mutate(time_label = factor(paste0(condition, " - Day ", time), levels = unique(paste0(
      condition, " - Day ", time
    ))))
  
  ## STEP 3: FILTER GENESETS BY COMMON DIFFERENTIAL EXPRESSION ##
  if (filter_commonDE == "global") {
    # GLOBAL COMMON DIFFERENTIAL EXPRESSION CASE
    
    # Filter by genesets which are commonly DE across a proportion of comparisons across selected times
    results_df_heatmap <- results_df_heatmap %>%
      group_by(gs.name) %>%
      mutate(frac_sig = sum(as.logical(significant)) / n()) %>%
      ungroup() %>%
      filter(frac_sig >= common_proportion) %>%
      select(-frac_sig)
    
  } else if (filter_commonDE == "withinTime") {
    # WITHIN-TIME COMMON DIFFERENTIAL EXPRESSION CASE
    
    # Compute within-time differential expression
    frac_tbl <- results_df_heatmap %>%
      group_by(gs.name, time) %>%
      summarise(frac_sig = sum(significant == "TRUE") / n(),
                .groups = "drop")
    
    # Extract genesets which pass the within-time thresholds for at least one selected time
    keep_sets <- frac_tbl %>%
      filter(frac_sig >= common_proportion) %>%
      pull(gs.name) %>%
      unique()
    
    # 3) filter your full data
    results_df_heatmap <- results_df_heatmap %>%
      filter(gs.name %in% keep_sets)
    
  } else if (filter_commonDE == "score") {
    # SHARING SCORES FILTRATION CASE
    
    # Sharing scores are the number of unique vaccines for which a geneset is DE in the same direction across all timepoints
    # and across selected vaccines
    sharing_score_df = results_df %>%
      filter(condition %in% conditions, method == method_name)
    
    calc_sharing_df_multi <- function(x, alpha) {
      map_dfr(adj_Pval_column_name, function(cur_col) {
        x %>%
          mutate(
            direction      = ifelse(.data[[filter_variable_column_name]] > 0, "up", "down"),
            is_significant = !is.na(.data[[cur_col]]) &
              .data[[cur_col]] < alpha
          ) %>%
          select(condition, gs.name, direction, is_significant) %>%
          distinct() %>%
          group_by(gs.name, direction) %>%
          summarise(n = sum(is_significant), .groups = "drop") %>%
          group_by(gs.name) %>%
          arrange(-n) %>%
          slice(1) %>%
          ungroup() %>%
          mutate(adjPval_method = cur_col)
      }) %>%
        pivot_wider(
          id_cols    = gs.name,
          names_from  = adjPval_method,
          values_from = n,
          values_fill = 0
        )
    }
    
    sharing_score_df <-
      sharing_score_df %>%
      calc_sharing_df_multi(alpha = p_threshold)
    
    keep_sets = sharing_score_df %>%
      filter(.data[[adj_Pval_column_name]] >= score_threshold) %>%
      pull(gs.name)
    
    results_df_heatmap <- results_df_heatmap %>%
      filter(gs.name %in% keep_sets)
  }
  
  
  if (nrow(results_df_heatmap) == 0) {
    stop("No genesets to plot under current parameters!")
  }
  
  ## STEP 4: AESTHETIC OPTIONS ##
  
  # Clip scores if desired
  # Find the threshold for score clipping
  clip_threshold  <- results_df_heatmap %>%
    pull(.data[[filter_variable_column_name]]) %>%
    abs() %>%
    quantile(quantile_scoreclip, na.rm = TRUE)
  
  # Clip the scores
  results_df_heatmap <- results_df_heatmap %>%
    mutate(
      !!filter_variable_column_name := case_when(
        .data[[filter_variable_column_name]] >  clip_threshold  ~  clip_threshold ,
        .data[[filter_variable_column_name]] < -clip_threshold  ~ -clip_threshold ,
        TRUE                                     ~  .data[[filter_variable_column_name]]
      )
    )
  
  
  # category_order = levels(results_df_heatmap$gs.aggregate)
  #
  # full_order <- levels(results_df_heatmap$condition)
  
  ## STEP 5: HEATMAP PLOTTING ##
  score_col <- scores
  results_df_heatmap <- results_df_heatmap %>% mutate(time = droplevels(time))
  time_order <- results_df_heatmap %>% pull(time) %>% levels()
  sep <- "___"
  
  df_wide <- results_df_heatmap %>%
    mutate(col_id = paste0(condition, sep, time)) %>%
    select(gs.name.description, col_id, !!rlang::sym(score_col)) %>%
    pivot_wider(names_from = col_id,
                values_from = !!rlang::sym(score_col))
  
  mat <- as.matrix(df_wide %>% select(-gs.name.description))
  rownames(mat) <- df_wide$gs.name.description
  col_full <- colnames(mat)
  col_time_raw <- sub(".*___", "", col_full)
  col_condition <- sub("___.*", "", col_full)
  
  time_to_daylabel <- function(t) {
    num <- stringr::str_extract(as.character(t), "\\d+")
    if (!is.na(num))
      paste0("Day ", num)
    else
      paste0("Day ", as.character(t))
  }
  time_display_levels <- vapply(time_order, time_to_daylabel, FUN.VALUE = character(1))
  col_time_display <- factor(vapply(col_time_raw, time_to_daylabel, FUN.VALUE = character(1)),
                             levels = time_display_levels)
  
  colnames(mat) <- col_condition
  all_aggregates <- levels(results_df$gs.aggregate)
  rows_meta <- results_df_heatmap %>%
    distinct(gs.name.description, gs.aggregate, gs.colour) %>%
    arrange(match(gs.name.description, rownames(mat)))
  row_groups_vec <- factor(rows_meta$gs.aggregate, levels = all_aggregates)
  names(row_groups_vec) <- rows_meta$gs.name.description
  row_group_cols <- results_df %>%
    distinct(gs.aggregate, gs.colour) %>%
    filter(gs.aggregate %in% all_aggregates) %>%
    {
      set_names(.$gs.colour, .$gs.aggregate)
    }
  
  row_ha <- rowAnnotation(
    `Geneset aggregate` = row_groups_vec,
    col = list(`Geneset aggregate` = row_group_cols),
    show_annotation_name = FALSE,
    width = unit(6, "mm")
  )
  
  cond_col_map <- results_df_heatmap %>% distinct(condition, condition.colour) %>% deframe()
  
  # --- ROBUST SIGNIFICANCE HANDLING FIX ---
  sig_wide <- results_df_heatmap %>%
    mutate(col_id = paste0(condition, sep, time)) %>%
    select(gs.name.description, col_id, significant) %>%
    pivot_wider(names_from = col_id, values_from = significant)
  
  sig_mat <- as.matrix(sig_wide %>% select(-gs.name.description))
  rownames(sig_mat) <- sig_wide$gs.name.description
  sig_mat <- matrix(
    as.logical(sig_mat),
    nrow = nrow(sig_mat),
    ncol = ncol(sig_mat),
    dimnames = dimnames(sig_mat)
  )
  
  sig_mat <- sig_mat[rownames(mat), col_full, drop = FALSE]
  colnames(sig_mat) <- col_condition
  sig_mat[is.na(sig_mat)] <- FALSE
  
  if (is.null(legend_max)) {
    vals <- mat[!is.na(mat)]
    max_abs <- max(abs(vals), na.rm = TRUE)
  } else {
    max_abs = legend_max
  }
  col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("blue", "white", "red"))
  
  cluster_rows_argument <- ifelse(y_order == "cluster", TRUE, FALSE)
  cluster_cols_argument <- ifelse(x_order == "cluster", TRUE, FALSE)
  
  if (!cluster_rows_argument) {
    ord <- match(row_groups_vec, levels(row_groups_vec))
    mat <- mat[order(ord), , drop = FALSE]
    row_groups_vec <- row_groups_vec[order(ord)]
    sig_mat <- sig_mat[order(ord), , drop = FALSE]
    show_row_dend <- FALSE
  } else
    show_row_dend <- TRUE
  
  row_ha <- rowAnnotation(
    `Geneset aggregate` = row_groups_vec,
    col = list(`Geneset aggregate` = row_group_cols),
    show_annotation_name = FALSE,
    width = unit(6, "mm")
  )
  
  top_ha <- HeatmapAnnotation(
    Vaccine = col_condition,
    col = list(Vaccine = cond_col_map),
    show_annotation_name = FALSE,
    annotation_height = unit(2, "mm")
  )
  
  ht <- Heatmap(
    mat,
    name = "Mean fold-change",
    col = col_fun,
    na_col = "grey95",
    cluster_rows = cluster_rows_argument,
    show_row_dend = show_row_dend,
    row_dend_side = "left",
    cluster_columns = cluster_cols_argument,
    column_split = col_time_display,
    cluster_column_slices = FALSE,
    show_column_dend = TRUE,
    column_dend_side = "top",
    column_dend_height = unit(8, "mm"),
    column_names_side = "bottom",
    top_annotation = top_ha,
    left_annotation = row_ha,
    show_column_names = TRUE,
    show_row_names = TRUE,
    heatmap_legend_param = list(title = "Mean fold-change"),
    row_names_max_width = fixed_row_names_width,
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (isTRUE(sig_mat[i, j])) {
        grid.text("*",
                  x = x,
                  y = y,
                  gp = gpar(fontsize = 10))
      }
    }
  )
  
  return(ht)
}


ht1 = plot_circos(
  method_name = "dearseq",
  conditions = levels(results_df$condition),
  times = c(1, 3, 7, 10),
  aggregates = levels(results_df$gs.aggregate)[-which(levels(results_df$gs.aggregate) == "NA")],
  fixed_row_names_width = unit(500, "mm"),
  p_correction = "BH",
  p_approach = "global",
  p_threshold = 0.05,
  scores = "fc.score",
  filter_mode = "none",
  user_threshold = 0.5,
  quantile_threshold = 0.5,
  y_order = "cluster",
  x_order = "cluster",
  filter_commonDE = "score",
  common_proportion = 0,
  score_threshold = 8,
  quantile_scoreclip = 0.995,
  legend_max = 2.2
)

ht2  = plot_circos(
  method_name = "qusage",
  conditions = levels(results_df$condition),
  times = c(1, 3, 7, 10),
  aggregates = levels(results_df$gs.aggregate)[-which(levels(results_df$gs.aggregate) == "NA")],
  fixed_row_names_width = unit(500, "mm"),
  p_correction = "BH",
  p_approach = "global",
  p_threshold = 0.05,
  scores = "fc.score",
  filter_mode = "none",
  user_threshold = 0.5,
  quantile_threshold = 0.5,
  y_order = "cluster",
  x_order = "cluster",
  filter_commonDE = "score",
  common_proportion = 0,
  score_threshold = 8,
  quantile_scoreclip = 0.995,
  legend_max = 2.2
)

figures_folder = fs::path("output", "figures", "dgsa")

# Open PDF
pdf(
  fs::path(figures_folder, "heatmap_stacked.pdf"),
  width = 20,
  height = 20
)

# Define grid layout: 2 rows, 1 column
pushViewport(viewport(layout = grid.layout(
  nrow = 2,
  ncol = 1,
  heights = unit.c(unit(1, "npc") * 0.5, unit(0.5, "npc") * 0.5)
)))

# Draw top heatmap
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(
  ht1,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  newpage = FALSE
)
upViewport()

# Draw bottom heatmap
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(
  ht2,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  newpage = FALSE
)
upViewport()

dev.off()
