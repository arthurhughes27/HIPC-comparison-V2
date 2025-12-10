# # File to find trial-level gene expression surrogates of influenza vaccination
# At each timepoint, we derive a signature from all but one study and test in on the left-out study

library(tidyverse)
library(SurrogateRank)

# Directory containing engineered / processed data files
processed_data_folder <- "data"

# Paths to processed gene-level data and gene-set objects
p_load_expr_all_norm <- fs::path(processed_data_folder, "hipc_merged_all_norm.rds")

# Load data objects
hipc_merged_all_norm <- readRDS(p_load_expr_all_norm)

# Timepoints of interest (numeric)
timepoints_to_keep <- c(0, 1, 3, 7)

# Filter to samples with non-missing immune response, Influenza vaccine,
# and collected at one of the specified timepoints.
hipc_merged_all_norm_filtered <- hipc_merged_all_norm %>%
  filter(
    !is.na(immResp_MFC_anyAssay_log2_MFC),
    vaccine_name == "Influenza (IN)",
    study_time_collected %in% timepoints_to_keep
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

timepoints_of_interest = c(1, 3, 7)

stability_selection_threshold = 0.5

# Initialise results vector
rise_studysignatures_results_list = vector("list", length(timepoints_of_interest))
names(rise_studysignatures_results_list) = paste0("Day ", timepoints_of_interest)

n_backup = 200

for (tp in timepoints_of_interest) {
  tp_string = paste0("Day ", tp)
  
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
  
  rise_studysignatures_results_list[[tp_string]] = vector("list", n_studies_tp)
  names(rise_studysignatures_results_list[[tp_string]]) = study_names_tp
  
  message(paste0(n_studies_tp, ' studies to analyse.'))
  
  # First derive fixed signatures for each study
  
  signature_list_train = vector("list", length(study_names_tp))
  names(signature_list_train) = study_names_tp
  
  for (sdy_train in study_names_tp) {
    # train on each of the other studies
    
    # Arrange data
    rise_df_train = hipc_merged_all_norm_filtered_tp %>%
      filter(study_accession == sdy_train) %>%
      dplyr::select(
        participant_id,
        study_time_collected,
        immResp_MFC_anyAssay_pre_value,
        immResp_MFC_anyAssay_post_value,
        all_of(gene_names)
      ) %>%
      arrange(participant_id)
    
    # Pre-vaccination response
    yzero_sdy_train = rise_df_train %>%
      filter(study_time_collected == 0) %>%
      pull(immResp_MFC_anyAssay_pre_value)
    
    # Post-vaccination response
    yone_sdy_train = rise_df_train %>%
      filter(study_time_collected == tp) %>%
      pull(immResp_MFC_anyAssay_post_value)
    
    # Pre-vaccination gene expression
    szero_sdy_train = rise_df_train %>%
      filter(study_time_collected == 0) %>%
      dplyr::select(all_of(gene_names))
    
    # Post-vaccination gene expression
    sone_sdy_train = rise_df_train %>%
      filter(study_time_collected == tp) %>%
      dplyr::select(all_of(gene_names))
    
    # Screen all markers
    rise_screen_result_train = rise.screen(
      yone = yone_sdy_train,
      yzero = yzero_sdy_train,
      sone = sone_sdy_train,
      szero = szero_sdy_train,
      alpha = 0.05,
      power.want.s = 0.8,
      # epsilon = 0.2,
      p.correction = "BH",
      n.cores = 1,
      alternative = "two.sided",
      paired = TRUE
    )
    
    if(length(rise_screen_result_train[["significant.markers"]]) == 0){
      
      threshold <- rise_screen_result_train[["screening.metrics"]] %>% 
        arrange(abs(delta)) %>% 
        slice_head(n = n_backup) %>%
        pull(delta) %>%
        abs() %>%
        max()
      
      sig_markers = rise_screen_result_train[["screening.metrics"]] %>% 
        filter(abs(delta) <= threshold) %>% 
        pull(marker)
      
    } else {
      sig_markers = rise_screen_result_train[["significant.markers"]]
    }
    
    signature_list_train[[sdy_train]] = sig_markers
    
  }
  
  for (sdy_eval in study_names_tp) {
    # for each study to be left out
    
    rise_studysignatures_results_list[[tp_string]][[sdy_eval]] = vector("list", 2)
    names(rise_studysignatures_results_list[[tp_string]][[sdy_eval]]) = c("signature", "evaluate")
    
    # Derive a signature not using this study based on stability selection
    
    signature_list_train_filtered = signature_list_train[-which(names(signature_list_train) == sdy_eval)]
    
    # derive signature
    signature_train <- enframe(signature_list_train_filtered,
                               name = "list_name",
                               value = "strings") %>%
      unnest(strings) %>%
      distinct(list_name, strings) %>%
      count(strings, name = "n_lists") %>%
      mutate(
        total_lists = length(signature_list_train_filtered),
        proportion = n_lists / total_lists
      ) %>%
      filter(proportion >= stability_selection_threshold) %>%
      pull(strings)
    
    if(length(signature_train) == 0){
      adaptive_threshold = enframe(signature_list_train_filtered,
                                   name = "list_name",
                                   value = "strings") %>%
        unnest(strings) %>%
        distinct(list_name, strings) %>%
        count(strings, name = "n_lists") %>%
        mutate(
          total_lists = length(signature_list_train_filtered),
          proportion = n_lists / total_lists
        ) %>% 
        pull(proportion) %>% 
        max()
      
      signature_train <- enframe(signature_list_train_filtered,
                                 name = "list_name",
                                 value = "strings") %>%
        unnest(strings) %>%
        distinct(list_name, strings) %>%
        count(strings, name = "n_lists") %>%
        mutate(
          total_lists = length(signature_list_train_filtered),
          proportion = n_lists / total_lists
        ) %>%
        filter(proportion >= adaptive_threshold) %>%
        pull(strings)
    }
    
    # Now evaluate that signature in the left out study
    
    # Arrange data
    rise_df_eval = hipc_merged_all_norm_filtered_tp %>%
      filter(study_accession == sdy_eval) %>%
      dplyr::select(
        participant_id,
        study_time_collected,
        immResp_MFC_anyAssay_pre_value,
        immResp_MFC_anyAssay_post_value,
        all_of(gene_names)
      ) %>%
      arrange(participant_id)
    
    # Pre-vaccination response
    yzero_sdy_eval = rise_df_eval %>%
      filter(study_time_collected == 0) %>%
      pull(immResp_MFC_anyAssay_pre_value)
    
    # Post-vaccination response
    yone_sdy_eval = rise_df_eval %>%
      filter(study_time_collected == tp) %>%
      pull(immResp_MFC_anyAssay_post_value)
    
    # Pre-vaccination gene expression
    szero_sdy_eval = rise_df_eval %>%
      filter(study_time_collected == 0) %>%
      dplyr::select(all_of(gene_names))
    
    # Post-vaccination gene expression
    sone_sdy_eval = rise_df_eval %>%
      filter(study_time_collected == tp) %>%
      dplyr::select(all_of(gene_names))
    
    # Extract the estimated effect on the response on the relevant samples
    yresult_sdy_eval = SurrogateRank::test.surrogate.extension(
      yone = yone_sdy_eval,
      yzero = yzero_sdy_eval,
      sone = yone_sdy_eval,
      szero = yzero_sdy_eval,
      power.want.s = 0.8,
      paired = T,
      alternative = "two.sided"
    )[["u.y"]]
    
    # Evaluate the signature
    rise_result_eval = rise.evaluate(
      yone = yone_sdy_eval,
      yzero = yzero_sdy_eval,
      sone = sone_sdy_eval,
      szero = szero_sdy_eval,
      alpha = 0.05,
      power.want.s = 0.8,
      # epsilon = 0.2,
      p.correction = "BH",
      n.cores = 1,
      alternative = "two.sided",
      paired = TRUE,
      markers = signature_train,
      evaluate.weights = F
    )
    
    rise_result_eval_df = data.frame(as.list(rise_result_eval[["gamma.s.evaluate"]]))
    
    # Store the treatment effects on each marker
    rise_result_eval_df$u_y = yresult_sdy_eval
    rise_result_eval_df$u_s = yresult_sdy_eval - rise_result_eval_df$delta
    
    # Store the number of observations
    rise_result_eval_df$n = 2 * length(yone_sdy_eval)
    
    rise_studysignatures_results_list[[tp_string]][[sdy_eval]][["signature"]] = signature_train
    rise_studysignatures_results_list[[tp_string]][[sdy_eval]][["evaluate"]] = rise_result_eval_df
  }
}

# result: a single data.frame with all evaluate rows stacked
all_evaluates_df <- imap_dfr(
  rise_studysignatures_results_list,
  function(studies_list, day_label) {
    # try to parse numeric day from "Day x" (NA if not parseable)
    day_numeric <- suppressWarnings(as.numeric(gsub("^Day\\s*", "", day_label)))
    
    # iterate studies inside this day
    imap_dfr(
      studies_list,
      function(study_elem, study_name) {
        ev <- study_elem[["evaluate"]]
        
        # skip if evaluate is NULL or has no rows
        if (is.null(ev) || (is.data.frame(ev) && nrow(ev) == 0)) return(NULL)
        
        # coerce to data.frame (handles tibbles)
        ev_df <- as.data.frame(ev, stringsAsFactors = FALSE)
        
        # ensure a study identifier column exists (use existing study_accession if present,
        # otherwise use the inner list name)
        if (!"study_accession" %in% names(ev_df)) {
          ev_df$study_accession <- study_name
        }
        
        # add time columns
        ev_df$timepoint_label <- day_label    # e.g. "Day 7"
        ev_df$timepoint_day   <- day_numeric  # e.g. 7 (NA if parse failed)
        
        ev_df
      }
    )
  }
) %>%
  select(study_accession, timepoint_label, timepoint_day, everything())


# 1) compute WMAE per timepoint (weights from n)
wmae_by_time <- all_evaluates_df %>%
  group_by(timepoint_label, timepoint_day) %>%
  summarise(
    total_n = sum(n[!is.na(delta) & !is.na(n)], na.rm = TRUE),
    wmae = if (total_n > 0) {
      valid <- !is.na(delta) & !is.na(n)
      weighted.mean(x = abs(delta[valid]), w = n[valid])
    } else {
      NA_real_
    },
    .groups = "drop"
  )

plots_by_day <- all_evaluates_df %>%
  # keep only rows that have at least one of u_y or u_s non-missing
  filter(!is.na(u_y) | !is.na(u_s)) %>%
  group_by(timepoint_label, timepoint_day) %>%
  group_split(.keep = TRUE) %>%
  set_names(map_chr(., ~ unique(.x$timepoint_label))) %>%
  map(function(df_day) {
    day_label <- unique(df_day$timepoint_label)
    day_num   <- unique(df_day$timepoint_day)
    
    # ensure numeric columns
    df_day <- df_day %>%
      mutate(
        u_y = as.numeric(u_y),
        u_s = as.numeric(u_s),
        n   = as.numeric(n)   # ensure numeric for size mapping
      )
    
    # compute axis minimum from data across both axes
    min_val <- min(c(df_day$u_y, df_day$u_s), na.rm = TRUE)
    if (!is.finite(min_val)) min_val <- 0
    
    # small margin based on range so points don't hug the border
    range_ref <- 1 - min_val
    if (!is.finite(range_ref) || range_ref <= 0) range_ref <- 0.1
    small_margin <- 0.1
    min_lim <- min_val - small_margin
    
    # increase upper limit slightly so points at 1 are visible
    upper_increase <- 0.015 * range_ref  # small absolute increase
    upper_lim <- 1 + upper_increase
    
    # sanity clamp: if min_lim >= upper_lim, set sensible defaults
    if (min_lim >= upper_lim) {
      min_lim <- 0
      upper_lim <- 1 + 0.02
    }
    
    # find WMAE for this day (if present)
    wmae_val <- wmae_by_time %>%
      filter(timepoint_label == day_label) %>%
      pull(wmae)
    wmae_val <- if (length(wmae_val) == 0) NA_real_ else wmae_val[[1]]
    
    # annotation position (upper-left)
    x_pos <- min_lim + 0.04 * (upper_lim - min_lim)
    y_pos <- upper_lim - 0.04 * (upper_lim - min_lim)
    
    # If some rows have NA for n, let them still plot with a small default size:
    # create a plotting-size variable so size legend remains meaningful.
    df_day <- df_day %>%
      mutate(n_plot = ifelse(is.na(n), 1, pmax(1, n))) # replace NA with 1, ensure >=1
    
    # determine legend breaks and limits from the actual n_plot values
    min_n <- min(df_day$n_plot, na.rm = TRUE)
    max_n <- max(df_day$n_plot, na.rm = TRUE)
    # use pretty to get 3-ish ticks, then ensure min/max are represented
    br <- pretty(c(min_n, max_n), n = 3)
    br <- unique(c(min_n, br[br >= min_n & br <= max_n], max_n))
    # if pretty returned only one value (rare), create sensible breaks
    if (length(br) == 1) br <- unique(round(c(min_n, (min_n + max_n)/2, max_n)))
    # ensure increasing numeric vector
    br <- sort(as.numeric(br))
    
    # plot
    p <- ggplot(df_day, aes(x = u_y, y = u_s)) +
      geom_point(aes(size = n_plot), alpha = 0.6) +
      geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 0.7) +
      coord_fixed(xlim = c(min_lim, upper_lim), ylim = c(min_lim, upper_lim), ratio = 1, expand = FALSE) +
      labs(
        title = paste0("Leave-One-Study-Out Cross Validation - ", day_label),
        # LaTeX-like axis labels using expression(): U_Y on x, U_S on y
        x = expression(U[Y]),
        y = expression(U[S])
      ) +
      # WMAE annotation: show NA nicely if missing
      annotate(
        "text",
        x = x_pos,
        y = y_pos,
        label = if (!is.na(wmae_val)) paste0("wMAE = ", formatC(wmae_val, digits = 3, format = "f")) else "wMAE = NA",
        hjust = 0,
        vjust = 1,
        color = "red",
        size = 5
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
  })


plots_by_day
