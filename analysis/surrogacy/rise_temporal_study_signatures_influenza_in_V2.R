# File to find trial-level gene expression surrogates of influenza vaccination
# At each timepoint, we derive a signature from the study with the largest sample size
# We then apply the signature to each study one by one and perform trial-level meta-analysis

library(tidyverse)
library(SurrogateRank)
library(patchwork)


# Directory containing engineered / processed data files
processed_data_folder <- "data"

# Paths to processed gene-level data and gene-set objects
p_load_expr_all_norm <- fs::path(processed_data_folder, "hipc_merged_all_norm.rds")

# Load data objects
hipc_merged_all_norm <- readRDS(p_load_expr_all_norm)

# Timepoints of interest (numeric)
timepoints_to_keep <- c(0, 1, 3, 7, 10, 14)

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

timepoints_of_interest = c(1, 7)

n_backup = 10

# Initialise results vector
rise_studysignatures_results_list = vector("list", length(timepoints_of_interest))
names(rise_studysignatures_results_list) = paste0("Day ", timepoints_of_interest)

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
  
  for (sdy_main in study_names_tp) {
    # For each of the available studies
    
    # Derive a signature based on this study with epsilon = 0.1
    
    # Filter the data
    hipc_merged_all_norm_filtered_tp_sdy_main = hipc_merged_all_norm_filtered_tp %>%
      filter(study_accession == sdy_main)
    
    n_sdy_main = hipc_merged_all_norm_filtered_tp_sdy_main %>%
      nrow()
    
    message(paste0(
      "Deriving signature from study ",
      sdy_main,
      " (N = ",
      n_sdy_main,
      ")"
    ))
    
    # Arrange data
    rise_df_screen_sdy_main = hipc_merged_all_norm_filtered_tp_sdy_main %>%
      dplyr::select(
        participant_id,
        study_time_collected,
        immResp_MFC_anyAssay_pre_value,
        immResp_MFC_anyAssay_post_value,
        all_of(gene_names)
      ) %>%
      arrange(participant_id)
    
    # Pre-vaccination response
    yzero_sdy_main = rise_df_screen_sdy_main %>%
      filter(study_time_collected == 0) %>%
      pull(immResp_MFC_anyAssay_pre_value)
    
    # Post-vaccination response
    yone_sdy_main = rise_df_screen_sdy_main %>%
      filter(study_time_collected == tp) %>%
      pull(immResp_MFC_anyAssay_post_value)
    
    # Pre-vaccination gene expression
    szero_sdy_main = rise_df_screen_sdy_main %>%
      filter(study_time_collected == 0) %>%
      dplyr::select(all_of(gene_names))
    
    # Post-vaccination gene expression
    sone_sdy_main = rise_df_screen_sdy_main %>%
      filter(study_time_collected == tp) %>%
      dplyr::select(all_of(gene_names))
    
    # Extract the estimated effect on the response on the relevant samples
    yresult_sdy_main = SurrogateRank::test.surrogate.extension(
      yone = yone_sdy_main,
      yzero = yzero_sdy_main,
      sone = yone_sdy_main,
      szero = yzero_sdy_main,
      power.want.s = 0.8,
      paired = T,
      alternative = "two.sided"
    )[["u.y"]]
    
    # Screen all markers
    rise_screen_result_sdy_main = rise.screen(
      yone = yone_sdy_main,
      yzero = yzero_sdy_main,
      sone = sone_sdy_main,
      szero = szero_sdy_main,
      alpha = 0.05,
      # power.want.s = 0.8,
      epsilon = 0.2,
      p.correction = "BH",
      n.cores = 1,
      alternative = "two.sided",
      paired = TRUE
    )
    
    rise_screen_result_sdy_main_df <- rise_screen_result_sdy_main[["screening.metrics"]] %>%
      arrange(abs(delta))
    
    # Number of significant markers
    n_signif_sdy_main = length(rise_screen_result_sdy_main[["significant.markers"]])
    
    if (n_signif_sdy_main > 0) {
      message(
        paste0(
          n_signif_sdy_main,
          " significant markers found. Evaluating this signature in all studies."
        )
      )
      
      signature_sdy_main = rise_screen_result_sdy_main[["significant.markers"]]
      
    } else {
      message(
        paste0(
          "No significant markers found. Taking top ",
          n_backup,
          " genes (+ ties) and evaluating this signature in all studies."
        )
      )
      
      threshold <- rise_screen_result_sdy_main_df %>%
        slice_head(n = n_backup) %>%
        pull(delta) %>%
        abs() %>%
        max()
      
      signature_sdy_main <- rise_screen_result_sdy_main_df %>%
        filter(abs(delta) <= threshold) %>%
        pull(marker)
    }
    
    # First, evaluate the signature in the current study
    
    min_delta = rise_screen_result_sdy_main_df %>%
      filter(marker %in% signature_sdy_main, abs(delta) > 0) %>%
      pull(delta) %>%
      min()
    
    screening_weights_sdy_main = rise_screen_result_sdy_main_df %>%
      filter(marker %in% signature_sdy_main) %>%
      mutate(weight = ifelse(delta != 0, 1 / abs(delta), 1 / abs(min_delta))) %>%
      select(marker, weight)
    
    
    if (is.infinite(min_delta)) {
      rise_evaluate_result_sdy_main = rise.evaluate(
        yone = yone_sdy_main,
        yzero = yzero_sdy_main,
        sone = sone_sdy_main,
        szero = szero_sdy_main,
        alpha = 0.05,
        # power.want.s = 0.8,
        epsilon = 0.2,
        p.correction = "BH",
        n.cores = 1,
        alternative = "two.sided",
        paired = TRUE,
        markers = signature_sdy_main,
        evaluate.weights = FALSE
      )
    } else {
      rise_evaluate_result_sdy_main = rise.evaluate(
        yone = yone_sdy_main,
        yzero = yzero_sdy_main,
        sone = sone_sdy_main,
        szero = szero_sdy_main,
        alpha = 0.05,
        # power.want.s = 0.8,
        epsilon = 0.2,
        p.correction = "BH",
        n.cores = 1,
        alternative = "two.sided",
        paired = TRUE,
        markers = signature_sdy_main,
        evaluate.weights = T,
        screening.weights = screening_weights_sdy_main
      )
    }
    
    rise_evaluate_result_sdy_main_df = data.frame(as.list(rise_evaluate_result_sdy_main[["gamma.s.evaluate"]]))
    
    # Store the treatment effects on each marker
    rise_evaluate_result_sdy_main_df$u_y = yresult_sdy_main
    rise_evaluate_result_sdy_main_df$u_s = yresult_sdy_main - rise_evaluate_result_sdy_main_df$delta
    
    # Store the number of observations
    rise_evaluate_result_sdy_main_df$n = 2 * length(yone_sdy_main)
    
    rise_evaluate_result_sdy_main_df$study_accession = sdy_main
    
    # Now evaluate all other studies
    for (sdy_alt in study_names_tp[study_names_tp != sdy_main]) {
      
      # Filter the data
      hipc_merged_all_norm_filtered_tp_sdy_alt = hipc_merged_all_norm_filtered_tp %>%
        filter(study_accession == sdy_alt)
      
      n_sdy_alt = hipc_merged_all_norm_filtered_tp_sdy_alt %>%
        nrow()
      
      message(paste0(
        "Deriving signature from study ",
        sdy_alt,
        " (N = ",
        n_sdy_alt,
        ")"
      ))
      
      # Arrange data
      rise_df_screen_sdy_alt = hipc_merged_all_norm_filtered_tp_sdy_alt %>%
        dplyr::select(
          participant_id,
          study_time_collected,
          immResp_MFC_anyAssay_pre_value,
          immResp_MFC_anyAssay_post_value,
          all_of(gene_names)
        ) %>%
        arrange(participant_id)
      
      # Pre-vaccination response
      yzero_sdy_alt = rise_df_screen_sdy_alt %>%
        filter(study_time_collected == 0) %>%
        pull(immResp_MFC_anyAssay_pre_value)
      
      # Post-vaccination response
      yone_sdy_alt = rise_df_screen_sdy_alt %>%
        filter(study_time_collected == tp) %>%
        pull(immResp_MFC_anyAssay_post_value)
      
      # Pre-vaccination gene expression
      szero_sdy_alt = rise_df_screen_sdy_alt %>%
        filter(study_time_collected == 0) %>%
        dplyr::select(all_of(gene_names))
      
      # Post-vaccination gene expression
      sone_sdy_alt = rise_df_screen_sdy_alt %>%
        filter(study_time_collected == tp) %>%
        dplyr::select(all_of(gene_names))
      
      # Extract the estimated effect on the response on the relevant samples
      yresult_sdy_alt = SurrogateRank::test.surrogate.extension(
        yone = yone_sdy_alt,
        yzero = yzero_sdy_alt,
        sone = yone_sdy_alt,
        szero = yzero_sdy_alt,
        power.want.s = 0.8,
        paired = T,
        alternative = "two.sided"
      )[["u.y"]]
      
      if(is.infinite(min_delta)){
        rise_evaluate_result_sdy_alt = rise.evaluate(
          yone = yone_sdy_alt,
          yzero = yzero_sdy_alt,
          sone = sone_sdy_alt,
          szero = szero_sdy_alt,
          alpha = 0.05,
          # power.want.s = 0.8,
          epsilon = 0.2,
          p.correction = "BH",
          n.cores = 1,
          alternative = "two.sided",
          paired = TRUE,
          markers = signature_sdy_main,
          evaluate.weights = F
        )
      } else {
        rise_evaluate_result_sdy_alt = rise.evaluate(
          yone = yone_sdy_alt,
          yzero = yzero_sdy_alt,
          sone = sone_sdy_alt,
          szero = szero_sdy_alt,
          alpha = 0.05,
          # power.want.s = 0.8,
          epsilon = 0.2,
          p.correction = "BH",
          n.cores = 1,
          alternative = "two.sided",
          paired = TRUE,
          markers = signature_sdy_main,
          evaluate.weights = T,
          screening.weights = screening_weights_sdy_main
        )
      }
    
      
      rise_evaluate_result_sdy_alt_df = data.frame(as.list(rise_evaluate_result_sdy_alt[["gamma.s.evaluate"]]))
      
      # Store the treatment effects on each marker
      rise_evaluate_result_sdy_alt_df$u_y = yresult_sdy_alt
      rise_evaluate_result_sdy_alt_df$u_s = yresult_sdy_alt - rise_evaluate_result_sdy_alt_df$delta
      
      # Store the number of observations
      rise_evaluate_result_sdy_alt_df$n = 2 * length(yone_sdy_alt)
      
      rise_evaluate_result_sdy_alt_df$study_accession = sdy_alt
      
      # bind to previous dataframe
      
      rise_evaluate_result_sdy_main_df = rbind(rise_evaluate_result_sdy_main_df,
                                               rise_evaluate_result_sdy_alt_df)
    }
    
    rise_studysignatures_results_list[[tp_string]][[sdy_main]] = rise_evaluate_result_sdy_main_df
  }
}

# Specify path for saving list of results
p_results = fs::path("output", "results", "surrogacy", "rise_studysignatures_results_list.rds")

# Save results
saveRDS(rise_studysignatures_results_list, file = p_results)

# rise_studysignatures_results_list = readRDS(p_results)

main_study = "SDY61"

combined_df  = rise_studysignatures_results_list[["Day 7"]][[main_study]]

main_label <- paste0("Training study")

combined_df <- combined_df %>%
  mutate(
    is_main_flag = ifelse(study_accession == main_study, main_label, "Evaluation studies"),
    is_main_flag = factor(is_main_flag, levels = c("Evaluation studies", main_label))
  )

plot_min <- 0

# n legend calculation (your original logic)
n_min <- min(combined_df$n, na.rm = TRUE)
n_max <- max(combined_df$n, na.rm = TRUE)
legend_min <- ceiling(n_min / 10) * 10
legend_max <- floor(n_max / 10) * 10
legend_ns <- seq(legend_min, legend_max, length.out = 4)
legend_ns[c(2, 3)] <- round(legend_ns[c(2, 3)] / 50) * 50

n_study_main = combined_df %>% 
  filter(study_accession == main_study) %>% 
  pull(n)

p <- ggplot(combined_df, aes(x = u_y, y = u_s)) +
  geom_point(
    aes(size = n, fill = is_main_flag, colour = is_main_flag),
    shape = 21,
    alpha = 0.5,
    stroke = 0.6
  ) +
  geom_abline(slope = 1, intercept = 0, color = "#FF2128", linetype = "dashed", size = 0.8) +
  scale_x_continuous(limits = c(plot_min, 1.01), expand = c(0, 0)) +
  scale_y_continuous(limits = c(plot_min, 1.01), expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  labs(
    x = expression(U[Y]),
    y = expression(U[S]),
    title = paste0(
      "Influenza (IN) : Treatment effects on response vs marker across trials (Day ",
      tp,
      ")"
    ),
    subtitle = paste0("Signature derived from study: ", main_study,
      
      " (N = ",
      n_study_main, ")"),
    size = "Trial N",
    fill = NULL
  ) +
  scale_size_continuous(range = c(1.5, 7), breaks = legend_ns, labels = as.character(legend_ns)) +
  scale_fill_manual(
    values = setNames(c("black", "#F60012"), c("Evaluation studies", main_label)),
    guide = guide_legend(
      override.aes = list(
        size = 5,
        shape = 21,
        colour = c("black", "#F60012"),
        fill = c("black", "#F60012")
      ),
      order = 1
    )
  ) +
  # color scale matches fill but we remove its legend to avoid duplication
  scale_colour_manual(
    values = setNames(c("black", "#F60012"), c("Evaluation studies", main_label)),
    guide = "none"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 20, hjust = 0.5),
    axis.title = element_text(size = 20)
  )

print(p)

# Compute an evaluation metric (weighted absolute delta across evaluation studies)

# --- 1) compute weighted metric per training study and time (wMAE) -----------
weighted_delta_df <- 
  names(rise_studysignatures_results_list) %>%        # iterate over times (e.g. "Day 7", "Day 1", ...)
  map_dfr(function(day_name) {
    day_list <- rise_studysignatures_results_list[[day_name]]
    
    # iterate over training studies for that day
    names(day_list) %>% 
      map_dfr(function(train_name) {
        df <- day_list[[train_name]]
        
        # keep evaluation studies only
        df_eval <- df %>%
          filter(study_accession != train_name) %>%
          filter(!is.na(delta), !is.na(n))
        
        total_w <- sum(df_eval$n, na.rm = TRUE)
        
        wmae_val <- if (total_w > 0) {
          sum(abs(df_eval$delta) * df_eval$n, na.rm = TRUE) / total_w
        } else {
          NA_real_
        }
        
        tibble(
          training_study = train_name,
          time = day_name,
          wMAE = wmae_val
        )
      })
  })

# --- 2) Build grid for a chosen day, adding wMAE into subtitles -------------
day_name <- "Day 1" 
day_number = 1 # change as needed
day_list <- rise_studysignatures_results_list[[day_name]]

# Global size-legend breaks across all studies for a consistent legend
all_n <- unlist(lapply(day_list, function(df) df$n))
n_min <- min(all_n, na.rm = TRUE)
n_max <- max(all_n, na.rm = TRUE)
legend_min <- ceiling(n_min / 10) * 10
legend_max <- floor(n_max / 10) * 10
legend_ns <- seq(legend_min, legend_max, length.out = 4)
legend_ns[c(2,3)] <- round(legend_ns[c(2,3)] / 50) * 50
legend_ns <- unique(legend_ns)

# Labels / colours
main_label_text <- "Training study"
eval_label_text <- "Evaluation studies"
fill_values <- setNames(c("black", "#F60012"), c(eval_label_text, main_label_text))
colour_values <- fill_values

# helper to fetch wMAE text for this training study/time (rounded)
get_wmae_text <- function(train_name, day_name, digits = 3) {
  val <- weighted_delta_df %>%
    filter(training_study == train_name, time == day_name) %>%
    pull(wMAE)
  if (length(val) == 0 || is.na(val)) {
    return("wMAE = NA")
  } else {
    return(paste0("wMAE = ", formatC(round(val, digits), format = "f", digits = digits)))
  }
}

# function to create a single plot given the name of the training study
make_plot_for_training <- function(main_study_name) {
  combined_df <- day_list[[main_study_name]] %>%
    mutate(
      is_main_flag = ifelse(study_accession == main_study_name, main_label_text, eval_label_text),
      is_main_flag = factor(is_main_flag, levels = c(eval_label_text, main_label_text))
    )
  
  # number (N) for the training study (may be NA if not present)
  n_study_main <- combined_df %>%
    filter(study_accession == main_study_name) %>%
    pull(n) %>%
    unique()
  if (length(n_study_main) == 0) n_study_main <- NA_integer_
  
  wmae_text <- get_wmae_text(main_study_name, day_name, digits = 3)
  
  subtitle_text <- if (!is.na(n_study_main)) {
    paste0("N = ", n_study_main, " — ", wmae_text)
  } else {
    wmae_text
  }
  
  p <- ggplot(combined_df, aes(x = u_y, y = u_s)) +
    geom_point(
      aes(size = n, fill = is_main_flag, colour = is_main_flag),
      shape = 21,
      alpha = 0.5,
      stroke = 0.6
    ) +
    geom_abline(slope = 1, intercept = 0, color = "#FF2128", linetype = "dashed", size = 0.8, alpha = 0.5) +
    scale_x_continuous(limits = c(0, 1.01), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.01), expand = c(0, 0)) +
    coord_fixed(ratio = 1) +
    labs(
      title = main_study_name,                        
      subtitle = subtitle_text,
      x = expression(U[Y]),
      y = expression(U[S]),
      size = "Trial N",
      fill = NULL
    ) +
    scale_size_continuous(range = c(1.5, 7), breaks = legend_ns, labels = as.character(legend_ns)) +
    scale_fill_manual(
      values = fill_values,
      guide = guide_legend(
        override.aes = list(size = 4, shape = 21, colour = unname(colour_values), fill = unname(fill_values)),
        order = 1
      )
    ) +
    scale_colour_manual(values = colour_values, guide = "none") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 10),
      legend.position = "right",
      # small margin around each panel to increase separation
      plot.margin = unit(c(4, 4, 4, 4), "mm")
    )
  
  p
}

# Create plots and assemble grid
plot_list <- lapply(names(day_list), make_plot_for_training)
n_plots <- length(plot_list)
ncol_grid <- ceiling(sqrt(n_plots))

# Increase overall spacing by adding plot margin at the assembly stage as well
grid_plot <- wrap_plots(plot_list, ncol = ncol_grid) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("Influenza (IN) — Treatment effects on response vs marker across trials (", day_name, ")"),
    subtitle = paste0("Each panel: signature trained in the named study; highlighted point = training study"),
    theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 18, hjust = 0.5))
  ) & theme(
    legend.position = "right",
    # larger margin around the assembled grid to visually separate panels
    plot.margin = unit(c(5, 5, 5, 5), "mm")
  )

# show plot
print(grid_plot)

# Save the assembled grid plot as a PDF
p_figure = fs::path("output", "figures", "surrogacy", paste0("influenzain_rise_study_grid_day", day_number, ".pdf"))

pdf(p_figure, width = 15, height = 15)
print(grid_plot)
dev.off()
