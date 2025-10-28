# Here, we select certain results from the prediction task and augment the figures for better representation

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggtext)
library(forcats)
library(scales)
library(patchwork)
library(grid)
library(gridExtra)

# Directory to store engineered data
processed_data_folder = "data"

# Directory to store figures
prediction_figures_folder = fs::path("output", "figures", "prediction")

# Assign each prediction set a colour code for later visualisation
prediction_set_colours = list(
  "clinical" = "#001219",
  "Day 0" = "#0a9396",
  "Day 1" = "#e9d8a6",
  "Day 3" = "#ee9b00",
  "Day 7" = "#bb3e03",
  "Day 10" = "#ae2012",
  "Day 14" = "#780000"
)

# Load the merged dataframe to extract the sample sizes for selected vaccines at selected timepoints
vaccine_list = c(
  "Meningococcus (PS)",
  "Meningococcus (CJ)" ,
  "Influenza (IN)"  ,
  "Hepatitis A/B (IN/RP)" ,
  "Yellow Fever (LV)"
)

time_list = c(0, 1, 3, 7, 14)

# Directory to store figures
descriptive_figures_folder = fs::path("output", "figures", "descriptive")

# Path to processed gene-level data
p_load_expr_all_norm <- fs::path(processed_data_folder, "hipc_merged_all_norm.rds")

# Load merged gene-level data
hipc_merged_all_norm = readRDS(p_load_expr_all_norm)

# Filter out observations which do not have an immune response value
hipc_merged_all_norm_response = hipc_merged_all_norm %>%
  filter(
    !is.na(immResp_MFC_anyAssay_log2_MFC),
    vaccine_name %in% vaccine_list,
    study_time_collected %in% time_list
  ) %>%
  filter(!(vaccine_name == "Yellow Fever (LV)" &
             study_time_collected == 1))

# Counts per vaccine x time
# Counts per vaccine x time
counts <- hipc_merged_all_norm_response %>%
  filter(!is.na(study_time_collected)) %>%
  group_by(vaccine_name, vaccine_colour, study_time_collected) %>%
  summarise(n = n(), .groups = "drop")

# Order the time points and make study_time_collected an ordered factor
time_levels <- counts %>%
  distinct(study_time_collected) %>%
  arrange(as.numeric(study_time_collected)) %>%
  pull(study_time_collected) %>%
  as.character()   # factor levels must be character

counts <- counts %>%
  mutate(study_time_collected = factor(
    as.character(study_time_collected),
    levels = time_levels,
    ordered = TRUE
  ))

# Make bubble size proportional to count/sqrt(count)
counts <- counts %>%
  mutate(size_var = n / sqrt(n))

# Choose a few representative count values for the legend
n_breaks <- c(10, 50, 100, 250, 800)       # auto selects ~4 visually pleasing numbers
size_breaks <- n_breaks / sqrt(n_breaks)  # convert to size-space

# Bubble plot
p1 <- ggplot(counts, aes(x = study_time_collected, y = vaccine_name)) +
  geom_point(
    aes(size = size_var, fill = vaccine_colour),
    shape = 21,
    colour = "black",
    alpha = 0.95
  ) +
  geom_text(
    aes(label = n),
    colour = "white",
    size = 3.5,
    vjust = 0.5
  ) +
  scale_fill_identity(guide = "none") +
  scale_size_area(
    max_size = 28,
    breaks = size_breaks,
    # positions in size-space
    labels = n_breaks,
    # what to display in the legend
    name = "Number of samples"
  ) +
  guides(size = guide_legend(override.aes = list(
    fill = "grey70",
    colour = "black",
    shape = 21
  ))) +
  labs(
    x = "Days post-vaccination",
    y = "Vaccine",
    title = "Number of participants with transcriptomic & immune response samples",
    subtitle = "Selected timepoints and vaccines"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(hjust = 1),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 15, hjust = 0.5),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )

print(p1)

ggsave(
  filename = "custom_bubble_plot_sequential.pdf",
  path = descriptive_figures_folder,
  plot = p1,
  width = 30,
  height = 17,
  units = "cm"
)

# Load the sequential, with clinical, without TBA results
# We are going to adjust the evaluation plot
p_prediction_results_all_sequential_withClinical_withoutTBA = fs::path(
  "output",
  "results",
  "prediction_results_all_sequential_withClinical_withoutTBA.rds"
)

prediction_results_all_sequential_withClinical_withoutTBA = readRDS(p_prediction_results_all_sequential_withClinical_withoutTBA)

# --- 1. Build vaccine x set grid and extract Rspearman + sRMSE in one pass ----
vaccines <- names(prediction_results_all_sequential_withClinical_withoutTBA)

# gather all predictor-set names that appear anywhere
all_sets <- c("clinical", "Day 0", "Day 1", "Day 3", "Day 7", "Day 14")

# expand grid of all combinations
grid <- expand.grid(vaccine = vaccines,
                    set     = all_sets,
                    stringsAsFactors = FALSE)

# pre-allocate numeric vectors for efficiency
n <- nrow(grid)
Rs_vec    <- rep(NA_real_, n)
sRMSE_vec <- rep(NA_real_, n)

# single-pass extraction (avoids rowwise / repeated list traversal overhead)
for (i in seq_len(n)) {
  vac  <- grid$vaccine[i]
  setn <- grid$set[i]
  # guard against missing vaccine or set
  if (is.null(prediction_results_all_sequential_withClinical_withoutTBA[[vac]]))
    next
  res_set <- prediction_results_all_sequential_withClinical_withoutTBA[[vac]][[setn]]
  if (is.null(res_set) || !is.list(res_set))
    next
  if (!("metrics" %in% names(res_set)))
    next
  metrics <- res_set[["metrics"]]
  # extract if present (coerce to numeric, keep NA if missing)
  if ("Rspearman" %in% names(metrics)) {
    Rs_vec[i] <- as.numeric(metrics[["Rspearman"]])
  }
  if ("sRMSE" %in% names(metrics)) {
    sRMSE_vec[i] <- as.numeric(metrics[["sRMSE"]])
  }
}

# combine into data frame
plot_df <- bind_cols(as_tibble(grid), tibble(Rs = Rs_vec, sRMSE = sRMSE_vec))

# formatted labels: empty string for NA (so no annotation on missing tiles)
plot_df <- plot_df %>%
  mutate(
    Rs_label    = ifelse(is.na(Rs), "", sprintf("%.2f", Rs)),
    sRMSE_label = ifelse(is.na(sRMSE), "", sprintf("%.2f", sRMSE))
  )

# factors for plotting order (adjust ordering if desired)
plot_df <- plot_df %>%
  mutate(
    set     = factor(set, levels = all_sets),
    vaccine = factor(vaccine, levels = vaccines)
  )

# --- 2. Spearman R heatmap (white = 0 -> blue = 1; NA -> grey) ----
heatmap_plot_R <- ggplot(plot_df, aes(x = set, y = vaccine, fill = Rs)) +
  geom_tile(colour = "white") +
  # annotate only non-NA values (labels are empty for NA)
  geom_text(aes(label = Rs_label),
            colour = "black",
            size = 5) +
  scale_fill_gradient(
    name = expression(rho),
    low = "white",
    high = "#0072B2",
    # blue
    na.value = "grey80",
    limits = c(0, 1),
    # anchor scale to 0..1
    oob = scales::squish,
    # values outside 0..1 will be squished to limits
    breaks = seq(0, 1, by = 0.25)
  ) +
  coord_fixed(ratio = 1) +   # square tiles
  labs(x = "Predictor set",
       y = "Vaccine",
       title = expression(paste(
         "Spearman ", rho, " by vaccine and predictor set"
       ))) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 25, hjust = 0.5),
  )
# Print Spearman R heatmap
print(heatmap_plot_R)

# --- 3. sRMSE heatmap (0 = green -> 1 = white; values > 1 shown as white; NA -> grey) ----
heatmap_plot_sRMSE <- ggplot(plot_df, aes(x = set, y = vaccine, fill = sRMSE)) +
  geom_tile(colour = "white") +
  # annotate only non-NA values
  geom_text(aes(label = sRMSE_label),
            colour = "black",
            size = 5) +
  scale_fill_gradient(
    name = "sRMSE",
    low = "#2ca02c",
    # green at 0
    high = "white",
    # white at 1 and above (squished)
    na.value = "grey80",
    limits = c(0, 1),
    oob = scales::squish,
    breaks = seq(0, 1, by = 0.25)
  ) +
  coord_fixed(ratio = 1) +   # square tiles
  labs(x = "Predictor set", y = "Vaccine", title = "Standardised RMSE by vaccine and predictor set") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 25, hjust = 0.5),
  )

# Print sRMSE heatmap
print(heatmap_plot_sRMSE)


# --- Modify heatmap plots ---
heatmap_plot_R_mod <- heatmap_plot_R +
  ggtitle(expression(paste("Spearman ", rho))) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 22, hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.margin = margin(
      t = 10,
      r = 2,
      b = 5,
      l = 10
    ),
    legend.position = "none"   # remove individual legend
  )

heatmap_plot_sRMSE_mod <- heatmap_plot_sRMSE +
  ggtitle("Standardised RMSE") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 22, hjust = 0.5),
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.margin = margin(
      t = 5,
      r = 2,
      b = 10,
      l = 10
    ),
    legend.position = "none"   # remove individual legend
  )

# --- Combine with shared legend and common title ---
combined <- heatmap_plot_R_mod /
  plot_spacer() /
  heatmap_plot_sRMSE_mod +
  plot_layout(
    ncol = 1,
    heights = c(1, 0.1, 1),
    # adjust 0.05 to control spacing
    guides = "collect"
  ) +
  plot_annotation(
    title = "Evaluation metrics of CV predictions across vaccines and timepoints",
    theme = theme(
      plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  ) &
  theme(
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.key.width = unit(1, "cm")
  )

# --- Display combined plot ---
print(combined)

# --- Save figure ---
ggsave(
  filename = "custom_evaluation_combined_sequential_withClinical_withoutTBA.pdf",
  path = prediction_figures_folder,
  plot = combined,
  width = 40,
  height = 25,
  units = "cm"
)

# Next we represent all CV predictions in a grid (vaccine x feature set)

# --- user inputs --------------------------------------------------------------
lst <- prediction_results_all_sequential_withClinical_withoutTBA
feature_set_set <- c("clinical", paste0("Day ", c(0, 3, 7)))
vaccines <- names(lst)   # expecting 5 vaccines as you mentioned

# --- blank placeholder (pure white, no border) --------------------------------
blank_placeholder <- function() {
  # an empty ggplot with no axes or content, pure white
  ggplot() +
    geom_blank() +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(6, 6, 6, 6)
    )
}

# --- build grob matrix --------------------------------------------------------
n_row <- length(vaccines)
n_col <- length(feature_set_set)

total_rows <- n_row + 1  # one extra top row for feature-set headers
total_cols <- n_col + 1  # one extra left column for vaccine labels

grob_list <- vector("list", total_rows * total_cols)

# top-left empty cell (no border)
grob_list[[1]] <- nullGrob()

# top row: feature set headers (larger font)
for (j in seq_len(n_col)) {
  grob_list[[1 + j]] <- textGrob(label = feature_set_set[j],
                                 gp = gpar(fontsize = 18),
                                 just = "center")
}

# left column: vaccine labels (larger) and body plots
for (i in seq_len(n_row)) {
  # left label position (row i+1, col 1)
  left_index <- (i + 1 - 1) * total_cols + 1
  grob_list[[left_index]] <- textGrob(label = vaccines[i],
                                      gp = gpar(fontsize = 16),
                                      just = "left")
  
  for (j in seq_len(n_col)) {
    idx <- (i + 1 - 1) * total_cols + (j + 1)
    feat <- feature_set_set[j]
    vac  <- vaccines[i]
    
    got_plot <- NULL
    if (!is.null(lst[[vac]]) && !is.null(lst[[vac]][[feat]])) {
      inner <- lst[[vac]][[feat]]
      if (!is.null(inner[["plots"]]) &&
          !is.null(inner[["plots"]][["cv_results"]])) {
        got_plot <- inner[["plots"]][["cv_results"]]
      }
    }
    
    if (!is.null(got_plot) && inherits(got_plot, "ggplot")) {
      # 1) remove title/subtitle
      p_mod <- got_plot + labs(title = NULL, subtitle = NULL)
      
      # 2) remove any geom_text / geom_label layers (annotation metrics)
      if (length(p_mod$layers) > 0) {
        p_mod$layers <- Filter(function(l) {
          !(inherits(l$geom, "GeomText") || inherits(l$geom, "GeomLabel"))
        }, p_mod$layers)
      }
      
      if (length(p_mod$layers) > 0) {
        for (k in seq_along(p_mod$layers)) {
          # GeomPoint layers
          if (inherits(p_mod$layers[[k]]$geom, "GeomPoint")) {
            # Override alpha; either fixed or mapped alpha
            p_mod$layers[[k]]$aes_params$alpha <- 0.45    # set fixed alpha
            p_mod$layers[[k]]$geom_params$alpha <- 0.45   # also set geom default alpha
          }
        }
      }
      
      # 3) tidy theme: remove axis titles, reduce axis text size, add panel border to the plot
      p_mod <- p_mod +
        theme(
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.margin = margin(6, 6, 6, 6),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          # this creates a border *inside* the plot (around the panel)
          panel.border = element_rect(
            colour = "grey70",
            fill = NA,
            size = 0.8
          )
        )
      
      # convert to grob (plot already contains the border)
      grob_list[[idx]] <- ggplotGrob(p_mod)
    } else {
      # missing -> blank white placeholder (no border)
      grob_list[[idx]] <- ggplotGrob(blank_placeholder())
    }
  }
}

# --- layout matrix ------------------------------------------------------------
layout_mat <- matrix(
  seq_len(total_rows * total_cols),
  nrow = total_rows,
  ncol = total_cols,
  byrow = TRUE
)

# relative widths/heights (tweak to taste)
col_widths <- c(4.5, rep(5, n_col))
row_heights <- c(1.5, rep(5, n_row))

# Draw the arranged grid with title
# Create a single grob object without drawing it immediately
final_grob <- arrangeGrob(
  grobs = grob_list,
  layout_matrix = layout_mat,
  widths = col_widths,
  heights = row_heights,
  top = textGrob(
    "CV predictions across feature sets and vaccines",
    gp = gpar(fontsize = 30, fontface = "bold")
  )
)
# You can also save as PDF
ggsave(
  filename = "custom_predictions_grid.pdf",
  path = prediction_figures_folder,
  plot = final_grob,
  width = 17,
  height = 12
)

# Now, just for influenza

# --- user inputs --------------------------------------------------------------
lst <- prediction_results_all_sequential_withClinical_withoutTBA
vac <- "Influenza (IN)"

feature_sets <- c("clinical", paste0("Day ", c(0, 1, 3, 7, 14)))

# --- annotation tuning (change these if you want smaller/larger annotations) ---
ann_size  <- 4   # geom_text/geom_label size (ggplot2 units) — smaller number => smaller text
ann_alpha <- 0.9   # transparency for annotation text

# --- blank placeholder (pure white, no border) --------------------------------
blank_placeholder <- function() {
  ggplot() + geom_blank() + theme_void() +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(6, 6, 6, 6)
    )
}

# --- build grobs: top row headers, bottom row plots ----------------------------
n_col <- length(feature_sets)
if (n_col == 0)
  stop("No feature sets found for vaccine: ", vac)

# We'll create 2 rows (headers + plots) and n_col columns
total_rows <- 2
total_cols <- n_col

grob_list <- vector("list", total_rows * total_cols)

# top row: feature set headers
for (j in seq_len(n_col)) {
  h <- textGrob(
    label = feature_sets[j],
    gp = gpar(fontsize = 16, fontface = "bold"),
    just = "center"
  )
  grob_list[[j]] <- h
}

# bottom row: the plots (with panel border inside each plot); blanks if missing
for (j in seq_len(n_col)) {
  feat <- feature_sets[j]
  got_plot <- NULL
  
  if (!is.null(lst[[vac]]) && !is.null(lst[[vac]][[feat]])) {
    inner <- lst[[vac]][[feat]]
    if (!is.null(inner[["plots"]]) &&
        !is.null(inner[["plots"]][["cv_results"]])) {
      got_plot <- inner[["plots"]][["cv_results"]]
    }
  }
  
  if (!is.null(got_plot) && inherits(got_plot, "ggplot")) {
    # Keep metric annotations but reduce their size/alpha
    p_mod <- got_plot + labs(title = NULL, subtitle = NULL)
    
    # shrink any geom_text / geom_label layers
    if (length(p_mod$layers) > 0) {
      for (k in seq_along(p_mod$layers)) {
        if (inherits(p_mod$layers[[k]]$geom, "GeomText") ||
            inherits(p_mod$layers[[k]]$geom, "GeomLabel")) {
          # set fixed text size & alpha
          p_mod$layers[[k]]$aes_params$size  <- ann_size
          p_mod$layers[[k]]$geom_params$size <- ann_size
          p_mod$layers[[k]]$aes_params$alpha  <- ann_alpha
          p_mod$layers[[k]]$geom_params$alpha <- ann_alpha
          
          # for labels, reduce padding/rounding so boxes are smaller
          if (inherits(p_mod$layers[[k]]$geom, "GeomLabel")) {
            p_mod$layers[[k]]$geom_params$label.padding <- grid::unit(0.05, "lines")
            p_mod$layers[[k]]$geom_params$label.r <- grid::unit(0.05, "lines")
          }
        }
      }
    }
    
    # apply compact theme and panel border
    p_mod <- p_mod +
      theme(
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.margin = margin(6, 6, 6, 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        panel.border = element_rect(
          colour = "grey70",
          fill = NA,
          size = 0.8
        ),
        panel.background = element_rect(fill = "white", colour = NA)
      )
    
    grob_list[[n_col + j]] <- ggplotGrob(p_mod) # second row index = n_col + j
  } else {
    # missing -> blank white placeholder (no border)
    grob_list[[n_col + j]] <- ggplotGrob(blank_placeholder())
  }
}

# --- layout matrix ------------------------------------------------------------
layout_mat <- matrix(
  seq_len(total_rows * total_cols),
  nrow = total_rows,
  ncol = total_cols,
  byrow = TRUE
)

# sizing: tweak these if the result looks cramped
col_widths <- rep(4, n_col)
row_heights <- c(1, 5)  # 1st row small for headers, second row big for plots

# Create the arranged grob
final_grob <- arrangeGrob(
  grobs = grob_list,
  layout_matrix = layout_mat,
  widths = col_widths,
  heights = row_heights,
  top = textGrob(paste0("CV results — ", vac), gp = gpar(
    fontsize = 18, fontface = "bold"
  ))
)

# Draw it to the screen (optional)
grid.newpage()
grid.draw(final_grob)

# Save to file
ggsave(
  filename = "custom_predictions_grid_influenza.pdf",
  plot = final_grob,
  # <-- pass the grob here
  path = prediction_figures_folder,
  width = 15,
  height = 4
)

# Now, just for influenza

# --- user inputs --------------------------------------------------------------
lst <- prediction_results_all_sequential_withClinical_withoutTBA
vac <- "Yellow Fever (LV)"

feature_sets <- c("clinical", paste0("Day ", c(0, 3, 7, 14)))

# --- annotation tuning (change these if you want smaller/larger annotations) ---
ann_size  <- 4   # geom_text/geom_label size (ggplot2 units) — smaller number => smaller text
ann_alpha <- 0.9   # transparency for annotation text

# --- blank placeholder (pure white, no border) --------------------------------
blank_placeholder <- function() {
  ggplot() + geom_blank() + theme_void() +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(6, 6, 6, 6)
    )
}

# --- build grobs: top row headers, bottom row plots ----------------------------
n_col <- length(feature_sets)
if (n_col == 0)
  stop("No feature sets found for vaccine: ", vac)

# We'll create 2 rows (headers + plots) and n_col columns
total_rows <- 2
total_cols <- n_col

grob_list <- vector("list", total_rows * total_cols)

# top row: feature set headers
for (j in seq_len(n_col)) {
  h <- textGrob(
    label = feature_sets[j],
    gp = gpar(fontsize = 16, fontface = "bold"),
    just = "center"
  )
  grob_list[[j]] <- h
}

# bottom row: the plots (with panel border inside each plot); blanks if missing
for (j in seq_len(n_col)) {
  feat <- feature_sets[j]
  got_plot <- NULL
  
  if (!is.null(lst[[vac]]) && !is.null(lst[[vac]][[feat]])) {
    inner <- lst[[vac]][[feat]]
    if (!is.null(inner[["plots"]]) &&
        !is.null(inner[["plots"]][["cv_results"]])) {
      got_plot <- inner[["plots"]][["cv_results"]]
    }
  }
  
  if (!is.null(got_plot) && inherits(got_plot, "ggplot")) {
    # Keep metric annotations but reduce their size/alpha
    p_mod <- got_plot + labs(title = NULL, subtitle = NULL)
    
    # shrink any geom_text / geom_label layers
    if (length(p_mod$layers) > 0) {
      for (k in seq_along(p_mod$layers)) {
        if (inherits(p_mod$layers[[k]]$geom, "GeomText") ||
            inherits(p_mod$layers[[k]]$geom, "GeomLabel")) {
          # set fixed text size & alpha
          p_mod$layers[[k]]$aes_params$size  <- ann_size
          p_mod$layers[[k]]$geom_params$size <- ann_size
          p_mod$layers[[k]]$aes_params$alpha  <- ann_alpha
          p_mod$layers[[k]]$geom_params$alpha <- ann_alpha
          
          # for labels, reduce padding/rounding so boxes are smaller
          if (inherits(p_mod$layers[[k]]$geom, "GeomLabel")) {
            p_mod$layers[[k]]$geom_params$label.padding <- grid::unit(0.05, "lines")
            p_mod$layers[[k]]$geom_params$label.r <- grid::unit(0.05, "lines")
          }
        }
      }
    }
    
    # apply compact theme and panel border
    p_mod <- p_mod +
      theme(
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.margin = margin(6, 6, 6, 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        panel.border = element_rect(
          colour = "grey70",
          fill = NA,
          size = 0.8
        ),
        panel.background = element_rect(fill = "white", colour = NA)
      )
    
    grob_list[[n_col + j]] <- ggplotGrob(p_mod) # second row index = n_col + j
  } else {
    # missing -> blank white placeholder (no border)
    grob_list[[n_col + j]] <- ggplotGrob(blank_placeholder())
  }
}

# --- layout matrix ------------------------------------------------------------
layout_mat <- matrix(
  seq_len(total_rows * total_cols),
  nrow = total_rows,
  ncol = total_cols,
  byrow = TRUE
)

# sizing: tweak these if the result looks cramped
col_widths <- rep(4, n_col)
row_heights <- c(1, 5)  # 1st row small for headers, second row big for plots

# Create the arranged grob
final_grob <- arrangeGrob(
  grobs = grob_list,
  layout_matrix = layout_mat,
  widths = col_widths,
  heights = row_heights,
  top = textGrob(paste0("CV results — ", vac), gp = gpar(
    fontsize = 18, fontface = "bold"
  ))
)

# Draw it to the screen (optional)
grid.newpage()
grid.draw(final_grob)

# Save to file
ggsave(
  filename = "custom_predictions_grid_yellowfeverlv.pdf",
  plot = final_grob,
  # <-- pass the grob here
  path = prediction_figures_folder,
  width = 15,
  height = 4
)

# For influenza, variable importance

# --- user inputs --------------------------------------------------------------
vac <- "Influenza (IN)"
feature_sets <- c(paste0("Day ", c(0, 1, 3, 7)))
topN <- 15

lst <- prediction_results_all_sequential_withClinical_withoutTBA

# --- blank placeholder (pure white, no border) --------------------------------
blank_placeholder <- function() {
  ggplot() + geom_blank() + theme_void() +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(6, 6, 6, 6)
    )
}

# --- helper function to build a single VI plot ---------------------------------
make_vi_plot <- function(vac, feature_set, topN, wrap_width = 20) {
  selected_results_vaccine <- lst[[vac]]
  selected_results_vaccine_feature_set <- selected_results_vaccine[[feature_set]]
  
  if (is.null(selected_results_vaccine_feature_set))
    return(blank_placeholder())
  
  metrics <- selected_results_vaccine_feature_set[["metrics"]]
  plot_obj <- selected_results_vaccine_feature_set[["var_imp"]]
  
  # Prepare topN data and wrap feature names
  plot_df <- plot_obj %>%
    slice_min(order_by = mean_rank, n = topN) %>%
    mutate(
      mean_imp = as.numeric(mean_imp),
      sd_imp   = as.numeric(sd_imp),
      feature  = forcats::fct_reorder(stringr::str_wrap(feature, width = wrap_width), mean_imp),
      xmin = pmax(0, mean_imp - sd_imp),
      xmax = pmin(100, mean_imp + sd_imp)
    )
  
  # colours
  colour_df <- plot_obj %>%
    dplyr::select(feature_group, feature_colour) %>%
    distinct()
  pred_cols <- colour_df$feature_colour
  names(pred_cols) <- colour_df$feature_group
  
  # Build plot
  p <- ggplot(plot_df, aes(x = feature, y = mean_imp)) +
    geom_errorbar(
      aes(
        ymin = xmin,
        ymax = xmax,
        colour = feature_group
      ),
      width = 0.25,
      size = 1,
      orientation = "x"
    ) +
    geom_point(aes(colour = feature_group), size = 2) +
    scale_colour_manual(values = pred_cols, na.value = "#777777") +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0.02, 0.12))) +
    labs(
      x = NULL,
      y = "Mean standardized VI",
      title = feature_set,
      subtitle = bquote(sRMSE == .(sprintf(
        "%.2f", metrics$sRMSE
      )) ~ ", " ~ rho == .(sprintf(
        "%.2f", metrics$Rspearman
      )))
    ) +
    coord_flip() +  # flip coordinates
    theme_minimal() +
    theme(
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      panel.border = element_rect(
        colour = "grey70",
        fill = NA,
        size = 0.8
      ),
      plot.margin = margin(4, 4, 4, 4)
    )
  
  return(p)
}


# --- build grobs for all feature sets ----------------------------------------
grob_list <- lapply(feature_sets, function(fs)
  make_vi_plot(vac, fs, topN))

# Assign the arranged grob to a variable
final_vi_grob <- arrangeGrob(
  grobs = grob_list,
  nrow = 1,
  top = textGrob(
    paste0("Variable Importance — ", vac),
    gp = gpar(fontsize = 18, fontface = "bold")
  )
)

# Save as PDF (or change to PNG/JPEG if desired)
ggsave(
  filename = "custom_variable_importance_influenza.pdf",
  plot = final_vi_grob,
  path = prediction_figures_folder,  # replace with your folder
  width = 15,   # adjust width to fit all plots nicely
  height = 6
)

# Optional: also display it
grid.newpage()
grid.draw(final_vi_grob)

# For influenza, variable importance

# --- user inputs --------------------------------------------------------------
vac <- "Yellow Fever (LV)"
feature_sets <- c(paste0("Day ", c(0, 3, 7, 14)))
topN <- 15

lst <- prediction_results_all_sequential_withClinical_withoutTBA

# --- blank placeholder (pure white, no border) --------------------------------
blank_placeholder <- function() {
  ggplot() + geom_blank() + theme_void() +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(6, 6, 6, 6)
    )
}

# --- helper function to build a single VI plot ---------------------------------
make_vi_plot <- function(vac, feature_set, topN, wrap_width = 20) {
  selected_results_vaccine <- lst[[vac]]
  selected_results_vaccine_feature_set <- selected_results_vaccine[[feature_set]]
  
  if (is.null(selected_results_vaccine_feature_set))
    return(blank_placeholder())
  
  metrics <- selected_results_vaccine_feature_set[["metrics"]]
  plot_obj <- selected_results_vaccine_feature_set[["var_imp"]]
  
  # Prepare topN data and wrap feature names
  plot_df <- plot_obj %>%
    slice_min(order_by = mean_rank, n = topN) %>%
    mutate(
      mean_imp = as.numeric(mean_imp),
      sd_imp   = as.numeric(sd_imp),
      feature  = forcats::fct_reorder(stringr::str_wrap(feature, width = wrap_width), mean_imp),
      xmin = pmax(0, mean_imp - sd_imp),
      xmax = pmin(100, mean_imp + sd_imp)
    )
  
  # colours
  colour_df <- plot_obj %>%
    dplyr::select(feature_group, feature_colour) %>%
    distinct()
  pred_cols <- colour_df$feature_colour
  names(pred_cols) <- colour_df$feature_group
  
  # Build plot
  p <- ggplot(plot_df, aes(x = feature, y = mean_imp)) +
    geom_errorbar(
      aes(
        ymin = xmin,
        ymax = xmax,
        colour = feature_group
      ),
      width = 0.25,
      size = 1,
      orientation = "x"
    ) +
    geom_point(aes(colour = feature_group), size = 2) +
    scale_colour_manual(values = pred_cols, na.value = "#777777") +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0.02, 0.12))) +
    labs(
      x = NULL,
      y = "Mean standardized VI",
      title = feature_set,
      subtitle = bquote(sRMSE == .(sprintf(
        "%.2f", metrics$sRMSE
      )) ~ ", " ~ rho == .(sprintf(
        "%.2f", metrics$Rspearman
      )))
    ) +
    coord_flip() +  # flip coordinates
    theme_minimal() +
    theme(
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      panel.border = element_rect(
        colour = "grey70",
        fill = NA,
        size = 0.8
      ),
      plot.margin = margin(4, 4, 4, 4)
    )
  
  return(p)
}


# --- build grobs for all feature sets ----------------------------------------
grob_list <- lapply(feature_sets, function(fs)
  make_vi_plot(vac, fs, topN))

# --- arrange horizontally ----------------------------------------------------
# Assign the arranged grob to a variable
final_vi_grob <- arrangeGrob(
  grobs = grob_list,
  nrow = 1,
  top = textGrob(
    paste0("Variable Importance — ", vac),
    gp = gpar(fontsize = 18, fontface = "bold")
  )
)

# Save as PDF (or change to PNG/JPEG if desired)
ggsave(
  filename = "custom_variable_importance_yellowfeverlv.pdf",
  plot = final_vi_grob,
  path = prediction_figures_folder,  # replace with your folder
  width = 15,   # adjust width to fit all plots nicely
  height = 7
)

# Optional: also display it
grid.newpage()
grid.draw(final_vi_grob)
