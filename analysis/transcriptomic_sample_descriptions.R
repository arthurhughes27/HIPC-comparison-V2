# R script to perform some descriptions of the available transcriptomic samples for those with valid immune response data
## In particular, we are interested in exploring how many samples are available at each timepoint and vaccine for those with immune response measurements
## We use this to justify a prediction strategy

# Packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

# Directory to store engineered data
processed_data_folder = "data"

# Directory to store figures
descriptive_figures_folder = fs::path("output", "figures", "descriptive")

# Path to processed gene-level data
p_load_expr_all_norm <- fs::path(processed_data_folder, "hipc_merged_all_norm.rds")

# Load merged gene-level data
hipc_merged_all_norm = readRDS(p_load_expr_all_norm)

# Filter out observations which do not have an immune response value
hipc_merged_all_norm_response = hipc_merged_all_norm %>%
  filter(!is.na(immResp_MFC_anyAssay_log2_MFC))

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

# Order the counts by the study time
counts <- counts %>%
  mutate(study_time_collected = factor(
    as.character(study_time_collected),
    levels = time_levels,
    ordered = TRUE
  ))

# Since the range of the number of samples is large, we make the bubble size proportional to the count/square root of count
counts <- counts %>%
  mutate(size_var = n / sqrt(n))

# Bubble plot
p1 <- ggplot(counts, aes(x = study_time_collected, y = vaccine_name)) +
  # points coloured by vaccine hex codes; shape 21 allows fill + black border
  geom_point(
    aes(size = size_var, fill = vaccine_colour),
    shape = 21,
    colour = "black",
    alpha = 0.95,
    show.legend = FALSE
  ) +
  # numeric labels inside bubbles; label shows raw counts n and is white
  geom_text(
    aes(label = n),
    colour = "white",
    size = 3.5,
    vjust = 0.5,
    show.legend = FALSE
  ) +
  scale_fill_identity(guide = "none") +     # use hex codes directly, no fill legend
  scale_size_area(max_size = 28, guide = "none") +
  labs(
    x = "Days post-vaccination",
    y = "Vaccine",
    title = "Number of participants with transcriptomic samples at each timepoint",
    subtitle = "All timepoints, only participants with valid immune response data included"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 15, hjust = 0.5)
  )

print(p1)

ggsave(
  filename = "bubble_plot_sequential.pdf",
  path = descriptive_figures_folder,
  plot = p1,
  width = 40,
  height = 20,
  units = "cm"
)

# Cumulative sample size bubble plot
## Define some timepoints. Samples only count if they have samples in all the previous timepoints.

# timepoints of interest (ordered)
timepoints_of_interest <- c(0, 1, 3, 7)

# Get vaccine -> colour map (in case some vaccines don't appear in filtered data)
vaccine_colours <- hipc_merged_all_norm_response %>%
  distinct(vaccine_name, vaccine_colour) %>%
  filter(!is.na(vaccine_name))

# Keep only rows for the timepoints of interest and deduplicate by participant/time
df_small <- hipc_merged_all_norm_response %>%
  filter(
    !is.na(participant_id),!is.na(study_time_collected),
    study_time_collected %in% timepoints_of_interest
  ) %>%
  distinct(vaccine_name,
           vaccine_colour,
           participant_id,
           study_time_collected)

# For each vaccine + participant collect the set of times they have
participant_times <- df_small %>%
  group_by(vaccine_name, vaccine_colour, participant_id) %>%
  summarise(times = list(sort(unique(
    as.numeric(study_time_collected)
  ))), .groups = "drop")

# Expand for each timepoint of interest, check whether the participant has ALL times up to that time
# (i.e., cumulative requirement).
time_grid <- tibble(study_time_collected = timepoints_of_interest)

participant_times_expanded <- participant_times %>%
  tidyr::crossing(time_grid) %>%
  mutate(study_time_collected = as.numeric(study_time_collected)) %>%
  rowwise() %>%
  mutate(has_all_up_to_t = all(timepoints_of_interest[timepoints_of_interest <= study_time_collected] %in% times)) %>%
  ungroup()

# Summarise: number of participants per vaccine x time that satisfy the cumulative condition
counts_cumulative <- participant_times_expanded %>%
  group_by(vaccine_name, vaccine_colour, study_time_collected) %>%
  summarise(n = sum(has_all_up_to_t, na.rm = TRUE),
            .groups = "drop")

# Ensure we have rows for all vaccine × timepoint_of_interest combinations (so zeros appear)
all_vaccines <- vaccine_colours %>% filter(!is.na(vaccine_name))
all_grid <- tidyr::crossing(all_vaccines, study_time_collected = timepoints_of_interest)

counts_cumulative <- all_grid %>%
  left_join(counts_cumulative,
            by = c("vaccine_name", "vaccine_colour", "study_time_collected")) %>%
  mutate(n = replace_na(n, 0))

# Make study_time_collected an ordered factor so spacing on x is equal
counts_cumulative <- counts_cumulative %>%
  mutate(study_time_collected = factor(
    as.character(study_time_collected),
    levels = as.character(timepoints_of_interest),
    ordered = TRUE
  ))

# Size variable proportional to n / sqrt(n) but handle n == 0 safely
counts_cumulative <- counts_cumulative %>%
  mutate(size_var = ifelse(n > 0, n / sqrt(n), 0))

# Plot — same visual styling as before: white labels, no legend, centered title, larger text
p2 <- ggplot(counts_cumulative,
             aes(x = study_time_collected, y = vaccine_name)) +
  geom_point(
    aes(size = size_var, fill = vaccine_colour),
    shape = 21,
    colour = "black",
    alpha = 0.95,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = n),
    colour = "white",
    size = 3.5,
    vjust = 0.5,
    show.legend = FALSE
  ) +
  scale_fill_identity(guide = "none") +
  scale_size_area(max_size = 28, guide = "none") +
  labs(
    x = "Days post-vaccination",
    y = "Vaccine",
    title = "Cumulative number of participants with transcriptomic samples at all prior timepoints",
    subtitle = "Selected timepoints, only participants with valid immune response data included"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 15, hjust = 0.5)
  )

print(p2)

ggsave(
  filename = "bubble_plot_cumulative.pdf",
  path = descriptive_figures_folder,
  plot = p2,
  width = 38,
  height = 20,
  units = "cm"
)

# Select specific timepoints for each vaccine to retain
## Plot : sequential bubble plot with non-specified timepoints greyed out

# Counts per vaccine x time
counts <- hipc_merged_all_norm_response %>%
  filter(!is.na(study_time_collected)) %>%
  group_by(vaccine_name, vaccine_colour, study_time_collected) %>%
  summarise(n = n(), .groups = "drop")

# Define the ordering
vaccine_levels <- unique(as.character(counts$vaccine_name))

# Specify the timepoints
specified_timepoints_list <- list(
  "Pneumococcus (PS)" = c(0, 1, 3, 7, 10, 14, 21, 28),
  "Meningococcus (PS)" = c(0, 3, 7),
  "Meningococcus (CJ)" = c(0, 3, 7),
  "Influenza (LV)" = c(0, 3, 7),
  "Influenza (IN)" = c(0, 1, 3, 7, 14, 28),
  "Hepatitis A/B (IN/RP)" = c(0, 7),
  "Yellow Fever (LV)" = c(0, 3, 7, 10, 14, 28),
  "Varicella Zoster (LV)" = c(0, 1, 3, 7),
  "Tuberculosis (RVV)" = c(0, 2, 7)
)

# Convert the named-list spec to a data.frame of allowed pairs
specified_pairs_df <- bind_rows(lapply(names(specified_timepoints_list), function(vac) {
  tibble(
    vaccine_name = vac,
    study_time_collected = as.character(specified_timepoints_list[[vac]])
  )
}))

# Prepare the counts data
counts_plot <- counts %>%
  mutate(study_time_collected = as.character(study_time_collected)) %>%
  # join to mark which pairs are specified
  left_join(
    specified_pairs_df %>% mutate(specified = TRUE),
    by = c("vaccine_name", "study_time_collected")
  ) %>%
  mutate(specified = replace_na(specified, FALSE))

# create visual-mapping columns (colors, alphas, label colors, border colors)
counts_plot <- counts_plot %>%
  mutate(
    fill_color   = ifelse(specified, vaccine_colour, "#D0D0D0"),
    # real hex or grey for non-specified
    border_color = ifelse(specified, "black", "#9A9A9A"),
    text_color   = ifelse(specified, "white", "#6B6B6B"),
    alpha_val    = ifelse(specified, 0.95, 0.35),
    # size_var as you used previously (n / sqrt(n) or log1p(n) etc.)
    size_var     = ifelse(n > 0, n / sqrt(n), 0)
  )

# ensure study_time_collected is an ordered factor with the timelevels you want
time_levels <- counts_plot %>%
  distinct(study_time_collected) %>%
  mutate(study_time_collected = as.numeric(study_time_collected)) %>%
  arrange(study_time_collected) %>%
  pull(study_time_collected) %>%
  as.character()

# Order the x labels by the timepoint ordering
counts_plot <- counts_plot %>%
  mutate(study_time_collected = factor(study_time_collected, levels = time_levels, ordered = TRUE))

# Order the y labels by the vaccine factor ordering
counts_plot <- counts_plot %>%
  mutate(vaccine_name = factor(
    as.character(vaccine_name),
    levels = vaccine_levels,
    ordered = TRUE
  ))

# Plot
p3 <- ggplot(counts_plot, aes(x = study_time_collected, y = vaccine_name)) +
  geom_point(
    aes(
      size = size_var,
      fill = fill_color,
      colour = border_color,
      alpha = alpha_val
    ),
    shape = 21,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = n, colour = text_color),
    size = 3.5,
    vjust = 0.5,
    show.legend = FALSE
  ) +
  # these identity scales make ggplot use the hex/literal colors we placed in columns
  scale_fill_identity(guide = "none") +
  scale_colour_identity(guide = "none") +
  scale_alpha_identity(guide = "none") +
  scale_size_area(max_size = 28, guide = "none") +
  labs(
    x = "Days post-vaccination",
    y = "Vaccine",
    title = "Number of participants with transcriptomic samples at selected timepoints",
    subtitle = "only participants with valid immune response data included"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 15, hjust = 0.5)
  )

print(p3)

ggsave(
  filename = "bubble_plot_sequential_specified.pdf",
  path = descriptive_figures_folder,
  plot = p3,
  width = 40,
  height = 20,
  units = "cm"
)

# sample identifier column
id_col <- "participant_id"


# build counts where a sample is only counted at time t if it has samples at all previous specified times
p4_counts_list <- lapply(names(specified_timepoints_list), function(vac) {
  times <- as.character(specified_timepoints_list[[vac]])
  df_v <- hipc_merged_all_norm_response %>%
    filter(vaccine_name == vac, !is.na(study_time_collected)) %>%
    mutate(study_time_collected = as.character(study_time_collected))
  
  if (nrow(df_v) == 0)
    return(NULL)
  
  # presence matrix per sample x specified-timepoint (1 if sample has an observation at that time)
  presence <- df_v %>%
    distinct(!!rlang::sym(id_col), study_time_collected) %>%
    filter(study_time_collected %in% times) %>%
    mutate(present = 1) %>%
    tidyr::complete(
      !!rlang::sym(id_col),
      study_time_collected = times,
      fill = list(present = 0)
    ) %>%
    arrange(!!rlang::sym(id_col), match(study_time_collected, times)) %>%
    group_by(!!rlang::sym(id_col)) %>%
    mutate(cum_all_prev = cumprod(present)) %>%   # 1 only when all up-to-and-including-this-timepoint are present
    ungroup()
  
  counts <- presence %>%
    group_by(study_time_collected) %>%
    summarise(n = sum(cum_all_prev), .groups = "drop") %>%
    mutate(vaccine_name = vac)
  
  # attach vaccine colour if present (take first if multiple)
  vaccine_colour <- df_v %>% distinct(vaccine_colour) %>% pull(vaccine_colour) %>% .[1]
  counts$vaccine_colour <- vaccine_colour
  
  counts
})

p4_counts <- bind_rows(p4_counts_list)

# make plotting columns (colors, sizes, labels)
# x-axis ordering: union of all specified timepoints, sorted numerically
time_levels <- sort(unique(as.numeric(unlist(
  specified_timepoints_list
)))) %>% as.character()

# Keep the vaccine factor ordering
vaccine_levels <- rev(unique(as.character(p4_counts$vaccine_name)))

counts_plot <- p4_counts %>%
  mutate(study_time_collected = as.character(study_time_collected)) %>%
  mutate(
    fill_color   = vaccine_colour,
    border_color = "black",
    text_color   = "white",
    alpha_val    = 0.95,
    size_var     = ifelse(n > 0, n / sqrt(n), 0)
  ) %>%
  mutate(
    study_time_collected = factor(study_time_collected, levels = time_levels, ordered = TRUE),
    vaccine_name = factor(vaccine_name, levels = vaccine_levels, ordered = TRUE)
  )

#  plot
p4 <- ggplot(counts_plot, aes(x = study_time_collected, y = vaccine_name)) +
  geom_point(
    aes(
      size = size_var,
      fill = fill_color,
      colour = border_color,
      alpha = alpha_val
    ),
    shape = 21,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = n, colour = text_color),
    size = 3.5,
    vjust = 0.5,
    show.legend = FALSE
  ) +
  scale_fill_identity(guide = "none") +
  scale_colour_identity(guide = "none") +
  scale_alpha_identity(guide = "none") +
  scale_size_area(max_size = 28, guide = "none") +
  labs(
    x = "Days post-vaccination",
    y = "Vaccine",
    title = "Cumulative number of participants with transcriptomic samples at all prior timepoints",
    subtitle = "Selected timepoints per-vaccine, only participants with valid immune response data included"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 15, hjust = 0.5)
  )

print(p4)

# Save (adjust path/filename as needed)
ggsave(
  filename = "bubble_plot_cumulative_speficied.pdf",
  path = descriptive_figures_folder,
  plot = p4,
  width = 38,
  height = 20,
  units = "cm"
)

# Based on exploration of antibody responses, we chose the following 5 vaccines on which to perform predictive analyses
# Yellow Fever (LV)
# Hepatitis A/B (IN/RP)
# Influenza (IN)
# Meningococcus (CJ)
# Meningococcus (PS)
# Below is the cumulative sample plot only for these 5 vaccines
p5 <- counts_plot %>% 
  filter(vaccine_name %in% c("Yellow Fever (LV)", "Hepatitis A/B (IN/RP)", "Influenza (IN)", "Meningococcus (CJ)", "Meningococcus (PS)")) %>% 
  ggplot(aes(x = study_time_collected, y = vaccine_name)) +
  geom_point(
    aes(
      size = size_var,
      fill = fill_color,
      colour = border_color,
      alpha = alpha_val
    ),
    shape = 21,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = n, colour = text_color),
    size = 3.5,
    vjust = 0.5,
    show.legend = FALSE
  ) +
  scale_fill_identity(guide = "none") +
  scale_colour_identity(guide = "none") +
  scale_alpha_identity(guide = "none") +
  scale_size_area(max_size = 28, guide = "none") +
  labs(
    x = "Days post-vaccination",
    y = "Vaccine",
    title = "Cumulative number of participants with transcriptomic samples at all prior timepoints",
    subtitle = "Selected vaccines and timepoints, only participants with valid immune response data included"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 15, hjust = 0.5)
  )

print(p5)

# Save (adjust path/filename as needed)
ggsave(
  filename = "bubble_plot_cumulative_speficied_selectedVaccines.pdf",
  path = descriptive_figures_folder,
  plot = p5,
  width = 38,
  height = 15,
  units = "cm"
)
