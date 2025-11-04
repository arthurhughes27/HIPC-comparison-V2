# R Script to describe the different studies present in the HIPC IS2 dataset

# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(tibble)

# Directory to store engineered data
processed_data_folder = "data"

# Directory to store figures
descriptive_figures_folder = fs::path("output", "figures", "descriptive")

# Path to processed gene-level data
p_load_expr_all_norm <- fs::path(processed_data_folder, "hipc_merged_all_norm.rds")

# Load merged gene-level data
hipc_merged_all_norm = readRDS(p_load_expr_all_norm)

# Counts per vaccine x time
counts <- hipc_merged_all_norm %>%
  filter(!is.na(time_post_last_vax)) %>%
  group_by(study_accession_unique, vaccine_colour, time_post_last_vax, vaccine_name) %>%
  summarise(n = n(), .groups = "drop")

# Order the time points and make time_post_last_vax an ordered factor
time_levels <- counts %>%
  distinct(time_post_last_vax) %>%
  arrange(as.numeric(time_post_last_vax)) %>%
  pull(time_post_last_vax) %>%
  as.character()   # factor levels must be character

# Order the counts by the study time
counts <- counts %>%
  mutate(time_post_last_vax = factor(
    as.character(time_post_last_vax),
    levels = time_levels,
    ordered = TRUE
  ))

# Since the range of the number of samples is large, we make the bubble size proportional to the count/square root of count
counts <- counts %>%
  mutate(size_var = n^{2/3})

# named vector of fill values for scale_fill_manual()
vaccine_map <- hipc_merged_all_norm %>%
  distinct(study_accession_unique, vaccine_name, vaccine_colour)

fill_values <- vaccine_map %>%
  distinct(vaccine_name, vaccine_colour) %>%
  { setNames(.$vaccine_colour, .$vaccine_name) }

# the raw counts you want to show in the size legend:
size_breaks_counts <- c(10, 50, 100, 200)
# convert them to the scale used in the plot (sqrt)
size_breaks <- (size_breaks_counts)^{2/3}

# Plot
p1 <- ggplot(counts, aes(x = time_post_last_vax, y = study_accession_unique)) +
  geom_point(
    aes(size = size_var, fill = vaccine_name),
    shape = 21, colour = "black", alpha = 0.95, show.legend = TRUE
  ) +
  geom_text(
    aes(label = n),
    colour = "white", size = 3.5, vjust = 0.5, show.legend = FALSE
  ) +
  # Fill legend: vaccine names with their hex colours
  scale_fill_manual(
    name = "Vaccine",
    values = fill_values,
    guide = guide_legend(override.aes = list(shape = 21, size = 6, colour = "black"))
  ) +
  # Size legend: show 10 / 50 / 100 / 200 as legend entries (we pass their sqrt values)
  scale_size_area(
    name = "Counts",
    max_size = 28,
    breaks = size_breaks,
    labels = size_breaks_counts,
    guide = guide_legend(override.aes = list(fill = "grey80", colour = "black"))
  ) +
  labs(
    x = "Days post-vaccination",
    y = "Vaccine",
    title = "Number of participants with transcriptomic samples at each timepoint per study",
    subtitle = "All timepoints"
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
  filename = "study_bubble_plot_sequential.pdf",
  path = descriptive_figures_folder,
  plot = p1,
  width = 40,
  height = 50,
  units = "cm"
)

# Based on exploration of antibody responses, we chose the following 5 vaccines on which to perform predictive analyses
# Yellow Fever (LV)
# Hepatitis A/B (IN/RP)
# Influenza (IN)
# Meningococcus (CJ)
# Meningococcus (PS)

specified_timepoints_list <- list(
  "Meningococcus (PS)" = c(0, 3, 7),
  "Meningococcus (CJ)" = c(0, 3, 7),
  "Influenza (IN)" = c(0, 1, 3, 7, 14),
  "Hepatitis A/B (IN/RP)" = c(0, 7),
  "Yellow Fever (LV)" = c(0, 3, 7, 10, 14, 28)
)

# counts2 per vaccine x time
counts2 <- hipc_merged_all_norm %>%
  filter(!is.na(time_post_last_vax),
         !is.na(immResp_MFC_anyAssay_log2_MFC)) %>%
  group_by(study_accession_unique, vaccine_colour, time_post_last_vax, vaccine_name) %>%
  summarise(n = n(), .groups = "drop")

# Order the time points and make time_post_last_vax an ordered factor
time_levels <- counts2 %>%
  distinct(time_post_last_vax) %>%
  arrange(as.numeric(time_post_last_vax)) %>%
  pull(time_post_last_vax) %>%
  as.character()   # factor levels must be character

# Order the counts2 by the study time
counts2 <- counts2 %>%
  mutate(time_post_last_vax = factor(
    as.character(time_post_last_vax),
    levels = time_levels,
    ordered = TRUE
  ))

# Since the range of the number of samples is large, we make the bubble size proportional to the count/square root of count
counts2 <- counts2 %>%
  mutate(size_var = n^{2/3})

# named vector of fill values for scale_fill_manual()
fill_values <- vaccine_map %>%
  distinct(vaccine_name, vaccine_colour) %>%
  { setNames(.$vaccine_colour, .$vaccine_name) }

# the raw counts2 you want to show in the size legend:
size_breaks_counts2 <- c(10, 50, 100, 200)
# convert them to the scale used in the plot (sqrt)
size_breaks <- (size_breaks_counts2)^{2/3}

# expand to long table (one row per vaccine × timepoint)
allowed_tbl <- enframe(specified_timepoints_list, name = "vaccine_name", value = "timepoints") %>%
  unnest_longer(timepoints) %>%
  mutate(timepoints = as.character(timepoints))

# join & filter counts22 so only allowed vaccine × time combos remain
counts2_sel <- counts2 %>%
  inner_join(allowed_tbl, by = "vaccine_name") %>%
  # keep only exact timepoint matches (time_post_last_vax is factor of characters)
  filter(as.character(time_post_last_vax) == timepoints) %>%
  select(-timepoints)

# if you want time levels ordered only by the selected times (useful for x-axis)
selected_time_levels <- counts2_sel %>%
  distinct(time_post_last_vax) %>%
  arrange(as.numeric(as.character(time_post_last_vax))) %>%
  pull(time_post_last_vax) %>%
  as.character()

counts2_sel <- counts2_sel %>%
  mutate(time_post_last_vax = factor(as.character(time_post_last_vax),
                                     levels = selected_time_levels,
                                     ordered = TRUE))

# recompute fill_values just for selected vaccines
fill_values <- counts2_sel %>%
  distinct(vaccine_name, vaccine_colour) %>%
  { setNames(.$vaccine_colour, .$vaccine_name) }

# size legend choices and mapping (you used exponent 2/3)
size_breaks_counts2 <- c(10, 50, 100, 200)
size_breaks <- (size_breaks_counts2)^(2/3)

# ---- Plot restricted to selected vaccines × timepoints ----
p2 <- ggplot(counts2_sel, aes(x = time_post_last_vax, y = study_accession_unique)) +
  geom_point(
    aes(size = size_var, fill = vaccine_name),
    shape = 21, colour = "black", alpha = 0.95, show.legend = TRUE
  ) +
  geom_text(
    aes(label = n),
    colour = "white", size = 3.5, vjust = 0.5, show.legend = FALSE
  ) +
  scale_fill_manual(
    name = "Vaccine",
    values = fill_values,
    guide = guide_legend(override.aes = list(shape = 21, size = 6, colour = "black"))
  ) +
  scale_size_area(
    name = "Count",
    max_size = 28,
    breaks = size_breaks,
    labels = size_breaks_counts2,
    guide = guide_legend(override.aes = list(fill = "grey80", colour = "black"))
  ) +
  labs(
    x = "Days post-vaccination",
    y = "Vaccine",
    title = "Number of participants with transcriptomic samples at each timepoint per study",
    subtitle = "Selected vaccines and timepoints"
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
  filename = "study_bubble_plot_sequential.pdf",
  path = descriptive_figures_folder,
  plot = p2,
  width = 40,
  height = 30,
  units = "cm"
)
