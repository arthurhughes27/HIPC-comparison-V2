# R script to perform some descriptions of the available immune response samples
## In particular, we are interested in exploring which vaccines have sufficient heterogeneity in antibody responses
## to justify performing a predictive analysis
## We will explore both individual trajectories as well as the distribution of the maximum fold-changes

# Packages
library(tidyverse)

# Directory to store engineered data
processed_data_folder = "data"

# Directory to store figures
descriptive_figures_folder = fs::path("output", "figures", "descriptive")

# Path to processed, merged data
p_load_expr_all_norm <- fs::path(processed_data_folder, "hipc_merged_all_norm.rds")

# Load merged data
hipc_merged_all_norm = readRDS(p_load_expr_all_norm)

# Filter out observations which do not have an immune response value
hipc_response = hipc_merged_all_norm %>%
  filter(!is.na(immResp_MFC_anyAssay_log2_MFC)) %>%
  select(
    participant_id,
    vaccine_name,
    vaccine_colour,
    immResp_MFC_anyAssay_pre_value,
    immResp_MFC_anyAssay_post_value,
    immResp_MFC_anyAssay_log2_MFC
  ) %>%
  distinct()

# Individual trajectories of each vaccine

# build a named vector of hex colours: names = vaccine_name, values = hex codes
vaccine_colours <- hipc_response %>%
  distinct(vaccine_name, vaccine_colour) %>%
  # in case of duplicates keep first occurrence
  group_by(vaccine_name) %>%
  slice(1) %>%
  ungroup() %>%
  with(setNames(vaccine_colour, vaccine_name))

# pivot to long format (pre / post)
long_hipc_response <- hipc_response %>%
  pivot_longer(
    cols = c(
      immResp_MFC_anyAssay_pre_value,
      immResp_MFC_anyAssay_post_value
    ),
    names_to = "timepoint",
    # capture pre / post from the column names
    names_pattern = ".*_(pre|post)_value",
    values_to = "antibody_value"
  ) %>%
  mutate(
    # make timepoint an ordered factor for plotting on x axis
    timepoint = factor(
      timepoint,
      levels = c("pre", "post"),
      labels = c("Baseline", "Day 28 +- 7 days")
    )
  ) %>%
  arrange(vaccine_name, participant_id, timepoint)

# Plot
p1 <- ggplot(
  long_hipc_response,
  aes(
    x = timepoint,
    y = log2(antibody_value + 0.001),
    group = participant_id,
    color = vaccine_name
  )
) +
  geom_line(alpha = 0.3, size = 0.6) +        # trajectory line per participant
  geom_point(size = 1.8, alpha = 0.5) +                    # points at baseline & post
  scale_color_manual(values = vaccine_colours, name = "Vaccine") +
  facet_wrap(~ vaccine_name, ncol = 3, scales = "free_y") +
  labs(x = "Timepoint", y = "log2(Antibody value)", title = "Individual trajectories of antibody response per vaccine") +
  # add a small expansion to avoid clipping at the bottom or top
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  theme_minimal(base_size = 14) +  # increase base font size
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    strip.text = element_text(face = "bold", size = 14),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    ),
    legend.position = "none"
  )

print(p1)

ggsave(
  filename = "individual_antibody_trajectories.pdf",
  path = descriptive_figures_folder,
  plot = p1,
  width = 40,
  height = 20,
  units = "cm"
)

p2 <- ggplot(
  long_hipc_response,
  aes(
    x = timepoint,
    y = log2(antibody_value + 0.001),
    group = participant_id,
    color = vaccine_name
  )
) +
  
  # individual points (transparent)
  geom_point(
    size = 1.6,
    alpha = 0.35,
    position = position_jitter(width = 0.05, height = 0)
  ) +
  
  # mean trajectory per vaccine (connect mean at each timepoint)
  stat_summary(
    fun = mean,
    geom = "line",
    aes(group = vaccine_name),
    # ensures one mean line per vaccine within each facet
    size = 1.2
  ) +
  
  # mean points
  stat_summary(
    fun = mean,
    geom = "point",
    aes(group = vaccine_name),
    size = 3.0,
    shape = 21,
    stroke = 0.75,
    fill = "white"
  ) +
  
  scale_color_manual(values = vaccine_colours, name = "Vaccine") +
  
  facet_wrap(~ vaccine_name, ncol = 3, scales = "free_y") +
  
  labs(x = "Timepoint", y = "log2(Antibody value)", title = "Mean trajectory of antibody response per vaccine") +
  
  # add a small expansion to avoid clipping at the bottom or top
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    strip.text = element_text(face = "bold", size = 14),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    ),
    legend.position = "none"
  )

print(p2)

ggsave(
  filename = "mean_antibody_trajectories.pdf",
  path = descriptive_figures_folder,
  plot = p2,
  width = 40,
  height = 20,
  units = "cm"
)

p3 <- ggplot(
  hipc_response,
  aes(x = vaccine_name, y = immResp_MFC_anyAssay_log2_MFC, fill = vaccine_name)
) +
  
  geom_boxplot(outlier.shape = 16,
               outlier.alpha = 0.4,
               width = 0.6) +
  
  scale_fill_manual(values = vaccine_colours, name = "Vaccine") +
  
  labs(x = "Vaccine", y = "log2(Max fold-change)", title = "Distribution of antibody responses per vaccine") +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12
    ),
    # show rotated vaccine names
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    ),
    legend.position = "none"
  )

print(p3)

# Based on this plot, we choose to make predictions for the following 5 vaccines which have sufficient heterogeneity : 
# Yellow Fever (LV)
# Hepatitis A/B (IN/RP)
# Influenza (IN)
# Meningococcus (CJ)
# Meningococcus (PS)

ggsave(
  filename = "MFC_boxplots.pdf",
  path = descriptive_figures_folder,
  plot = p3,
  width = 40,
  height = 20,
  units = "cm"
)
