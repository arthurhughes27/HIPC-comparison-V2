# In this script, we explore potential study effects 
# In particular, we are interested in differences in distributions in baseline covariates between studies
# as well as differences in antibody distributions between studies
# For continuous variables we do histograms for each study 
# For categorical variables we do bar charts for each category within each study

# Libraries
library(dplyr)
library(forcats)
library(ggplot2)
library(ggtext)
library(scales)

# Directory to store engineered data
processed_data_folder = "data"

# Directory to store figures
descriptive_figures_folder = fs::path("output", "figures", "descriptive")

# Path to processed gene-level data and processed genesets
p_load_expr_all_norm <- fs::path(processed_data_folder, "hipc_merged_all_norm.rds")

# Load merged gene-level data
hipc_merged_all_norm = readRDS(p_load_expr_all_norm)

# Vaccines to analyse
vaccine_set = c("Yellow Fever (LV)", "Influenza (IN)", "Hepatitis A/B (IN/RP)", "Meningococcus (CJ)", "Meningococcus (PS)")

# ANTIBODIES #

# Select the relevant columns for this task
df_clinical_antibodies <- hipc_merged_all_norm %>%
  dplyr::select(
    participant_id,
    study_accession,
    vaccine_name,
    vaccine_colour,
    age_imputed,
    gender,
    race,
    ethnicity,
    immResp_MFC_anyAssay_log2_MFC
  ) %>%
  filter(!is.na(immResp_MFC_anyAssay_log2_MFC),
         vaccine_name %in% vaccine_set) %>%
  distinct() %>%
  mutate(
    study_accession = as.factor(study_accession),
    vaccine_name = as.factor(vaccine_name)
  )

# prepare a plotting dataframe
df_plot <- df_clinical_antibodies %>%
  mutate(
    study_accession = as.character(study_accession),   # ensure character for paste
    vaccine_name = as.character(vaccine_name),
    # one x value per study-vaccine combination
    study_vaccine = paste(study_accession, vaccine_name, sep = " - ")
  ) %>%
  # decide ordering so that boxes are grouped by vaccine_name (vaccine groups follow each other)
  arrange(vaccine_name, study_accession) %>%
  mutate(
    # set factor levels in the arranged order
    study_vaccine = factor(study_vaccine, levels = unique(study_vaccine)),
    vaccine_name = factor(vaccine_name)  # factor for consistent legend order
  )

# build a named vector vaccine_name -> vaccine_colour (ensure unique)
palette_df <- df_plot %>%
  distinct(vaccine_name, vaccine_colour) %>%
  arrange(vaccine_name)

vaccine_colors <- setNames(palette_df$vaccine_colour, palette_df$vaccine_name)

# plot
p1 <- ggplot(df_plot, aes(x = study_vaccine, y = immResp_MFC_anyAssay_log2_MFC, fill = vaccine_name)) +
  geom_boxplot(outlier.size = 1) +
  # use the provided colours mapped to vaccine_name
  scale_fill_manual(values = vaccine_colors, name = "Vaccine") +
  # label x ticks with study_accession only (text before the " - " we created)
  scale_x_discrete(labels = function(x) sub(" - .*", "", x)) +
  labs(
    x = "Study",
    y = "Log2(MFC)",
    title = "Distribution of Max Fold-Change per study"
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 25)
  )

print(p1)

ggsave(
  filename = "study_effects_MFC.pdf",
  path = descriptive_figures_folder,
  plot = p1,
  width = 40,
  height = 20,
  units = "cm"
)

# AGE #
p2 <- ggplot(df_plot, aes(x = study_vaccine, y = age_imputed, fill = vaccine_name)) +
  geom_boxplot(outlier.size = 1) +
  # use the provided colours mapped to vaccine_name
  scale_fill_manual(values = vaccine_colors, name = "Vaccine") +
  # label x ticks with study_accession only (text before the " - " we created)
  scale_x_discrete(labels = function(x) sub(" - .*", "", x)) +
  labs(
    x = "Study",
    y = "Age",
    title = "Distribution of age per study"
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 25)
  )

print(p2)

ggsave(
  filename = "study_effects_age.pdf",
  path = descriptive_figures_folder,
  plot = p2,
  width = 40,
  height = 20,
  units = "cm"
)

# GENDER #

sv_levels <- levels(df_plot$study_vaccine)

labels_html <- vapply(sv_levels, function(sv) {
  study_label <- sub(" - .*", "", sv)
  vaccine_label <- sub(".* - ", "", sv)
  col <- vaccine_colors[[vaccine_label]]
  if (is.null(col) || is.na(col) || col == "") col <- "#000000"
  paste0("<span style='color:", col, ";font-weight:700;'>", study_label, "</span>")
}, FUN.VALUE = character(1), USE.NAMES = FALSE)
names(labels_html) <- sv_levels

df_legend_points <- df_plot %>%
  distinct(study_vaccine, vaccine_name) %>%
  group_by(vaccine_name) %>%
  slice(1) %>%
  ungroup()

# ---- corrected plot: prevent geom_point from inheriting fill=gender ----
p3 <- ggplot(df_plot %>% filter(!is.na(gender)), aes(x = study_vaccine, fill = gender)) +
  geom_bar(position = "fill", width = 0.9) +
  # IMPORTANT: set inherit.aes = FALSE so this layer doesn't expect 'gender'
  geom_point(
    data = df_legend_points,
    mapping = aes(x = study_vaccine, y = 0, color = vaccine_name),
    size = 3,
    alpha = 0,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_discrete(labels = labels_html) +
  scale_color_manual(values = vaccine_colors, name = "Vaccine") +
  labs(
    x = "Study",
    y = "Proportion",
    title = "Gender composition per Study",
    fill = "Gender",
    color = "Vaccine"
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1, size = 12),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 18)
  ) +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2, override.aes = list(alpha = 1, size = 5))
  )

print(p3)

ggsave(
  filename = "study_gender_proportions.pdf",
  path = descriptive_figures_folder,
  plot = p3,
  width = 40,
  height = 20,
  units = "cm"
)

## RACE ##

p4 <- ggplot(df_plot %>% filter(!is.na(race)), aes(x = study_vaccine, fill = race)) +
  geom_bar(position = "fill", width = 0.9) +
  # IMPORTANT: set inherit.aes = FALSE so this layer doesn't expect 'race'
  geom_point(
    data = df_legend_points,
    mapping = aes(x = study_vaccine, y = 0, color = vaccine_name),
    size = 3,
    alpha = 0,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_discrete(labels = labels_html) +
  scale_color_manual(values = vaccine_colors, name = "Vaccine") +
  labs(
    x = "Study",
    y = "Proportion",
    title = "Race composition per Study",
    fill = "race",
    color = "Vaccine"
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1, size = 12),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 18)
  ) +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2, override.aes = list(alpha = 1, size = 5))
  )

print(p4)

ggsave(
  filename = "study_race_proportions.pdf",
  path = descriptive_figures_folder,
  plot = p4,
  width = 40,
  height = 20,
  units = "cm"
)

## ETHNICITY ##

p5 <- ggplot(df_plot %>% filter(!is.na(ethnicity)), aes(x = study_vaccine, fill = ethnicity)) +
  geom_bar(position = "fill", width = 0.9) +
  # IMPORTANT: set inherit.aes = FALSE so this layer doesn't expect 'ethnicity'
  geom_point(
    data = df_legend_points,
    mapping = aes(x = study_vaccine, y = 0, color = vaccine_name),
    size = 3,
    alpha = 0,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_discrete(labels = labels_html) +
  scale_color_manual(values = vaccine_colors, name = "Vaccine") +
  labs(
    x = "Study",
    y = "Proportion",
    title = "Ethnicity composition per Study",
    fill = "ethnicity",
    color = "Vaccine"
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1, size = 12),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 18)
  ) +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2, override.aes = list(alpha = 1, size = 5))
  )

print(p5)

ggsave(
  filename = "study_ethnicity_proportions.pdf",
  path = descriptive_figures_folder,
  plot = p5,
  width = 40,
  height = 20,
  units = "cm"
)

