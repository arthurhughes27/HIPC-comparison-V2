# File to pre-process the HIPC IS2 clinical characteristics data for use
## This includes some data engineering on columns (e.g. collapse "unknown" and "Not Specified" gender to same value)
## and also the creation of some data for downstream analysis purposes (e.g. the associating a colour with each vaccine)

# Load necessary packages
library(fs)
library(Biobase)
library(dplyr)
library(forcats)
library(stringr)

# Specify folder within folder root where the raw data lives
raw_data_folder = "data-raw"

# Use fs::path() to specify the data path robustly
p_load <- fs::path(raw_data_folder, "all_norm_eset.rds")

# Read in the rds file
all_norm_eset <- readRDS(p_load)

# Extract the clinical data as a dataframe
hipc_clinical = all_norm_eset@phenoData@data %>% 
  as.data.frame()

# Data engineering
# "collapsing" the values from certain important columns
# For gender and race, "Not Specified" and "Unknown" represent the same thing
# For ethnicity, the values "Not Hispanic or Latino" and "Other" represent the same thing (the only other values are
# "Not Specified" and "Hispanic or Latino").
# In addition, we round numerical time columns to avoid long recurring numbers
hipc_clinical <- hipc_clinical %>%
  mutate(
    gender    = fct_collapse(factor(gender), Unknown = c("Not Specified", "Unknown")),
    race      = fct_collapse(factor(race), Unknown = c("Not Specified", "Unknown")),
    ethnicity = fct_collapse(factor(ethnicity), Other   = c("Not Hispanic or Latino", "Other")),
    study_time_collected = round(study_time_collected, 2),
    time_post_last_vax = round(time_post_last_vax, 2)
  )

# Abbreviate the vaccine type column
hipc_clinical$vaccine_type <- hipc_clinical$vaccine_type %>%
  as.factor() %>%
  recode_factor(
    "Conjugate" = "CJ",
    "Inactivated" = "IN",
    "Inactivated/Recombinant protein" = "IN/RP",
    "Live virus" = "LV",
    "Polysaccharide" = "PS",
    "Recombinant Viral Vector" = "RVV",
    "Recombinant protein" = "RP"
  )

# Create vaccine name column by combining pathogen and vaccine type
hipc_clinical <- hipc_clinical %>%
  mutate(vaccine_name = str_c(pathogen, " (", vaccine_type, ")"))

# Define an ordering for the vaccines (this is for later to make figures consistent)
conditions_order <- c(
  "Tuberculosis (RVV)",
  "Varicella Zoster (LV)",
  "Yellow Fever (LV)",
  "Ebola (RVV)",
  "Hepatitis A/B (IN/RP)",
  "HIV (RVV)",
  "Influenza (IN)",
  "Influenza (LV)",
  "Malaria (RP)",
  "Meningococcus (CJ)",
  "Meningococcus (PS)",
  "Pneumococcus (PS)",
  "Smallpox (LV)"
)


# Assign this order to the vaccine names
hipc_clinical <- hipc_clinical %>% 
  mutate(vaccine_name = factor(vaccine_name, levels = conditions_order))

# Define a colour for each vaccine
## This colour palette was chosen to maximise visual distinctiveness for 13 vaccines
## Using the "iwanthue" tool (https://medialab.github.io/iwanthue/)
color_palette = c(
  "#b94a73",
  "#c6aa3c",
  "#6f71d9",
  "#64c46a",
  "#be62c2",
  "#7d973c",
  "#563382",
  "#4ea76e",
  "#bc69b0",
  "#33d4d1",
  "#bb4c41",
  "#6a87d3",
  "#b57736"
)

# Write a helper function to assign the colours
assign_color <- function(vaccine_name) {
  return(color_palette[match(hipc_clinical$vaccine_name,
                             levels(hipc_clinical$vaccine_name))])
}

# Assign the colours to the vaccine names
hipc_clinical$vaccine_colour <-
  assign_color(hipc_clinical$vaccine_name)

# Final dataframe to be saved has samples as rows and variables as columns
dim(hipc_clinical)

# Save processed dataframe

# Specify folder within folder root where the processed data lives
processed_data_folder = "data"

# Use fs::path() to specify the data path robustly
p_save <- fs::path(processed_data_folder, "hipc_clinical.rds")

# Save dataframe
saveRDS(hipc_clinical, file = p_save)

rm(list = ls())
