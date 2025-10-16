# File for preprocessing the raw immune response data files
## In particular, we are interested in deriving "maximum fold-change in immune response" (MFC) values in the first 28 (+-7) days after vaccination
## We use the fold-change to take into account baseline immunity
## We use the maximum of fold-change values to choose one measurement across multiple analytes, timepoints, or assays
## In this script, the final result will be an immune response dataframe with one row per unique participant identifier
## We will have one column for the MFC values for each of three assay types (neutralising antibodies, ELISA, and HAI)
## We will also create a column giving the maximum of any of the assay types

# Load packages
library(fs)
library(dplyr)
library(stringr)
library(janitor)
library(readxl)
library(purrr)


# Specify folder within folder root where the raw data lives
raw_data_folder = "data-raw"
processed_data_folder = "data"

# Use fs::path() to specify the data paths robustly
p_load_elisa <- fs::path(raw_data_folder, "elisa_2025-01-10_01-14-21.xlsx")
p_load_nAb <- fs::path(raw_data_folder, "neut_ab_titer_2025-01-10_01-13-22.xlsx")
p_load_hai <- fs::path(raw_data_folder, "hai_2025-01-10_01-13-41.xlsx")
p_load_clinical <- fs::path(processed_data_folder, "hipc_clinical.rds")

# HAI response data
response_hai = read_excel(p_load_hai) %>%
  clean_names()

# nAb response data
response_nAb = read_excel(p_load_nAb) %>%
  clean_names()

# ELISA response data
response_elisa = read_excel(p_load_elisa) %>%
  clean_names()

# Clinical data for filtration of participants
hipc_clinical = readRDS(p_load_clinical)

# First, filter each dataframe to only contain information on participants for which we have gene expression data
# Find identifiers of participants with gene expression measurements
participants = hipc_clinical %>%
  pull(participant_id) %>%
  unique()

# Filter immune response data by these participants
response_hai = response_hai %>%
  filter(participant_id %in% participants)

response_nAb = response_nAb %>%
  filter(participant_id %in% participants)

response_elisa = response_elisa %>%
  filter(participant_id %in% participants)

# Depending on the assay, there may be multiple viral strains or analytes measured.
# We rename some columns and add a column for the assay name, so that we can merge these data together.

response_hai = response_hai %>%
  rename(response_strain_analyte = virus) %>%
  mutate(assay = "hai")

response_nAb = response_nAb %>%
  rename(response_strain_analyte = virus) %>%
  mutate(assay = "nAb")

# Filter ELISA data for IgG and Hepatitis B antibodies
response_elisa = response_elisa %>%
  filter(grepl("IgG", analyte) |
           analyte == "Hepatitis B Virus Surface Antibody") %>%
  rename(response_strain_analyte = analyte) %>%
  mutate(assay = "elisa")

# Now merge all the raw response data
response_raw_merged = bind_rows(response_nAb, response_elisa, response_hai) %>%
  arrange(participant_id) %>%
  distinct()

# First get the studies and other clinical data corresponding to each participant id
hipc_studies = hipc_clinical %>%
  select(participant_id,
         study_accession,
         vaccine,
         vaccine_type,
         pathogen) %>%
  distinct()

# Merge the study names into the immune response data (it is not directly given)
response_raw_merged_studies = merge(x = response_raw_merged,
                                    y = hipc_studies,
                                    by = "participant_id",
                                    all = F)

# There are a series of study-specific errors which have been inferred from the "immuneResponseCallGeneration.R" script in the ImmuneSignatures2 GitHub page.
# We fix these here.
# First, "SDY1276" is apparently scaled in log4.

response_raw_merged_studies = response_raw_merged_studies %>%
  mutate(value_preferred = if_else(
    study_accession == "SDY1276",
    4^value_preferred,
    value_preferred
  ))

# Next, "SDY1289" has some nAb baseline values at exactly 0, these are set to 1.

response_raw_merged_studies = response_raw_merged_studies %>%
  mutate(value_preferred = if_else(
    (
      study_accession == "SDY1289" &
        value_preferred == 0 &
        as.numeric(study_time_collected) == 0
    ) ,
    1,
    value_preferred
  ))

# "SDY1264" lacks baseline data. The investigators create this row, where the baseline value is set to 1 for everyone.

day_zero <- response_raw_merged_studies %>%
  filter(study_accession == "SDY1264") %>%        # pull out SDY1264
  distinct(participant_id, .keep_all = TRUE) %>%   # keep only one row per participant
  mutate(
    study_time_collected = 0,
    # reset time to “0”
    value_preferred      = 1                     # set preferred value to 1
  )

# Bind them back on to the full dataset
response_raw_merged_studies <- bind_rows(response_raw_merged_studies, day_zero)

# "SDY1328" has an incorrect time label. Apparently, the data which is claimed to be taken at day 7 should actually be at day 30.

response_raw_merged_studies = response_raw_merged_studies %>%
  mutate(study_time_collected = if_else(
    (study_accession == "SDY1328" &
       study_time_collected == 7) ,
    30,
    study_time_collected
  ))

# Additionally, this study has some anomalies at day 0, where data suggests these participants have antibody responses prior to vaccination.
# The investigators set baseline values as well as values of 2.5 to 1.

response_raw_merged_studies = response_raw_merged_studies %>%
  mutate(value_preferred = if_else((
    study_accession == "SDY1328" &
      (value_preferred == 2.5 | study_time_collected == 0)
  ), 1, value_preferred))

# Furthermore, apparently these values need to be transformed with log2. This also applies to SDY984.

response_raw_merged_studies = response_raw_merged_studies %>%
  mutate(value_preferred = if_else(((study_accession == "SDY1328" |
                                       study_accession == "SDY984") &
                                      assay == "elisa"
  ), log2(value_preferred), value_preferred))

# For "SDY1260", the values need to be de-transformed from log2.

response_raw_merged_studies = response_raw_merged_studies %>%
  mutate(value_preferred = if_else(
    (study_accession == "SDY1260" &
       assay == "elisa"),
    2^value_preferred,
    value_preferred
  ))

# Furthermore, the IgG serotype values need to be summed within each participant, timepoint, and vaccine.

response_raw_merged_studies <- response_raw_merged_studies %>%
  group_by(participant_id,
           study_time_collected,
           vaccine,
           vaccine_type,
           pathogen) %>%
  mutate(value_preferred = if_else(
    (study_accession == "SDY1260" & assay == "elisa"),
    # sum only the SDY1260 values within this group, then log2
    log2(sum(value_preferred[study_accession == "SDY1260"])),
    # otherwise leave unchanged
    value_preferred
  )) %>%
  ungroup()

# We intend to predict on immune response data at day 28 (plus or minus 7 days), so we filter to get the immune data at these timepoints or pre-vaccination.

response_raw_merged_studies = response_raw_merged_studies %>%
  filter((study_time_collected < 36 &
            study_time_collected > 20) | study_time_collected <= 0)

response_raw_merged = response_raw_merged_studies %>%
  select(-c(study_accession, vaccine, vaccine_type, pathogen))

# Now derive MFC values, taking most recent measurement as baseline in the case where there are two pre-vaccination measurements

pre_vax <- response_raw_merged %>%
  filter(study_time_collected <= 0) %>%
  group_by(participant_id, assay, response_strain_analyte) %>%
  slice_max(study_time_collected, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(response_MFC_pre_value = value_preferred,
         response_MFC_pre_time = study_time_collected)

post_vax <- response_raw_merged %>%
  filter(study_time_collected > 0) %>%
  rename(response_MFC_post_value = value_preferred,
         response_MFC_post_time = study_time_collected)

# Now we merge pre- and post- values to get a dataframe in wide format with one line per unique combination of participant/assay/strain-analyte
merged_vax <- full_join(
  pre_vax,
  post_vax,
  by = c("participant_id", "response_strain_analyte", "assay"),
  relationship = "many-to-many"
) %>%
  select(
    participant_id,
    assay,
    response_strain_analyte,
    response_MFC_pre_time,
    response_MFC_post_time,
    response_MFC_pre_value,
    response_MFC_post_value
  )

# Now we can exclude any rows with NA values in either pre- or post-.
merged_vax = merged_vax %>%
  filter(!is.na(response_MFC_pre_value),
         !is.na(response_MFC_post_value))

# Derive the log2 fold changes.
response_MFC_df <- merged_vax %>%
  mutate(
    # Define fold change
    response_MFC = response_MFC_post_value / response_MFC_pre_value,
    # If pre-vaccination value is 0, set fold change to post-vaccination value
    response_MFC = ifelse(
      is.infinite(response_MFC),
      response_MFC_post_value,
      response_MFC
    ),
    # Set fold change to 1 if both response_MFC_pre_value and response_MFC_post_value are 0
    response_MFC = ifelse(
      response_MFC_pre_value == 0 &
        response_MFC_post_value == 0,
      1,
      response_MFC
    )
  ) %>%
  # Calculate log2 fold change
  mutate(
    response_log2_MFC = log2(response_MFC),
    response_log2_MFC = ifelse(is.infinite(response_log2_MFC), NA, response_log2_MFC)
  ) %>%
  # Select relevant columns
  select(
    participant_id,
    response_strain_analyte,
    assay,
    response_MFC_pre_time,
    response_MFC_post_time,
    response_MFC_pre_value,
    response_MFC_post_value,
    response_MFC,
    response_log2_MFC
  )

# Now, for each participant, we can derive the maximum of these fold changes.
## Across all assays
max_response_MFC_df_any <- response_MFC_df %>%
  group_by(participant_id) %>%
  slice_max(response_MFC, n = 1, with_ties = F) %>%
  ungroup() %>%
  rename(
    immResp_MFC_anyAssay_response_strain_analyte = response_strain_analyte,
    immResp_MFC_anyAssay_assay = assay,
    immResp_MFC_anyAssay_pre_time = response_MFC_pre_time,
    immResp_MFC_anyAssay_post_time = response_MFC_post_time,
    immResp_MFC_anyAssay_pre_value = response_MFC_pre_value,
    immResp_MFC_anyAssay_post_value = response_MFC_post_value,
    immResp_MFC_anyAssay_MFC = response_MFC,
    immResp_MFC_anyAssay_log2_MFC = response_log2_MFC
  )

## And within each assay
max_response_MFC_df_each <- response_MFC_df %>%
  group_by(participant_id, assay) %>%
  slice_max(response_MFC, n = 1, with_ties = F) %>%
  ungroup()

max_response_MFC_df_nAb = max_response_MFC_df_each %>%
  filter(assay == "nAb") %>%
  rename(
    immResp_MFC_nAb_response_strain_analyte = response_strain_analyte,
    immResp_MFC_nAb_pre_time = response_MFC_pre_time,
    immResp_MFC_nAb_post_time = response_MFC_post_time,
    immResp_MFC_nAb_pre_value = response_MFC_pre_value,
    immResp_MFC_nAb_post_value = response_MFC_post_value,
    immResp_MFC_nAb_MFC = response_MFC,
    immResp_MFC_nAb_log2_MFC = response_log2_MFC
  ) %>%
  select(-assay)


max_response_MFC_df_elisa = max_response_MFC_df_each %>%
  filter(assay == "elisa") %>%
  rename(
    immResp_MFC_elisa_response_strain_analyte = response_strain_analyte,
    immResp_MFC_elisa_pre_time = response_MFC_pre_time,
    immResp_MFC_elisa_post_time = response_MFC_post_time,
    immResp_MFC_elisa_pre_value = response_MFC_pre_value,
    immResp_MFC_elisa_post_value = response_MFC_post_value,
    immResp_MFC_elisa_MFC = response_MFC,
    immResp_MFC_elisa_log2_MFC = response_log2_MFC
  ) %>%
  select(-assay)

max_response_MFC_df_hai = max_response_MFC_df_each %>%
  filter(assay == "hai") %>%
  rename(
    immResp_MFC_hai_response_strain_analyte = response_strain_analyte,
    immResp_MFC_hai_pre_time = response_MFC_pre_time,
    immResp_MFC_hai_post_time = response_MFC_post_time,
    immResp_MFC_hai_pre_value = response_MFC_pre_value,
    immResp_MFC_hai_post_value = response_MFC_post_value,
    immResp_MFC_hai_MFC = response_MFC,
    immResp_MFC_hai_log2_MFC = response_log2_MFC
  ) %>%
  select(-assay)

# Now merge those dataframes together to get the final immune response data
hipc_immResp <- list(
  max_response_MFC_df_any,
  max_response_MFC_df_nAb,
  max_response_MFC_df_elisa,
  max_response_MFC_df_hai
) %>%
  reduce(full_join, by = "participant_id") %>%
  arrange(participant_id)

# Use fs::path() to specify the data path robustly
p_save <- fs::path(processed_data_folder, "hipc_immResp.rds")

# Save dataframe
saveRDS(hipc_immResp, file = p_save)

rm(list = ls())