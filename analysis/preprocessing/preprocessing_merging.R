# R script to merge the clinical, immune response and expression data together

# Packages
library(fs)
library(dplyr)

# Specify folder within folder root where the raw data lives
processed_data_folder = "data"

# Use fs::path() to specify the data paths robustly
p_load_expr_young_noNorm <- fs::path(processed_data_folder, "young_noNorm_expr.rds")
p_load_expr_all_norm <- fs::path(processed_data_folder, "all_norm_expr.rds")
p_load_clinical <- fs::path(processed_data_folder, "hipc_clinical.rds")
p_load_immResp <- fs::path(processed_data_folder, "hipc_immResp.rds")

# Read in the files
expr_young_noNorm <- readRDS(p_load_expr_young_noNorm)
expr_all_norm <- readRDS(p_load_expr_all_norm)
hipc_clinical <- readRDS(p_load_clinical)
hipc_immResp <- readRDS(p_load_immResp)

# Merge together the clinical and immune response dataframes
merged_hipc_clinical_immresp = full_join(x = hipc_clinical, y = hipc_immResp, by = "participant_id")

# We are going to merge this to the expression dataframe by the participant id and the the timepoint.
# Let's check that these uniquely identify samples
count_unique = merged_hipc_clinical_immresp %>%
  summarise(count = n(),
            .by = c(participant_id, study_time_collected)) %>%
  arrange(desc(count))

head(count_unique)

# We see that there is one participant with 2 samples at day 0.
# We can choose to arbitrarily remove one of the rows of this individual
# (since it is only one individual from a large study of YF17D, this shouldn't make any difference)

merged_hipc_clinical_immresp  = merged_hipc_clinical_immresp %>%
  distinct(participant_id, study_time_collected, .keep_all = T)

# Check this worked
count_unique = merged_hipc_clinical_immresp %>%
  summarise(count = n(),
            .by = c(participant_id, study_time_collected)) %>%
  arrange(desc(count))

head(count_unique)

# Let's do the same check for the expression data
count_unique_expr = expr_young_noNorm %>%
  summarise(count = n(),
            .by = c(participant_id, study_time_collected)) %>%
  arrange(desc(count))

head(count_unique_expr)

# Again, one participant has multiple measurements for a given timepoint

expr_young_noNorm  = expr_young_noNorm %>%
  distinct(participant_id, study_time_collected, .keep_all = T)

expr_all_norm  = expr_all_norm %>%
  distinct(participant_id, study_time_collected, .keep_all = T)

# Now we can merge these dataframes together by pid and study time
hipc_merged_young_noNorm = right_join(
  x = merged_hipc_clinical_immresp,
  y = expr_young_noNorm,
  by = c("participant_id", "study_time_collected")
)

hipc_merged_all_norm = full_join(
  x = merged_hipc_clinical_immresp,
  y = expr_all_norm,
  by = c("participant_id", "study_time_collected")
)

# For participants with multiple baseline measurements, take only the most recent. 
# Rename all baseline sample timepoints to time = 0
# Work on a copy

hipc_merged_young_noNorm <- hipc_merged_young_noNorm %>%
  # ensure times are numeric (if they are character/factor)
  mutate(study_time_collected = as.numeric(as.character(study_time_collected)),
         .orig_row = row_number()) %>%
  
  # separate baseline rows and non-baseline rows
  {
    d <- .
    baseline_keep <- d %>%
      filter(!is.na(study_time_collected) &
               study_time_collected <= 0) %>%
      # pick the single "most recent" baseline per participant
      group_by(participant_id) %>%
      slice_max(order_by = study_time_collected,
                n = 1,
                with_ties = FALSE) %>%
      ungroup() %>%
      # if the selected baseline is negative, rename it to 0 (only for the kept row)
      mutate(study_time_collected = ifelse(study_time_collected < 0, 0, study_time_collected))
    
    non_baseline <- d %>%
      # keep all positive-time rows and rows with NA times (these are not considered baseline here)
      filter(is.na(study_time_collected) | study_time_collected > 0)
    
    # combine back
    bind_rows(non_baseline, baseline_keep) %>%
      # optional: restore original ordering (or change as desired)
      arrange(participant_id, .orig_row) %>%
      select(-.orig_row)
  }

hipc_merged_all_norm <- hipc_merged_all_norm %>%
  # ensure times are numeric (if they are character/factor)
  mutate(study_time_collected = as.numeric(as.character(study_time_collected)),
         .orig_row = row_number()) %>%
  
  # separate baseline rows and non-baseline rows
  {
    d <- .
    baseline_keep <- d %>%
      filter(!is.na(study_time_collected) &
               study_time_collected <= 0) %>%
      # pick the single "most recent" baseline per participant
      group_by(participant_id) %>%
      slice_max(order_by = study_time_collected,
                n = 1,
                with_ties = FALSE) %>%
      ungroup() %>%
      # if the selected baseline is negative, rename it to 0 (only for the kept row)
      mutate(study_time_collected = ifelse(study_time_collected < 0, 0, study_time_collected))
    
    non_baseline <- d %>%
      # keep all positive-time rows and rows with NA times (these are not considered baseline here)
      filter(is.na(study_time_collected) | study_time_collected > 0)
    
    # combine back
    bind_rows(non_baseline, baseline_keep) %>%
      # optional: restore original ordering (or change as desired)
      arrange(participant_id, .orig_row) %>%
      select(-.orig_row)
  }

# Save the merged dataframes

# Specify folder within folder root where the processed data lives
processed_data_folder = "data"

# Use fs::path() to specify the data path robustly
p_save_all_norm <- fs::path(processed_data_folder, "hipc_merged_all_norm.rds")
p_save_young_noNorm <- fs::path(processed_data_folder, "hipc_merged_young_noNorm.rds")

# Save dataframe
saveRDS(hipc_merged_all_norm, file = p_save_all_norm)
saveRDS(hipc_merged_young_noNorm, file = p_save_young_noNorm)

rm(list = ls())
