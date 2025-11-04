# In this script, we extract the variable importance results for selected prediction tasks and extract the variable names
# in order to perform an interpretation

# Libraries
library(clipr)

# Load the results
# Sequential approach, with clinical variables, without TBA
p_prediction_results_all_sequential_withClinical_withoutTBA = fs::path(
  "output",
  "results",
  "prediction",
  "prediction_results_all_sequential_withClinical_withoutTBA.rds"
)

prediction_results_all_sequential_withClinical_withoutTBA = readRDS(p_prediction_results_all_sequential_withClinical_withoutTBA)

vaccine = "Influenza (IN)"
feature_sets = c("Day 0", "Day 1", "Day 3", "Day 7", "Day 14")

feature_sets_vaccine = list(
  "Yellow Fever (LV)" = c("Day 0", "Day 3", "Day 7", "Day 14"),
  "Influenza (IN)" = c("Day 0", "Day 1", "Day 3", "Day 7", "Day 14"),
  "Meningococcus (CJ)" = c("Day 0", "Day 7"),
  "Meningococcus (PS)" = c("Day 0", "Day 7")
)

vaccine_names = names(feature_sets_vaccine)

varimp_by_vaccine <- setNames(vector("list", length(vaccine_names)), vaccine_names)

for (v in vaccine_names) {
  # feature set identifiers for this vaccine:
  fs_obj <- feature_sets_vaccine[[v]]
  # handle either a named list or a character vector
  fs_names <- as.character(fs_obj)
  
  dfs <- list()
  for (fs in fs_names) {
    # safe extraction of the var_imp object
    varimp_obj <- prediction_results_all_sequential_withClinical_withoutTBA[[v]][[fs]][["var_imp"]]
    
    
    # store using the feature-set name as the list name so bind_rows creates a column
    dfs[[fs]] <- varimp_obj
  }
  
  # bind and add a column "feature_set" with the list names
  varimp_by_vaccine[[v]] <- bind_rows(dfs, .id = "feature_set")
}

get_top_vars_string <- function(vaccine_name,
                                n_top = 10,
                                varimp_list = varimp_by_vaccine) {
  # extract the combined var_imp dataframe for the vaccine
  df <- varimp_list[[vaccine_name]]
  
  # order by mean_rank ascending and pick top n
  top_vars <- df %>%
    filter(feature_group != "clinical", mean_rank <= n_top) %>%
    pull(feature)
  
  # collapse into a comma-separated string
  paste(top_vars, collapse = "; ")
}

# Example usage:
top20_vars_YellowFeverLV = get_top_vars_string("Yellow Fever (LV)", n_top = 20)
writeLines(
  top20_vars_YellowFeverLV ,
  fs::path(
    "output",
    "results",
    "prediction",
    "top20_vars_YellowFeverLV.txt"
  )
)

top20_vars_InfluenzaIN = get_top_vars_string("Influenza (IN)", n_top = 20)
writeLines(
  top20_vars_InfluenzaIN ,
  fs::path(
    "output",
    "results",
    "prediction",
    "top20_vars_InfluenzaIN.txt"
  )
)

top20_vars_MeningococcusCJ = get_top_vars_string("Meningococcus (CJ)", n_top = 20)
writeLines(
  top20_vars_MeningococcusCJ ,
  fs::path(
    "output",
    "results",
    "prediction",
    "top20_vars_MeningococcusCJ.txt"
  )
)

top20_vars_MeningococcusPS = get_top_vars_string("Meningococcus (PS)", n_top = 20)
writeLines(
  top20_vars_MeningococcusPS ,
  fs::path(
    "output",
    "results",
    "prediction",
    "top20_vars_MeningococcusPS.txt"
  )
)