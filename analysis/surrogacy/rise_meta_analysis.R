# R script to perform meta-analysis with rise

library(tidyverse)
library(SurrogateRank)
library(grid)    # for unit()
library(scales)  # pretty formatting

tp = 7

# Directory containing engineered / processed data files
processed_data_folder <- "data"

# Paths to processed gene-level data and gene-set objects
p_load_expr_all_norm <- fs::path(processed_data_folder, "hipc_merged_all_norm.rds")

# Load data objects
hipc_merged_all_norm <- readRDS(p_load_expr_all_norm)

# Extract the gene names
gene_names = hipc_merged_all_norm %>%
  dplyr::select(a1cf:zzz3) %>%
  colnames()

# Timepoints of interest (numeric)
timepoints_to_keep <- c(0,tp)

# Filter to samples with non-missing immune response, Influenza vaccine,
# and collected at one of the specified timepoints.
hipc_merged_all_norm_filtered <- hipc_merged_all_norm %>%
  filter(
    !is.na(immResp_MFC_anyAssay_log2_MFC),
    vaccine_name == "Influenza (IN)",
    study_time_collected %in% timepoints_to_keep
  ) %>% 
  select(participant_id, 
         study_accession,
         study_time_collected,
         immResp_MFC_anyAssay_pre_value,
         immResp_MFC_anyAssay_post_value,
         all_of(gene_names))

# Extract the study names
study_names = hipc_merged_all_norm_filtered %>%
  filter(study_time_collected > 0) %>%
  pull(study_accession) %>%
  unique()

timepoints_of_interest = c(tp)

epsilon = 0.2

for (sdy in study_names) {
  
  hipc_merged_all_norm_filtered_sdy = hipc_merged_all_norm_filtered %>%
    filter(study_time_collected %in% c(0, tp),
           study_accession == sdy) %>%
    group_by(participant_id) %>%
    filter(all(c(0, tp) %in% study_time_collected)) %>%
    ungroup()
  
  # Pre-vaccination response
  yzero_sdy = hipc_merged_all_norm_filtered_sdy %>%
    filter(study_time_collected == 0) %>%
    pull(immResp_MFC_anyAssay_pre_value)
  
  # Post-vaccination response
  yone_sdy = hipc_merged_all_norm_filtered_sdy %>%
    filter(study_time_collected == tp) %>%
    pull(immResp_MFC_anyAssay_post_value)
  
  # Pre-vaccination gene expression
  szero_sdy = hipc_merged_all_norm_filtered_sdy %>%
    filter(study_time_collected == 0) %>%
    select(all_of(gene_names))
  
  # Post-vaccination gene expression
  sone_sdy = hipc_merged_all_norm_filtered_sdy %>%
    filter(study_time_collected == tp) %>%
    select(all_of(gene_names))
  
  # Extract the estimated effect on the response on the relevant samples
  yresult_sdy_main = SurrogateRank::test.surrogate.extension(
    yone = yone_sdy,
    yzero = yzero_sdy,
    sone = yone_sdy,
    szero = yzero_sdy,
    power.want.s = 0.8,
    paired = T,
    alternative = "two.sided"
  )[["u.y"]]
  
  test_result_sdy = SurrogateRank::rise.screen(
    yone = yone_sdy,
    yzero = yzero_sdy,
    sone = sone_sdy,
    szero = szero_sdy,
    epsilon = epsilon,
    alternative = "two.sided",
    paired = T,alpha = 0.05,
    p.correction = "none"
  )
  
  
  test_result_df = test_result_sdy[["screening.metrics"]]
  test_result_df$study_accession = sdy
  test_result_df$n = 2*length(yone_sdy)
  test_result_df$u.y = yresult_sdy_main
  test_result_df$u.s = yresult_sdy_main - test_result_df$delta
  
  test_result_df = test_result_df %>%
    select(
      marker,
      study_accession,
      n,
      epsilon,
      u.y,
      u.s,
      delta,
      ci_upper,
      ci_lower,
      sd,
      p_unadjusted
    )
  
  if (sdy == study_names[1]) {
    all_results = test_result_df
  } else {
    all_results = rbind(all_results, test_result_df)
  }
}

# ---- compute per-marker fixed-effect pooled results -----------------------
summary_df <- all_results %>%
  # keep only required columns (tolerant to extra cols)
  transmute(marker,
            study_accession,
            n,
            delta,
            ci_upper,
            ci_lower,
            sd,
            p_unadjusted) %>%
  
  # group by marker and compute pooled FE; exclude rows with sd == 0 or NA for pooling
  group_by(marker) %>%
  group_modify(~{
    dat <- .x
    
    # exclude studies with sd == 0 or NA for the pooling step
    dat_use <- dat %>% filter(!is.na(sd) & sd > 0)
    
    if(nrow(dat_use) == 0) {
      # no usable studies -> return NA row
      return(tibble(
        pooled_delta = NA_real_,
        se_pooled = NA_real_,
        var_pooled = NA_real_,
        ci_pooled_lower = NA_real_,
        ci_pooled_upper = NA_real_,
        p_two_sided = NA_real_,
        p1_pooled = NA_real_,
        p2_pooled = NA_real_,
        p_tost_pooled = NA_real_,
        k_used = 0L,
        total_N = sum(dat$n, na.rm = TRUE),
        W_sum = NA_real_
      ))
    }
    
    # study-level variances and weights
    dat_use <- dat_use %>%
      mutate(var = sd^2,
             weight = 1 / var)
    
    W_sum <- sum(dat_use$weight, na.rm = TRUE)
    
    # pooled fixed-effect estimate (inverse-variance)
    pooled_delta <- sum(dat_use$weight * dat_use$delta, na.rm = TRUE) / W_sum
    var_pooled <- 1 / W_sum
    se_pooled <- sqrt(var_pooled)
    
    # 95% CI (normal approximation)
    ci_pooled_lower <- pooled_delta - qnorm(0.975) * se_pooled
    ci_pooled_upper <- pooled_delta + qnorm(0.975) * se_pooled
    
    # two-sided test vs 0 (for reporting)
    p_two_sided <- 2 * (1 - pnorm(abs(pooled_delta / se_pooled)))
    
    # pooled one-sided p-values for TOST (normal approx)
    Z1 <- (pooled_delta - epsilon) / se_pooled   # test delta < epsilon
    Z2 <- (pooled_delta + epsilon) / se_pooled   # test delta > -epsilon
    
    p1_pooled <- pnorm(Z1)           # one-sided p for delta < epsilon
    p2_pooled <- 1 - pnorm(Z2)       # one-sided p for delta > -epsilon
    
    p_tost_pooled <- max(p1_pooled, p2_pooled)  # TOST pooled p-value
    
    # return a single-row tibble for this marker
    tibble(
      pooled_delta = pooled_delta,
      se_pooled = se_pooled,
      var_pooled = var_pooled,
      ci_pooled_lower = ci_pooled_lower,
      ci_pooled_upper = ci_pooled_upper,
      p_two_sided = p_two_sided,
      p1_pooled = p1_pooled,
      p2_pooled = p2_pooled,
      p_tost_pooled = p_tost_pooled,
      k_used = nrow(dat_use),
      total_N = sum(dat$n, na.rm = TRUE),
      W_sum = W_sum
    )
  }, .keep = TRUE) %>%
  ungroup()

# ---- optional: reorder columns and view -----------------------------------
summary_df <- summary_df %>%
  mutate(p_adjusted = p.adjust(p_tost_pooled, method = "BH")) %>% 
  select(marker,
         k_used,
         total_N,
         W_sum,
         pooled_delta, se_pooled, var_pooled,
         ci_pooled_lower, ci_pooled_upper,
         p_two_sided,
         p1_pooled, p2_pooled, p_tost_pooled, p_adjusted)

signature_genes = summary_df %>% 
  filter(p_adjusted < 0.05) %>% 
  pull(marker)

all_results_filtered = all_results %>% 
  filter(marker == "jchain")

m.gen = metagen(
  data = all_results_filtered,
  TE = delta,
  seTE = sd,
  studlab = study_accession,
  sm = "RD",
  common = TRUE,
  random = TRUE,
  method.tau = "REML",
  method.random.ci = "HK"
)

# signature_genes = summary_df %>% 
#   arrange(p_tost_pooled) %>% 
#   slice_head(n = 10) %>% 
#   pull(marker)

weight_df = summary_df %>% 
  filter(marker %in% signature_genes) %>% 
  mutate(weight = 1/abs(pooled_delta)) %>% 
  select(marker, weight)

for (sdy in study_names) {
  hipc_merged_all_norm_filtered_sdy_eval = hipc_merged_all_norm_filtered %>%
    filter(study_time_collected %in% c(0, tp), study_accession == sdy) %>%
    group_by(participant_id) %>%
    filter(all(c(0, tp) %in% study_time_collected)) %>%
    ungroup()
  
  # Pre-vaccination response
  yzero_sdy = hipc_merged_all_norm_filtered_sdy_eval %>%
    filter(study_time_collected == 0) %>%
    pull(immResp_MFC_anyAssay_pre_value)
  
  # Post-vaccination response
  yone_sdy = hipc_merged_all_norm_filtered_sdy_eval %>%
    filter(study_time_collected == tp) %>%
    pull(immResp_MFC_anyAssay_post_value)
  
  # Pre-vaccination gene expression
  szero_sdy = hipc_merged_all_norm_filtered_sdy_eval %>%
    filter(study_time_collected == 0) %>%
    select(all_of(signature_genes))
  
  # Post-vaccination gene expression
  sone_sdy = hipc_merged_all_norm_filtered_sdy_eval %>%
    filter(study_time_collected == tp) %>%
    select(all_of(signature_genes))
  
  # Extract the estimated effect on the response on the relevant samples
  yresult_sdy_main = SurrogateRank::test.surrogate.extension(
    yone = yone_sdy,
    yzero = yzero_sdy,
    sone = yone_sdy,
    szero = yzero_sdy,
    power.want.s = 0.8,
    paired = T,
    alternative = "two.sided"
  )[["u.y"]]
  
  evaluate_result_sdy = SurrogateRank::rise.evaluate(
    yone = yone_sdy,
    yzero = yzero_sdy,
    sone = sone_sdy,
    szero = szero_sdy,
    epsilon = epsilon,
    alternative = "two.sided",
    paired = T,
    alpha = 0.05,
    p.correction = "none",
    evaluate.weights = T,
    screening.weights = weight_df,
    markers = signature_genes
  )
  
  
  evaluate_result_df = as.data.frame(as.list(evaluate_result_sdy[["gamma.s.evaluate"]]))
  evaluate_result_df$study_accession = sdy
  evaluate_result_df$n = 2 * length(yone_sdy)
  evaluate_result_df$u.y = yresult_sdy_main
  evaluate_result_df$u.s = yresult_sdy_main - evaluate_result_df$delta
  
  evaluate_result_df = evaluate_result_df %>%
    select(study_accession,
           n,
           epsilon,
           u.y,
           u.s,
           delta,
           ci_upper,
           ci_lower,
           sd,
           p_unadjusted)
  
  if (sdy == study_names[1]) {
    all_results_eval = evaluate_result_df
  } else {
    all_results_eval = rbind(all_results_eval, evaluate_result_df)
  }
}

# ----------------- 0. Data prep and pooled estimates ------------------------
# Use your all_results dataframe
df <- all_results_eval %>%
  filter(sd != 0)

# Ensure epsilon is single
eps_vals <- unique(df$epsilon)

# Recover sd from symmetric CI if needed
df <- df %>%
  mutate(
    sd = ifelse(is.na(sd) & !is.na(ci_upper) & !is.na(delta),
                      (ci_upper - delta) / qnorm(0.975),
                      sd),
    sd = ifelse(sd == 0, NA, sd) # avoid zero sd
  )

if(any(is.na(df$sd))) {
  stop("sd missing for some studies and cannot be recovered from CIs. Please supply sd or symmetric CIs for all studies.")
}

# compute study-level variance and weight
df <- df %>%
  mutate(var = sd^2, weight = 1 / var)

# fixed-effect pooled estimate
W_sum <- sum(df$weight, na.rm = TRUE)
pooled_delta <- sum(df$weight * df$delta, na.rm = TRUE) / W_sum
var_pooled <- 1 / W_sum
se_pooled <- sqrt(var_pooled)
ci_pooled_lower <- pooled_delta - qnorm(0.975) * se_pooled
ci_pooled_upper <- pooled_delta + qnorm(0.975) * se_pooled
p_two_sided <- 2 * (1 - pnorm(abs(pooled_delta / se_pooled)))   # test vs 0 (kept for reporting)

# pooled TOST p-value (match per-study TOST = max of the two one-sided p-values)
Z1_pooled <- (pooled_delta - epsilon) / se_pooled
Z2_pooled <- (pooled_delta + epsilon) / se_pooled
p1_pooled <- pnorm(Z1_pooled)        # p for delta < epsilon
p2_pooled <- 1 - pnorm(Z2_pooled)    # p for delta > -epsilon
p_tost_pooled <- max(p1_pooled, p2_pooled)

# console summary
cat("\nPooled (FE) summary:\n")
cat(sprintf(" pooled delta = %0.6f\n", pooled_delta))
cat(sprintf(" SE pooled    = %0.6f\n", se_pooled))
cat(sprintf(" 95%% CI       = [%0.6f, %0.6f]\n", ci_pooled_lower, ci_pooled_upper))
cat(sprintf(" two-sided p  = %0.6g\n", p_two_sided))
cat(sprintf(" pooled TOST p = %0.6g (p1=%0.6g, p2=%0.6g)  (epsilon=%0.6g)\n\n",
            p_tost_pooled, p1_pooled, p2_pooled, epsilon))

# ------ prepare per-study labels and numeric columns ------------------------
df <- df %>%
  mutate(
    pct_weight = 100 * weight / W_sum,
    label_pval = formatC(p_unadjusted, format = "f", digits = 3),     # per-study p_unadjusted (assumed TOST)
    label_wgt = formatC(pct_weight, format = "f", digits = 1),
    label_n = ifelse(is.na(n), "", as.character(n))
  )

# prepare order and y positions
k <- nrow(df)
df <- df %>% mutate(study_order = row_number(), y = k - study_order + 1)

# Summary row (bottom). If you prefer total N, replace label_n = "" with sum(df$n, na.rm = TRUE)
summary_row <- data.frame(
  study_accession = "Summary",
  n = NA,
  delta = pooled_delta,
  ci_lower = ci_pooled_lower,
  ci_upper = ci_pooled_upper,
  sd = se_pooled,
  p_unadjusted = p_tost_pooled,   # store pooled TOST p-value
  epsilon = epsilon,
  var = var_pooled,
  weight = NA,
  pct_weight = 100,
  label_pval = formatC(p_tost_pooled, format = "f", digits = 3),
  label_wgt = formatC(100, format = "f", digits = 1),
  label_n = "",   # or as.character(sum(df$n, na.rm=TRUE))
  study_order = k + 1,
  y = 0,
  stringsAsFactors = FALSE
)

plot_df <- bind_rows(df, summary_row) %>%
  mutate(study_label = study_accession,
         is_summary = study_label == "Summary")

# ----- plotting parameters --------------------------------------------------
base_text_size <- 14
# vertical extents: allow space below y=0 for epsilon labels
y_min <- -1
y_max <- k + 1

# left/middle/right relative widths (tune as needed)
rel_w_left  <- 0.45
rel_w_mid   <- 0.95
rel_w_right <- 0.50

# ----------------- Left panel: study labels ---------------------------------
left_labels <- ggplot(plot_df, aes(y = y)) +
  geom_text(aes(x = 0, label = study_label), hjust = 0, size = 4.2) +
  # fix a narrow x-range so text cannot overflow into the middle panel
  scale_x_continuous(limits = c(-0.1, 0.95), expand = c(0, 0)) +
  scale_y_continuous(breaks = plot_df$y, limits = c(y_min, y_max), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = unit(c(1.2, 0.6, 1, 1.2), "lines")  # extra right padding
  )

# ----------------- Middle panel: forest plot (strict x-limits) --------------
forest_mid <- ggplot(plot_df, aes(x = delta, y = y)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.15, size = 0.6) +
  geom_point(data = filter(plot_df, !is_summary), shape = 16, size = 3) +
  geom_point(data = filter(plot_df, is_summary), shape = 18, size = 5) +
  scale_y_continuous(breaks = plot_df$y, limits = c(y_min, y_max), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-1.05, 1.05), ylim = c(y_min, y_max), expand = FALSE, clip = "off") +
  labs(x = expression(delta), y = NULL, title = paste0("Fixed-effect inverse-variance meta-analysis")) +
  theme_minimal(base_size = base_text_size) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = base_text_size + 4),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = base_text_size + 10),
    axis.text.x = element_text(size = base_text_size - 1),
    # add left margin so left panel never intrudes
    plot.margin = unit(c(1.2, 0.6, 1, 1.0), "lines")
  ) +
  geom_vline(xintercept = c(-epsilon, epsilon), linetype = "dashed", color = "red", size = 0.6) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey80") 

# ----------------- Right panel: p-value / weight / N ------------------------
# use dummy x positions 1,2,3 for the three columns; align y with middle and left
right_table <- ggplot(plot_df, aes(y = y)) +
  # column headers
  annotate("text", x = 1, y = k + 0.9, label = "p-value", fontface = "bold", hjust = 0, size = 4.5) +
  annotate("text", x = 2, y = k + 0.9, label = "Weight", fontface = "bold", hjust = 0, size = 4.5) +
  annotate("text", x = 3, y = k + 0.9, label = "N", fontface = "bold", hjust = 0, size = 4.5) +
  geom_text(aes(x = 1, label = label_pval), hjust = 0, size = 4) +
  geom_text(aes(x = 2, label = paste0(label_wgt, "%")), hjust = 0, size = 4) +
  geom_text(aes(x = 3, label = label_n), hjust = 0, size = 4) +
  scale_y_continuous(breaks = plot_df$y, limits = c(y_min, y_max), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0.8, 3.6), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = unit(c(1.2, 1.2, 1, 0.6), "lines"))

# ----------------- Combine panels -------------------------------------------
combined <- plot_grid(left_labels, forest_mid, right_table, nrow = 1,
                      rel_widths = c(rel_w_left, rel_w_mid, rel_w_right), align = "h")

# show on screen
print(combined)



