# R script to produce figures showing specification heatmaps summarising DGSA analysis results specification at the level of modules

# Libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(ggtext)
library(ggh4x)  
library(stringr)
library(scales)
library(patchwork)
library(tibble)
library(ggforce)   # for facet_wrap2
library(stringr)

p_dearseq_dgsa_results_processed = fs::path("output",
                                            "results",
                                            "dearseq",
                                            "dearseq_dgsa_results_processed.rds")
p_qusage_dgsa_results_processed = fs::path("output",
                                           "results",
                                           "qusage",
                                           "qusage_dgsa_results_processed.rds")

results_dearseq = readRDS(p_dearseq_dgsa_results_processed)
results_qusage = readRDS(p_qusage_dgsa_results_processed)

results_df = bind_rows(results_dearseq, results_qusage)

# results_df$condition_short <- recode(
#   results_df$condition,
#   "Ebola (RVV)"           = "Ebola (RVV)",
#   "Yellow Fever (LV)"     = "Y.F. (LV)",
#   "Smallpox (LV)"         = "Smallpox (LV)",
#   "Tuberculosis (RVV)"    = "T.B. (RVV)",
#   "Hepatitis A/B (IN/RP)" = "Hep. (IN/RP)",
#   "Meningococcus (CJ)"    = "Men. (CJ)",
#   "Meningococcus (PS)"    = "Men. (PS)",
#   "Malaria (RP)"          = "Malaria (RP)",
#   "Influenza (IN)"        = "Inf. (IN)",
#   "Influenza (LV)"        = "Inf. (LV)",
#   "Pneumococcus (PS)"     = "Pneumo. (PS)",
#   "Varicella Zoster (LV)" = "Varicella (LV)",
#   "HIV (RVV)"             = "HIV (RVV)"
# )

results_df$condition_short <- recode(
  results_df$condition,
  "Ebola (RVV)"           = "Ebola (RVV)",
  "Yellow Fever (LV)"     = "Yellow Fever (LV)",
  "Smallpox (LV)"         = "Smallpox (LV)",
  "Tuberculosis (RVV)"    = "Tuberculosis (RVV)",
  "Hepatitis A/B (IN/RP)" = "Hepatitis (IN/RP)",
  "Meningococcus (CJ)"    = "Meningococcus (CJ)",
  "Meningococcus (PS)"    = "Meningococcus (PS)",
  "Malaria (RP)"          = "Malaria (RP)",
  "Influenza (IN)"        = "Influenza (IN)",
  "Influenza (LV)"        = "Influenza (LV)",
  "Pneumococcus (PS)"     = "Pneumococcus (PS)",
  "Varicella Zoster (LV)" = "Varicella (LV)",
  "HIV (RVV)"             = "HIV (RVV)"
)

rm(results_dearseq, results_qusage)

methods = c("dearseq", "qusage")
conditions   = levels(results_df$condition_short)
aggregates   = levels(results_df$gs.aggregate)
timepoint    = 1
p_approach = c("global", "withinTime", "withinComparison")
p_correction = c("BH", "BY", "bonferroni", "hochberg", "holm", "hommel")
p_threshold = c(0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1)
fc_threshold = seq(0, 2, 0.1)
order = c("available")
legend_scale_exponent = 0.6
short_names = TRUE
strip_text_size = 11
condition_text_size = 13


plot_specification_heatmap_modules = function(
  methods = c("dearseq", "qusage"),
  conditions = NULL,
  timepoint = NULL,
  aggregates = NULL,
  p_approach = c("global", "withinTime", "withinComparison"),
  p_correction = c("BH", "BY", "bonferroni", "hochberg", "holm", "hommel"),
  p_threshold = c(0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1),
  fc_threshold = seq(0, 2, 0.1),
  order = c("available", "set"),
  legend_scale_exponent = 0.6,
  short_names = TRUE,
  strip_text_size = 17,
  condition_text_size = 13,
  topN = 16
) {
    
    # Prepare input table and rename column for clarity
    results_df2 = results_df
    if (short_names) {
      results_df2$condition = results_df2$condition_short
    }
    
    DT = as.data.table(
      results_df2 %>%
        filter(
          condition    %in% conditions,
          time         == timepoint,
          gs.aggregate %in% aggregates,
          method       %in% methods
        )
    )
    setnames(DT, "method", "analysis_method")
    rm(results_df2); gc()  # free memory
    
    # Build list of adjusted-p value column names to melt
    combos = expand.grid(
      approach = p_approach,
      method   = p_correction,
      stringsAsFactors = FALSE
    )
    adj_cols = with(combos, paste0(approach, ".adjPval_", method))
    rm(combos); gc()
    
    # Melt p-value columns to long format
    dt_long = data.table::melt(
      DT,
      id.vars      = c(
        "time",
        "gs.name.description",
        "condition",
        "comparison",
        "fc.score",
        "gs.colour",
        "gs.aggregate",
        "analysis_method"
      ),
      measure.vars = adj_cols,
      variable.name = "pval_type",
      value.name   = "adjp"
    )
    rm(adj_cols); gc()
    
    # Split pval_type into approach and correction method
    dt_long[, c("p_approach_spec", "p_method_spec") :=
              tstrsplit(pval_type, ".adjPval_", fixed = TRUE)]
    
    # Construct the grid of specifications (the multiverse)
    spec_grid = CJ(
      p_spec           = as.numeric(p_threshold),
      p_method_spec    = p_correction,
      p_approach_spec  = p_approach,
      filtration_spec  = as.numeric(fc_threshold),
      analysis_method  = methods
    )
    n_specifications = nrow(spec_grid)
    rm(p_threshold, fc_threshold); gc()
    
    # Precompute per-geneset metadata
    gs_meta = DT[, .(
      n.comparisons = uniqueN(comparison),
      gs.colour     = gs.colour[1],
      gs.aggregate  = gs.aggregate[1]
    ), by = gs.name.description]
    
    # Full cartesian grid of genesets × specifications
    full_grid = CJ(
      gs.name.description = unique(DT$gs.name.description),
      spec_id = seq_len(nrow(spec_grid))
    )[, spec_id := NULL][, cbind(.SD, spec_grid), by = gs.name.description]
    
    # Count observed DE by specification and geneset
    obs_counts = spec_grid[
      dt_long,
      on = .(p_approach_spec, p_method_spec, analysis_method),
      allow.cartesian = TRUE
    ][abs(fc.score) >= filtration_spec & adjp <= p_spec,
      .(n.DE = uniqueN(condition)),
      by = .(
        gs.name.description,
        p_spec,
        p_method_spec,
        p_approach_spec,
        filtration_spec,
        analysis_method
      )
    ]
    gc()
    
    # Join metadata and counts, compute percent.DE per geneset × spec
    setkey(gs_meta, gs.name.description)
    setkey(full_grid, gs.name.description)
    setkey(obs_counts, gs.name.description)
    
    res = full_grid[gs_meta]
    res[, n.DE := 0L]
    rm(full_grid); gc()
    
    key_cols = c(
      "gs.name.description",
      "p_spec",
      "p_method_spec",
      "p_approach_spec",
      "filtration_spec",
      "analysis_method"
    )
    setkeyv(res, key_cols)
    setkeyv(obs_counts, key_cols)
    res[obs_counts, n.DE := i.n.DE]
    rm(obs_counts); gc()
    
    specification_gs_results = res[, percent.DE := fifelse(n.DE > 0, 100L * n.DE / n.comparisons, 0L)][, .(
      gs.name.description,
      analysis_method,
      p_spec,
      p_method_spec,
      p_approach_spec,
      filtration_spec,
      n.DE,
      n.comparisons,
      percent.DE,
      gs.colour,
      gs.aggregate
    )]
    rm(res); gc()
    
    # Aggregate mean percent.DE across geneset aggregates and rank specifications
    summary_df = specification_gs_results %>%
      group_by(
        analysis_method,
        p_method_spec,
        p_approach_spec,
        filtration_spec,
        p_spec,
        gs.aggregate
      ) %>%
      summarise(percent.DE = mean(percent.DE, na.rm = TRUE), .groups = "drop")
    
    ranked_specs = summary_df %>%
      group_by(
        analysis_method,
        p_spec,
        p_method_spec,
        p_approach_spec,
        filtration_spec
      ) %>%
      summarise(avg_percent_DE = mean(percent.DE, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(avg_percent_DE), desc(filtration_spec)) %>%
      mutate(rank = row_number()) %>%
      mutate(
        spec_id = paste0(
          "p=",
          p_spec,
          "|m=",
          p_method_spec,
          "|a=",
          p_approach_spec,
          "|f=",
          filtration_spec,
          "|d=",
          analysis_method
        )
      )
    rm(summary_df); gc()
    
    # Unique conditions and geneset names
    conds    = DT[, unique(condition)]
    gs_names = DT[, unique(gs.name.description)]
    rm(DT); gc()
    
    # Grid of condition × geneset × specifications
    full_grid2 = CJ(
      condition      = conds,
      gs.name.description = gs_names,
      spec_id        = seq_len(nrow(spec_grid)),
      unique         = TRUE
    )[, spec_id := NULL][, cbind(.SD, spec_grid), by = .(condition, gs.name.description)]
    rm(conds, gs_names); gc()
    
    # Flag DE per condition × geneset × specification
    # 1) Merge dt_long with spec_grid to attach all specification parameters
    dt_long_specs <- dt_long[
      spec_grid,
      on = .(p_approach_spec, p_method_spec, analysis_method),
      allow.cartesian = TRUE
    ]
    
    # 2) Flag DE directly in dt_long
    de_flags <- dt_long_specs[
      abs(fc.score) >= filtration_spec & adjp <= p_spec,
      .(DE = 1L),
      by = .(
        condition,
        gs.name.description,
        p_spec,
        p_method_spec,
        p_approach_spec,
        filtration_spec,
        analysis_method
      )
    ]
    
    # 3) Remove duplicates (common when multiple comparisons per condition)
    de_flags <- unique(de_flags)
    
    # 4) Join with full_grid2 (all condition × gs.name × spec combinations) and fill missing DE
    setkeyv(full_grid2, c("condition", "gs.name.description", "p_spec",
                          "p_method_spec", "p_approach_spec", "filtration_spec", "analysis_method"))
    setkeyv(de_flags, key(full_grid2))
    
    full_de <- de_flags[full_grid2]      # left join
    full_de[is.na(DE), DE := 0L]         # fill missing
    
    rm(specification_gs_results); gc()
    
    # Merge flags into full grid; missing → DE = 0
    setkeyv(
      full_grid2,
      c(
        "condition",
        "gs.name.description",
        "p_spec",
        "p_method_spec",
        "p_approach_spec",
        "filtration_spec",
        "analysis_method"
      )
    )
    setkeyv(de_flags, key(full_grid2))
    
    full_de = de_flags[full_grid2]
    full_de[is.na(DE), DE := 0L]
    rm(de_flags, full_grid2); gc()
    
    # Add aggregate labels to the geneset-level table
    agg_map = gs_meta[, .(gs.name.description, gs.aggregate)]  # use gs_meta, not spec_grid
    setkeyv(agg_map, "gs.name.description")
    setkeyv(full_de, "gs.name.description")
    full_de = agg_map[full_de]  # join to add gs.aggregate
    rm(agg_map, spec_grid); gc()
  

  # Create geneset-level summary table with spec identifiers and ranks
  gs_de_summary = full_de %>%
    mutate(
      spec_id = paste0(
        "p=",
        p_spec,
        "|m=",
        p_method_spec,
        "|a=",
        p_approach_spec,
        "|f=",
        filtration_spec,
        "|d=",
        analysis_method
      )
    ) %>%
    left_join(ranked_specs %>% select(spec_id, rank), by = "spec_id") %>%
    mutate(spec_id = factor(spec_id, levels = ranked_specs$spec_id[order(ranked_specs$rank)])) %>%
    left_join(as.data.frame(gs_meta) %>% select(gs.name.description, gs.colour, gs.aggregate),
              by = "gs.name.description") %>%
    select(
      condition,
      p_spec,
      p_method_spec,
      p_approach_spec,
      filtration_spec,
      analysis_method,
      gs.name.description,
      DE,
      spec_id,
      rank,
      gs.colour
    )
  gs_de_summary$DE = as.integer(gs_de_summary$DE)

  # Determine ordering of conditions for plotting
  if (order == "available") {
    condition_order = gs_de_summary %>%
      group_by(condition) %>%
      summarise(mean_prop_DE = mean(DE, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_prop_DE)) %>%
      pull(condition)
  } else {
    condition_order = conditions
  }

  gs_de_summary = gs_de_summary %>%
    mutate(condition = factor(condition, levels = rev(condition_order)))

  # Total number of unique specifications
  total_specs <- gs_de_summary %>%
    pull(spec_id) %>%
    unique() %>%
    length()

  # Compute robustness per geneset × condition
  geneset_robustness_df <- gs_de_summary %>%
    group_by(condition, gs.name.description, gs.colour) %>%
    summarise(
      n_specifications = n(),
      n_DE             = sum(DE, na.rm = TRUE),
      robustness       = ifelse(n_specifications > 0, n_DE / n_specifications, NA_real_),
      robustness_pct   = 100 * robustness,
      .groups = "drop"
    ) %>%
    mutate(expected_total_specs = total_specs) %>%
    arrange(condition, desc(robustness))

  # Summarise robustness across conditions for each geneset
  geneset_robustness_across_conditions_df <- geneset_robustness_df %>%
    group_by(gs.name.description, gs.colour) %>%
    summarise(
      mean_robustness = mean(robustness, na.rm = TRUE),
      sd_robustness   = sd(robustness, na.rm = TRUE),
      min_robustness  = min(robustness, na.rm = TRUE),
      max_robustness  = max(robustness, na.rm = TRUE),
      n_conditions    = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_robustness))

  # Select top N genesets by mean robustness
  N = 16
  top_N_genesets <- geneset_robustness_across_conditions_df %>%
    arrange(desc(mean_robustness)) %>%
    slice(1:N) %>%
    pull(gs.name.description)

  # Prepare facet strip colours for top genesets
  strip_colours = gs_de_summary %>%
    filter(gs.name.description %in% top_N_genesets) %>%
    distinct(gs.name.description, gs.colour) %>%
    arrange(gs.name.description) %>%
    deframe()

  facet_levels = top_N_genesets

  # Build strip theme elements from colours
  strip_background_elements = lapply(facet_levels, function(f) element_rect(fill = strip_colours[[f]]))
  strip_text_elements       = lapply(facet_levels, function(f) element_text(face = "bold", color = "white"))

  # Plotting grid for top genesets
  full_plot_grid = CJ(
    condition = condition_order,
    spec_id   = levels(gs_de_summary$spec_id),
    gs.name.description = facet_levels,
    unique    = TRUE
  )

  # Merge geneset-level DE flags into plotting grid
  heatmap_df = full_plot_grid %>%
    left_join(
      gs_de_summary %>% select(condition, spec_id, gs.name.description, DE),
      by = c("condition", "spec_id", "gs.name.description")
    ) %>%
    mutate(
      condition = factor(condition, levels = rev(condition_order)),
      spec_id   = factor(spec_id, levels = rev(levels(gs_de_summary$spec_id))),
      gs.name.description = factor(gs.name.description, levels = facet_levels)
    )

  # Identify conditions with all NA values
  na_conditions = heatmap_df %>%
    group_by(condition) %>%
    summarize(all_na = all(is.na(DE)), .groups = "drop") %>%
    filter(all_na) %>%
    pull(condition)

  label_map = setNames(
    vapply(levels(heatmap_df$condition), function(lvl) {
      if (lvl %in% na_conditions) {
        lvl
      } else {
        lvl
      }
    }, FUN.VALUE = character(1)),
    levels(heatmap_df$condition)
  )

  color_transform <- function(x, exp = legend_scale_exponent) x^exp

  # Build the heatmap: binary fill (DE vs Not DE) with legend outlines
  heatmap_plot <- ggplot(
    heatmap_df,
    aes(x = spec_id, y = condition, fill = factor(DE))
  ) +
    geom_tile() +
    scale_y_discrete(labels = label_map) +
    scale_fill_manual(
      values = c("0" = "#FFFFFF", "1" = "#3FAF00"),
      na.value = "#FFFFFF",
      name = "Diff. Expressed",
      labels = c("0" = "Not DE", "1" = "DE"),
      guide = guide_legend(
        title.hjust = 0.5,
        label.theme = element_text(size = 16),
        override.aes = list(color = "grey50")
      )
    ) +
    labs(
      x     = "Specifications",
      y     = "Vaccines",
      title = paste0("Day ", timepoint)
    ) +
    facet_wrap2(
      ~ gs.name.description,
      scales = "free_y",
      axes   = "y",
      strip  = strip_themed(background_x = strip_background_elements, text_x = strip_text_elements)
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 28, face = "bold"),
      axis.title   = element_text(size = 30),
      strip.text.x = element_text(
        size = strip_text_size,
        face = "bold",
        color = "black"
      ),
      axis.text.x  = element_blank(),
      axis.text.y = ggtext::element_markdown(size = condition_text_size)
    )

  heatmap_plot

  return(heatmap_plot)
}



# ---------------------------
# Inputs
# ---------------------------
p1 = plot_specification_heatmap_modules(
  methods      = c("dearseq", "qusage"),
  conditions   = levels(results_df$condition_short),
  aggregates   = levels(results_df$gs.aggregate)[-16],
  timepoint    = 1,
  #p_approach   = c("global", "withinTime", "withinComparison"),
  p_approach   = c("withinTime"),
  p_correction = c("BH", "BY", "bonferroni", "holm", "hommel", "hochberg"),
  p_threshold  = c(0.0001, 0.001, 0.01, 0.05, 0.1),
  fc_threshold = seq(0, 1, 0.5),
  order = "available",
  legend_scale_exponent = 0.5,
  short_names = T,
  strip_text_size = 17,
  condition_text_size = 18
)

p2 = plot_specification_heatmap_modules(
  methods      = c("dearseq", "qusage"),
  conditions   = levels(results_df$condition_short),
  aggregates   = levels(results_df$gs.aggregate)[-16],
  timepoint    = 3,
  # p_approach   = c("global", "withinTime", "withinComparison"),
  p_approach   = c("withinTime"),
  p_correction = c("BH", "BY", "bonferroni", "holm", "hommel", "hochberg"),
  p_threshold  = c(0.0001, 0.001, 0.01, 0.05, 0.1),
  fc_threshold = seq(0, 1, 0.5),
  order = "available",
  legend_scale_exponent = 0.5,
  short_names = T,
  strip_text_size = 17,
  condition_text_size = 18
)

p3 = plot_specification_heatmap_modules(
  methods      = c("dearseq", "qusage"),
  conditions   = levels(results_df$condition_short),
  aggregates   = levels(results_df$gs.aggregate)[-16],
  timepoint    = 7,
  # p_approach   = c("global", "withinTime", "withinComparison"),
  p_approach   = c("withinTime"),
  p_correction = c("BH", "BY", "bonferroni", "holm", "hommel", "hochberg"),
  p_threshold  = c(0.0001, 0.001, 0.01, 0.05, 0.1),
  fc_threshold = seq(0, 1, 0.5),
  order = "available",
  legend_scale_exponent = 0.5,
  short_names = T,
  strip_text_size = 17,
  condition_text_size = 18
)

p1 <- p1 +
  labs(title = "a) Day 1") +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 50,
      hjust = -0.08
    ),
    axis.title = element_text(face = "plain", size = 50),
    plot.margin = margin(
      t = 60,
      r = 10,
      b = 40,
      l = 20
    )  # increased top margin
  ) 

p2 <- p2 +
  labs(title = "Day 3") +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 50,
      hjust = -0.08
    ),
    axis.title = element_text(face = "plain", size = 50),
    plot.margin = margin(
      t = 40,
      r = 10,
      b = 40,
      l = 20
    )
  )

p3 <- p3 +
  labs(title = "b) Day 7") +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 50,
      hjust = -0.08
    ),
    axis.title = element_text(face = "plain", size = 50),
    plot.margin = margin(
      t = 40,
      r = 10,
      b = 20,
      l = 20
    )
  )

# Combine vertically with shared legend and main title
combined <- (p1 / p2/ p3) +
  plot_layout(ncol = 1, heights = c(9/12, 1, 14/12), guides = "collect") +
  plot_annotation(
    title = "Specification Heatmaps Across Time",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 55,
        hjust = 0.5,
        margin = margin(b = 30, l = 20)
      )
    )
  ) &
  theme(
    legend.position = "right",
  )

# Display
print(combined)


figures_folder = fs::path("output", "figures", "dgsa")

# Save
ggsave(
  fs::path(
    figures_folder,
    "specification_heatmap_module_combined.pdf"
  ),
  combined,
  width = 35,
  height = 45,
  dpi = 300
)

ggsave(
  fs::path(
    figures_folder,
    "specification_heatmap_day1.pdf"
  ),
  p1,
  width = 22,
  height = 12,
  dpi = 300
)

ggsave(
  fs::path(
    figures_folder,
    "specification_heatmap_day3.pdf"
  ),
  p2,
  width = 33,
  height = 17,
  dpi = 300
)

ggsave(
  fs::path(
    figures_folder,
    "specification_heatmap_day7.pdf"
  ),
  p3,
  width = 22,
  height = 12,
  dpi = 300
)

