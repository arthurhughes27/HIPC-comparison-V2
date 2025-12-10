# R script to produce figures showing sensitivity heatmaps summarising DGSA analysis results sensitivity

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

results_df$condition_short <- recode(
  results_df$condition,
  "Ebola (RVV)"           = "Ebola (RVV)",
  "Yellow Fever (LV)"     = "Y.F. (LV)",
  "Smallpox (LV)"         = "Smallpox (LV)",
  "Tuberculosis (RVV)"    = "T.B. (RVV)",
  "Hepatitis A/B (IN/RP)" = "Hep. (IN/RP)",
  "Meningococcus (CJ)"    = "Men. (CJ)",
  "Meningococcus (PS)"    = "Men. (PS)",
  "Malaria (RP)"          = "Malaria (RP)",
  "Influenza (IN)"        = "Inf. (IN)",
  "Influenza (LV)"        = "Inf. (LV)",
  "Pneumococcus (PS)"     = "Pneumo. (PS)",
  "Varicella Zoster (LV)" = "Varicella (LV)",
  "HIV (RVV)"             = "HIV (RVV)"
)

rm(results_dearseq, results_qusage)


plot_sensitivity_heatmap = function(methods = c("dearseq", "qusage"),
                                    conditions = NULL,
                                    timepoint = NULL,
                                    aggregates = NULL,
                                    p_approach = c("global", "withinTime", "withinComparison"),
                                    p_correction = c("BH", "BY", "bonferroni", "hochberg", "holm", "hommel"),
                                    p_threshold = c(0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1),
                                    fc_threshold = seq(0, 2, 0.1),
                                    order = c("available", "set"),
                                    legend_scale_exponent = 0.6,
                                    short_names = TRUE) {
  
  
  # Convert results to data.table and rename 'method' to 'analysis_method'
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
  
  # Extract all desired adjusted p-value columns (e.g., "global.adjPval_BH", etc.)
  combos = expand.grid(approach = p_approach,
                        method   = p_correction,
                        stringsAsFactors = FALSE)
  adj_cols = with(combos, paste0(approach, ".adjPval_", method))
  
  # Reshape from wide to long format for each p-value approach × correction method
  dt_long = data.table::melt(
    DT,
    id.vars      = c(
      "time",
      "gs.name",
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
  dt_long[, c("p_approach_spec", "p_method_spec") :=
            tstrsplit(pval_type, ".adjPval_", fixed = TRUE)]
  
  # Create a grid of reasonable parameter specifications, including analysis_method
  spec_grid = CJ(
    p_spec           = as.numeric(p_threshold),
    p_method_spec    = p_correction,
    p_approach_spec  = p_approach,
    filtration_spec  = as.numeric(fc_threshold),
    analysis_method  = methods
  )
  n_specifications = nrow(spec_grid)
  
  # 6) Precompute per‑gene‑set metadata
  gs_meta = DT[, .(
    n.comparisons = uniqueN(comparison),
    gs.colour     = gs.colour[1],
    gs.aggregate  = gs.aggregate[1]
  ), by = gs.name]
  
  # 7) Build the full “gs.name × spec” cartesian grid
  full_grid = CJ(gs.name = unique(DT$gs.name), spec_id = seq_len(nrow(spec_grid)))[, spec_id := NULL][, cbind(.SD, spec_grid), by = gs.name]
  
  # 8) Count observed DE per spec × gs.name
  obs_counts = spec_grid[dt_long, on = .(p_approach_spec, p_method_spec, analysis_method), allow.cartesian = TRUE][abs(fc.score) >= filtration_spec &
                                                                                                                      adjp <= p_spec, .(n.DE = uniqueN(condition)), by = .(gs.name,
                                                                                                                                                                           p_spec,
                                                                                                                                                                           p_method_spec,
                                                                                                                                                                           p_approach_spec,
                                                                                                                                                                           filtration_spec,
                                                                                                                                                                           analysis_method)]
  
  
  # 9) Stitch everything together and compute percentages
  setkey(gs_meta, gs.name)
  setkey(full_grid, gs.name)
  setkey(obs_counts, gs.name)
  
  res = full_grid[gs_meta]
  res[, n.DE := 0L]
  
  # join in observed counts
  key_cols = c(
    "gs.name",
    "p_spec",
    "p_method_spec",
    "p_approach_spec",
    "filtration_spec",
    "analysis_method"
  )
  setkeyv(res, key_cols)
  setkeyv(obs_counts, key_cols)
  res[obs_counts, n.DE := i.n.DE]
  
  # finalize sensitivity_gs_results
  sensitivity_gs_results = res[, percent.DE := fifelse(n.DE > 0, 100L * n.DE / n.comparisons, 0L)][, .(
    gs.name,
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
  
  # Summarise and rank specifications
  summary_df = sensitivity_gs_results %>%
    group_by(
      analysis_method,
      p_method_spec,
      p_approach_spec,
      filtration_spec,
      p_spec,
      gs.aggregate
    ) %>%
    summarise(percent.DE = mean(percent.DE, na.rm = TRUE),
              .groups = "drop")
  
  ranked_specs = summary_df %>%
    group_by(analysis_method,
             p_spec,
             p_method_spec,
             p_approach_spec,
             filtration_spec) %>%
    summarise(avg_percent_DE = mean(percent.DE, na.rm = TRUE),
              .groups = "drop") %>%
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
  
  # 1) Unique conditions & gene‐sets
  conds    = DT[, unique(condition)]
  gs_names = DT[, unique(gs.name)]
  
  # 2) Build full grid of (condition × gs.name × specs)
  full_grid2 = CJ(
    condition      = conds,
    gs.name        = gs_names,
    spec_id        = seq_len(nrow(spec_grid)),
    unique         = TRUE
  )[, spec_id := NULL][, cbind(.SD, spec_grid), by = .(condition, gs.name)]
  
  # 3) Flag DE in the long table per spec
  de_flags = spec_grid[dt_long, on = .(p_approach_spec, p_method_spec, analysis_method), allow.cartesian = TRUE][abs(fc.score) >= filtration_spec &
                                                                                                                    adjp <= p_spec, .(DE = 1L), by = .(
                                                                                                                      condition,
                                                                                                                      gs.name,
                                                                                                                      p_spec,
                                                                                                                      p_method_spec,
                                                                                                                      p_approach_spec,
                                                                                                                      filtration_spec,
                                                                                                                      analysis_method
                                                                                                                    )]
  
  # 4) Left‑join flags into full_grid2, fill missing → DE = 0
  setkeyv(
    full_grid2,
    c(
      "condition",
      "gs.name",
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
  
  # 5) Bring in aggregate labels
  agg_map = gs_meta[, .(gs.name, gs.aggregate)]
  setkeyv(agg_map, "gs.name")
  setkeyv(full_de, "gs.name")
  full_de = agg_map[full_de]
  
  # 6) Compute proportion DE per spec × condition × aggregate
  agg_de_summary = full_de[, .(prop.DE = sum(DE) / .N), by = .(
    condition,
    p_spec,
    p_method_spec,
    p_approach_spec,
    filtration_spec,
    analysis_method,
    gs.aggregate
  )] %>%
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
    left_join(DT %>% distinct(gs.aggregate, gs.colour), by = "gs.aggregate")
  
  # Determine condition ordering
  if (order == "available") {
    condition_order = agg_de_summary %>%
      group_by(condition) %>%
      summarise(mean_prop_DE = mean(prop.DE, na.rm = TRUE),
                .groups = "drop") %>%
      arrange(desc(mean_prop_DE)) %>%
      pull(condition)
  } else {
    condition_order = conditions
  }
  
  agg_de_summary = agg_de_summary %>%
    mutate(condition = factor(condition, levels = rev(condition_order)))
  
  # Prepare facet strip colours
  strip_colours = agg_de_summary %>%
    distinct(gs.aggregate, gs.colour) %>%
    arrange(gs.aggregate) %>%
    deframe()
  facet_levels = aggregates
  strip_background_elements = lapply(facet_levels, function(f)
    element_rect(fill = strip_colours[[f]]))
  strip_text_elements       = lapply(facet_levels, function(f)
    element_text(face = "bold", color = "white"))
  
  # Build plotting grid and merge in prop.DE
  full_plot_grid = CJ(
    condition    = condition_order,
    spec_id      = levels(agg_de_summary$spec_id),
    gs.aggregate = aggregates,
    unique       = TRUE
  )
  heatmap_df = full_plot_grid %>%
    left_join(
      agg_de_summary %>% select(condition, spec_id, gs.aggregate, prop.DE),
      by = c("condition", "spec_id", "gs.aggregate")
    ) %>%
    mutate(
      condition    = factor(condition, levels = rev(condition_order)),
      spec_id      = factor(spec_id, levels = rev(levels(
        agg_de_summary$spec_id
      ))),
      gs.aggregate = factor(gs.aggregate, levels = aggregates)
    )
  
  # Identify conditions with all NA
  na_conditions = heatmap_df %>%
    group_by(condition) %>%
    summarize(all_na = all(is.na(prop.DE)), .groups = "drop") %>%
    filter(all_na) %>%
    pull(condition)
  label_map = setNames(vapply(levels(heatmap_df$condition), function(lvl) {
    if (lvl %in% na_conditions) {
      # sprintf('<span style="color:#ced4da">%s</span>', lvl)
      lvl
    } else {
      lvl
    }
  }, FUN.VALUE = character(1)),
  levels(heatmap_df$condition))
  
  color_transform <- function(x, exp = legend_scale_exponent) x^exp
  
  heatmap_plot <- ggplot(
    heatmap_df,
    aes(x = spec_id, y = condition, fill = color_transform(prop.DE))
  ) +
    geom_tile() +
    scale_y_discrete(labels = label_map) +
    scale_fill_gradient(
      low = "#FFFFFF",
      high = "#3FAF00",
      na.value = "#FFFFFF",
      limits = c(0, 1),  # legend still 0–1
      name = str_wrap("Proportion of DE genesets within aggregate", 20),
      guide = guide_colorbar(
        title.position = "top",
        title.hjust    = 0.5,
        title.theme    = element_text(size = 20, margin = margin(b = 30)),
        label.theme    = element_text(size = 16),
        barwidth       = unit(4.5, "lines"),
        barheight      = unit(10, "lines")
      ),
      breaks = seq(0, 1, 0.25),   # linear tick marks
      labels = seq(0, 1, 0.25)    # keep linear legend labels
    ) +
    labs(
      x     = "Specifications",
      y     = "Vaccines",
      title = paste0(
        "Day ",
        timepoint,
        ""
      )
    ) +
    facet_wrap2(
      ~ gs.aggregate,
      scales = "free_y",
      axes   = "y",
      strip  = strip_themed(background_x = strip_background_elements, text_x       = strip_text_elements)
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 28, face = "bold"),
      axis.title   = element_text(size = 30),
      strip.text.x = element_text(
        size = 12,
        face = "bold",
        color = "black"
      ),
      axis.text.x  = element_blank(),
      axis.text.y = ggtext::element_markdown(size = 10),
    )

  heatmap_plot
  
  return(heatmap_plot)
}


# ---------------------------
# Inputs
# ---------------------------
p1 = plot_sensitivity_heatmap(
  methods      = c("dearseq", "qusage"),
  conditions   = levels(results_df$condition_short),
  aggregates   = levels(results_df$gs.aggregate),
  timepoint    = 1,
  #p_approach   = c("global", "withinTime", "withinComparison"),
  p_approach   = c("withinTime"),
  p_correction = c("BH", "BY", "bonferroni", "holm", "hommel"),
  p_threshold  = c(0.0001, 0.001, 0.01, 0.05, 0.1),
  fc_threshold = seq(0, 1.5, 0.5),
  order = "available",
  legend_scale_exponent = 0.5
)

p2 = plot_sensitivity_heatmap(
  methods      = c("dearseq", "qusage"),
  conditions   = levels(results_df$condition_short),
  aggregates   = levels(results_df$gs.aggregate),
  timepoint    = 3,
  # p_approach   = c("global", "withinTime", "withinComparison"),
  p_approach   = c("withinTime"),
  p_correction = c("BH", "BY", "bonferroni", "holm", "hommel"),
  p_threshold  = c(0.0001, 0.001, 0.01, 0.05, 0.1),
  fc_threshold = seq(0, 1.5, 0.5),
  order = "available",
  legend_scale_exponent = 0.5
)

p3 = plot_sensitivity_heatmap(
  methods      = c("dearseq", "qusage"),
  conditions   = levels(results_df$condition_short),
  aggregates   = levels(results_df$gs.aggregate),
  timepoint    = 7,
  # p_approach   = c("global", "withinTime", "withinComparison"),
  p_approach   = c("withinTime"),
  p_correction = c("BH", "BY", "bonferroni", "holm", "hommel"),
  p_threshold  = c(0.0001, 0.001, 0.01, 0.05, 0.1),
  fc_threshold = seq(0, 1.5, 0.5),
  order = "available",
  legend_scale_exponent = 0.5
)

p1 <- p1 +
  labs(title = "a) Day 1") +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 50,
      hjust = -0.08
    ),
    axis.title.x = element_text(face = "plain", size = 35),
    plot.margin = margin(
      t = 60,
      r = 10,
      b = 40,
      l = 20
    )  # increased top margin
  ) +
  scale_fill_gradient(
    low = "#FFFFFF",
    high = "#3FAF00",
    na.value = "#FFFFFF",
    limits = c(0, 1),  # legend still 0–1
    name = str_wrap("Proportion of DE genesets within aggregate", 20),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      title.theme    = element_text(size = 35, margin = margin(b = 30, l = 30)),
      label.theme    = element_text(size = 28),
      barwidth       = unit(3.5, "lines"),
      barheight      = unit(17, "lines")
    ),
    breaks = seq(0, 1, 0.25),   # linear tick marks
    labels = seq(0, 1, 0.25)    # keep linear legend labels
  )

p2 <- p2 +
  labs(title = "b) Day 3") +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 50,
      hjust = -0.08
    ),
    axis.title.x = element_text(face = "plain", size = 35),
    plot.margin = margin(
      t = 40,
      r = 10,
      b = 40,
      l = 20
    )
  ) +
  scale_fill_gradient(
    low = "#FFFFFF",
    high = "#3FAF00",
    na.value = "#FFFFFF",
    limits = c(0, 1),  # legend still 0–1
    name = str_wrap("Proportion of DE genesets within aggregate", 20),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      title.theme    = element_text(size = 35, margin = margin(b = 30, l = 30)),
      label.theme    = element_text(size = 28),
      barwidth       = unit(3.5, "lines"),
      barheight      = unit(17, "lines")
    ),
    breaks = seq(0, 1, 0.25),   # linear tick marks
    labels = seq(0, 1, 0.25)    # keep linear legend labels
  )

p3 <- p3 +
  labs(title = "c) Day 7") +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 50,
      hjust = -0.08
    ),
    axis.title.x = element_text(face = "plain", size = 35),
    plot.margin = margin(
      t = 40,
      r = 10,
      b = 20,
      l = 20
    )
  ) +
  scale_fill_gradient(
    low = "#FFFFFF",
    high = "#3FAF00",
    na.value = "#FFFFFF",
    limits = c(0, 1),  # legend still 0–1
    name = str_wrap("Proportion of DE genesets within aggregate", 20),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      title.theme    = element_text(size = 35, margin = margin(b = 30, l = 30)),
      label.theme    = element_text(size = 28),
      barwidth       = unit(3.5, "lines"),
      barheight      = unit(17, "lines")
    ),
    breaks = seq(0, 1, 0.25),   # linear tick marks
    labels = seq(0, 1, 0.25)    # keep linear legend labels
  )

# Combine vertically with shared legend and main title
combined <- (p1 / p2 / p3) +
  plot_layout(ncol = 1, heights = c(7/12, 11/12, 1), guides = "collect") +
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
    "sensitivity_heatmap_combined.pdf"
  ),
  combined,
  width = 22,
  height = 29.5,
  dpi = 300
)

ggsave(
  fs::path(
    figures_folder,
    "sensitivity_heatmap_day1.pdf"
  ),
  p1,
  width = 22,
  height = 12,
  dpi = 300
)

ggsave(
  fs::path(
    figures_folder,
    "sensitivity_heatmap_day3.pdf"
  ),
  p2,
  width = 22,
  height = 12,
  dpi = 300
)

ggsave(
  fs::path(
    figures_folder,
    "sensitivity_heatmap_day7.pdf"
  ),
  p3,
  width = 22,
  height = 12,
  dpi = 300
)

