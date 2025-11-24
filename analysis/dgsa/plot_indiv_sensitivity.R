# R script to produce figures showing sensitivity results to changes in singular hyperparameters

# Libraries
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(glue)
library(purrr)
library(stringr)
library(patchwork)

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

rm(results_dearseq, results_qusage)

indiv_sensitivity_plot = function(conditions = NULL,
                                  timepoint = NULL,
                                  aggregates = NULL,
                                  parameter = c("method",
                                                "p_correction",
                                                "p_approach",
                                                "p_threshold",
                                                "fc_threshold"),
                                  groupby = c("byparameter", "byvaccine"),
                                  p_threshold_range = c(0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1),
                                  p_threshold_fix = 0.05,
                                  method_range = c("dearseq", "qusage"),
                                  method_fix = c("dearseq"),
                                  p_correction_range = c("BH", "BY", "bonferroni", "hochberg", "holm", "hommel"),
                                  p_correction_fix = c("BH"),
                                  p_approach_range = c("global", "withinTime", "withinComparison"),
                                  p_approach_fix = "global",
                                  fc_threshold_range = seq(0, 2, 0.1),
                                  fc_threshold_fix = 0) {
  # --- Helpers ---------------------------------------------------------------
  base_theme <- function(angle_x = 0) {
    theme_minimal(base_size = 14) +
      theme(
        plot.title            = element_text(
          hjust = 0.5,
          face = "bold",
          size = 40
        ),
        axis.title            = element_text(face = "bold", size = 30),
        axis.text             = element_text(size = 25),
        axis.text.x           = element_text(angle = angle_x, hjust = ifelse(angle_x == 0, 0.5, 1)),
        legend.title          = element_text(face = "bold", size = 25),
        legend.text           = element_text(size = 20),
        plot.caption          = element_text(
          size = 14,
          hjust = 0,
          lineheight = 1.2,
          margin = margin(t = 15)
        ),
        plot.caption.position = "plot"
      )
  }
  
  make_caption <- function(fixed_msg_lines) {
    paste0("Fixed parameters:\n",
           paste0("- ", fixed_msg_lines, collapse = "\n"))
  }
  
  # convert wide -> long for the adjPval columns
  dt_long <- results_df %>%
    filter(
      time == timepoint,
      condition %in% conditions,
      gs.aggregate %in% aggregates,!is.na(rawPval)
    ) %>%
    pivot_longer(
      cols = matches("\\.adjPval_"),
      names_to = c("approach", "p_method"),
      names_pattern = "(.*)\\.adjPval_(.*)",
      values_to = "adj_pval"
    )
  
  # small utility to compute prop_signif for given filtered dt_long and grouping variables
  compute_prop <- function(df, group_cols) {
    df %>%
      mutate(significant = (
        adj_pval <= p_threshold_fix &
          abs(fc.score) >= fc_threshold_fix
      )) %>%
      group_by(!!!rlang::syms(group_cols), .add = FALSE) %>%
      summarise(prop_signif = mean(significant), .groups = "drop")
  }
  
  # Start main switch on parameter
  param <- parameter
  
  if (param == "method") {
    # keep the same filters as original
    dtf <- dt_long %>%
      filter(approach == p_approach_fix,
             p_method == p_correction_fix,
             method %in% method_range)
    df_prop_de <- compute_prop(dtf, c("condition", "method", "condition.colour")) %>%
      mutate(condition = fct_reorder(condition, prop_signif, .fun = mean))
    
    if (groupby == "byparameter") {
      plot <- df_prop_de %>%
        ggplot(aes(x = method, y = prop_signif, fill = condition)) +
        scale_fill_manual(
          name = "Vaccine",
          values = set_names(df_prop_de$condition.colour, df_prop_de$condition),
          breaks = levels(df_prop_de$condition),
          labels = levels(df_prop_de$condition)
        ) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        labs(
          title = glue(
            "Proportion of DE Genesets by DGSA Method\n(Day {timepoint})"
          ),
          x = "DGSA Method",
          y = "Proportion of DE genesets",
          caption = make_caption(c(
            paste0(
              "p value correction: ",
              p_correction_fix,
              "     p value subset: ",
              p_approach_fix
            ),
            paste0(
              "significance level: ",
              p_threshold_fix,
              "     abs log2 FC threshold: ",
              fc_threshold_fix
            )
          ))
        ) +
        base_theme() + ylim(0, 1)
      
    } else {
      # byvaccine
      plot <- df_prop_de %>%
        ggplot(aes(x = condition, y = prop_signif, fill = method)) +
        scale_fill_manual(
          name = "DGSA Method",
          values = set_names(c("#1b9e77", "#d95f02"), unique(df_prop_de$method)),
          breaks = unique(df_prop_de$method),
          labels = unique(df_prop_de$method)
        ) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        labs(
          title = glue(
            "Proportion of DE Genesets by DGSA Method\n(Day {timepoint})"
          ),
          x = "Vaccine",
          y = "Proportion of DE genesets",
          caption = make_caption(c(
            paste0(
              "p value correction: ",
              p_correction_fix,
              "    p value subset: ",
              p_approach_fix
            ),
            paste0(
              "significance level: ",
              p_threshold_fix,
              "    abs log2 FC threshold: ",
              fc_threshold_fix
            )
          ))
        ) +
        base_theme(angle_x = 45) + ylim(0, 1)
    }
    
  } else if (param == "p_correction") {
    dtf <- dt_long %>%
      filter(approach == p_approach_fix, method == method_fix)  # preserved original variable usage
    df_prop_de <- compute_prop(dtf, c("condition", "p_method", "condition.colour")) %>%
      mutate(condition = fct_reorder(condition, prop_signif, .fun = mean))
    
    if (groupby == "byparameter") {
      plot <- df_prop_de %>%
        ggplot(aes(x = p_method, y = prop_signif, fill = condition)) +
        scale_fill_manual(
          name = "Vaccine",
          values = set_names(df_prop_de$condition.colour, df_prop_de$condition),
          breaks = levels(df_prop_de$condition),
          labels = levels(df_prop_de$condition)
        ) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        labs(
          title = paste0(
            "Proportion of DE Genesets by p-value Correction Method \n(Day ",
            timepoint,
            ")"
          ),
          x = "p-value Correction Method",
          y = "Proportion of DE genesets",
          caption = make_caption(c(
            paste0(
              "DGSA Method: ",
              method_fix,
              "     p value subset: ",
              p_approach_fix
            ),
            paste0(
              "significance level: ",
              p_threshold_fix,
              "     abs log2 FC threshold: ",
              fc_threshold_fix
            )
          ))
        ) +
        base_theme() + ylim(0, 1)
      
    } else {
      plot <- df_prop_de %>%
        ggplot(aes(x = condition, y = prop_signif, fill = p_method)) +
        scale_fill_manual(
          name = "p-value Correction Method",
          values = set_names(
            c(
              "#1b9e77",
              "#d95f02",
              "#7570b3",
              "#e7298a",
              "#66a61e",
              "#e6ab02"
            ),
            unique(df_prop_de$p_method)
          ),
          breaks = unique(df_prop_de$p_method),
          labels = unique(df_prop_de$p_method)
        ) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        labs(
          title = glue(
            "Proportion of DE Genesets by Correction Method \n(Day {timepoint})"
          ),
          x = "Vaccine",
          y = "Proportion of DE genesets",
          caption = make_caption(c(
            paste0(
              "DGSA Method: ",
              method_fix,
              "    p-value subset: ",
              p_approach_fix
            ),
            paste0(
              "significance level: ",
              p_threshold_fix,
              "    abs log2-FC threshold: ",
              fc_threshold_fix
            )
          ))
        ) +
        base_theme(angle_x = 45) + ylim(0, 1)
    }
    
  } else if (param == "p_approach") {
    dtf <- dt_long %>%
      filter(p_method == p_correction_fix, method == method_fix)
    df_prop_de <- compute_prop(dtf, c("condition", "approach", "condition.colour")) %>%
      mutate(condition = fct_reorder(condition, prop_signif, .fun = mean))
    
    if (groupby == "byparameter") {
      plot <- df_prop_de %>%
        ggplot(aes(x = approach, y = prop_signif, fill = condition)) +
        scale_fill_manual(
          name = "Vaccine",
          values = set_names(df_prop_de$condition.colour, df_prop_de$condition),
          breaks = levels(df_prop_de$condition),
          labels = levels(df_prop_de$condition)
        ) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        labs(
          title = paste0(
            "Proportion of DE Genesets by p-value adjustment subset \n(Day ",
            timepoint,
            ")"
          ),
          x = "p-value adjustment subset",
          y = "Proportion of DE genesets",
          caption = make_caption(c(
            paste0(
              "DGSA Method: ",
              method_fix,
              "     p value correction: ",
              p_correction_fix
            ),
            paste0(
              "significance level: ",
              p_threshold_fix,
              "     abs log2 FC threshold: ",
              fc_threshold_fix
            )
          ))
        ) +
        base_theme() + ylim(0, 1)
      
    } else {
      plot <- df_prop_de %>%
        ggplot(aes(x = condition, y = prop_signif, fill = approach)) +
        scale_fill_manual(
          name = "p-value Subset",
          values = set_names(
            c("#1b9e77", "#d95f02", "#7570b3"),
            unique(df_prop_de$approach)
          ),
          breaks = unique(df_prop_de$approach),
          labels = unique(df_prop_de$approach)
        ) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        labs(
          title = glue(
            "Proportion of DE Genesets by p-value adjustment subset\n(Day {timepoint})"
          ),
          x = "Vaccine",
          y = "Proportion of DE genesets",
          caption = make_caption(c(
            paste0(
              "DGSA Method: ",
              method_fix,
              "    p value correction: ",
              p_correction_fix
            ),
            paste0(
              "significance level: ",
              p_threshold_fix,
              "    abs log2-FC threshold: ",
              fc_threshold_fix
            )
          ))
        ) +
        base_theme(angle_x = 45) + ylim(0, 1)
    }
    
  } else if (param == "p_threshold") {
    # preserved original behaviour: map over p_threshold_range (as in original script)
    dtf <- dt_long %>%
      filter(p_method == p_correction_fix, method == method_fix)
    
    dtf2 <- bind_cols(dtf, map_dfc(p_threshold_range, function(thresh_char) {
      thresh <- as.numeric(thresh_char)
      indicator <- with(dtf, abs(fc.score) >= fc_threshold_fix &
                          adj_pval <= thresh)
      colname <- paste0("significant_threshold_", gsub("\\.", "_", thresh_char))
      tibble(!!colname := indicator)
    }))
    
    df_prop_de <- dtf2 %>%
      pivot_longer(
        cols = starts_with("significant_threshold_"),
        names_to = "threshold",
        names_prefix = "significant_threshold_",
        values_to = "significant"
      ) %>%
      mutate(threshold = as.numeric(str_replace(threshold, "_", "."))) %>%
      group_by(condition, condition.colour, threshold) %>%
      summarise(prop_signif = mean(significant), .groups = "drop") %>%
      mutate(condition = fct_reorder(condition, prop_signif, .fun = mean))
    
    # formatting helper copied from original
    custom_pct <- function(x) {
      pct <- x * 100
      vapply(pct, function(p) {
        if (p < 1) {
          txt <- formatC(p, format = "f", digits = 5)
          txt <- sub("\\.?0+$", "", txt)
        } else {
          txt <- as.character(round(p))
        }
        paste0(txt, "%")
      }, FUN.VALUE = character(1))
    }
    
    plot <- df_prop_de %>%
      ggplot(aes(
        x = threshold,
        y = prop_signif,
        colour = condition,
        group = condition
      )) +
      scale_color_manual(
        name = "Vaccine",
        values = set_names(df_prop_de$condition.colour, df_prop_de$condition),
        breaks = levels(df_prop_de$condition),
        labels = levels(df_prop_de$condition)
      ) +
      geom_point(size = 2) + geom_line(linewidth = 2, alpha = 0.7) +
      scale_x_log10(breaks = as.numeric(p_threshold_range), labels = custom_pct) +
      labs(
        title = paste0(
          "Proportion of DE Genesets by significance level\n(Day ",
          timepoint,
          ")"
        ),
        x = "Significance level",
        y = "Proportion of DE genesets",
        caption = make_caption(c(
          paste0(
            "DGSA Method: ",
            method_fix,
            "     p value correction: ",
            p_correction_fix
          ),
          paste0(
            "p-value subset: ",
            p_approach_fix,
            "     abs log2 FC threshold: ",
            fc_threshold_fix
          )
        ))
      ) +
      base_theme() + ylim(0, 1)
    
  } else if (param == "fc_threshold") {
    # preserved original behaviour: map over fc_threshold_range
    dtf <- dt_long %>%
      filter(p_method == p_correction_fix, method == method_fix)
    
    dtf2 <- bind_cols(dtf, map_dfc(fc_threshold_range, function(thresh) {
      indicator <- with(dtf, abs(fc.score) >= thresh &
                          adj_pval <= p_threshold_fix)
      colname <- paste0("significant_threshold_", gsub("\\.", "_", format(thresh)))
      tibble(!!colname := indicator)
    }))
    
    df_prop_de <- dtf2 %>%
      pivot_longer(
        cols = starts_with("significant_threshold_"),
        names_to = "threshold",
        names_prefix = "significant_threshold_",
        values_to = "significant"
      ) %>%
      mutate(threshold = as.numeric(str_replace(threshold, "_", "."))) %>%
      group_by(condition, condition.colour, threshold) %>%
      summarise(prop_signif = mean(significant), .groups = "drop")
    
    plot <- df_prop_de %>%
      ggplot(aes(
        x = threshold,
        y = prop_signif,
        colour = condition,
        group = condition
      )) +
      scale_color_manual(
        name = "Vaccine",
        values = set_names(df_prop_de$condition.colour, df_prop_de$condition),
        breaks = levels(df_prop_de$condition),
        labels = levels(df_prop_de$condition)
      ) +
      geom_point(size = 2) + geom_line(linewidth = 2, alpha = 0.7) +
      labs(
        title = paste0(
          "Proportion of DE Genesets by Absolute log2-fold-change threshold\n(Day ",
          timepoint,
          ")"
        ),
        x = "Absolute log2-FC threshold",
        y = "Proportion of DE genesets",
        caption = make_caption(c(
          paste0(
            "DGSA Method: ",
            method_fix,
            "     p value correction: ",
            p_correction_fix
          ),
          paste0(
            "p-value subset: ",
            p_approach_fix,
            "     significance level: ",
            p_threshold_fix
          )
        ))
      ) +
      base_theme() + ylim(0, 1)
  }
  
  # final object (same as original)
  return(plot)
  
}

conditions = results_df %>%
  pull(condition) %>%
  unique()

aggregates = results_df %>%
  pull(gs.aggregate) %>%
  unique()

day = 1

p1 = indiv_sensitivity_plot(
  conditions = conditions,
  timepoint = day,
  aggregates = aggregates,
  parameter = "method",
  groupby = "byvaccine",
  method_range = c("dearseq", "qusage"),
  p_threshold_range = c(0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1),
  p_correction_range = c("BH", "BY", "bonferroni", "hochberg", "holm", "hommel"),
  p_approach_range = c("global", "withinTime", "withinComparison"),
  fc_threshold_range = seq(0, 2, 0.1),
  method_fix = "dearseq",
  p_threshold_fix = 0.05,
  p_correction_fix = "BH",
  p_approach_fix = "global",
  fc_threshold_fix = 0
)

p1

p2 = indiv_sensitivity_plot(
  conditions = conditions,
  timepoint = day,
  aggregates = aggregates,
  parameter = "p_threshold",
  groupby = "byvaccine",
  method_range = c("dearseq", "qusage"),
  p_threshold_range = c(0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1),
  p_correction_range = c("BH", "BY", "bonferroni", "hochberg", "holm", "hommel"),
  p_approach_range = c("global", "withinTime", "withinComparison"),
  fc_threshold_range = seq(0, 2, 0.1),
  method_fix = "dearseq",
  p_threshold_fix = 0.05,
  p_correction_fix = "BH",
  p_approach_fix = "global",
  fc_threshold_fix = 0
)

p2

p3 = indiv_sensitivity_plot(
  conditions = conditions,
  timepoint = day,
  aggregates = aggregates,
  parameter = "fc_threshold",
  groupby = "byvaccine",
  method_range = c("dearseq", "qusage"),
  p_threshold_range = c(0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1),
  p_correction_range = c("BH", "BY", "bonferroni", "hochberg", "holm", "hommel"),
  p_approach_range = c("global", "withinTime", "withinComparison"),
  fc_threshold_range = seq(0, 2, 0.1),
  method_fix = "dearseq",
  p_threshold_fix = 0.05,
  p_correction_fix = "BH",
  p_approach_fix = "global",
  fc_threshold_fix = 0
)

p3

# Add vertical spacing and left-align titles
p1 <- p1 +
  labs(title = "a) DGSA method") +
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
  )

p2 <- p2 +
  labs(title = "b) Significance level") +
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
  )

p3 <- p3 +
  labs(title = "c) Fold-change threshold") +
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
  )

# Combine vertically with shared legend and main title
combined <- (p1 / p2 / p3) +
  plot_layout(ncol = 1, heights = c(1, 1, 1)) +
  plot_annotation(
    title = paste0(
      "Sensitivity of results to individual hyperparameter changes (day ",
      day,
      ")"
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 55,
        hjust = 0.5,
        margin = margin(b = 1)
      )  # extra bottom margin
    )
  ) &
  theme(legend.position = "right")

# Display
print(combined)

figures_folder = fs::path("output", "figures", "dgsa")

# Save
ggsave(
  fs::path(
    figures_folder,
    paste0("indiv_sensitivity_combined_day", day, ".pdf")
  ),
  combined,
  width = 25,
  height = 35,
  dpi = 300
)

# Plot individual figures

ggsave(
  fs::path(
    figures_folder,
    paste0("indiv_sensitivity_dgsamethod_day", day, ".pdf")
  ),
  p1,
  width = 25,
  height = 15,
  dpi = 300
)

ggsave(
  fs::path(
    figures_folder,
    paste0("indiv_sensitivity_pthreshold_day", day, ".pdf")
  ),
  p2,
  width = 25,
  height = 15,
  dpi = 300
)

ggsave(
  fs::path(
    figures_folder,
    paste0("indiv_sensitivity_fcthreshold_day", day, ".pdf")
  ),
  p3,
  width = 25,
  height = 15,
  dpi = 300
)
