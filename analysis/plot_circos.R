# R function to plot a circos plot given results from a multiple-condition differential gene set analysis

# Libraries
suppressPackageStartupMessages({
  library(tidyselect)
  library(circlize)
})

# Load data
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

plot_circos = function(method_name,
                       conditions,
                       timepoint,
                       aggregates_name,
                       p_correction = c("BH", "bonferroni", "holm", "hommel", "hochberg", "BY"),
                       p_approach = c("global", "withinTime", "withinComparison"),
                       p_threshold = 0.05,
                       filter_variable = c("none", "fc.score", "activation.score"),
                       filter_mode = c("user", "data"),
                       user_threshold = 0.5,
                       quantile_threshold = 0.5,
                       scores = c("fc.score", "activation.score"),
                       correlation = c("corr.mean", "mean.corr"),
                       arc = c("any", "direction", "positive", "negative"),
                       ring = c("all", "expression", "none"),
                       order = c("set_all", "set_available", "cluster"),
                       quantile_scoreclip = 0.995,
                       legend = TRUE,
                       plot_title = NULL,
                       title_size = 2.5,
                       title_line = 0) {
  # Subset the results by the user selected method, conditions and aggregates
  results_df_circos = results_df %>%
    filter(gs.aggregate %in%  aggregates_name,
           condition %in% conditions,
           method == method_name) %>%
    mutate(gs.aggregate = droplevels(gs.aggregate))
  
  # Based on the remaining genesets, define a global ordering
  global_order <- results_df_circos %>%
    select(gs.name, gs.aggregate) %>%
    distinct() %>%
    arrange(gs.aggregate) %>%
    pull(gs.name)
  
  # Apply this global gene set ordering to the dataframe
  results_df_circos <- results_df_circos %>%
    mutate(gs.name = factor(gs.name, levels = global_order))
  
  # Extract a list of conditions to plot
  ## If user selects all to appear (circos_order == set_all), they should all appear regardless of if there are significant results
  if (order == "set_available" | order == "cluster") {
    condition_names = results_df_circos %>%
      filter(time == timepoint) %>%
      pull(condition) %>%
      unique()
  } else if (order == "set_all") {
    condition_names <- unique(results_df$condition)
  }
  
  # Drop levels which are not to be included
  condition_order = levels(results_df$condition)
  condition_order = condition_order[condition_order %in% condition_names]
  
  # Create a dataframe containing the metadata for the condition_names
  # This contains the ordering of the conditions in the circos plot, the availability of data, and the availability of the response correlation
  
  if (order == "set_all" | order == "set_available") {
    # If the user selects a fixed ordering :
    
    condition_circos_metadata <- data.frame(
      order = match(condition_names, condition_order),
      condition = factor(condition_names, levels = condition_order),
      stringsAsFactors = FALSE
    ) %>%
      arrange(condition)
    
  } else if (order == "cluster") {
    # Else if the user wants the order to be determined by hierarchical clustering of scores
    sc <- scores
    mat_df <- results_df_circos %>%
      filter(time == timepoint) %>%
      select(condition, gs.name, score = .data[[sc]]) %>%
      pivot_wider(names_from  = gs.name, values_from = score) %>%
      # now every missing goes to NA; replace them with 0
      mutate(across(-condition, ~ tidyr::replace_na(.x, 0)))
    
    # 2) Turn into a numeric matrix with condition names as row names
    mat <- mat_df %>%
      column_to_rownames("condition") %>%
      as.matrix()
    
    # 3) Compute distances & cluster
    #    Here we use Euclidean; you could also do:
    #      d <- as.dist(1 - cor(t(mat)))    # 1 - Pearson correlation
    d <- dist(mat, method = "euclidean")
    hc <- hclust(d, method = "ward.D2")   # or "complete", "ward.D2", ...
    
    # 4) Extract the condition order
    ordered_conditions <- hc$labels[hc$order]
    
    condition_circos_metadata <- data.frame(
      order = match(condition_names, ordered_conditions),
      condition = factor(condition_names, levels = ordered_conditions),
      stringsAsFactors = FALSE
    ) %>%
      arrange(condition)
  }
  
  # Create columns in the condition metadata to describe the condition availability at each timepoint (gene expression + immune response)
  unique_times <- unique(results_df_circos$time)
  
  condition_colors <- results_df %>%
    select(condition, condition.colour) %>%
    distinct() %>%
    deframe()
  
  for (t in unique_times) {
    # Get condition_names for the expression condition
    conditions_expr_at_time <- results_df_circos %>%
      filter(time == t) %>%
      pull(condition) %>%
      unique()
    
    # Construct dynamic column names for expression availability and color
    expr_col_name <- paste0("condition_expr_available_", t)
    color_col <- paste0("condition_color_", t)
    text_col <- paste0("text_color_", t)
    
    # Get condition_names for the response condition (with non-NA cor_MFC_mean)
    conditions_resp_at_time <- results_df %>%
      filter(time == t, !is.na(corr.mean)) %>%
      pull(condition) %>%
      unique()
    
    # Construct the dynamic column name for response availability
    resp_col_name <- paste0("condition_response_available_", t)
    
    # Update condition_circos_metadata with the new columns
    condition_circos_metadata <- condition_circos_metadata %>%
      mutate(
        !!expr_col_name := condition %in% conditions_expr_at_time,!!resp_col_name := condition %in% conditions_resp_at_time,!!color_col := if_else(.data[[expr_col_name]], condition_colors[as.character(condition)], "#D3D3D3"),!!text_col  := if_else(.data[[expr_col_name]], "black", "#9999a1")
      )
  }
  
  # Set circos options
  condition_circos_metadata$xmin <- 0
  condition_circos_metadata$xmax <- length(global_order)
  
  # Map gene sets to fixed positions based on global order
  gene_set_positions <- setNames(seq(0.5, length(global_order) - 0.5, length.out = length(global_order)), global_order)
  
  # Define colours for each aggregate
  aggregate_colors <- results_df %>%
    select(gs.aggregate, gs.colour) %>%
    distinct() %>%
    {
      setNames(.$gs.colour, .$gs.aggregate)
    }
  
  # Compute metadata for the positions of the gene set aggregates on the segments
  aggregate_distinct <- results_df_circos %>%
    select(gs.name, gs.aggregate) %>%
    distinct()
  
  pt <- prop.table(table(aggregate_distinct$gs.aggregate))
  
  aggregate_proportions <- data.frame(
    gs.aggregate = names(pt),
    prop         = as.vector(pt),
    stringsAsFactors = FALSE
  )
  aggregate_proportions$colour = aggregate_colors[aggregate_proportions$gs.aggregate]
  
  # Compute cumulative proportions for each aggregate
  aggregate_proportions$cum_start <- c(0, head(cumsum(aggregate_proportions$prop), -1))
  aggregate_proportions$cum_end <- cumsum(aggregate_proportions$prop)
  
  # FILTRATION : filter significant results based on user inputs
  if (p_correction == "none") {
    results_df_significant = results_df_circos %>%
      filter(rawPval < p_threshold, time == timepoint)
  } else {
    results_df_significant = results_df_circos %>%
      filter(!!sym(paste0(
        p_approach, ".adjPval_", p_correction
      )) < p_threshold, time == timepoint)
    
  }
  
  # Filter results by effect size if desired
  if (filter_variable != "none" & filter_mode == "data") {
    # Global threshold for filtration
    if (p_approach == "global") {
      effect_size_threshold <- results_df_circos %>%
        pull(filter_variable) %>%
        abs() %>%
        quantile(quantile_threshold)
      results_df_significant = results_df_significant %>%
        filter(abs(.data[[filter_variable]]) > effect_size_threshold)
      
    } else if (p_approach == "withinTime") {
      effect_size_threshold <- results_df_circos %>%
        filter(time == timepoint) %>%
        pull(filter_variable) %>%
        abs() %>%
        quantile(quantile_threshold)
      
      results_df_significant = results_df_significant %>%
        filter(abs(.data[[filter_variable]]) > effect_size_threshold)
      
    } else if (p_approach == "withinCondition") {
      results_df_significant = results_df_circos %>%
        filter(time == timepoint) %>%
        group_by(condition) %>%
        mutate(abs_score = abs(.data[[filter_variable]])) %>%
        filter(abs_score >= quantile(abs_score, probs = quantile_threshold, na.rm = TRUE)) %>%
        ungroup()
    }
  } else if (filter_variable != "none" & filter_mode == "user") {
    results_df_significant = results_df_significant %>%
      filter(abs(.data[[filter_variable]]) > user_threshold)
  }
  
  
  # Identify the shared differentially expressed modules
  if (arc == "any") {
    shared_modules <- results_df_significant %>%
      group_by(gs.name) %>%
      filter(n() > 1) %>%
      summarise(condition_names = list(unique(condition)), .groups = "drop")
    
    any_links = (nrow(shared_modules) != 0)
    
    if (any_links) {
      # Prepare the data in the format to plot
      shared_links <- shared_modules %>%
        mutate(pairs = purrr::map(condition_names, ~ combn(.x, 2, simplify = FALSE))) %>%
        select(gs.name, pairs) %>%
        unnest(pairs) %>%
        unnest_wider(pairs, names_sep = "_") %>%
        rename(condition1 = pairs_1, condition2 = pairs_2) %>%
        mutate(position = gene_set_positions[gs.name])  # Add position for gene
    }
  } else if (arc == "direction") {
    column_name1 <- paste0(scores, 1)
    column_name2 <- paste0(scores, 2)
    
    # 1. Identify gene‐sets with >1 significant hit and nest their data
    shared_links <- results_df_significant %>%
      # keep only the columns we need
      select(gs.name, condition, .data[[scores]]) %>%
      
      # self-join on gs.name to get all condition pairs
      inner_join(
        results_df_significant %>% select(gs.name, condition, .data[[scores]]),
        by = "gs.name",
        suffix = c("1", "2"),
        relationship = "many-to-many"
      ) %>%
      mutate(!!column_name1 := as.numeric(.data[[column_name1]]),!!column_name2 := as.numeric(.data[[column_name2]])) %>%
      
      # avoid self‐pairs and duplicate pairs (A–B vs B–A)
      filter(condition1 != condition2) %>%
      
      # require same sign of activation_score (product > 0)
      filter(.data[[column_name1]] * .data[[column_name2]] > 0) %>%
      
      # select / rename for plotting
      transmute(
        gs.name,
        condition1 = condition1,
        condition2 = condition2,
        position = gene_set_positions[gs.name]
      ) %>%
      distinct()
    
    any_links = (nrow(shared_links) != 0)
  } else if (arc == "positive") {
    shared_modules <- results_df_significant %>%
      filter(.data[[scores]] >= 0) %>%
      group_by(gs.name) %>%
      filter(n() > 1) %>%
      summarise(condition_names = list(unique(condition)), .groups = "drop")
    
    any_links = (nrow(shared_modules) != 0)
    
    if (any_links) {
      # Prepare the data in the format to plot
      shared_links <- shared_modules %>%
        mutate(pairs = purrr::map(condition_names, ~ combn(.x, 2, simplify = FALSE))) %>%
        select(gs.name, pairs) %>%
        unnest(pairs) %>%
        unnest_wider(pairs, names_sep = "_") %>%
        rename(condition1 = pairs_1, condition2 = pairs_2) %>%
        mutate(position = gene_set_positions[gs.name])  # Add position for gene
    }
    
  } else if (arc == "negative") {
    shared_modules <- results_df_significant %>%
      filter(.data[[scores]] <= 0) %>%
      group_by(gs.name) %>%
      filter(n() > 1) %>%
      summarise(condition_names = list(unique(condition)), .groups = "drop")
    
    any_links = (nrow(shared_modules) != 0)
    
    if (any_links) {
      # Prepare the data in the format to plot
      shared_links <- shared_modules %>%
        mutate(pairs = purrr::map(condition_names, ~ combn(.x, 2, simplify = FALSE))) %>%
        select(gs.name, pairs) %>%
        unnest(pairs) %>%
        unnest_wider(pairs, names_sep = "_") %>%
        rename(condition1 = pairs_1, condition2 = pairs_2) %>%
        mutate(position = gene_set_positions[gs.name])  # Add position for gene
    }
    
  }
  
  # Initialize circos plot
  # par(mar = rep(0, 4))
  
  # Set layout: 1 row, 2 columns (75% for circos, 25% for legend)
  # layout(matrix(1:2, ncol = 2), widths = c(0.75, 0.25))
  
  circos.clear()
  # Set circos dimensions
  # Before you start plotting
  circos.par(
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0, 0.01),
    start.degree = 81,
    gap.degree = 2
  )
  circos.par("canvas.xlim" = c(-1.2, 1.2),
             "canvas.ylim" = c(-1.2, 1.2))
  # Initialise
  circos.initialize(
    factors = condition_circos_metadata$condition,
    xlim = cbind(
      condition_circos_metadata$xmin,
      condition_circos_metadata$xmax
    )
  )
  
  colour_variable_name = paste0("condition_color_", timepoint)
  text_variable_name   = paste0("text_color_", timepoint)
  expr_variable_name = paste0("condition_expr_available_", timepoint)
  
  if (!is.null(plot_title)) {
    title(plot_title, line = title_line, cex.main = title_size)  # Adjust line and cex.main as needed
  }
  
  
  ### LAYER 1 : SEGMENTS LABELLED BY CONDITION ###
  suppressMessages({
    # Plot condition sectors
    circos.trackPlotRegion(
      ylim = c(0, 1),
      factors = condition_circos_metadata$condition,
      track.height = 0.1,
      bg.lwd = 1,
      panel.fun = function(x, y) {
        name <- get.cell.meta.data("sector.index")
        i <- get.cell.meta.data("sector.numeric.index")
        xlim <- get.cell.meta.data("xlim")
        ylim <- get.cell.meta.data("ylim")
        theta <- circlize(mean(xlim), 1.3)[1, 1] %% 360
        dd <- ifelse(theta > 200 &&
                       theta < 340, "outside", "inside")
        circos.text(
          x = mean(xlim),
          y = ylim[2] + 0.4,
          # Adjusted y position
          labels = name,
          facing = dd,
          cex = 1.5,
          adj = c(0.5, 0.5),
          niceFacing = T,
          col = condition_circos_metadata[[text_variable_name]][i]
        )
        
        # Draw segment with black border
        circos.rect(
          xleft = xlim[1],
          ybottom = ylim[1],
          xright = xlim[2],
          ytop = ylim[2],
          col = condition_circos_metadata[[colour_variable_name]][i],
          border = "black"
        )  # Black border
        
        circos.axis(labels = FALSE, major.tick = FALSE)
      }
    )
  }) # End suppressMessages
  
  # Layer 2 : Correlation with response
  
  if (ring == "all") {
    response_variable_name = paste0("condition_response_available_", timepoint)
    suppressMessages({
      cor_data <- results_df_circos %>%
        filter(time == timepoint, !is.na(.data[[correlation]])) %>%
        group_by(condition) %>%
        summarise(
          cor_scores = list(.data[[correlation]]),
          gs.names = list(gs.name),
          .groups = "drop"
        ) %>%
        split(.$condition)
      
      cor_data <- cor_data[sapply(cor_data, nrow) > 0]
      
      # Draw the response layer track using the pre-calculated cor_data
      circos.trackPlotRegion(
        ylim = c(-1, 1),
        factors = condition_circos_metadata$condition,
        track.height = 0.1,
        bg.lwd = 1,
        bg.border = condition_circos_metadata[[response_variable_name]] %>% ifelse("black", "white"),
        panel.fun = function(x, y) {
          name <- get.cell.meta.data("sector.index")
          i <- get.cell.meta.data("sector.numeric.index")
          xlim <- get.cell.meta.data("xlim")
          ylim <- get.cell.meta.data("ylim")
          
          # Use pre-calculated data for the current sector
          if (name %in% names(cor_data)) {
            sec_data <- cor_data[[name]]
            # Extract correlation scores and gene set names
            correlation_scores <- sec_data$cor_scores[[1]]
            gs.names <- sec_data$gs.names[[1]]
            
            # Get positions for the gene sets relevant to this condition
            positions <-
              gene_set_positions[as.character(gs.names)]
            
            # Determine bar colors vectorized over correlation_scores
            bar_colors <-
              ifelse(correlation_scores > 0, "purple", "orange")
            
            # Normalize score lengths using vectorized computation
            bar_lengths <- correlation_scores
            
            # Draw all bars in one go if circos.rect supports vectorized arguments
            circos.rect(
              xleft = positions - 0.5,
              ybottom = 0,
              xright = positions + 0.5,
              ytop = bar_lengths,
              col = bar_colors,
              border = NA
            )
          }
        }
      )
      
    }) # End suppressMessages
  }
  
  ## LAYER 3 - ACTIVATION SCORES
  if (ring %in% c("all", "expression")) {
    all_scores <- results_df_circos %>%
      filter(time == timepoint) %>%
      pull(!!sym(scores)) %>%
      as.numeric()
    
    threshold <- quantile(abs(all_scores), quantile_scoreclip)
    
    # Clip the activation scores according to the threshold
    max_clipped_score <- max(abs(pmin(
      pmax(all_scores, -threshold), threshold
    )))
    
    # Aggregate the activation scores per condition and gene set.
    # For example, here we use the mean of duplicate entries.
    act_data <- results_df_circos %>%
      filter(time == timepoint) %>%
      group_by(condition, gs.name) %>%
      summarise(avg_score = mean(!!sym(scores), na.rm = TRUE),
                .groups = "drop") %>%
      group_by(condition) %>%
      summarise(
        gs.names = list(gs.name),
        raw_scores = list(avg_score),
        .groups = "drop"
      ) %>%
      mutate(clipped_scores = purrr::map(raw_scores, ~ pmin(pmax(.x, -threshold), threshold))) %>%
      split(.$condition)
    
    act_data <- act_data[sapply(act_data, nrow) > 0]
    
    suppressMessages({
      circos.trackPlotRegion(
        ylim = c(-1, 1),
        factors = condition_circos_metadata$condition,
        track.height = 0.1,
        bg.lwd = 1,
        bg.border = condition_circos_metadata[[paste0("condition_expr_available_", timepoint)]] %>% ifelse("black", "white"),
        panel.fun = function(x, y) {
          # Get current sector (condition)
          sector_name <- get.cell.meta.data("sector.index")
          xlim <- get.cell.meta.data("xlim")
          ylim <- get.cell.meta.data("ylim")
          
          if (sector_name %in% names(act_data)) {
            sector_act <- act_data[[sector_name]]
            
            # Extract pre-calculated data vectors
            clipped_scores <- sector_act$clipped_scores[[1]]
            gs.names <- sector_act$gs.names[[1]]
            
            # Use match to get positions so that vector lengths are consistent.
            positions <- gene_set_positions[match(gs.names, names(gene_set_positions))]
            
            # Remove any missing positions (if some gene sets are not in the global order)
            keep <- !is.na(positions)
            positions <- positions[keep]
            clipped_scores <- clipped_scores[keep]
            
            # Define colors for each score: red for positive, blue for negative
            bar_colors <- ifelse(clipped_scores > 0, "red", "blue")
            
            # Normalize score lengths to scale within the track
            bar_lengths <- ifelse(
              clipped_scores < 0,
              clipped_scores * -(ylim[1] / max_clipped_score),
              clipped_scores * (ylim[2] / max_clipped_score)
            )
            
            # Draw all bars in a single, vectorized call
            circos.rect(
              xleft = positions - 0.5,
              ybottom = 0,
              xright = positions + 0.5,
              ytop = bar_lengths,
              col = bar_colors,
              border = NA
            )
          }
        }
      )
    }) # End SuppressMessages
  }
  
  # LAYER 4 : GENE SET AGGREGATE POSITIONS
  suppressMessages({
    circos.trackPlotRegion(
      ylim = c(0, 1),
      factors = condition_circos_metadata$condition,
      track.height = 0.05,
      bg.lwd = 1,
      bg.border = condition_circos_metadata[[expr_variable_name]] %>% ifelse("black", "white"),
      panel.fun = function(x, y) {
        # Get sector and drawing parameters once per panel call
        sector_index <- get.cell.meta.data("sector.index")
        xlim <- get.cell.meta.data("xlim")
        
        # Directly extract condition availability using indexing
        condition_available <- condition_circos_metadata[condition_circos_metadata$condition == sector_index, expr_variable_name]
        
        # Vectorize the calculation of start and end positions
        x_start <- xlim[1] + (xlim[2] - xlim[1]) * aggregate_proportions$cum_start
        x_end   <- xlim[1] + (xlim[2] - xlim[1]) * aggregate_proportions$cum_end
        
        # Determine colours based on condition availability
        if (condition_available) {
          col_fill <- aggregate_proportions$colour
          border_col <- rep(NA, length(col_fill))
        } else {
          col_fill <- rep("white", length(x_start))
          border_col <- rep("white", length(x_start))
        }
        
        # Draw all segments in a single, vectorized call
        circos.rect(
          xleft = x_start,
          ybottom = 0,
          xright = x_end,
          ytop = 1,
          col = col_fill,
          border = border_col
        )
      }
    )
  }) # End suppressMessages
  
  # PLOT ARCS
  suppressMessages({
    if (any_links) {
      # Pre-calculate the aggregate values for each gene set (if not already present)
      gs.aggregates <- unique(results_df_circos[, c("gs.name", "gs.aggregate")])
      shared_links <- merge(shared_links, gs.aggregates, by = "gs.name")
      
      # make sure gs.aggregate is character
      shared_links$gs.aggregate <- as.character(shared_links$gs.aggregate)
      
      # now name-based lookup “just works”
      shared_links$link_colour <- aggregate_colors[shared_links$gs.aggregate]
      
      
      # Use mapply to vectorize the drawing of links
      mapply(
        function(condition1,
                 condition2,
                 position,
                 link_colour) {
          circos.link(
            sector.index1 = condition1,
            point1 = c(position, position),
            sector.index2 = condition2,
            point2 = c(position, position),
            col = adjustcolor(link_colour, alpha.f = 0.3)
          )
        },
        shared_links$condition1,
        shared_links$condition2,
        shared_links$position,
        shared_links$link_colour
      )
    }
  }) # End suppressMessages
  
  # PLOT LEGENDS
  if (legend) {
    if (ring == "all") {
      legend_labels <- c(
        expression(bold("Gene Set Aggregate")),
        aggregates_name,
        expression(bold("Direction of Regulation")),
        "Upregulated",
        "Downregulated",
        expression(bold("Correlation with Ab Response")),
        "Positive",
        "Negative"
      )
      
      legend_colours <- c(
        NA,
        aggregate_colors[aggregates_name],
        # Gene Set Aggregate
        NA,
        "red",
        "blue",
        # Gene Set Activity
        NA,
        "purple",
        "orange" # Correlation with MFC
      )
      
      # Plot the combined legend with optimized spacing and manually larger titles
      legend(
        "right",
        legend = legend_labels,
        fill = legend_colours,
        border = "white",
        cex = c(
          1.5,
          rep(1.2, length(aggregates_name)),
          # Gene Set Aggregate
          1.5,
          1.2,
          1.2,
          # Gene Set Activity
          1.5,
          1.2,
          1.2
        ),
        # Correlation with MFC
        text.width = 1.2 * max(strwidth(legend_labels, cex = 1.5)),
        # Ensures alignment
        y.intersp = c(
          1,
          rep(0.8, length(aggregates_name)),
          # Gene Set Aggregate
          1.5,
          0.8,
          0.8,
          # Gene Set Activity
          1.5,
          0.8,
          0.8
        ),
        # Correlation with MFC
        inset = c(0, 0.01),
        ncol = 1,
        bty = "o",
        xjust = 0,
        text.col = rep("black", length(legend_labels)),
        title = NULL
      )  # Disable default title
    } else if (ring == "expression") {
      legend_labels <- c(
        expression(bold("Gene Set Aggregate")),
        aggregates_name,
        expression(bold("Direction of Regulation")),
        "Upregulated",
        "Downregulated"
      )
      
      legend_colours <- c(
        NA,
        aggregate_colors[aggregates_name],
        # Gene Set Aggregate
        NA,
        "red",
        "blue"     # Gene Set Activity
      )
      
      # Plot the combined legend with optimized spacing and manually larger titles
      legend(
        "right",
        legend = legend_labels,
        fill = legend_colours,
        border = "white",
        cex = c(
          1.5,
          rep(1.2, length(aggregates_name)),
          # Gene Set Aggregate
          1.5,
          1.2,
          1.2
        ),
        # Correlation with MFC
        text.width = 1.2 * max(strwidth(legend_labels, cex = 1.5)),
        # Ensures alignment
        y.intersp = c(
          1,
          rep(0.8, length(aggregates_name)),
          # Gene Set Aggregate
          1.5,
          0.8,
          0.8
        ),
        # Correlation with MFC
        inset = c(0, 0.01),
        ncol = 1,
        bty = "o",
        xjust = 0,
        text.col = rep("black", length(legend_labels)),
        title = NULL
      )  # Disable default title
      
      
    } else if (ring == "none") {
      legend_labels <- c(expression(bold("Gene Set Aggregate")), aggregates_name)
      
      legend_colours <- c(NA, aggregate_colors[aggregates_name])
      
      # Plot the combined legend with optimized spacing and manually larger titles
      legend(
        "right",
        legend = legend_labels,
        fill = legend_colours,
        border = "white",
        cex = c(1.5, rep(1.2, length(
          aggregates_name
        ))),
        # Correlation with MFC
        text.width = 1.2 * max(strwidth(legend_labels, cex = 1.5)),
        # Ensures alignment
        y.intersp = c(1, rep(0.8, length(
          aggregates_name
        ))),
        # Correlation with MFC
        inset = c(0, 0.01),
        ncol = 1,
        bty = "o",
        xjust = 0,
        text.col = rep("black", length(legend_labels)),
        title = NULL
      )  # Disable default title
      
    }
  }
}

# Function to draw the legend
# Function to draw circos legend separately
draw_circos_legend <- function(aggregates_name,
                               ring = c("all", "expression", "none"),
                               placeholder = FALSE,
                               scale = 1) {
  ring <- match.arg(ring)
  
  # Extract colours automatically, just like in plot_circos
  aggregate_colors <- results_df %>%
    select(gs.aggregate, gs.colour) %>%
    distinct() %>%
    {
      setNames(.$gs.colour, .$gs.aggregate)
    }
  
  # Define legend content based on ring type
  if (ring == "all") {
    legend_labels <- c(
      expression(bold("Gene Set Aggregate")),
      aggregates_name,
      expression(bold("Direction of Regulation")),
      "Upregulated",
      "Downregulated",
      expression(bold("Correlation with Ab Response")),
      "Positive",
      "Negative"
    )
    legend_colours <- c(NA,
                        aggregate_colors[aggregates_name],
                        NA,
                        "red",
                        "blue",
                        NA,
                        "purple",
                        "orange")
    legend_cex <- c(1.5, rep(1.2, length(aggregates_name)), 1.5, 1.2, 1.2, 1.5, 1.2, 1.2) * scale
    legend_yintersp <- c(1, rep(0.8, length(aggregates_name)), 1.5, 0.8, 0.8, 1.5, 0.8, 0.8) * scale
  } else if (ring == "expression") {
    legend_labels <- c(
      expression(bold("Gene Set Aggregate")),
      aggregates_name,
      expression(bold("Direction of Regulation")),
      "Upregulated",
      "Downregulated"
    )
    legend_colours <- c(NA, aggregate_colors[aggregates_name], NA, "red", "blue")
    legend_cex <- c(1.5, rep(1.2, length(aggregates_name)), 1.5, 1.2, 1.2) * scale
    legend_yintersp <- c(1, rep(0.8, length(aggregates_name)), 1.5, 0.8, 0.8) * scale
  } else if (ring == "none") {
    legend_labels <- c(expression(bold("Gene Set Aggregate")), aggregates_name)
    legend_colours <- c(NA, aggregate_colors[aggregates_name])
    legend_cex <- c(1.5, rep(1.2, length(aggregates_name))) * scale
    legend_yintersp <- c(1, rep(0.8, length(aggregates_name))) * scale
  }
  
  max_cex <- max(legend_cex)
  
  if (!placeholder) {
    legend(
      "center",
      legend = legend_labels,
      fill = legend_colours,
      border = "white",
      cex = legend_cex,
      y.intersp = legend_yintersp,
      text.width = 1.2 * max(strwidth(legend_labels, cex = max_cex)) * scale,
      ncol = 1,
      bty = "o",
      xjust = 0,
      text.col = rep("black", length(legend_labels)),
      title = NULL
    )
  } else {
    blank_labels <- rep(" ", length(legend_labels))
    
    legend(
      "center",
      legend = blank_labels,
      fill = rep(NA, length(blank_labels)),
      border = rep(NA, length(blank_labels)),
      cex = legend_cex,
      y.intersp = legend_yintersp,
      text.width = 1.2 * max(strwidth(legend_labels, cex = max_cex)) * scale,
      ncol = 1,
      bty = "n",
      xjust = 0,
      text.col = rep(NA, length(blank_labels)),
      title = NULL
    )
  }
}


# Annotate the row
plot_row_annotation <- function(label = NULL,
                                height = 0.05,
                                size = 5,
                                placeholder = FALSE) {
  par(mar = c(0, 0, 0, 0))
  plot(
    0,
    0,
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    xlim = c(0, 1),
    ylim = c(0, 1)
  )
  
  if (!placeholder && !is.null(label)) {
    text(0.5, 0.5, label, cex = size, font = 2)  # centered, bold
  }
  # if placeholder = TRUE, nothing is drawn, but space is reserved
}



figures_folder = fs::path("output", "figures", "dgsa")

pdf(
  fs::path(figures_folder, "circos_side_by_side.pdf"),
  width = 24,
  height = 25
)

# Define layout: 3 columns: left plot, right plot, legend
layout(
  matrix(1:12, nrow = 3, byrow = TRUE),
  # 3 rows × 4 columns
  widths  = c(0.1, 0.35, 0.35, 0.25),
  # column widths
  heights = c(0.33, 0.33, 0.33)
)         # row heights
par(mar = rep(0, 4))  # no margins

# DAY 1 #
plot_row_annotation("Day 1", size = 6, placeholder = FALSE)

circos.clear()# make sure canvas is clean

plot_circos(
  method_name = "qusage",
  conditions = levels(results_df$condition),
  timepoint = 1,
  aggregates_name = c(
    "Antigen Presentation",
    "Inflammatory/TLR/Chemokines",
    "Interferon/Antiviral Sensing",
    "Monocytes",
    "DC Activation",
    "Neutrophils",
    "NK Cells",
    "Signal Transduction",
    "ECM And Migration",
    "Energy Metabolism",
    "Cell Cycle",
    "Platelets",
    "T Cells",
    "B Cells",
    "Plasma Cells"
  ),
  p_correction = "BH",
  p_approach = "global",
  p_threshold = 0.05,
  filter_variable = "none",
  filter_mode = "user",
  user_threshold = 0.5,
  quantile_threshold = 0.5,
  scores = "fc.score",
  correlation = "corr.mean",
  arc = "positive",
  ring = "expression",
  order = "set_all",
  quantile_scoreclip = 0.995,
  legend = FALSE,
  plot_title = "qusage",
  title_size = 6,
  title_line = -4
)

circos.clear()
plot_circos(
  method_name = "dearseq",
  conditions = levels(results_df$condition),
  timepoint = 1,
  aggregates_name = c(
    "Antigen Presentation",
    "Inflammatory/TLR/Chemokines",
    "Interferon/Antiviral Sensing",
    "Monocytes",
    "DC Activation",
    "Neutrophils",
    "NK Cells",
    "Signal Transduction",
    "ECM And Migration",
    "Energy Metabolism",
    "Cell Cycle",
    "Platelets",
    "T Cells",
    "B Cells",
    "Plasma Cells"
  ),
  p_correction = "BH",
  p_approach = "global",
  p_threshold = 0.05,
  filter_variable = "none",
  filter_mode = "user",
  user_threshold = 0.5,
  quantile_threshold = 0.5,
  scores = "fc.score",
  correlation = "corr.mean",
  arc = "positive",
  ring = "expression",
  order = "set_all",
  quantile_scoreclip = 0.995,
  legend = FALSE,
  plot_title = "dearseq",
  title_size = 6,
  title_line = -4
)

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1)) # empty plot
draw_circos_legend(
  aggregates_name = c(
    "Antigen Presentation",
    "Inflammatory/TLR/Chemokines",
    "Interferon/Antiviral Sensing",
    "Monocytes",
    "DC Activation",
    "Neutrophils",
    "NK Cells",
    "Signal Transduction",
    "ECM And Migration",
    "Energy Metabolism",
    "Cell Cycle",
    "Platelets",
    "T Cells",
    "B Cells",
    "Plasma Cells"
  ),
  ring = "expression",
  placeholder = TRUE,
  scale = 1.3
)
# END DAY 1 #

# DAY 3 #
plot_row_annotation("Day 3", size = 6, placeholder = FALSE)

circos.clear()# make sure canvas is clean

plot_circos(
  method_name = "qusage",
  conditions = levels(results_df$condition),
  timepoint = 3,
  aggregates_name = c(
    "Antigen Presentation",
    "Inflammatory/TLR/Chemokines",
    "Interferon/Antiviral Sensing",
    "Monocytes",
    "DC Activation",
    "Neutrophils",
    "NK Cells",
    "Signal Transduction",
    "ECM And Migration",
    "Energy Metabolism",
    "Cell Cycle",
    "Platelets",
    "T Cells",
    "B Cells",
    "Plasma Cells"
  ),
  p_correction = "BH",
  p_approach = "global",
  p_threshold = 0.05,
  filter_variable = "none",
  filter_mode = "user",
  user_threshold = 0.5,
  quantile_threshold = 0.5,
  scores = "fc.score",
  correlation = "corr.mean",
  arc = "positive",
  ring = "expression",
  order = "set_all",
  quantile_scoreclip = 0.995,
  legend = FALSE,
  plot_title = NULL,
  title_size = 6,
  title_line = -4
)

circos.clear()
plot_circos(
  method_name = "dearseq",
  conditions = levels(results_df$condition),
  timepoint = 3,
  aggregates_name = c(
    "Antigen Presentation",
    "Inflammatory/TLR/Chemokines",
    "Interferon/Antiviral Sensing",
    "Monocytes",
    "DC Activation",
    "Neutrophils",
    "NK Cells",
    "Signal Transduction",
    "ECM And Migration",
    "Energy Metabolism",
    "Cell Cycle",
    "Platelets",
    "T Cells",
    "B Cells",
    "Plasma Cells"
  ),
  p_correction = "BH",
  p_approach = "global",
  p_threshold = 0.05,
  filter_variable = "none",
  filter_mode = "user",
  user_threshold = 0.5,
  quantile_threshold = 0.5,
  scores = "fc.score",
  correlation = "corr.mean",
  arc = "positive",
  ring = "expression",
  order = "set_all",
  quantile_scoreclip = 0.995,
  legend = FALSE,
  plot_title = NULL,
  title_size = 6,
  title_line = -4
)

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1)) # empty plot
draw_circos_legend(
  aggregates_name = c(
    "Antigen Presentation",
    "Inflammatory/TLR/Chemokines",
    "Interferon/Antiviral Sensing",
    "Monocytes",
    "DC Activation",
    "Neutrophils",
    "NK Cells",
    "Signal Transduction",
    "ECM And Migration",
    "Energy Metabolism",
    "Cell Cycle",
    "Platelets",
    "T Cells",
    "B Cells",
    "Plasma Cells"
  ),
  ring = "expression",
  placeholder = FALSE,
  scale = 1.3
)
# END DAY 3 #

# DAY 1 #
plot_row_annotation("Day 7", size = 6, placeholder = FALSE)

circos.clear()# make sure canvas is clean

plot_circos(
  method_name = "qusage",
  conditions = levels(results_df$condition),
  timepoint = 7,
  aggregates_name = c(
    "Antigen Presentation",
    "Inflammatory/TLR/Chemokines",
    "Interferon/Antiviral Sensing",
    "Monocytes",
    "DC Activation",
    "Neutrophils",
    "NK Cells",
    "Signal Transduction",
    "ECM And Migration",
    "Energy Metabolism",
    "Cell Cycle",
    "Platelets",
    "T Cells",
    "B Cells",
    "Plasma Cells"
  ),
  p_correction = "BH",
  p_approach = "global",
  p_threshold = 0.05,
  filter_variable = "none",
  filter_mode = "user",
  user_threshold = 0.5,
  quantile_threshold = 0.5,
  scores = "fc.score",
  correlation = "corr.mean",
  arc = "positive",
  ring = "expression",
  order = "set_all",
  quantile_scoreclip = 0.995,
  legend = FALSE,
  plot_title = NULL,
  title_size = 6,
  title_line = -4
)

circos.clear()
plot_circos(
  method_name = "dearseq",
  conditions = levels(results_df$condition),
  timepoint = 7,
  aggregates_name = c(
    "Antigen Presentation",
    "Inflammatory/TLR/Chemokines",
    "Interferon/Antiviral Sensing",
    "Monocytes",
    "DC Activation",
    "Neutrophils",
    "NK Cells",
    "Signal Transduction",
    "ECM And Migration",
    "Energy Metabolism",
    "Cell Cycle",
    "Platelets",
    "T Cells",
    "B Cells",
    "Plasma Cells"
  ),
  p_correction = "BH",
  p_approach = "global",
  p_threshold = 0.05,
  filter_variable = "none",
  filter_mode = "user",
  user_threshold = 0.5,
  quantile_threshold = 0.5,
  scores = "fc.score",
  correlation = "corr.mean",
  arc = "positive",
  ring = "expression",
  order = "set_all",
  quantile_scoreclip = 0.995,
  legend = FALSE,
  plot_title = NULL,
  title_size = 6,
  title_line = -4
)

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1)) # empty plot
draw_circos_legend(
  aggregates_name = c(
    "Antigen Presentation",
    "Inflammatory/TLR/Chemokines",
    "Interferon/Antiviral Sensing",
    "Monocytes",
    "DC Activation",
    "Neutrophils",
    "NK Cells",
    "Signal Transduction",
    "ECM And Migration",
    "Energy Metabolism",
    "Cell Cycle",
    "Platelets",
    "T Cells",
    "B Cells",
    "Plasma Cells"
  ),
  ring = "expression",
  placeholder = TRUE,
  scale = 1.3
)
# END DAY 7 #

dev.off()
