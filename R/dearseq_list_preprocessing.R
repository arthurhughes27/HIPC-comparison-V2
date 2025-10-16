#'Pre-processing of many gene-set differential analyses for comparison
#'
#'A function which takes lists of results from multiple gene-set analyses with
#'the \code{dgsa_seq()} function and preprocesses it for visualisation
#'
#'@param raw_pvals_list a list of raw-values from a \code{dgsa_seq()} analysis.
#' Each element of the list is the output of one differential analysis i.e. a vector of raw p-values
#' where each entry of the vector corresponds to one geneset.
#'
#'@param scores_list a list of outputs from the \code{compute_scores()} function.
#' Each list entry should be a list with at least two elements \code{activation.scores} and \code{fc.scores}.
#' The ordering and names of the elements should match exactly that of the \code{dgsa_seq_list} argument.
#'
#'@param correlations_list optional list giving correlations between each geneset and a given downstream
#' response to the condition under investigation, given by the \code{calculate_gs_correlation()} function.
#' Each list entry should be a list with at least two elements \code{mean.corr} and \code{corr.mean}.
#'
#'@param conditions a vector of strings giving the active condition being investigated in each element.
#' If not given, assume all elements are from distinct conditions and will give generic names to each.
#'
#'@param condition_colors an optional vector of strings containing hex codes corresponding to each
#'of the \code{conditions}. If not given, will auto-generate colours for each condition with
#'\code{Rcolorbrewer}.
#'
#'@param times a vector of numeric times (in days post-perturbation) corresponding to each element of
#'the results list. If not given, assume all elements are from the same timepoint.
#'
#'@param GSA a \code{GSA} object containing the following named elements \itemize{
#' \item \code{genesets} : a list where each element gives the names of the genes within a set
#' \item \code{geneset.names} : a vector of strings containing the names of the genesets
#' \item \code{geneset.descriptions} : a vector of strings containing descriptions for the genesets
#' \item \code{geneset.aggregates} : a vector of strings containing aggregate groups for the genesets.
#' If this is given as a factor variable, the ordering will be preserved for visualisation. Else,
#' the default ordering will be alphabetical.
#'}
#'
#'@param aggregate_colors an optional vector of strings containing hex codes corresponding to each of
#'the \code{geneset.aggregates}. If given, the ordering that these colours will be assigned is given
#'by the factor ordering of \code{geneset.aggregates}. If the factor ordering is not given, this
#'defaults to alphabetical. If not given, colours will be assigned automatically for each aggregate
#'with \code{Rcolorbrewer}.
#'
#'@returns comparison_dataframe : a dataframe of results ready to be passed for visualisation.
#'
#'@import dplyr, Rcolorbrewer, purrr
#'@author Arthur Hughes
#'

dearseq_list_preprocessing = function(raw_pvals_list,
                                      scores_list,
                                      correlations_list = NULL,
                                      conditions = NULL,
                                      condition_colors = NULL,
                                      times = NULL,
                                      GSA,
                                      aggregate_colors = NULL) {
  n_gene_sets = length(raw_pvals_list[[1]])
  n_comparisons = length(raw_pvals_list)
  
  ### INPUT TRANSFORMATION ###
  # If no conditions are specified, assume each element is a disinct condition and assign generic names
  if (is.null(conditions)) {
    conditions = paste0('Condition ', 1:n_comparisons)
  }
  
  # If no times are specified, assume each element is from the same time
  if (is.null(times)) {
    times = rep("NA", n_comparisons)
  }
  
  condition_times = paste0(conditions, " - Day ", times)
  
  ### VALIDITY CHECKS ###
  # Make sure everything has the same number of elements
  
  if (length(scores_list[[1]][["activation.scores"]]) != n_gene_sets |
      length(correlations_list[[1]][["mean.corr"]]) != n_gene_sets |
      length(GSA[["geneset.names"]]) != n_gene_sets) {
    stop("Number of genesets does not match across inputs!")
  }
  
  if (length(conditions) != n_comparisons |
      length(times) != n_comparisons |
      length(scores_list) != n_comparisons |
      length(correlations_list) != n_comparisons) {
    stop("Number of comparisons does not match across inputs!")
  }
  
  # Transform GSA object into geneset metadata
  gs_metadata <- list(
    names       = GSA$geneset.names,
    descriptions = GSA$geneset.descriptions,
    aggregates  = GSA$geneset.aggregates,
    sets        = GSA$genesets
  )
  
  # Define a colour for each distinct condition
  # If conditions is given as a factor, conserve the ordering
  # Else, order alphabetically
  if (is.factor(conditions)) {
    unique_conditions <- levels(conditions)
  } else {
    unique_conditions <- sort(unique(conditions))
    conditions = factor(conditions, levels = unique_conditions)
  }
  
  n_conditions = length(unique_conditions)
  
  # Assign colours to conditions
  ## If colours given, use these
  if (!is.null(condition_colors)) {
    condition_colours = setNames(condition_colors, unique_conditions)
  } else {
    # Else generate with Rcolorbrewer
    colours <-
      brewer.pal(min(n_conditions, brewer.pal.info["Paired", "maxcolors"]), "Paired")
    
    # If more unique conditions than available colors, you can interpolate
    if (n_conditions > length(colours)) {
      colours <- colorRampPalette(colours)(n_conditions)
    }
    
    # Assign colours to condition names
    condition_colours <- setNames(colours, unique_conditions)
  }
  
  # Now define colours for geneset aggregates
  # If genesets is given as a factor, conserve the ordering
  # Else, order alphabetically
  if (is.factor(gs_metadata$aggregates)) {
    unique_aggregates <- levels(gs_metadata$aggregates)
  } else {
    unique_aggregates <- sort(unique(gs_metadata$aggregates))
    
    gs_metadata$aggregates = factor(gs_metadata$aggregates, levels = unique_aggregates)
  }
  n_aggregates = length(unique_aggregates)
  
  # If colours are provided, use them
  if (!is.null(aggregate_colors)) {
    aggregate_colours <- setNames(aggregate_colors, unique_aggregates)
  } else {
    # Else generate with Rcolorbrewer
    colours <-
      brewer.pal(min(n_aggregates, brewer.pal.info["Spectral", "maxcolors"]), "Spectral")
    
    # If more unique aggregates than available colors, you can interpolate
    if (n_aggregates > length(colours)) {
      colours <- colorRampPalette(colours)(n_aggregates)
    }
    # Assign colours to condition names
    aggregate_colours <- setNames(colours, unique_aggregates)
  }
  
  
  # Assemble results into a single data frame
  if (is.null(correlations_list)) {
    results_df <-
      purrr::map_dfr(seq_along(results_list), function(idx) {
        # Build base tibble
        tibble(
          comparison          = condition_times[idx],
          condition           = conditions[idx],
          time                = times[idx],
          condition.colour    = condition_colours[conditions][idx],
          gs.name             = gs_metadata$names,
          gs.description      = gs_metadata$descriptions,
          gs.name.description = paste0(gs_metadata$names, " - ", gs_metadata$descriptions),
          gs.aggregate        = gs_metadata$aggregates,
          gs.colour           = aggregate_colours[gs_metadata$aggregates],
          activation.score    = scores_list[[idx]][["activation.scores"]] %>% unlist(),
          fc.score            = scores_list[[idx]][["fc.scores"]] %>% unlist(),
          corr.mean           = NA_real_,
          mean.corr           = NA_real_,
          rawPval             = raw_pvals_list[[idx]]
        )
      })
  } else {
    results_df <-
      purrr::map_dfr(seq_along(results_list), function(idx) {
        # Build base tibble
        tibble(
          comparison          = condition_times[idx],
          condition           = conditions[idx],
          time                = times[idx],
          condition.colour    = condition_colours[conditions][idx],
          gs.name             = gs_metadata$names,
          gs.description      = gs_metadata$descriptions,
          gs.name.description = paste0(gs_metadata$names, " - ", gs_metadata$descriptions),
          gs.aggregate        = gs_metadata$aggregates,
          gs.colour           = aggregate_colours[gs_metadata$aggregates],
          activation.score    = scores_list[[idx]][["activation.scores"]] %>% unlist(),
          fc.score            = scores_list[[idx]][["fc.scores"]] %>% unlist(),
          corr.mean           = correlations_list[[idx]][["corr.mean"]],
          mean.corr           = correlations_list[[idx]][["mean.corr"]],
          rawPval             = raw_pvals_list[[idx]]
        )
      })
  }
  
  # Now perform different p-value correction methods/approaches
  correction_methods <-
    c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY")
  for (method in correction_methods) {
    # Global correction across all gene sets
    results_df <- results_df %>%
      mutate(!!paste0("global.adjPval_", method) := p.adjust(rawPval, method = method))
    
    # Timepoint-wise correction
    results_df <- results_df %>%
      group_by(time) %>%
      mutate(!!paste0("withinTime.adjPval_", method) := p.adjust(rawPval, method = method)) %>%
      ungroup()
    
    # Comparison-wise correction
    results_df <- results_df %>%
      group_by(comparison) %>%
      mutate(
        !!paste0("withinComparison.adjPval_", method) := p.adjust(rawPval, method = method)
      ) %>%
      ungroup()
  }
  
  
  return(results_df)
}