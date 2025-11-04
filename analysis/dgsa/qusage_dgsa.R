# File to perform differential analysis using the code provided by Hagan et al. (2022), using QuSAGE meta-analysis.

# Libraries

suppressPackageStartupMessages({
  library(package = "Biobase")
  library(package = "tidyverse")
  library(package = "qusage")
  library(package = "boot")
  library(package = "ggthemes")
  library(package = "pheatmap")
  library(package = "matrixStats")
  library(package = "RColorBrewer")
  library(package = "reshape2")
  library(package = "lme4")
  #library(package = "metap")
  library(package = "circlize")
  library(package = "corrplot")
  library(package = "ComplexHeatmap")
  library(package = "ggpubr")
  library(package = "GSA")
  library(package = "stringr")
  library(package = "purrr")
  library(package = "data.table")
  library(package = "Biobase")
  library(package = "GSVA")
  library(package = "limma")
  library(package = "caret")
  library(package = "caretEnsemble")
  library(package = "WeightedROC")
  library(package = "nlme")
  library(package = "emmeans")
  library(package = "plotROC")
})

# Define some helper functions
# For every gene set, find out in how many vaccines it is significantly up at
# any one postvax time point.
# The same but for down at any postvax time point
calc_sharing_df <- function(x) {
  x %>%
    dplyr::mutate(
      direction = ifelse(
        # Directionality is defined as positive of negative activity score
        test = activity_score > 0,
        yes = "up",
        no = "down"
      ),
      is_significant = !is.na(FDR) &
        FDR < 0.05
    ) %>% # significance threshold is always 0.05 with B-H correction
    dplyr::select(pathogen, vaccine_type, geneset, direction, is_significant) %>%
    dplyr::distinct() %>%
    dplyr::group_by(geneset, direction) %>% # They only consider shared DE if same direction
    dplyr::summarise(n = sum(is_significant)) %>%
    dplyr::group_by(geneset) %>%
    dplyr::arrange(-n) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
}

# Qusage (version 2.18.0) runs very slow on certain datasets, because it uses
# the convolve() function from stats which itself uses the fft() function. The
# fft() function runs very slow on certain inputs (e.g. 3 minutes instead of
# under 1 second with other fft implementations for the same input).
#
# Found a report of speed issues with the convolve() function at:
# https://stat.ethz.ch/pipermail/r-help/2007-December/148777.html
#
# Found a solution for the slow convolve() (actually the fft() it uses) at:
# https://stat.ethz.ch/pipermail/r-help/2007-December/148779.html
#
# The solution below requires the fftw library, which on ubuntu required the
# installation of the fftw-dev package (sudo apt install fftw3-dev) before
# install.packages("fftw") was successful.
#
# Here we override the convolve function to use a different (faster) fft
# function.
convolve <- function (x,
                      y,
                      conj = TRUE,
                      type = c("circular", "open", "filter")) {
  type <- match.arg(type)
  n <- length(x)
  ny <- length(y)
  Real <- is.numeric(x) && is.numeric(y)
  if (type == "circular") {
    if (ny != n)
      stop("length mismatch in convolution")
  }
  else {
    n1 <- ny - 1
    x <- c(rep.int(0, n1), x)
    n <- length(y <- c(y, rep.int(0, n - 1)))
  }
  # Original fft lines:
  #x <- fft(fft(x) * (if (conj)
  #  Conj(fft(y))
  #  else fft(y)), inverse = TRUE)
  
  # Use fftw instead of fft
  x <- fftw::FFT(fftw::FFT(x) * (if (conj)
    Conj(fftw::FFT(y))
    else
      fftw::FFT(y)), inverse = TRUE)
  
  if (type == "filter")
    (if (Real)
      Re(x)
     else
       x)[-c(1L:n1, (n - n1 + 1L):n)] / n
  else
    (if (Real)
      Re(x)
     else
       x) / n
}

# Helper function for DE thresholds
create_signif_label <- function(x) {
  cut(x, c(-Inf, .001, .01, .05, Inf), labels = c("***", "**", "*", "ns"))
}

#' Keep only specified symbols in a given list of gene sets
#'
#' @param genesets.symbol list of character vectors
#' @param symbols2keep character vector
#' @param min_size minimum number of symbols to keep a gene set
keep_symbols <- function(genesets.symbol, symbols2keep, min_size = 2) {
  stopifnot(is.character(symbols2keep))
  stopifnot(is.list(genesets.symbol))
  common_symbols <- intersect(Reduce(union, genesets.symbol), symbols2keep)
  purrr::map(.x = genesets.symbol, ~ .x[.x %in% common_symbols]) %>%
    purrr::discard( ~ length(.x) < min_size)
}

# NOTE: function below was copied from
# https://github.com/RGLab/ImmuneSignatures2/blob/master/R/geneExpressionPreProcessing.R
#
#' Remove incomplete rows
#'
#' @param eset expressionSet
#' @export
#'
removeAllNArows <- function(eset) {
  em <- Biobase::exprs(eset)
  allNARows <- apply(em, 1, function(x) {
    all(is.na(x))
  })
  eset <- eset[!allNARows, ]
}
remove_any_NA_rows <- function(eset) {
  any_NA_rows <- apply(
    X = Biobase::exprs(eset),
    MARGIN = 1,
    FUN = function(x)
      any(is.na(x))
  )
  eset[!any_NA_rows, ]
}

#' Keep most recent pre-vaccination time point for each participant
#'
#' @param eset expressionSet
#' @param drop_postvax if TRUE, postvax time points will be discarded
#' @export
#'
keep_most_recent_prevax_time_point <- function(eset, drop_postvax) {
  prevax_uids <- Biobase::pData(eset) %>%
    dplyr::filter(time_post_last_vax <= 0) %>%
    dplyr::arrange(-time_post_last_vax) %>%
    dplyr::group_by(participant_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::pull(uid)
  postvax_uids <- Biobase::pData(eset) %>%
    dplyr::filter(time_post_last_vax > 0) %>%
    dplyr::pull(uid)
  if (drop_postvax)
    uids <- prevax_uids
  else
    uids <- union(prevax_uids, postvax_uids)
  eset[, eset$uid %in% uids]
}

#' Add 'response_class' column containing maxRBA_p30 values for Influenza and
#' MFC_p30 values for all other vaccines
#'
#' @param eset expressionSet
#' @export
#'
add_response_class_column <- function(eset) {
  # Note: maxRBA measure is used to determine high and low responders for
  # Influenza, because many participants are unlikely to be naive for this virus
  # (resulting in high baseline titers and consequently lower fold changes which
  # maxRBA is designed to take into account)
  Biobase::pData(eset)$response_class <-
    as.character(Biobase::pData(eset)$MFC_p30)
  Biobase::pData(eset)$response_class[Biobase::pData(eset)$pathogen == "Influenza"] <-
    as.character(Biobase::pData(eset)$maxRBA_p30[Biobase::pData(eset)$pathogen == "Influenza"])
  eset
}
#' Keep participants that are high- or low-responders
#'
#' @param eset expressionSet
#' @param col name of column containing high or low responders
#' @export
#'
keep_only_high_and_low_responders <- function(eset, col) {
  eset[, eset[[col]] %in% c("lowResponder", "highResponder")]
}
#' Perform qusage meta analysis (combinePDFs) on data frame
#'
#' @param x data frame with list column containing QSarray objects
#' @param group_cols character vector with names of columns to perform grouping
#' on for meta analysis
#' @return x with additional column containing output from meta analysis
meta_qusage <- function(x, group_cols) {
  require(qusage)
  if (find("convolve")[1] == "package:stats") {
    stop("stats::convolve is slow, please load convolve.R")
  }
  
  x %>%
    dplyr::group_by_at(group_cols) %>%
    tidyr::nest() %>%
    dplyr::ungroup() %>%
    dplyr::mutate(meta_qusage_output = purrr::map(.x = data, ~ {
      QSarray_list <- .x$qusage_output %>%
        purrr::keep( ~ class(.x) == "QSarray")
      if (length(QSarray_list) == 0)
        return(NULL)
      if (length(QSarray_list) > 1) {
        cat("Running combinePDFs\n")
        return(qusage::combinePDFs(QSarrayList = QSarray_list))
      } else {
        return(QSarray_list[[1]])
      }
    })) %>%
    dplyr::select(-data)
}
#' Replace QSarray objects with a few informative columns
#'
#' @param x data frame containing list column with QSarray objects
#' @param col string with name of list column
#' @return x with col replaced with "p.value", "geneset", and "activity_score"
#' columns
tidy_qusage_output <- function(x, col) {
  require(qusage)
  x %>%
    dplyr::filter(purrr::map_lgl(.x = !!sym(col), ~ class(.x) == "QSarray")) %>%
    dplyr::mutate(tidied =
                    purrr::map(
                      .x = !!sym(col),
                      ~ tibble(
                        p.value = pdf.pVal(.x),
                        geneset = colnames(.x$path.PDF),
                        activity_score = .x$path.mean
                      )
                    )) %>%
    dplyr::select(-!!sym(col)) %>%
    tidyr::unnest(tidied)
}
#' Get confidence intervals for every variable in a bootstrap.
#'
#' @param boot_result an object of class "boot" containing the output of a
#' bootstrap calculation
#' @param type character string representing the type of interval required
tidy_boot_ci <- function(boot_result, type) {
  stopifnot(class(boot_result) == "boot")
  stopifnot(is.character(type))
  purrr::map_df(.x = 1:ncol(boot_result$t), ~ {
    ci <- boot.ci(boot.out = boot_result,
                  type = type,
                  index = .x)
    # When using "perc" as type, the name in the output list ci is "percent"
    # Using grep to handle this inconsistency
    i <- grep(pattern = type, names(ci))[1]
    tibble(
      name = names(ci$t0),
      statistic = ci$t0,
      ci_low = ci[[i]][4],
      ci_high = ci[[i]][5]
    )
  }) %>%
    dplyr::mutate(bootstrapped_mean_of_statistic = colMeans(boot_result$t))
}
# Description: functions to extend pheatmap functionality
pheatmap_add_labels <- function(p,
                                annot_labels,
                                gp = gpar(col = "#000000", fontsize = 10),
                                angle_col = 0,
                                angle_row = 270) {
  require(pheatmap)
  require(grid)
  require(gtable)
  # Get the column and row annotation rectangles (called "track" in pheatmap)
  tracks <- list()
  annot_names <- vector()
  for (annot_dir in c("col", "row")) {
    track_grob_name <- paste0(annot_dir, "_annotation")
    track_index <- which(p$gtable$layout$name == track_grob_name)
    
    if (length(track_index) == 0) {
      next
    }
    track_grob <- p$gtable$grobs[[track_index]]
    tracks[[annot_dir]] <- list("index" = track_index, "grob" = track_grob)
    annot_names <- c(annot_names, colnames(track_grob$gp$fill))
  }
  
  # Get legend grob
  legend_index <- which(p$gtable$layout$name == "annotation_legend")
  legend_grob <- p$gtable$grobs[[legend_index]]
  
  # Get mapping between annotation labels and annotation colors
  elem2color <- list()
  color2label <- list()
  for (annot_name in annot_names) {
    # Retrieve from the legend the filled rectangles grob belonging to current annotation
    cur_grob <- getGrob(gTree = legend_grob, gPath = childNames(legend_grob)[[paste(annot_name, "r")]])
    elem2color[[annot_name]] <- cur_grob$gp$fill
    cur_colors <- as.character(elem2color[[annot_name]])
    cur_elems <- names(elem2color[[annot_name]])
    if (annot_name %in% names(annot_labels)) {
      cur_labels <- as.character(annot_labels[[annot_name]][cur_elems])
    } else {
      cur_labels <- rep("", length.out = length(cur_elems))
    }
    color2label[[annot_name]] <- cur_labels
    names(color2label[[annot_name]]) <- cur_colors
  }
  
  # Add labels to annotation legend
  for (annot_name in annot_names) {
    target_child_name <- paste(annot_name, "r")
    cur_grob <- getGrob(gTree = legend_grob, gPath = childNames(legend_grob)[[target_child_name]])
    legend_labels <- color2label[[annot_name]][as.character(cur_grob$gp$fill)]
    
    res <- gList()
    res[["rect"]] <- cur_grob
    res[["text"]] <- textGrob(
      x = cur_grob$x + 0.5 * cur_grob$width,
      y = cur_grob$y - 0.5 * cur_grob$height,
      label = legend_labels,
      gp = gp
    )
    legend_grob$children[[target_child_name]] <- gTree(children = res)
  }
  p$gtable$grobs[[legend_index]] <- gTree(children = legend_grob$children)
  
  # Add labels to annotation tracks
  for (track_dir in names(tracks)) {
    track_grob <- tracks[[track_dir]]$grob
    track_index <- tracks[[track_dir]]$index
    
    track_labels <- purrr::map(.x = colnames(track_grob$gp$fill), ~ color2label[[.x]][track_grob$gp$fill[, .x]]) %>%
      do.call(what = c, args = .)
    
    # Add labels to annotation track
    res <- gList()
    res[["rect"]] = track_grob
    res[["text"]] = textGrob(
      x = track_grob$x,
      y = track_grob$y,
      label = track_labels,
      rot = ifelse(track_dir == "col", yes = angle_col, no = angle_row),
      gp = gp
    )
    p$gtable$grobs[[track_index]] <- gTree(children = res)
  }
  return(p)
}
# Function based on (and slightly modified, removing the option to make
# vertically positioned dendrograms)
# https://github.com/cran/pheatmap/blob/master/R/pheatmap.r
# on 2019_02_14
draw_dendrogram <- function(hc) {
  h = hc$height / max(hc$height) / 1.05
  m = hc$merge
  o = hc$order
  n = length(o)
  
  m[m > 0] = n + m[m > 0]
  m[m < 0] = abs(m[m < 0])
  
  dist = matrix(0,
                nrow = 2 * n - 1,
                ncol = 2,
                dimnames = list(NULL, c("x", "y")))
  dist[1:n, 1] = 1 / n / 2 + (1 / n) * (match(1:n, o) - 1)
  
  for (i in 1:nrow(m)) {
    dist[n + i, 1] = (dist[m[i, 1], 1] + dist[m[i, 2], 1]) / 2
    dist[n + i, 2] = h[i]
  }
  
  draw_connection = function(x1, x2, y1, y2, y) {
    res = list(x = c(x1, x1, x2, x2), y = c(y1, y, y, y2))
    
    return(res)
  }
  
  x = rep(NA, nrow(m) * 4)
  y = rep(NA, nrow(m) * 4)
  id = rep(1:nrow(m), rep(4, nrow(m)))
  
  for (i in 1:nrow(m)) {
    c = draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
    k = (i - 1) * 4 + 1
    x[k:(k + 3)] = c$x
    y[k:(k + 3)] = c$y
  }
  
  x = unit(x, "npc")
  y = unit(y, "npc")
  
  return(polylineGrob(x = x, y = y, id = id))
}
# pheatmap variant enabling clustering of groups of columns.
col_grouped_pheatmap <- function(mat,
                                 col_groups,
                                 col_clust_fun = "hclust",
                                 col_dist_fun = "dist",
                                 method = "euclidean",
                                 treeheight_col = 50,
                                 display_numbers = FALSE,
                                 silent = TRUE,
                                 ...) {
  require(pheatmap)
  require(grid)
  require(gtable)
  stopifnot(is.matrix(mat))
  stopifnot(any(class(display_numbers) %in% c("logical", "matrix")))
  
  col_clust_fun <- match.fun(col_clust_fun)
  col_dist_fun <- match.fun(col_dist_fun)
  
  # Cluster the columns and draw the dendrograms
  dendro_grobs <- list()
  for (cur_group in unique(col_groups)) {
    col_indexes.subset <- which(col_groups == cur_group)
    #cl <- hclust(d = dist(x = t(mat[, col_indexes.subset]),
    #                      method = method))
    cl <- col_clust_fun(d = col_dist_fun(x = t(mat[, col_indexes.subset]), method = method))
    
    col_indexes <- 1:ncol(mat)
    col_indexes[col_indexes.subset] <- NA
    col_indexes[is.na(col_indexes)] <- col_indexes.subset[cl$order]
    mat <- mat[, col_indexes]
    
    dendro_grobs[[length(dendro_grobs) + 1]] <- draw_dendrogram(cl)
  }
  
  # Get column numbers after which a gap occurs
  gaps_col <- which(purrr::map_lgl(2:length(col_groups), ~ {
    col_groups[.x] != col_groups[.x - 1]
  }))
  
  gap_width <- unit("4", "bigpts")
  matrix_cell_width <- (1 / ncol(mat)) * (unit(1, "npc") - length(gaps_col) * gap_width)
  dendro_widths <- base::diff(c(0, gaps_col, ncol(mat))) * matrix_cell_width
  
  gap_grob <- rectGrob(gp = gpar(fill = NA, col = NA), width = gap_width)
  
  my_grobs <- list(dendro_grobs[[1]])
  my_widths <- dendro_widths[1]
  for (i in 1:length(gaps_col)) {
    my_grobs[[2 * i]] <- gap_grob
    my_grobs[[2 * i + 1]] <- dendro_grobs[[i + 1]]
    my_widths <- unit.c(my_widths, gap_width, dendro_widths[i + 1])
  }
  
  combined_dendros_grob <- grobTree(
    gtable_matrix(
      name = "my_dendrograms",
      grobs = matrix(my_grobs, ncol = length(my_grobs)),
      widths = my_widths,
      heights = unit(1, "null")
    )
  )
  
  if (is.matrix(display_numbers)) {
    display_numbers <- display_numbers[rownames(mat), colnames(mat)]
  }
  
  p <- pheatmap::pheatmap(
    mat = mat,
    display_numbers = display_numbers,
    gaps_col = gaps_col,
    treeheight_col = treeheight_col,
    cluster_cols = FALSE,
    filename = NA,
    silent = silent,
    ...
  )
  
  p$gtable <- gtable_add_grob(
    x = p$gtable,
    grobs = combined_dendros_grob,
    t = 2,
    l = 3,
    name = "col_tree"
  )
  p$gtable$heights[2] <- p$gtable$heights[2] + unit(treeheight_col, "bigpts")
  p$gtable$heights[4] <- p$gtable$heights[4] - unit(treeheight_col, "bigpts")
  return(p)
}
#Principal Variant Component Analysis (PVCA) function
pvca <- function (abatch, expInfo) {
  theDataMatrix = abatch
  dataRowN <- nrow(theDataMatrix)
  dataColN <- ncol(theDataMatrix)
  
  ########## Center the data (center rows) ##########
  theDataMatrixCentered <- matrix(data = 0,
                                  nrow = dataRowN,
                                  ncol = dataColN)
  theDataMatrixCentered_transposed = apply(theDataMatrix, 1, scale, center = TRUE, scale = FALSE)
  theDataMatrixCentered = t(theDataMatrixCentered_transposed)
  
  ########## Compute correlation matrix &  Obtain eigenvalues ##########
  
  theDataCor <- cor(theDataMatrixCentered)
  eigenData <- eigen(theDataCor)
  eigenValues = eigenData$values
  ev_n <- length(eigenValues)
  eigenVectorsMatrix = eigenData$vectors
  eigenValuesSum = sum(eigenValues)
  percents_PCs = eigenValues / eigenValuesSum
  
  ##===========================================
  ##  Getting the experimental information
  ##===========================================
  
  exp_design <- as.data.frame(expInfo)
  expDesignRowN <- nrow(exp_design)
  expDesignColN <- ncol(exp_design)
  
  ########## Merge experimental file and eigenvectors for n components ##########
  
  ## pc_n is the number of principal components to model
  
  ## Use fixed pc_n
  pc_n = 5
  
  pc_data_matrix <- matrix(data = 0,
                           nrow = (expDesignRowN * pc_n),
                           ncol = 1)
  mycounter = 0
  for (i in 1:pc_n) {
    for (j in 1:expDesignRowN) {
      mycounter <- mycounter + 1
      pc_data_matrix[mycounter, 1] = eigenVectorsMatrix[j, i]
    }
  }
  
  AAA <- exp_design[rep(1:expDesignRowN, pc_n), ]
  Data <- cbind(AAA, pc_data_matrix)
  
  
  ####### Edit these variables according to your factors #######
  
  variables <- c (colnames(exp_design))
  for (i in 1:length(variables))
  {
    Data$variables[i] <- as.factor(Data$variables[i])
  }
  
  
  ########## Mixed linear model ##########
  op <- options(warn = (-1))
  effects_n = expDesignColN + 1 #effects size without interaction terms
  randomEffectsMatrix <- matrix(data = 0,
                                nrow = pc_n,
                                ncol = effects_n)
  
  ##============================#
  ##  Get model functions
  ##============================#
  model.func <- c()
  index <- 1
  
  ##  level-1
  for (i in 1:length(variables))
  {
    mod = paste("(1|", variables[i], ")", sep = "")
    model.func[index] = mod
    index = index + 1
  }
  
  function.mods <- paste (model.func , collapse = " + ")
  
  ##============================#
  ##  Get random effects  #
  ##============================#
  
  for (i in 1:pc_n) {
    y = (((i - 1) * expDesignRowN) + 1)
    funct <- paste ("pc_data_matrix", function.mods, sep = " ~ ")
    Rm1ML <- lmer(
      funct ,
      Data[y:(((i - 1) * expDesignRowN) + expDesignRowN), ],
      REML = TRUE,
      verbose = FALSE,
      na.action = na.omit
    )
    randomEffects <- Rm1ML
    randomEffectsMatrix[i, ] <- c(unlist(VarCorr(Rm1ML)), resid = sigma(Rm1ML)^2)
  }
  effectsNames <- c(names(getME(Rm1ML, "cnms")), "resid")
  
  ########## Standardize Variance ##########
  randomEffectsMatrixStdze <- matrix(data = 0,
                                     nrow = pc_n,
                                     ncol = effects_n)
  for (i in 1:pc_n) {
    mySum = sum(randomEffectsMatrix[i, ])
    for (j in 1:effects_n) {
      randomEffectsMatrixStdze[i, j] = randomEffectsMatrix[i, j] / mySum
    }
  }
  
  ########## Compute Weighted Proportions ##########
  
  randomEffectsMatrixWtProp <- matrix(data = 0,
                                      nrow = pc_n,
                                      ncol = effects_n)
  for (i in 1:pc_n) {
    weight = eigenValues[i] / eigenValuesSum
    for (j in 1:effects_n) {
      randomEffectsMatrixWtProp[i, j] = randomEffectsMatrixStdze[i, j] * weight
    }
  }
  
  ########## Compute Weighted Ave Proportions ##########
  
  randomEffectsSums <- matrix(data = 0,
                              nrow = 1,
                              ncol = effects_n)
  randomEffectsSums <- colSums(randomEffectsMatrixWtProp)
  totalSum = sum(randomEffectsSums)
  randomEffectsMatrixWtAveProp <- matrix(data = 0,
                                         nrow = 1,
                                         ncol = effects_n)
  
  for (j in 1:effects_n) {
    randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j] / totalSum
  }
  return(list(dat = randomEffectsMatrixWtAveProp, label = effectsNames))
}

### Parameter definitions
FDR_THRESHOLD <- 0.05   # Threshold for significant gene set up/down regulation
NUM_PERMUTATIONS <- 10   # Number of permutations for downstream analyses

### Data directory
# Modify this to point to your local download folder
data_dir <- "./data-raw/"
extra_data_dir <- "./output/results/qusage"

### Load normalized expression sets (no normalization, no response)
young.noNorm.noResponse    <- readRDS(file.path(data_dir, "young_noNorm_eset.rds"))
old.noNorm.noResponse      <- readRDS(file.path(data_dir, "old_noNorm_eset.rds"))
young.noNorm.withResponse  <- readRDS(file.path(data_dir, "young_noNorm_withResponse_eset.rds"))
old.noNorm.withResponse    <- readRDS(file.path(data_dir, "old_noNorm_withResponse_eset.rds"))

### Merge pre/post-vax samples for SDY1529 into same matrix annotation
# This ensures both pre- and post-vaccine samples are labeled consistently
young.noNorm.noResponse$matrix[young.noNorm.noResponse$study_accession == 'SDY1529']   <- 'SDY1529_WholeBlood_HealthyAdults'
young.noNorm.withResponse$matrix[young.noNorm.withResponse$study_accession == 'SDY1529'] <- 'SDY1529_WholeBlood_HealthyAdults'

### Add adjuvant information for Hepatitis A/B vaccines
# All Hep A/B studies assumed to use Alum adjuvant
datasets <- list(young.noNorm.noResponse, old.noNorm.noResponse,
                 young.noNorm.withResponse, old.noNorm.withResponse)
for (i in seq_along(datasets)) {
  ds <- datasets[[i]]
  ds$adjuvant[ds$pathogen == 'Hepatitis A/B'] <- 'Alum'
  datasets[[i]] <- ds
}
# Restore updated adjuvants
young.noNorm.noResponse$adjuvant    <- datasets[[1]]$adjuvant
old.noNorm.noResponse$adjuvant      <- datasets[[2]]$adjuvant
young.noNorm.withResponse$adjuvant  <- datasets[[3]]$adjuvant
old.noNorm.withResponse$adjuvant    <- datasets[[4]]$adjuvant

### Abbreviate vaccine types consistently
# Keep full names in _vaccine_type_full, then recode to abbreviations
vaccine_abrev <- c(
  'Live virus'                    = 'LV',
  'Recombinant protein'           = 'RP',
  'Recombinant Viral Vector'      = 'RVV',
  'Polysaccharide'                = 'PS',
  'Conjugate'                     = 'CJ',
  'Inactivated'                   = 'IN',
  'Inactivated/Recombinant protein' = 'IN/RP'
)

# Backup full names
young.noNorm.noResponse$vaccine_type_full   <- young.noNorm.noResponse$vaccine_type
old.noNorm.noResponse$vaccine_type_full     <- old.noNorm.noResponse$vaccine_type
young.noNorm.withResponse$vaccine_type_full <- young.noNorm.withResponse$vaccine_type
old.noNorm.withResponse$vaccine_type_full   <- old.noNorm.withResponse$vaccine_type

# Apply abbreviations
young.noNorm.noResponse$vaccine_type   <- str_replace_all(young.noNorm.noResponse$vaccine_type,   vaccine_abrev)
old.noNorm.noResponse$vaccine_type     <- str_replace_all(old.noNorm.noResponse$vaccine_type,     vaccine_abrev)
young.noNorm.withResponse$vaccine_type <- str_replace_all(young.noNorm.withResponse$vaccine_type, vaccine_abrev)
old.noNorm.withResponse$vaccine_type   <- str_replace_all(old.noNorm.withResponse$vaccine_type,   vaccine_abrev)

# Colors and abbreviations
#load(file.path(data_dir, "adj_path_vt_colors_abbv_v3.RData"))
# This file has not been made available, so we will have to reconstruct it based on later code. I guess it is a list containing vaccine_type, pathogen and adjuvant each associated with a hex code.
meta_levels <- young.noNorm.noResponse@phenoData@data %>%
  select(vaccine_type, pathogen, adjuvant) %>%
  distinct()

vt_levels <- unique(meta_levels$vaccine_type)
pt_levels <- unique(meta_levels$pathogen)
ad_levels <- unique(meta_levels$adjuvant)

# Generate palettes matching the number of unique levels
cols_vt <- colorRampPalette(brewer.pal(8, "Paired"))(length(vt_levels))
cols_pt <- colorRampPalette(brewer.pal(8, "Set2"))(length(pt_levels))
cols_ad <- colorRampPalette(brewer.pal(8, "Dark2"))(length(ad_levels))

# Build the named‐vector mapping (duplicates now collapse to the same code)
adj_path_vt_colors <- list(
  vaccine_type = setNames(cols_vt, vt_levels),
  pathogen     = setNames(cols_pt, pt_levels),
  adjuvant     = setNames(cols_ad, ad_levels)
)


# adj_path_vt_abbv = adj_path_vt_colors
# names(adj_path_vt_abbv$vaccine_type)=str_replace_all(names(adj_path_vt_abbv$vaccine_type), vaccine_abrev)

names(adj_path_vt_colors$vaccine_type)=str_replace_all(names(adj_path_vt_colors$vaccine_type), vaccine_abrev)


### Combine young and old cohorts (with response)
# Find common genes between both datasets
symbols <- intersect(rownames(young.noNorm.withResponse),
                     rownames(old.noNorm.withResponse))

# Build phenotype data frame with age_group annotation
df_young <- Biobase::pData(young.noNorm.withResponse) %>%
  dplyr::mutate(age_group = 'young')
df_old   <- Biobase::pData(old.noNorm.withResponse) %>%
  dplyr::mutate(age_group = 'old')
pd <- rbind(df_young, df_old)

# Combine expression matrices and align to phenotype order
ge <- cbind(
  Biobase::exprs(young.noNorm.withResponse)[symbols, ],
  Biobase::exprs(old.noNorm.withResponse)[symbols, ]
)[, pd$uid]

# Set sample IDs as rownames for phenotype data
rownames(pd) <- pd$uid

# Create merged ExpressionSet object
young_and_old.noNorm.withResponse <- new(
  'ExpressionSet',
  exprs     = ge,
  phenoData = new('AnnotatedDataFrame', pd)
)

### Combine young and old cohorts (no response)
symbols <- intersect(rownames(young.noNorm.noResponse),
                     rownames(old.noNorm.noResponse))
df_young <- Biobase::pData(young.noNorm.noResponse) %>% dplyr::mutate(age_group = 'young')
df_old   <- Biobase::pData(old.noNorm.noResponse)   %>% dplyr::mutate(age_group = 'old')
pd <- rbind(df_young, df_old)

ge <- cbind(
  Biobase::exprs(young.noNorm.noResponse)[symbols, ],
  Biobase::exprs(old.noNorm.noResponse)[symbols, ]
)[, pd$uid]
rownames(pd) <- pd$uid
young_and_old.noNorm <- new(
  'ExpressionSet',
  exprs     = ge,
  phenoData = new('AnnotatedDataFrame', pd)
)

### Load Blood Transcription Module (BTM) gene sets via GSA
# Have to guess its structure
GSA <- GSA.read.gmt(file.path(data_dir, 'BTM.gmt'))
btm_list <- GSA[['genesets']]
names(btm_list) <- GSA[['geneset.names']]

### Prepare fold-change dataset from normalized eset
eset <- readRDS(file.path(data_dir, 'young_norm_eset.rds'))
colnames(eset) <- eset$uid

# Same pre/post-vax merge and adjuvant addition
eset$matrix[eset$study_accession == 'SDY1529']  <- 'SDY1529_WholeBlood_HealthyAdults'
eset$adjuvant[eset$pathogen == 'Hepatitis A/B'] <- 'Alum'

# Abbreviate vaccine types in fold-change eset
eset$vaccine_type_full <- eset$vaccine_type
eset$vaccine_type      <- str_replace_all(eset$vaccine_type, vaccine_abrev)

# Combine pathogen and vaccine type labels
eset$pt    <- paste(eset$pathogen, ' (', eset$vaccine_type, ')', sep='')
eset$SDY_pt <- paste(eset$study, eset$pt)

# Identify and retain timepoints >= 0
tp_int <- sort(unique(eset$time_post_last_vax[eset$time_post_last_vax >= 0]))
ind    <- lapply(tp_int, function(x) which(eset$time_post_last_vax == x))
ind_all <- Reduce(union, ind)
eset    <- eset[, ind_all]

# Recompute indices after filtering
ind <- lapply(tp_int, function(x) which(eset$time_post_last_vax == x))

# Remove studies with fewer than sample_cutoff samples at any timepoint
sample_cutoff <- 3
matrix_uni_tp <- lapply(ind, function(x) unique(eset[, x]$matrix))
matrix_ind <- lapply(seq_along(ind), function(i) {
  lapply(matrix_uni_tp[[i]], function(m) which(eset[, ind[[i]]]$matrix == m))
})

# Identify samples to drop
ind_cut_all <- c()
for (i in seq_along(matrix_ind)) {
  small_sets <- which(sapply(matrix_ind[[i]], length) < sample_cutoff)
  if (length(small_sets) > 0) {
    for (j in small_sets) {
      ind_cut_all <- c(ind_cut_all, ind[[i]][matrix_ind[[i]][[j]]])
    }
  }
}
if (length(ind_cut_all) > 0) {
  eset <- eset[, -ind_cut_all]
}

# Remove matrices with only Day 0 samples
matrix_uni <- unique(eset$matrix)
matrix_tp  <- lapply(matrix_uni, function(m) unique(eset$time_post_last_vax[eset$matrix == m]))
matrix_d0_only <- matrix_uni[sapply(matrix_tp, function(x) identical(x, 0))]
eset <- eset[, !(eset$matrix %in% matrix_d0_only)]

# Remove genes with NA expression values
eset <- eset[complete.cases(exprs(eset)), ]

# Recompute timepoints after all filtering
tp_int <- sort(unique(eset$time_post_last_vax[eset$time_post_last_vax >= 0]))
ind    <- lapply(tp_int, function(x) which(eset$time_post_last_vax == x))

# Filter BTMs to those with matching genes
btm_list <- btm_list[!sapply(btm_list, function(x) is_empty(intersect(x, rownames(eset))))]

# Collapse gene-level to BTM-level by arithmetic mean
eset_BTM <- do.call(rbind, lapply(btm_list, function(geneset) {
  matched <- na.omit(match(geneset, rownames(eset)))
  colMeans(exprs(eset[matched, ]), na.rm = TRUE)
}))
eset_BTM <- ExpressionSet(as.matrix(eset_BTM), eset@phenoData)

### Compute D0-normalized fold change for genes and BTMs
ind_D0  <- which(eset$time_post_last_vax == 0)
common   <- lapply(2:length(ind), function(i) intersect(eset$participant_id[ind[[i]]], eset$participant_id[ind_D0]))
ia       <- lapply(2:length(ind), function(i) na.omit(match(common[[i-1]], eset$participant_id[ind[[i]]])))
ib       <- lapply(2:length(ind), function(i) na.omit(match(common[[i-1]], eset$participant_id[ind_D0])))

# Gene-level fold change objects
# Calculates fold changes by subtracting baseline measurements from post-vaccination measurements for all samples. Same for gene sets but first summarises expression within-gene-set by the mean.
exp_FC     <- lapply(2:length(ind), function(i) eset[, ind[[i]][ia[[i-1]]]])
exp_FC     <- lapply(2:length(ind), function(i) { exprs(exp_FC[[i-1]]) <- exprs(exp_FC[[i-1]]) - exprs(eset[, ind_D0[ib[[i-1]]]]); exp_FC[[i-1]] })

# BTM-level fold change objects
exp_FC_BTM = lapply(2:length(ind),function(x) eset_BTM[,ind[[x]][ia[[x-1]]]])

exp_FC_BTM = lapply(2:length(ind),function(x) {exprs(exp_FC_BTM[[x-1]])=exprs(exp_FC_BTM[[x-1]])-exprs(eset_BTM[,ind_D0[ib[[x-1]]]]); exp_FC_BTM[[x-1]]})

tp_FC=tp_int[-1] #Store timepoints for FC data (remove Day 0)

# Path to store/read the QuSAGE results for young cohort
dgsa_results_directory = "./output/results/qusage"

longitudinal_qusage_young_df_path <- file.path(
  dgsa_results_directory,
  "longitudinal_qusage_young_df.rds"
)

if (file.exists(longitudinal_qusage_young_df_path)) {
  # If already computed, just load it
  longitudinal_qusage_young_df <- readRDS(longitudinal_qusage_young_df_path)
  
} else {
  # Filter to the most recent pre-vax sample per participant (keep post-vax too),
  # then remove any genes/samples with NA rows
  young.noNorm.noResponse.filtered <- young.noNorm.noResponse %>%
    keep_most_recent_prevax_time_point(drop_postvax = FALSE) %>%
    remove_any_NA_rows()
  
  # Ensure no participant has duplicate measurements at the same time
  stopifnot(
    Biobase::pData(young.noNorm.noResponse.filtered) %>%
      dplyr::group_by(participant_id, time_post_last_vax) %>%
      dplyr::filter(n() > 1) %>%
      nrow() == 0
  )
  
  # Ensure everyone has a baseline (time ≤ 0) sample
  stopifnot(
    Biobase::pData(young.noNorm.noResponse.filtered) %>%
      dplyr::group_by(participant_id) %>%
      dplyr::filter(!any(time_post_last_vax <= 0)) %>%
      nrow() == 0
  )
  
  # Restrict BTMs to those present in the filtered dataset
  btm_list.filtered <- keep_symbols(
    genesets.symbol = btm_list,
    symbols2keep     = rownames(young.noNorm.noResponse.filtered)
  )
  
  # Build the analysis design table
  analysis_df <- Biobase::pData(young.noNorm.noResponse.filtered) %>%
    dplyr::filter(time_post_last_vax > 0) %>%
    dplyr::select(
      pathogen,
      vaccine_type,
      time_post_last_vax,
      study_accession
    ) %>%
    dplyr::arrange(
      pathogen,
      vaccine_type,
      time_post_last_vax,
      study_accession
    ) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      n_participants  = NA,
      qusage_contrast = NA,
      qusage_output   = list(NULL),
      error           = list(NULL)
    ) %>%
    dplyr::as_tibble()
  
  # Wrap QuSAGE in a safely() call to catch errors without stopping the loop
  safe_qusage <- purrr::safely(
    .f    = qusage::qusage,
    quiet = TRUE,
    otherwise = NULL
  )
  
  # Loop through each row of the design table
  for (row_nr in seq_len(nrow(analysis_df))) {
    row <- analysis_df[row_nr, , drop = FALSE]
    
    # Print progress
    cat(
      sprintf(
        "Row %d/%d: %s %s Day%s %s\n",
        row_nr,
        nrow(analysis_df),
        row$pathogen,
        row$vaccine_type,
        row$time_post_last_vax,
        row$study_accession
      )
    )
    
    # Subset the ExpressionSet to baseline + current post-vax day
    cur_eset <- young.noNorm.noResponse.filtered[,
                                                 young.noNorm.noResponse.filtered$pathogen        == row$pathogen &
                                                   young.noNorm.noResponse.filtered$vaccine_type    == row$vaccine_type &
                                                   young.noNorm.noResponse.filtered$study_accession == row$study_accession &
                                                   (
                                                     young.noNorm.noResponse.filtered$time_post_last_vax <= 0 |
                                                       young.noNorm.noResponse.filtered$time_post_last_vax == row$time_post_last_vax
                                                   )
    ]
    
    # Define labels/contrast strings
    prevax_string    <- "baseline"
    postvax_string   <- paste0("Day", row$time_post_last_vax)
    qusage_contrast  <- paste0(postvax_string, "-", prevax_string)
    qusage_labels    <- ifelse(
      test = cur_eset$time_post_last_vax == row$time_post_last_vax,
      yes  = postvax_string,
      no   = prevax_string
    )
    
    # Run (safely) and collect results
    result <- safe_qusage(
      eset       = cur_eset,
      labels     = qusage_labels,
      contrast   = qusage_contrast,
      geneSets   = btm_list.filtered,
      pairVector = cur_eset$participant_id
    )
    
    # Record number of participants and contrast
    analysis_df$n_participants[row_nr]    <- length(unique(cur_eset$participant_id))
    analysis_df$qusage_contrast[row_nr]   <- qusage_contrast
    
    # Save output or error
    if (!is.null(result$result)) analysis_df$qusage_output[[row_nr]] <- result$result
    if (!is.null(result$error))  analysis_df$error[[row_nr]]        <- result$error
  }
  
  # Persist results to disk and load into the workspace
  saveRDS(analysis_df, longitudinal_qusage_young_df_path)
  longitudinal_qusage_young_df <- analysis_df
}

dgsa_results_directory = "./output/results/qusage"
# Path for saving/loading the meta‐analysis results
longitudinal_qusage_young_meta_df_path <- file.path(
  dgsa_results_directory,
  "longitudinal_qusage_young_meta_df.rds"
)

if (file.exists(longitudinal_qusage_young_meta_df_path)) {
  # If results already exist, load them
  longitudinal_qusage_young_meta_df <- readRDS(longitudinal_qusage_young_meta_df_path)
  
} else {
  # Otherwise run the meta‐analysis across pathogen, vaccine, and time groups
  longitudinal_qusage_young_meta_df <- meta_qusage(
    longitudinal_qusage_young_df,
    group_cols = c(
      "pathogen",
      "vaccine_type",
      "time_post_last_vax"
    )
  )
  
  # Save the newly computed results for future use
  saveRDS(
    longitudinal_qusage_young_meta_df,
    longitudinal_qusage_young_meta_df_path
  )
}

# Tidy up the meta‐QuSAGE output and annotate with FDR
tidy_longitudinal_qusage_young_meta_df <- longitudinal_qusage_young_meta_df %>%
  tidy_qusage_output(col = "meta_qusage_output") %>%
  dplyr::mutate(
    FDR              = p.adjust(p.value, method = "BH"),
    FDR_signif_label = create_signif_label(FDR)
  )

corrMethod='pearson'
tp_int=c(1,3,7)
ab_response_col='MFC'

#Load processed qusage results
qusage_df=tidy_longitudinal_qusage_young_meta_df

#Load BTM groups
BTM_groups <- read.delim(
  "./data-raw/BTM_functional_groups.txt",
  comment.char = "#"
)
#Optional:Remove BTMs without subgroup
BTM_groups=BTM_groups[sapply(BTM_groups$SUBGROUP, function(x) nchar(x)>0),]
#Set BTM group as factor and manually set order
BTM_groups$SUBGROUP=factor(BTM_groups$SUBGROUP,
                           levels=c('ANTIGEN PRESENTATION', 'INFLAMMATORY/TLR/CHEMOKINES','INTERFERON/ANTIVIRAL SENSING',
                                    'MONOCYTES', 'DC ACTIVATION', 'NEUTROPHILS', 'SIGNAL TRANSDUCTION', 'ECM AND MIGRATION',
                                    'ENERGY METABOLISM', 'CELL CYCLE', 'NK CELLS', 'T CELLS', 'B CELLS', 'PLASMA CELLS'))
#Sort BTM list by subgroup
BTM_groups=BTM_groups[order(BTM_groups$SUBGROUP),]
rownames(BTM_groups)=BTM_groups$NAME

#Match qusage BTM names with list
common_BTMs=intersect(sub('.*\\(S','S',sub('.*\\(M','M',BTM_groups$NAME)), sub('.*\\(S','S',sub('.*\\(M','M',unique(qusage_df$geneset))))
#Remove qusage BTMs without group
qusage_df=qusage_df[sub('.*\\(S','S',sub('.*\\(M','M',qusage_df$geneset)) %in% common_BTMs,]
#Rename to match BTM_groups
qusage_df$geneset=BTM_groups$NAME[match(sub('.*\\(S','S',sub('.*\\(M','M',qusage_df$geneset)),sub('.*\\(S','S',sub('.*\\(M','M',BTM_groups$NAME)))]
#Remove BTMs missing from qusage results
BTM_groups=BTM_groups[sub('.*\\(S','S',sub('.*\\(M','M',BTM_groups$NAME)) %in% common_BTMs,]

#Add pathogen_type (vaccine) column
qusage_df$pt=paste(qusage_df$pathogen," (",qusage_df$vaccine_type,")", sep='')

#Subset data to timepoints of interest
qusage_df=qusage_df[qusage_df$time_post_last_vax %in% tp_int,]
pt_uni=unique(qusage_df$pt)
pt_uni_tp=lapply(tp_int, function(x) unique(qusage_df$pt[qusage_df$time_post_last_vax==x])) #Get vaccines with data at each timepoint
#Split by timepoint
qusage_df=lapply(tp_int, function(x) qusage_df[qusage_df$time_post_last_vax==x,])
#Convert to wide format
pt_FC_NES=lapply(qusage_df, function(x) data.frame(spread(x[,c('pt','geneset','activity_score')], pt, activity_score)))
pt_FC_NES=lapply(pt_FC_NES, function(x) {rownames(x)=x$geneset; x[,-1]})
pt_FC_NES=lapply(1:length(pt_FC_NES), function(x) {colnames(pt_FC_NES[[x]])=pt_uni_tp[[x]]; pt_FC_NES[[x]]})
pt_FC_p=lapply(qusage_df, function(x) data.frame(spread(x[,c('pt','geneset','FDR')], pt, FDR)))
pt_FC_p=lapply(pt_FC_p, function(x) {rownames(x)=x$geneset; x[,-1]})
pt_FC_p=lapply(1:length(pt_FC_p), function(x) {colnames(pt_FC_p[[x]])=pt_uni_tp[[x]]; pt_FC_p[[x]]})
#Reorder qusage results to match BTM_groups
pt_FC_NES=lapply(pt_FC_NES, function(x) x[match(BTM_groups$NAME,rownames(x)),])
pt_FC_p=lapply(pt_FC_p, function(x) x[match(BTM_groups$NAME,rownames(x)),])

# pt_FC_NES and pt_FC_p contains a list of 3 data frames (times 1,3,7)
# Each row refers to a geneset and each column a vaccine. the entries are either
# P values or FC

#Generate Ab response correlation data
#Load GE
eset_wAb=young.noNorm.withResponse
colnames(eset_wAb)=eset_wAb$uid

#Create combined vaccine type/pathogen column
eset_wAb$pt=paste(eset_wAb$pathogen," (",eset_wAb$vaccine_type,")", sep='')
#Label columns by uid
colnames(eset_wAb)=eset_wAb$uid

#Find samples from timepoints of interest
tp_int=c(0,tp_int)
ind=lapply(tp_int, function(x) which(eset_wAb$study_time_collected==x))
#Combine indices of all timepoints of interest
ind_all=Reduce(union,ind)
#Retain only samples from timepoints of interest
eset_wAb=eset_wAb[,ind_all]
#Recompute timepoint indices after removing extraneous timepoints
ind=lapply(tp_int, function(x) which(eset_wAb$study_time_collected==x))
#Remove samples from a single study with fewer than sample_cutoff samples at any timepoint
sample_cutoff=3
matrix_uni_tp=lapply(ind,function(x) unique(eset_wAb[,x]$matrix))
matrix_ind=lapply(1:length(ind),function(x)
  lapply(1:length(matrix_uni_tp[[x]]), function(y) which(matrix_uni_tp[[x]][[y]]==eset_wAb[,ind[[x]]]$matrix)))
ind_cut_all=vector()
for (i in 1:length(matrix_ind)) {
  ind_cut=which(sapply(matrix_ind[[i]],length)<sample_cutoff)
  if (is_empty(ind_cut)==FALSE) {
    for (j in 1:length(ind_cut)){
      ind_cut_all=c(ind_cut_all,ind[[i]][matrix_ind[[i]][[ind_cut[j]]]])
    }
  }
}
if (is_empty(ind_cut_all)==FALSE) {
  eset_wAb=eset_wAb[,-ind_cut_all]
}
#Recompute timepoint indices after removing samples
tp_int=unique(eset_wAb$study_time_collected[which(eset_wAb$study_time_collected>=0)])
ind=lapply(tp_int, function(x) which(eset_wAb$study_time_collected==x))

#Remove genes with NA
eset_wAb = eset_wAb[complete.cases(exprs(eset_wAb)),]

# The above expression set contains expression and MFC

#Create unique list of studies
matrix_uni=unique(eset_wAb$matrix)

#Collapse to BTM level
#Load BTM list
GSA <- GSA.read.gmt(file.path("./data-raw/", 'BTM.gmt'))
BTM_list <- GSA[['genesets']]
names(BTM_list) <- GSA[['geneset.names']]
#match with BTM groups/qusage data
BTM_list=BTM_list[na.omit(match(sub('.*\\(S','S',sub('.*\\(M','M',BTM_groups$NAME)),sub('.*\\(S','S',sub('.*\\(M','M',names(BTM_list)))))]

#Collapse - find arithmetic mean of genes in each BTM for each sample within each study
exp_BTM=do.call(rbind, lapply(BTM_list, function(x) colMeans(exprs(eset_wAb[na.omit(match(x,rownames(eset_wAb))),]),na.rm=TRUE)))
#Create BTM eset_wAb
# Binds this to the MFC
eset_wAb_BTM=ExpressionSet(exp_BTM, eset_wAb@phenoData)
#Replace eset_wAb with BTM eset_wAb
eset_wAb=eset_wAb_BTM
rm(eset_wAb_BTM)

#eset_wAb is now study+time+geneset level FC + MFC values

#Compute D0 normalized FC
ind_D0=which(0==eset_wAb$study_time_collected)
common=lapply(2:length(ind),function(x) intersect(eset_wAb$participant_id[ind[[x]]],eset_wAb$participant_id[ind_D0]))
ia=lapply(2:length(ind),function(x) na.omit(match(common[[x-1]],eset_wAb$participant_id[ind[[x]]])))
ib=lapply(2:length(ind),function(x) na.omit(match(common[[x-1]],eset_wAb$participant_id[ind_D0])))
exp_FC_wAb=lapply(2:length(ind),function(x) eset_wAb[,ind[[x]][ia[[x-1]]]])
exp_FC_wAb=lapply(2:length(ind),function(x) {exprs(exp_FC_wAb[[x-1]])=exprs(exp_FC_wAb[[x-1]])-exprs(eset_wAb[,ind_D0[ib[[x-1]]]]); exp_FC_wAb[[x-1]]})

# exp_FC_wAb is now a list of expression sets (1 for each timepoint)
# Inside, the data is per-study FC values normalised by baseline

#For each study, average expression across all subjects
matrix_uni_tp=lapply(exp_FC_wAb,function(x) x$matrix[!duplicated(x$matrix)])
#Store study metadata
matrix_uni_tp_metaData=lapply(exp_FC_wAb,function(x) pData(x)[!duplicated(x$matrix),] %>%
                                rownames_to_column() %>%
                                dplyr::select(matrix, vaccine, vaccine_type, pathogen, adjuvant, pt) %>%
                                column_to_rownames(var = "matrix"))
#Indices for each study
matrix_ind=lapply(1:length(exp_FC_wAb),
                  function(x) lapply(1:length(matrix_uni_tp[[x]]), function(y) which(matrix_uni_tp[[x]][[y]]==exp_FC_wAb[[x]]$matrix)))
#Study size for weighting
exp_FC_wAb_n=lapply(1:length(exp_FC_wAb),
                    function(x) sapply(1:length(matrix_uni_tp[[x]]), function(y) length(matrix_ind[[x]][[y]])))

#Correlation with Ab response
# This is within-study/timepoint correlation of MFC with d0 normalised FC.
exp_FC_wAb_corr_stat =
  lapply(1:length(exp_FC_wAb),  #for each timepoint
         function(x) sapply(1:length(matrix_uni_tp[[x]]),  #for each study
                            function(y) apply(exprs(exp_FC_wAb[[x]][,matrix_ind[[x]][[y]]]), 1,  #for each geneset
                                              function(z) cor.test(z, pData(exp_FC_wAb[[x]][,matrix_ind[[x]][[y]]])[[ab_response_col]], method = corrMethod)$estimate)))

#Merge mean FC, correlation/t-test test statistics by vaccine (weighted average by study size)

#Find unique pathogen/vaccine types by timepoint
pt_uni_tp=lapply(matrix_uni_tp_metaData,function(x) x$pt[!duplicated(x$pt)])
pt_metaData=lapply(1:length(matrix_uni_tp_metaData),
                   function(x) dplyr::select(
                     matrix_uni_tp_metaData[[x]][!duplicated(matrix_uni_tp_metaData[[x]]$pt),],
                     vaccine_type, pathogen, adjuvant))

#Sample indices per vaccine
pt_ind=lapply(1:length(exp_FC_wAb),
              function(x) sapply(1:length(pt_uni_tp[[x]]),
                                 function(y) which(pt_uni_tp[[x]][[y]]==matrix_uni_tp_metaData[[x]]$pt)))
pt_corr_stat=
  lapply(1:length(exp_FC_wAb),
         function(x) sapply(1:length(pt_uni_tp[[x]]),
                            function(y) if(length(pt_ind[[x]][[y]])>1) rowWeightedMeans(exp_FC_wAb_corr_stat[[x]][,pt_ind[[x]][[y]]], w=exp_FC_wAb_n[[x]][pt_ind[[x]][[y]]], na.rm=TRUE) else exp_FC_wAb_corr_stat[[x]][,pt_ind[[x]][[y]]]))

pt_corr_stat=lapply(1:length(pt_corr_stat), function(x) {rownames(pt_corr_stat[[x]])=BTM_groups$NAME; colnames(pt_corr_stat[[x]])=pt_uni_tp[[x]]; pt_corr_stat[[x]]})

#Add vaccines with missing data at each timepoint
pt_FC_NES_full=vector("list",length(pt_FC_NES))
pt_FC_p_full=vector("list",length(pt_FC_NES))
pt_corr_stat_full=vector("list",length(pt_corr_stat))


for (i in 1:length(pt_corr_stat)) {
  pt_FC_NES_full[[i]]=matrix(NA, ncol=length(pt_uni), nrow=nrow(pt_FC_NES[[i]]), dimnames = list(rownames(pt_FC_NES[[i]]), pt_uni))
  
  pt_FC_NES_full[[i]][rownames(pt_FC_NES[[i]]), colnames(pt_FC_NES[[i]])]=unlist(pt_FC_NES[[i]])
  
  pt_FC_p_full[[i]]=matrix(NA, ncol=length(pt_uni), nrow=nrow(pt_FC_p[[i]]), dimnames = list(rownames(pt_FC_p[[i]]), pt_uni))
  
  pt_FC_p_full[[i]][rownames(pt_FC_p[[i]]), colnames(pt_FC_p[[i]])]=unlist(pt_FC_p[[i]])
  
  pt_corr_stat_full[[i]]=matrix(NA, ncol=length(pt_uni), nrow=nrow(pt_corr_stat[[i]]), dimnames = list(rownames(pt_corr_stat[[i]]), pt_uni))
  
  pt_corr_stat_full[[i]][rownames(pt_corr_stat[[i]]), colnames(pt_corr_stat[[i]])]= unlist(pt_corr_stat[[i]])
}
pt_FC_NES=pt_FC_NES_full
pt_FC_p=pt_FC_p_full
pt_corr_stat=pt_corr_stat_full
rm(pt_corr_stat_full, pt_FC_NES_full, pt_FC_p_full)
#Update vaccine per timepoint variable
pt_uni_tp=lapply(1:length(pt_FC_NES), function(x) pt_uni)

tp_int=c(1,3,7)
names(pt_corr_stat) = tp_int
names(pt_FC_NES) = tp_int

# Load aggregates
BTM_functional_groups <- read.delim("./data-raw/BTM_functional_groups.txt", comment.char="#")

# Aggregates preprocessing 
# Assign geneset aggregates
# Convert to Lowercase
BTM_functional_groups$SUBGROUP <- str_to_title(tolower(BTM_functional_groups$SUBGROUP))
# Replace specific terms to stay in uppercase
BTM_functional_groups$SUBGROUP <- gsub("\\bEcm\\b", "ECM", BTM_functional_groups$SUBGROUP)
BTM_functional_groups$SUBGROUP <- gsub("\\bDc\\b", "DC", BTM_functional_groups$SUBGROUP)
BTM_functional_groups$SUBGROUP <- gsub("\\bTlr\\b", "TLR", BTM_functional_groups$SUBGROUP)
BTM_functional_groups$SUBGROUP <- gsub("\\bNk\\b", "NK", BTM_functional_groups$SUBGROUP)


# Replace empty mappings with NA
BTM_functional_groups$SUBGROUP <- ifelse(is.na(BTM_functional_groups$SUBGROUP) | BTM_functional_groups$SUBGROUP == "", "NA", BTM_functional_groups$SUBGROUP)


results_df_qusage = tidy_longitudinal_qusage_young_meta_df %>%
  # keep original geneset text
  mutate(gs.name.description = geneset) %>%
  # extract the final "(M...)" or "(S...)" into gs.name
  mutate(
    gs.name = str_extract(geneset, "(?<=\\()([MS][^)]*)(?=\\)$)"),
    # drop that final "(M...)" or "(S...)" from the description
    gs.description = str_remove(geneset, " \\([MS][^)]*\\)$")
  ) %>%
  # build the rest of your columns
  mutate(
    condition        = paste0(pathogen, " (", vaccine_type, ")"),
    time             = round(time_post_last_vax,2),
    comparison       = paste0(condition, " - Day ", time),
    activation.score = activity_score,
    rawPval          = p.value
  ) %>%
  # 4) join in the SUBGROUP lookup and rename it
  left_join(
    BTM_functional_groups %>% select(BTM, SUBGROUP),
    by = c("gs.name" = "BTM")
  ) %>%
  rename(gs.aggregate = SUBGROUP) %>%
  select(
    condition,
    time,
    comparison,
    gs.name,
    gs.description,
    gs.name.description,
    gs.aggregate,
    activation.score,
    rawPval
  ) 

# Now add metadata (colours, aggregates etc.)
conditions_order <- c("Tuberculosis (RVV)",
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
                      "Smallpox (LV)")

results_df_qusage$condition = factor(results_df_qusage$condition, levels = conditions_order)

times <- as.numeric(results_df_qusage$time) %>% unique() %>% sort()

results_df_qusage$time = results_df_qusage$time %>% 
  factor(levels = times)

condition_colors = c(
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

# name the vector so we can index by the condition string
names(condition_colors) <- conditions_order

# now create the new column by lookup
results_df_qusage$condition.colour <- condition_colors[ as.character(results_df_qusage$condition) ]

category_order <- c(
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
  "Plasma Cells",
  "NA"
)

results_df_qusage$gs.aggregate = factor(results_df_qusage$gs.aggregate, levels = category_order)

aggregate_colors <- c(
  "#7c5fcd",
  "#57c39d",
  "#c165b9",
  "#55c463",
  "#7082ca",
  "#64a332",
  "#45aecf",
  "#df9545",
  "#b7b238",
  "#a6b36c",
  "#667328",
  "#662d2e",
  "#ff8fa3",
  "#c05299",
  "#8f2d56",
  "#adb5bd"
)

# name the vector so we can look up by the aggregate string
names(aggregate_colors) <- category_order

# create the new colour column
results_df_qusage$gs.colour <- 
  aggregate_colors[ as.character(results_df_qusage$gs.aggregate) ]

# Create multiple adjusted p-value columns 
p_methods <- c("BH", "bonferroni", "holm", "BY", "hommel", "hochberg")
results_df_qusage <- results_df_qusage %>%
  
  # 1) global corrections
  mutate(
    !!! setNames(
      lapply(p_methods, function(m) {
        # for each method, p.adjust all rawPval
        rlang::expr(p.adjust(rawPval, method = !!m))
      }),
      # name them global.adjPval_<method>
      paste0("global.adjPval_", p_methods)
    )
  ) %>%
  
  # 2) within‐time corrections
  group_by(time) %>%
  mutate(
    !!! setNames(
      lapply(p_methods, function(m) {
        rlang::expr(p.adjust(rawPval, method = !!m))
      }),
      paste0("withinTime.adjPval_", p_methods)
    )
  ) %>%
  ungroup() %>%
  
  # 3) within‐comparison corrections
  group_by(comparison) %>%
  mutate(
    !!! setNames(
      lapply(p_methods, function(m) {
        rlang::expr(p.adjust(rawPval, method = !!m))
      }),
      paste0("withinComparison.adjPval_", p_methods)
    )
  ) %>%
  ungroup()

# Finally add in fold change and correlation statistics
# Build a list of long data.tables, one per time-point
dt_list <- lapply(names(pt_corr_stat), function(tm) {
  mat <- pt_corr_stat[[tm]]
  # coerce to matrix if it’s a tibble or data.frame (so rownames exist)
  mat <- as.matrix(mat)
  # turn into a data.table, keeping the rownames column
  dt <- as.data.table(mat, keep.rownames = "gs.name.description")
  # add the time column and then melt to long form
  melt(
    dt[, time := tm],
    id.vars = c("time", "gs.name.description"),
    variable.name = "condition",
    value.name    = "corr.mean"
  )
})

# Stack them all into one big long table
pt_corr_long <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# 2) join back; unmatched get NA automatically
results_df_qusage <- results_df_qusage %>%
  left_join(
    pt_corr_long,
    by = c("time", "condition", "gs.name.description")
  )

# results_df_qusage <- results_df_qusage %>%
#   group_by(time) %>%
#   mutate(
#     corr.mean = case_when(
#       # if all the same (max==min), set to 0
#       max(corr.mean, na.rm = TRUE) == min(corr.mean, na.rm = TRUE) ~ 0,
#       TRUE ~ 2 * (corr.mean - min(corr.mean, na.rm = TRUE)) /
#              (max(corr.mean, na.rm = TRUE) - min(corr.mean, na.rm = TRUE)) - 1
#     )
#   ) %>%
#   ungroup()

results_df_qusage$corr.mean = results_df_qusage$corr.mean %>% as.numeric()

# Now FC
# Build a list of long data.tables, one per time-point
dt_list <- lapply(names(pt_FC_NES), function(tm) {
  mat <- pt_FC_NES[[tm]]
  # coerce to matrix if it’s a tibble or data.frame (so rownames exist)
  mat <- as.matrix(mat)
  # turn into a data.table, keeping the rownames column
  dt <- as.data.table(mat, keep.rownames = "gs.name.description")
  # add the time column and then melt to long form
  melt(
    dt[, time := tm],
    id.vars = c("time", "gs.name.description"),
    variable.name = "condition",
    value.name    = "fc.score"
  )
})

# Stack them all into one big long table
pt_FC_long <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# 2) join back; unmatched get NA automatically
results_df_qusage <- results_df_qusage %>%
  left_join(
    pt_FC_long,
    by = c("time", "condition", "gs.name.description")
  )

results_df_qusage$fc.score = as.numeric(results_df_qusage$fc.score)


results_df_qusage = results_df_qusage %>% 
  mutate(fc.score = activation.score,
         mean.corr = corr.mean)
# Reorder columns
results_df_qusage = results_df_qusage %>% 
  select(comparison, condition, time, condition.colour,
         gs.name, gs.description, gs.name.description, gs.aggregate, gs.colour,
         activation.score, fc.score, corr.mean, mean.corr,
         rawPval, 
         global.adjPval_holm, withinTime.adjPval_holm, withinComparison.adjPval_holm,
         global.adjPval_hommel, withinTime.adjPval_hommel, withinComparison.adjPval_hommel,
         global.adjPval_bonferroni, withinTime.adjPval_bonferroni, withinComparison.adjPval_bonferroni,
         global.adjPval_BH, withinTime.adjPval_BH, withinComparison.adjPval_BH,
         global.adjPval_hochberg, withinTime.adjPval_hochberg, withinComparison.adjPval_hochberg,
         global.adjPval_BY, withinTime.adjPval_BY, withinComparison.adjPval_BY)

results_df_qusage$time = results_df_qusage$time %>% 
  factor(levels = times)

results_df_qusage$method = "qusage"

dgsa_results_directory = "./output/results/qusage/"

saveRDS(results_df_qusage, file = paste0(dgsa_results_directory, "qusage_dgsa_results_processed.rds"))
