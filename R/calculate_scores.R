#'Calculate individual, gene-level, and set-level activation scores for binary conditions
#'
#'These scores can be used in downstream analysis.
#'
#'By default, gene-set level fold changes are calculated as follows: first, average expressions
#'are computed for each level of phi. A gene-level fold-change is calculated as the difference in
#'these averages. A gene-set level fold change is the mean of the gene-level fold-changes in the set.
#'
#'If data is paired, fold-changes may be calculated by first computing an individual gene-level fold-change by
#'subtracting the inactive-level expression from the active-level expression for each individual,
#'then averaging these over individuals and genes within each set.
#'
#'@param y a numeric matrix of size \code{G x n} containing the raw RNA-seq
#'counts or preprocessed expression from \code{n} samples for \code{G} genes.
#'
#'@param x a numeric matrix of size \code{n x p} containing the model
#'covariate(s) from \code{n} samples (design matrix).
#'
#'@param phi a numeric design matrix of size \code{n x K} containing the K
#'variable(s) of interest( e.g. bases of time).
#'
#'@param GSA a GSA object containing information on the genesets. This should contain
#'at minimum \itemize{
#'\item \code{genesets} a list where each element is a vector
#'containing names of rows from \code{y} for a given geneset.}
#'
#'@param use_phi a logical flag indicating whether conditional means should be
#'conditioned on \code{phi} and on covariate(s) \code{x}, or on \code{x} alone.
#'Default is \code{TRUE} in which case conditional means are estimated
#'conditionally on both \code{x} and \code{phi}.
#'
#'@param preprocessed a logical flag indicating whether the expression data have
#'already been preprocessed (e.g. log2 transformed). Default is \code{FALSE}, in
#'which case \code{y} is assumed to contain raw counts and is normalized into
#'log(counts) per million.
#'
#'@param gene_based a logical flag indicating whether to estimate weights at the
#'gene-level. Default is \code{FALSE}, when weights will be estimated at the
#'observation-level.
#'
#'@param bw a character string indicating the smoothing bandwidth selection
#'method to use. See \code{\link[stats]{bandwidth}} for details. Possible values
#'are \code{'ucv'}, \code{'SJ'}, \code{'bcv'}, \code{'nrd'} or \code{'nrd0'}.
#'Default is \code{'nrd'}.
#'
#'@param kernel a character string indicating which kernel should be used.
#'Possibilities are \code{'gaussian'}, \code{'epanechnikov'},
#'\code{'rectangular'}, \code{'triangular'}, \code{'biweight'},
#'\code{'tricube'}, \code{'cosine'}, \code{'optcosine'}. Default is
#'\code{'gaussian'} (NB: \code{'tricube'} kernel corresponds to the loess
#'method).
#'
#'@param transform a logical flag indicating whether values should be
#'transformed to uniform for the purpose of local linear smoothing. This may be
#'helpful if tail observations are sparse and the specified bandwidth gives
#'suboptimal performance there. Default is \code{TRUE}.
#'
#'@param verbose a logical flag indicating whether informative messages are
#'printed during the computation. Default is \code{TRUE}.
#'
#'@param na.rm logical: should missing values (including \code{NA} and
#'\code{NaN}) be omitted from the calculations? Default is \code{FALSE}.
#'
#'@param active_level the level of binary \code{phi} which should be considered the "active" group for the
#'purpose of calculating the scores. This determines the "directionality" of the scores, i.e.
#'reversing the active levels will simply reverse the signs of the all the scores.
#'
#'@param sample_group a vector of length \code{n} indicating whether the samples
#'are grouped in pairs. If data is paired, and all individuals under study have pairs,individual gene-level
#'fold-changes will be calculated by subtracting the non-active-level values from the active-level values.
#'Then, these are averaged across individuals and across genes in each geneset. If this is given,
#'but each individual does not have a pair (i.e. one value corresponding to \code{phi} =\code{active_level}
#'and another corresponding to \code{phi} !=\code{active_level}), fold-changes are calculated in the
#'default manner. Default is \code{NULL} in which case no grouping is performed.
#'
#'@return a list containing the following components:\itemize{
#'\item \code{individual.scores} a list of matrices giving the individual-level scores for each set.
#' These represent covariate-centered, standardised gene expression.
#'\item \code{gene.scores} a list of vectors giving the gene-level scores for each set. These
#' are simply the difference in the mean individual scores under the active group compared to the
#' the control group, for each gene.
#' \item \code{activation.scores} these are the set-level scores which are the mean of the gene scores within each set.
#' }
#'}
#'
#'@author Arthur Hughes

calculate_scores = function(y,
                            x = NULL,
                            phi,
                            GSA,
                            use_phi = TRUE,
                            preprocessed = TRUE,
                            gene_based = FALSE,
                            bw = c("nrd", "ucv", "SJ", "nrd0", "bcv"),
                            kernel = c(
                              "gaussian",
                              "epanechnikov",
                              "rectangular",
                              "triangular",
                              "biweight",
                              "tricube",
                              "cosine",
                              "optcosine"
                            ),
                            transform = TRUE,
                            verbose = TRUE,
                            na.rm = TRUE,
                            active_level,
                            sample_group = NULL) {
  ## dimensions & validity checks
  
  stopifnot(is.matrix(y))
  stopifnot(is.matrix(x))
  stopifnot(is.null(phi) | is.matrix(phi))
  
  g <- nrow(y)  # the number of genes measured
  n <- ncol(y)  # the number of samples measured
  qq <- ncol(x)  # the number of covariates
  stopifnot(nrow(x) == n)
  if (use_phi) {
    stopifnot(nrow(phi) == n)
  }
  
  # If sample_group given, check that data is exactly paired
  if (!is.null(sample_group)) {
    # Tabulate phi and sample_group to check levels for each individual
    tab <- table(sample_group, phi)
    # All elements must be 1 in this table if data is truly paired
    is_paired <- all(tab == 1)
    if (!is_paired) {
      message(
        "You have provided individuals but the data is not paired.
              Ignoring this when calculating the fold-change and proceeding in the default manner."
      )
    }
  }
  
  # removing genes never observed:
  observed <- which(rowSums(y, na.rm = TRUE) != 0)
  nb_g_sum0 <- length(observed) - g
  if (nb_g_sum0 > 0) {
    warning(nb_g_sum0,
            " y rows sum to 0 (i.e. are never observed)",
            "and have been removed")
  }
  
  # kernel <- match.arg(kernel)
  if (preprocessed) {
    y_lcpm <- t(y[observed, ])
  } else {
    # transforming raw counts to log-counts per million (lcpm)
    y_lcpm <- t(apply(y[observed, ], MARGIN = 2, function(v) {
      log2((v + 0.5) / (sum(v) + 1) * 10^6)
    }))
  }
  N <- length(y_lcpm)
  p <- ncol(y_lcpm)
  
  
  # fitting OLS to the lcpm
  if (na.rm) {
    y_lcpm0 <- y_lcpm
    y_lcpm0[is.na(y_lcpm0)] <- 0
    B_ols <- solve(crossprod(x)) %*% t(x) %*% y_lcpm0
  } else {
    B_ols <- solve(crossprod(x)) %*% t(x) %*% y_lcpm
  }
  mu <- x %*% B_ols
  
  
  if (gene_based) {
    mu_avg <- colMeans(mu, na.rm = na.rm)
    mu_x <- mu_avg
  } else {
    mu_x <- mu
    mu_x[is.na(y_lcpm)] <- NA
  }
  
  y_T = t(y)
  yt_mu = y_T - mu_x
  
  # calculate weights
  weights = sp_weights(y,
                       x,
                       phi,
                       use_phi,
                       preprocessed,
                       gene_based,
                       bw,
                       kernel,
                       transform,
                       verbose,
                       na.rm)
  
  weights_all = weights$weights
  
  yt_mu_standardised = yt_mu * t(weights_all)
  
  expr_mat <- as.matrix(yt_mu_standardised)
  
  # 2) Your gene‐sets, with names
  genesets <-
    GSA[["genesets"]]                   # list of character vectors
  names(genesets) <- GSA[["geneset.names"]]
  
  # 3) Precompute, for each set, the column‐indices into expr_mat
  geneset_idx <- lapply(genesets, function(g) {
    which(colnames(expr_mat) %in% g)
  })
  
  phi_inactive <- phi != active_level
  phi_active <- phi == active_level
  
  gene_scores_all <-
    colMeans(expr_mat[phi_active, , drop = FALSE], na.rm = TRUE) -
    colMeans(expr_mat[phi_inactive, , drop = FALSE], na.rm = TRUE)
  
  
  # Calculate gene-level FC scores
  if (is_paired) {
    #If paired, first calculate individual-gene-level scores
    active_idx   <- which(phi == active_level)
    inactive_idx <- which(phi != active_level)
    
    ids_active   <- sample_group[active_idx]
    ids_inactive <- sample_group[inactive_idx]
    
    # 2. order each so that rows line up by individual
    ord_active   <- order(ids_active)
    ord_inactive <- order(ids_inactive)
    
    Y_act <- y_T[active_idx[ord_active], , drop = FALSE]
    Y_in <- y_T[inactive_idx[ord_inactive], , drop = FALSE]
    
    # sanity check: same individual order
    stopifnot(all(ids_active[ord_active] == ids_inactive[ord_inactive]))
    
    
    # 3. per‐individual fold‐change matrix (rows = individuals, cols = predictors)
    ind_fc_mat <- Y_act - Y_in
    
    # 4. predictor‐level score = average across individuals
    fc_scores_all <- colMeans(ind_fc_mat)
    
  } else {
    #If not paired, take difference in gene averages between binary groups
    fc_scores_all = colMeans(y_T[phi_active, , drop = FALSE], na.rm = TRUE) -
      colMeans(y_T[phi_inactive, , drop = FALSE], na.rm = TRUE)
  }
  
  
  # === Extract per‐geneset results ===
  
  # 1) individual.scores: list of sub‐matrices (n_samples × |geneset|)
  individual.scores <- lapply(geneset_idx, function(idxs) {
    expr_mat[, idxs, drop = FALSE]
  })
  
  # 2) gene.scores: list of per‐gene score vectors
  gene.scores <- lapply(geneset_idx, function(idxs) {
    gene_scores_all[idxs]
  })
  
  fc_gene.scores <- lapply(geneset_idx, function(idxs) {
    fc_scores_all[idxs]
  })
  
  # 3) activation.scores: list of single values = mean of gene.scores
  activation.scores <- lapply(gene.scores, function(gs) {
    mean(gs, na.rm = TRUE)
  })
  
  fc.scores <- lapply(fc_gene.scores, function(gs) {
    mean(gs, na.rm = TRUE)
  })
  
  result <- list(
    individual.scores = individual.scores,
    gene.scores       = gene.scores,
    activation.scores = activation.scores,
    fc.scores = fc.scores
  )
}