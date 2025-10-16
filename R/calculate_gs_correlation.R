#' Computes different measures of correlation between geneset expression and a response
#'
#' The measures computed are the mean of spearman correlations within a set, and
#' the spearman correlation of mean standardised expression within a set.
#'
#'@param y a numeric matrix of size \code{G x n} containing the raw RNA-seq
#'counts or preprocessed expression from \code{n} samples for \code{G} genes.
#'
#'@param response a numeric vector of size \code{n x 1} containing the response of interest.
#'
#'@param GSA a GSA object containing information on the genesets. This should contain
#'at minimum \itemize{
#'\item \code{genesets} a list where each element is a vector
#'containing names of rows from \code{y} for a given geneset.}
#'
#'@author Arthur Hughes
#'
#'@returns a list containing \itemize{
#'\item \code{mean.corr} : a vector containing the mean of the correlations within each set
#'\item \code{cor.mean} : a vector containing the correlation of the mean expression within each set}
calculate_gs_correlation = function(y, response, GSA) {
  ### VALIDITY CHECKS ###
  
  response = response %>% as.matrix()
  
  # Check that number of samples are the same in y and response
  n = ncol(y)
  if (length(response) != n) {
    stop("Number of response values given does not match number of samples in y")
  }
  
  mean.corr <- sapply(seq_along(GSA[["genesets"]]), function(i) {
    gs <- GSA[["genesets"]][[i]]
    # keep only genes actually in y
    keep <- intersect(gs, rownames(y))
    if (length(keep) == 0) {
      return(NA_real_)  # or 0, or whatever you’d prefer
    }
    # compute Spearman cor for each gene
    cors <- apply(y[keep, , drop = FALSE], 1, function(gene_expr)
      cor(gene_expr, response, method = "spearman", use = "pairwise.complete.obs"))
    # average them
    mean(cors, na.rm = TRUE)
  })
  
  # assign names
  names(mean.corr) <- GSA[["geneset.names"]]
  
  corr.mean <- sapply(seq_along(GSA[["genesets"]]), function(i) {
    # 1. pick your genes
    gs    <- GSA[["genesets"]][[i]]
    keep  <- intersect(gs, rownames(y))
    if (length(keep) == 0)
      return(NA_real_)
    
    # 2. subset and standardize each row (gene) across samples
    y_sub <- y[keep, , drop = FALSE]
    #    t(scale(t(y_sub))) gives genes×samples Z-scores
    zsub  <- t(scale(t(y_sub), center = TRUE, scale = TRUE))
    
    # 3. meta-gene: mean Z per sample
    meta  <- colMeans(zsub, na.rm = TRUE)
    
    # 4. Spearman correlation with response
    cor(meta, response, method = "spearman", use = "pairwise.complete.obs")
  })
  
  # assign names
  names(corr.mean) <- GSA[["geneset.names"]]
  
  return(list(mean.corr = mean.corr, corr.mean = corr.mean))
}