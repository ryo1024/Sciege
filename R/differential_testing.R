#' Differential isoform expression testing by cell type
#'
#' @param isoform_matrix A matrix of isoform expression (samples × isoforms).
#' @param celltype_props A matrix of cell type proportions (samples × cell types).
#' @param metadata A data.frame with sample metadata, rownames matching samples.
#' @param group_col Column in metadata defining two groups for comparison.
#' @param n_boot Number of bootstrap iterations for standard error estimation.
#' @param out_dir Optional directory to save CSV results.
#'
#' @return A data.frame with cell type, isoform, logFC, se, Wald statistic, and p-value.
#' @export
differential_isoform_expression <- function(isoform_matrix,
                                            celltype_props,
                                            metadata,
                                            group_col,
                                            n_boot = 1000,
                                            out_dir = NULL) {
  if (!requireNamespace("nnls", quietly = TRUE)) stop("Package 'nnls' is required.")
  if (!requireNamespace("boot", quietly = TRUE)) stop("Package 'boot' is required.")

  samples <- intersect(intersect(rownames(isoform_matrix), rownames(celltype_props)), rownames(metadata))
  isoform_matrix <- isoform_matrix[samples, , drop = FALSE]
  celltype_props <- celltype_props[samples, , drop = FALSE]
  metadata <- metadata[samples, , drop = FALSE]
  
  groups <- unique(metadata[[group_col]])
  if (length(groups) != 2) stop("Only two-group comparison is supported.")

  group_a <- groups[1]
  group_b <- groups[2]

  get_group_result <- function(group_samples) {
    estimate_isoform_expression_by_celltype(
      isoform_matrix = isoform_matrix[group_samples, , drop = FALSE],
      celltype_props = celltype_props[group_samples, , drop = FALSE],
      n_boot = n_boot
    )
  }

  # Estimate for each group
  est_a <- get_group_result(rownames(metadata)[metadata[[group_col]] == group_a])
  est_b <- get_group_result(rownames(metadata)[metadata[[group_col]] == group_b])

  # Compute logFC, Wald statistic, and p-values
  mean_a <- est_a$mean_expression
  mean_b <- est_b$mean_expression
  se_a <- est_a$se_expression
  se_b <- est_b$se_expression

  all_celltypes <- rownames(mean_a)
  all_isoforms <- colnames(mean_a)

  results <- do.call(rbind, lapply(all_celltypes, function(celltype) {
    sapply(all_isoforms, function(iso) {
      mu1 <- mean_a[celltype, iso] + 1e-6  # add small value to avoid log(0)
      mu2 <- mean_b[celltype, iso] + 1e-6
      logFC <- log2(mu2 / mu1)
      diff <- mu2 - mu1
      se_total <- sqrt(se_a[celltype, iso]^2 + se_b[celltype, iso]^2)
      wald_stat <- diff / se_total
      pval <- 2 * pnorm(-abs(wald_stat))
      c(logFC = logFC, diff = diff, se = se_total, stat = wald_stat, pval = pval)
    }, simplify = "matrix")
  }))

#   print(results)
#   print(all_celltypes)
#   print(all_isoforms)

  result_df <- as.data.frame(results)
  result_df$cell_type <- rep(all_celltypes, each = length(c('logFC', 'diff', 'se', 'stat', 'pval') ) )
  result_df$stat_type <- rep(c('logFC', 'diff', 'se', 'stat', 'pval'), times = length(all_celltypes))
  #result_df$isoform <- rep(all_isoforms, times = length(all_celltypes))
  #result_df <- result_df[, c("cell_type", "isoform", "logFC", "diff", "se", "stat", "pval")]
  #rownames(result_df) <- NULL

  #print(result_df)

  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    write.csv(result_df, file.path(out_dir, "differential_isoform_expression.csv"), row.names = FALSE)
  }

  return(result_df)
}