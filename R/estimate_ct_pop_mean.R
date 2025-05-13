#' Estimate mean isoform expression per cell type using NNLS and bootstrap
#'
#' @param isoform_matrix A matrix of isoform expression (samples × isoforms).
#' @param celltype_props A matrix of cell type proportions (samples × cell types).
#' @param n_boot Number of bootstrap replicates to estimate standard errors.
#' @param out_dir Optional directory to write output CSVs.
#'
#' @return A list with:
#'   - mean_expression: matrix of estimated expression (cell types × isoforms)
#'   - se_expression: matrix of bootstrapped standard errors (cell types × isoforms)
#' @export
estimate_isoform_expression_by_celltype <- function(isoform_matrix,
                                                    celltype_props,
                                                    n_boot = 1000,
                                                    out_dir = NULL) {
  if (!requireNamespace("nnls", quietly = TRUE)) stop("Package 'nnls' is required.")
  if (!requireNamespace("boot", quietly = TRUE)) stop("Package 'boot' is required.")

  # Ensure samples match
  common_samples <- intersect(rownames(isoform_matrix), rownames(celltype_props))
  isoform_matrix <- isoform_matrix[common_samples, , drop = FALSE]
  celltype_props <- celltype_props[common_samples, , drop = FALSE]

  cell_types <- colnames(celltype_props)
  isoforms <- colnames(isoform_matrix)

  # Step 1: Estimate mean per cell type using NNLS for each isoform
  mean_expr <- sapply(seq_len(ncol(isoform_matrix)), function(j) {
    y <- isoform_matrix[, j]
    fit <- nnls::nnls(as.matrix(celltype_props), y)
    return(fit$x)
  })
  rownames(mean_expr) <- cell_types
  colnames(mean_expr) <- isoforms

  # Step 2: Bootstrap SE
  boot_fn <- function(data, indices, isoform_idx) {
    y <- data$isoform[indices]
    X <- data$props[indices, , drop = FALSE]
    fit <- nnls::nnls(X, y)
    return(fit$x)
  }

  se_matrix <- matrix(NA, nrow = ncol(celltype_props), ncol = ncol(isoform_matrix))
  rownames(se_matrix) <- cell_types
  colnames(se_matrix) <- isoforms

  for (j in seq_len(ncol(isoform_matrix))) {
    dat <- list(
      isoform = isoform_matrix[, j],
      props = as.matrix(celltype_props)
    )
    boot_out <- boot::boot(data = dat, statistic = function(d, i) boot_fn(d, i, j), R = n_boot)
    se_matrix[, j] <- apply(boot_out$t, 2, sd)
  }

  # Optional output
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(mean_expr, file.path(out_dir, "estimated_isoform_expression.csv"))
    write.csv(se_matrix, file.path(out_dir, "estimated_isoform_se.csv"))
  }

  return(list(mean_expression = mean_expr, se_expression = se_matrix))
}