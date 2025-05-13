#' Cell type decomposition using BisqueRNA
#'
#' @param expr_matrix Matrix of bulk expression (genes x samples)
#' @param reference_sce SingleCellExperiment or Seurat object as reference
#' @param reference_colname Column in reference metadata containing cell type labels
#' @return A data.frame with estimated cell type proportions (samples x cell types)
#' @export
decompose_cell_types_bisque <- function(expr_matrix,
                                        reference_sce,
                                        out_dir,
                                        reference_colname = "cell_type") {
  if (!requireNamespace("BisqueRNA", quietly = TRUE)) {
    stop("Please install the 'BisqueRNA' package.")
  }

  # Convert expression matrix to ExpressionSet (as required by Bisque)
  bulk_eset <- Biobase::ExpressionSet(assayData = as.matrix(expr_matrix))
  #reference_sce <- Biobase::ExpressionSet(assayData = as.matrix(reference_sce))

  # Ensure reference has cell type column
  if (!reference_colname %in% colnames(Biobase::pData(reference_sce))) {
    stop(paste("Reference object must contain column:", reference_colname))
  }

  # Run decomposition
  result <- BisqueRNA::ReferenceBasedDecomposition(
    bulk.eset = bulk_eset,
    sc.eset = reference_sce,
    use.overlap = FALSE,
    markers = NULL,
    cell.types = reference_colname
  )

  return(as.data.frame(result$bulk.props) )  # Data.frame of proportions
}
