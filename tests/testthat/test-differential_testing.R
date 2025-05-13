test_that("differential_isoform_expression works with simulated data", {
  skip_if_not_installed("nnls")
  skip_if_not_installed("boot")

  set.seed(123)

  n_samples <- 30
  n_isoforms <- 8
  n_celltypes <- 3

  # Simulate isoform expression (samples × isoforms)
  isoform_matrix <- matrix(
    abs(rnorm(n_samples * n_isoforms, mean = 20, sd = 5)),
    nrow = n_samples, ncol = n_isoforms
  )
  rownames(isoform_matrix) <- paste0("Sample", seq_len(n_samples))
  colnames(isoform_matrix) <- paste0("Isoform", seq_len(n_isoforms))

  # Simulate cell type proportions (rows sum to 1)
  props_raw <- matrix(runif(n_samples * n_celltypes), nrow = n_samples)
  celltype_props <- props_raw / rowSums(props_raw)
  colnames(celltype_props) <- paste0("CellType", seq_len(n_celltypes))
  rownames(celltype_props) <- rownames(isoform_matrix)

  # Simulate metadata for two groups
  metadata <- data.frame(
    group = rep(c("Control", "Treatment"), each = n_samples / 2),
    row.names = rownames(isoform_matrix)
  )

  # Run differential expression
  result <- differential_isoform_expression(
    isoform_matrix = isoform_matrix,
    celltype_props = celltype_props,
    metadata = metadata,
    group_col = "group",
    n_boot = 50  # keep small for test speed
  )

  expect_s3_class(result, "data.frame")
  #expect_true(all(c("cell_type", "isoform", "logFC", "diff", "se", "stat", "pval") %in% colnames(result)))
  expect_equal(nrow(result), n_celltypes * length(c("logFC", "diff", "se", "stat", "pval")))
  expect_false(any(is.na(result$logFC)))
  expect_false(any(is.na(result$pval)))
  expect_true(all(result$pval >= 0 & result$pval <= 1))
})
