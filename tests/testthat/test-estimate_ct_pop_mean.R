test_that("estimate_isoform_expression_by_celltype works on simulated data", {
  skip_if_not_installed("nnls")
  skip_if_not_installed("boot")

  set.seed(123)

  # Simulate isoform expression matrix (samples × isoforms)
  n_samples <- 20
  n_isoforms <- 10
  n_celltypes <- 3

  isoform_matrix <- matrix(
    abs(rnorm(n_samples * n_isoforms, mean = 10, sd = 3)),
    nrow = n_samples,
    ncol = n_isoforms
  )
  colnames(isoform_matrix) <- paste0("Isoform", 1:n_isoforms)
  rownames(isoform_matrix) <- paste0("Sample", 1:n_samples)

  # Simulate cell type proportions (rows sum to 1)
  raw_props <- matrix(runif(n_samples * n_celltypes), nrow = n_samples)
  celltype_props <- raw_props / rowSums(raw_props)
  colnames(celltype_props) <- paste0("CellType", 1:n_celltypes)
  rownames(celltype_props) <- paste0("Sample", 1:n_samples)

  # Run function with small bootstrap for test speed
  result <- estimate_isoform_expression_by_celltype(
    isoform_matrix = isoform_matrix,
    celltype_props = celltype_props,
    n_boot = 50
  )

  expect_type(result, "list")
  expect_named(result, c("mean_expression", "se_expression"))
  expect_equal(dim(result$mean_expression), c(n_celltypes, n_isoforms))
  expect_equal(dim(result$se_expression), c(n_celltypes, n_isoforms))

  # All values should be non-negative
  expect_true(all(result$mean_expression >= 0))
  expect_true(all(result$se_expression >= 0))
})