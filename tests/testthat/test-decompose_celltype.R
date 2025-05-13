test_that("decompose_cell_types_bisque works with simulated data", {
  skip_if_not_installed("BisqueRNA")
  skip_if_not_installed("Biobase")
  skip_if_not_installed("SummarizedExperiment")

  library(BisqueRNA)
  library(Biobase)

  output_dir <- tempfile("bambu_ref_test_")
  dir.create(output_dir)

  # Simulate reference and bulk expression
  cell.types <- c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial Cells")
  avg.props <- c(.5, .2, .2, .07, .03)
  sim <- SimulateData(n.ind=3, n.genes = 100, n.cells = 50,
          cell.types=cell.types, avg.props=avg.props)

  bulk_expr <- sim$bulk
  sc_ref <- sim$sc.eset

  #writeEset(bulk_expr, paste0(output_dir, '/bulk_expr.mat') )

  bulk_mat <- exprs(bulk_expr)
  sc_mat <- exprs(sc_ref)

  #print(sim)

  # print(bulk_mat)


  # print(bulk_expr)
  # #print(as.matrix(bulk_expr))
  # print(sc_mat)

  # Add fake cell type labels to colData (Bisque expects this)
  #colData(sc_ref)$cell_type <- sample(c("TypeA", "TypeB", "TypeC"), ncol(sc_ref), replace = TRUE)

  Biobase::pData(sc_ref)$cell_type <- sample(
    c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial Cells"),
    ncol(sc_ref),
    replace = TRUE
  )

  # Run your decomposition function
  props <- decompose_cell_types_bisque(expr_matrix = bulk_mat,
                                       reference_sce = sc_ref,
                                       out_dir = out_dir,
                                       reference_colname = "cell_type")

  # print(props)
  # print(ncol(props))

  # print(bulk_mat)
  # print(ncol(bulk_mat))

  # Check structure of result
  expect_true(is.data.frame(props))
  expect_equal(ncol(props), ncol(bulk_mat))  # samples
  expect_gt(ncol(props), 0)  # at least one cell type estimated

  # Check values are between 0 and 1
  expect_true(all(props >= 0))
  expect_true(all(props <= 1))
})