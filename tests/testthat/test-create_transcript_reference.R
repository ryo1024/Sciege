test_that("create_reference works with real Bambu example data", {
  skip_if_not_installed("bambu")
  skip_if(Sys.which("gffread") == "", "gffread is not installed or not in PATH")

  # Use example files from bambu package
  test.bam <- system.file("extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam", package = "bambu")
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", package = "bambu")
  gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")

  output_dir <- tempfile("bambu_ref_test_")
  dir.create(output_dir)

  # Run the real function
  result <- create_transcript_reference(
    bam_files = test.bam,
    annotation_gtf = gtf.file,
    genome_seq = fa.file,
    output_dir = output_dir
  )

  # Assertions
  expect_type(result, "list")
  expect_named(result, c("se", "transcript_fasta"))
  expect_s4_class(result$se, "SummarizedExperiment")
  expect_true(file.exists(result$transcript_fasta))
  expect_gt(file.info(result$transcript_fasta)$size, 100)  # Should not be empty
})
