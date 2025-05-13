test_that("quantify_isoforms_kallisto returns valid outputs", {
  skip_if_not(Sys.which("kallisto") != "", "kallisto not installed")

  fastq_files <- list(
    sample1 = c("data/s1_chr9_bambu_simulated_R1.fastq", "data/s1_chr9_bambu_simulated_R2.fastq"),
    sample2 = c("data/s2_chr9_bambu_simulated_R1.fastq", "data/s2_chr9_bambu_simulated_R2.fastq"),
    sample3 = c("data/s3_chr9_bambu_simulated_R1.fastq", "data/s3_chr9_bambu_simulated_R2.fastq")
  )

  output_dir <- tempfile("kallisto_quant_test_")
  dir.create(output_dir)

  result <- quantify_isoforms_kallisto(
    fastq_files = fastq_files,
    transcript_fasta = "data/bambu_transcript.fasta",
    output_dir = output_dir,
    threads = 1
  )

  expect_type(result, "list")
  expect_true("abundance_matrix" %in% names(result))
  expect_gt(nrow(result$abundance_matrix), 0)
  expect_equal(ncol(result$abundance_matrix), length(fastq_files))
})