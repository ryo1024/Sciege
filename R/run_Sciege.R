#' Run the full or partial Sciege pipeline
#'
#' @param steps Character vector of steps to run, or "all". Valid steps are:
#'        "create_transcript_reference", "estimate_isoforms", "decompose_celltypes",
#'        "estimate_ct_pop_mean", "differential_testing".
#' @param custom_inputs A named list of intermediate results to use instead of running prior steps.
#'        Names should match expected outputs: e.g., "reference", "isoform_matrix", "celltype_matrix".
#' @param config A list of configuration parameters for each step. Each step should have its own sublist.
#' @return A list of outputs from all executed steps.
#' @export
run_Sciege <- function(
    steps = "all",
    custom_inputs = list(),
    config = list()
) {
  valid_steps <- c("create_transcript_reference", "quantify_isoforms_kallisto", "decompose_celltype",
                   "estimate_ct_pop_mean", "differential_testing")

  if (identical(steps, "all")) steps <- valid_steps

  outputs <- list()

  # Step 1: Create Reference
  if ("create_transcript_reference" %in% steps) {
    message("Running: Create Reference")
    outputs$reference <- create_transcript_reference(
      fastq_dir = config$create_transcript_reference$fastq_dir,
      gtf_path = config$create_transcript_reference$gtf_path,
      output_dir = config$create_transcript_reference$output_dir
    )
  } else {
    outputs$reference <- custom_inputs$reference
  }

  # Step 2: Quantify Isoforms
  if ("quantify_isoforms_kallisto" %in% steps) {
    message("Running: Quantify Isoforms")
    outputs$isoform_matrix <- quantify_isoforms_kallisto(
      reference = outputs$reference,
      fastq_dir = config$quantify_isoforms_kallisto$fastq_dir,
      known_gtf = config$quantify_isoforms_kallisto$known_gtf,
      output_dir = config$quantify_isoforms_kallisto$output_dir
    )
  } else {
    outputs$isoform_matrix <- custom_inputs$isoform_matrix
  }

  # Step 3: Cell Type Decomposition
  if ("decompose_celltype" %in% steps) {
    message("Running: Decompose Cell Types")
    outputs$celltype_matrix <- decompose_celltypes_bisque(
      expression_matrix = outputs$isoform_matrix,
      reference_scrna = config$decompose_celltypes$reference_scrna,
      output_dir = config$decompose_celltypes$output_dir
    )
  } else {
    outputs$celltype_matrix <- custom_inputs$celltype_matrix
  }

  # Step 4: Estimate Population Mean Abundance
  if ("estimate_ct_pop_mean" %in% steps) {
    message("Running: Estimate Isoform Abundance")
    outputs$abundance_estimates <- estimate_isoform_expression_by_celltype(
      isoform_matrix = outputs$isoform_matrix,
      celltype_matrix = outputs$celltype_matrix,
      boot_iters = config$estimate_ct_pop_mean$boot_iters,
      output_dir = config$estimate_ct_pop_mean$output_dir
    )
  } else {
    outputs$abundance_estimates <- custom_inputs$abundance_estimates
  }

  # Step 5: Differential Testing (optional)
  if ("differential_testing" %in% steps) {
    message("Running: Differential Testing")
    outputs$differential_results <- run_differential_testing(
      abundance_estimates = outputs$abundance_estimates,
      condition_info = config$differential_testing$condition_info,
      output_dir = config$differential_testing$output_dir
    )
  }

  message("Pipeline complete.")
  return(outputs)
}
