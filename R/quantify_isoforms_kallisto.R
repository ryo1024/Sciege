#' Quantify isoform expression using kallisto
#'
#' @param fastq_files A named list or named character vector of FASTQ file paths (paired or single-end)
#' @param transcript_fasta Path to transcript FASTA file
#' @param output_dir Directory where kallisto output should be stored
#' @param kallisto_index Path to kallisto index (optional). If NULL, one will be created.
#' @param threads Number of threads to use
#' @param fragment_length Optional fragment length (for single-end reads)
#' @param sd Optional fragment length standard deviation (for single-end reads)
#'
#' @return A list with: abundance_matrix (TPM), and sample_ids
#' @export
quantify_isoforms_kallisto <- function(
  fastq_files,
  transcript_fasta,
  output_dir,
  kallisto_index = NULL,
  threads = 4,
  fragment_length = NULL,
  sd = NULL
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (is.null(kallisto_index)) {
    message("Building kallisto index...")
    kallisto_index <- file.path(output_dir, "transcripts.idx")
    system2("kallisto", c("index", "-i", kallisto_index, transcript_fasta))
  }

  results <- list()
  abundance_list <- list()

  for (sample_name in names(fastq_files)) {
    sample_fastqs <- fastq_files[[sample_name]]
    sample_dir <- file.path(output_dir, sample_name)
    if (!dir.exists(sample_dir)) dir.create(sample_dir)

    args <- c("quant", "-i", kallisto_index, "-o", sample_dir, "-t", threads)

    # Paired or single-end
    if (length(sample_fastqs) == 2) {
      args <- c(args, sample_fastqs)
    } else if (length(sample_fastqs) == 1 && !is.null(fragment_length) && !is.null(sd)) {
      args <- c(args, "--single", "-l", fragment_length, "-s", sd, sample_fastqs)
    } else {
      stop("Specify either paired FASTQ files or single-end with fragment_length and sd.")
    }

    message("Running kallisto for: ", sample_name)
    system2("kallisto", args)

    # Read TPM estimates
    abundance_file <- file.path(sample_dir, "abundance.tsv")
    tab <- read.delim(abundance_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    abundance_list[[sample_name]] <- setNames(tab$tpm, tab$target_id)
  }

  # Combine into a matrix (transcript x sample)
  all_ids <- unique(unlist(lapply(abundance_list, names)))
  abundance_matrix <- do.call(cbind, lapply(abundance_list, function(x) x[all_ids]))
  rownames(abundance_matrix) <- all_ids
  colnames(abundance_matrix) <- names(abundance_list)

  write.table(abundance_matrix, paste0(output_dir, '/all_abundances.tsv'))

  message("Saved results matricies to :", output_dir)

  return(list(
    abundance_matrix = abundance_matrix,
    sample_ids = names(abundance_list),
    kallisto_index = kallisto_index
  ))
}