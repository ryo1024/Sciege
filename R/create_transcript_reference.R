#' Create transcriptome reference using Bambu
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomicRanges mcols
#' @importFrom GenomeInfoDb seqlengths
#' @param bam_files Vector of paths to BAM files (aligned reads)
#' @param annotation_gtf Path to known GTF file (e.g., GENCODE)
#' @param genome_seq Path to reference genome FASTA (e.g., GRCh38.fa)
#' @param output_dir Directory to save outputs
#' @return A SummarizedExperiment object containing the new transcript annotations
#' @export
create_transcript_reference <- function(bam_files, annotation_gtf, genome_seq, output_dir) {
  if (!requireNamespace("bambu", quietly = TRUE)) {
    stop("The 'bambu' package is required but not installed.")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  message("Running Bambu on input BAMs...")
  se <- bambu::bambu(
    reads = bam_files,
    annotations = annotation_gtf,
    genome = genome_seq,
    discovery = TRUE,
    NDR = 0.95,
  )

  #print(class(se) )

  bambu::writeBambuOutput(se, path = output_dir)
  message("Transcript reference saved to: ", output_dir)

  gtf_path <- file.path(output_dir, "extended_annotations.gtf")
  combined_path <- file.path(output_dir, "merged")

  cmd_cat <- paste( 
    "gffcompare -o ",shQuote(combined_path),
     shQuote(annotation_gtf), 
    shQuote(gtf_path)
  )

  system(cmd_cat)

  fasta_path <- file.path(output_dir, "bambu_transcript.fasta")
  combined_gtf_path <- file.path(output_dir, "merged.combined.gtf")

  cmd <- paste(
    "gffread", shQuote(combined_gtf_path),
    "-g", shQuote(genome_seq),
    "-w", shQuote(fasta_path),
    "-F"
  )
  message("Running gffread to generate transcriptome FASTA...")
  system(cmd)

  # output_path <- file.path(output_dir, "discovered_transcripts.rds")
  # saveRDS(se, output_path)


  return(list(se=se, transcript_fasta=fasta_path))
}
