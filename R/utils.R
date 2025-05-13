#' Simulate reads using wgsim
#'
#' @param fasta Path to input transcript FASTA file
#' @param out_prefix Prefix for output files
#' @param n_reads Number of reads to simulate
#' @param read_len Length of each read
#' @param seed Random seed (default 11)
#' @param wgsim_path Path to wgsim binary (default: "wgsim")
#' @return Named list of output FASTQ paths
#' @export
simulate_reads_wgsim <- function(fasta, out_prefix, n_reads = 10000, read_len = 100, seed = 11, wgsim_path = "wgsim") {
  r1 <- paste0(out_prefix, "_R1.fastq")
  r2 <- paste0(out_prefix, "_R2.fastq")

  cmd <- sprintf("%s -e 0 -1 %d -2 %d -N %d -r 0 -R 0 -X 0 -A 0 -S %d %s %s %s",
                 wgsim_path, read_len, read_len, n_reads, seed, fasta, r1, r2)

  system(cmd, intern = TRUE)

  return(list(R1 = r1, R2 = r2))
}