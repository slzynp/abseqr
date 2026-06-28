#' Read a FASTA file into a tidy tibble
#'
#' Parses a FASTA file and returns a tibble with one row per sequence.
#' Multi-line sequences are collapsed automatically. Empty sequences and
#' comment lines (starting with \code{;}) are silently dropped.
#'
#' @param path Path to a FASTA file (\code{.fa}, \code{.fasta}, or \code{.fna}).
#' @param clean_headers Logical. If \code{TRUE} (default), strips leading \code{>}
#'   and trims whitespace from header lines.
#'
#' @return A \code{\link[tibble]{tibble}} with columns:
#'   \describe{
#'     \item{header}{Character. Sequence identifier (without \code{>}).}
#'     \item{sequence}{Character. Amino acid or nucleotide sequence in uppercase.}
#'     \item{length}{Integer. Number of residues/bases.}
#'   }
#'
#' @export
#'
#' @examples
#' fasta_file <- system.file("extdata", "example_nanobodies.fasta",
#'                            package = "abseqr")
#' seqs <- read_fasta(fasta_file)
#' seqs
read_fasta <- function(path, clean_headers = TRUE) {
  if (!file.exists(path)) {
    rlang::abort(paste0("File not found: ", path))
  }

  lines <- readLines(path, warn = FALSE)

  # Drop comment lines and blank lines
  lines <- lines[!grepl("^;", lines)]
  lines <- lines[nzchar(trimws(lines))]

  if (length(lines) == 0) {
    return(tibble::tibble(header = character(), sequence = character(),
                          length = integer()))
  }

  header_idx <- which(grepl("^>", lines))

  if (length(header_idx) == 0) {
    rlang::abort("No FASTA headers (lines starting with '>') found in file.")
  }

  headers <- lines[header_idx]
  if (clean_headers) {
    headers <- trimws(sub("^>", "", headers))
  }

  # Extract sequence blocks between headers
  seq_end <- c(header_idx[-1] - 1L, length(lines))
  sequences <- mapply(
    function(start, end) {
      paste(toupper(lines[(start + 1L):end]), collapse = "")
    },
    header_idx, seq_end,
    SIMPLIFY = TRUE
  )

  tibble::tibble(
    header   = headers,
    sequence = as.character(sequences),
    length   = nchar(sequences)
  )
}
