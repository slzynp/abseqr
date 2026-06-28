#' Extract CDR3 regions from antibody or nanobody sequences
#'
#' Identifies and extracts CDR3 regions using conserved flanking motifs.
#' For nanobodies and VH domains the CDR3 is bounded by a conserved
#' cysteine (C) on the N-terminal side and a conserved tryptophan (W) or
#' phenylalanine (F) on the C-terminal side, following the Kabat/IMGT
#' convention.
#'
#' The default anchor pattern \code{C} ... \code{WGxG} is the canonical
#' nanobody CDR3 boundary. Use \code{read_fasta()} to load sequences first.
#'
#' @param sequences A character vector of amino acid sequences, or a tibble
#'   as returned by \code{read_fasta()} (must contain a \code{sequence} column).
#' @param c_anchor Character. Regex for the conserved N-terminal residue
#'   preceding CDR3. Default \code{"C"}.
#' @param w_anchor Character. Regex for the conserved C-terminal motif
#'   following CDR3. Default \code{"WG.G"} (matches WGQG, WGAG, etc.).
#' @param min_length Minimum CDR3 length to keep. Default \code{5L}.
#' @param max_length Maximum CDR3 length to keep. Default \code{30L}.
#'
#' @return If \code{sequences} is a tibble, returns the same tibble with two
#'   added columns:
#'   \describe{
#'     \item{cdr3}{Character. Extracted CDR3 sequence, or \code{NA} if not found.}
#'     \item{cdr3_length}{Integer. Length of CDR3, or \code{NA}.}
#'   }
#'   If \code{sequences} is a character vector, returns a tibble with columns
#'   \code{sequence}, \code{cdr3}, and \code{cdr3_length}.
#'
#' @export
#'
#' @examples
#' fasta_file <- system.file("extdata", "example_nanobodies.fasta",
#'                            package = "abseqr")
#' seqs <- read_fasta(fasta_file)
#' extract_cdr3(seqs)
#'
#' # Works on a plain character vector too
#' extract_cdr3(c("QVQLVESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGDYYCSSTSCPIWYDYWGQGTQVTVSS"))
extract_cdr3 <- function(sequences,
                         c_anchor   = "C",
                         w_anchor   = "WG.G",
                         min_length = 5L,
                         max_length = 30L) {

  is_tbl <- inherits(sequences, "data.frame")

  if (is_tbl) {
    if (!"sequence" %in% names(sequences)) {
      rlang::abort("Input tibble must contain a 'sequence' column.")
    }
    seqvec <- sequences[["sequence"]]
  } else {
    seqvec <- as.character(sequences)
  }

  pattern <- paste0(c_anchor, "([A-Z]{", min_length, ",", max_length, "})", w_anchor)

  cdr3 <- stringr::str_match(seqvec, pattern)[, 2]

  result <- tibble::tibble(
    cdr3        = cdr3,
    cdr3_length = nchar(cdr3)
  )

  if (is_tbl) {
    return(dplyr::bind_cols(sequences, result))
  }

  dplyr::bind_cols(tibble::tibble(sequence = seqvec), result)
}
