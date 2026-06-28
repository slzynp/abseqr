#' Compute per-position amino acid composition
#'
#' Calculates the relative frequency of each amino acid at each position
#' across a set of aligned or length-filtered sequences. Useful for
#' identifying position-specific preferences in CDR3 regions.
#' Use \code{extract_cdr3()} to obtain CDR3 sequences first.
#'
#' @param data A tibble with a \code{cdr3} column (as returned by
#'   \code{extract_cdr3()}), or a character vector of sequences. All sequences
#'   must be the same length, or \code{length_filter} must be set.
#' @param col Column name (quoted or unquoted) containing the sequences to
#'   analyse. Default \code{"cdr3"}.
#' @param length_filter Integer or \code{NULL}. If set, only sequences of
#'   exactly this length are included. Useful for selecting the most common
#'   CDR3 length before computing composition.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{position}{Integer. 1-based position within the sequence.}
#'     \item{amino_acid}{Character. Single-letter amino acid code.}
#'     \item{n}{Integer. Count of sequences with this amino acid at this position.}
#'     \item{frequency}{Double. Relative frequency (0-1).}
#'   }
#'
#' @export
#'
#' @examples
#' fasta_file <- system.file("extdata", "example_nanobodies.fasta",
#'                            package = "abseqr")
#' seqs <- read_fasta(fasta_file) |> extract_cdr3()
#' aa_composition(seqs, length_filter = 12L)
aa_composition <- function(data, col = "cdr3", length_filter = NULL) {

  if (inherits(data, "data.frame")) {
    col <- rlang::as_name(rlang::ensym(col))
    if (!col %in% names(data)) {
      rlang::abort(paste0("Column '", col, "' not found in data."))
    }
    seqs <- data[[col]]
  } else {
    seqs <- as.character(data)
  }

  seqs <- seqs[!is.na(seqs) & nzchar(seqs)]

  if (!is.null(length_filter)) {
    seqs <- seqs[nchar(seqs) == length_filter]
    if (length(seqs) == 0) {
      rlang::abort(paste0("No sequences of length ", length_filter, " found."))
    }
  }

  lengths <- nchar(seqs)
  if (length(unique(lengths)) > 1) {
    rlang::abort(
      paste0("Sequences have different lengths (",
             paste(sort(unique(lengths)), collapse = ", "),
             "). Set `length_filter` to select a single length.")
    )
  }

  L <- lengths[1]
  n_seqs <- length(seqs)

  char_mat <- do.call(rbind, strsplit(seqs, ""))

  result <- lapply(seq_len(L), function(pos) {
    counts <- table(char_mat[, pos])
    tibble::tibble(
      position   = pos,
      amino_acid = names(counts),
      n          = as.integer(counts),
      frequency  = as.integer(counts) / n_seqs
    )
  })

  dplyr::bind_rows(result)
}


#' Plot per-position amino acid composition as a heatmap
#'
#' Visualises the output of \code{aa_composition()} as a heatmap where rows
#' are amino acids, columns are positions, and fill intensity represents
#' frequency. Amino acids are ordered by physicochemical property group.
#'
#' @param composition A tibble as returned by \code{aa_composition()}.
#' @param title Character. Plot title. Default \code{"Amino Acid Composition"}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @export
#'
#' @examples
#' fasta_file <- system.file("extdata", "example_nanobodies.fasta",
#'                            package = "abseqr")
#' seqs <- read_fasta(fasta_file) |> extract_cdr3()
#' comp <- aa_composition(seqs, length_filter = 12L)
#' aa_composition_plot(comp)
aa_composition_plot <- function(composition,
                                title = "Amino Acid Composition") {

  required <- c("position", "amino_acid", "frequency")
  missing  <- setdiff(required, names(composition))
  if (length(missing) > 0) {
    rlang::abort(paste0("Missing columns: ", paste(missing, collapse = ", "),
                        ". Run aa_composition() first."))
  }

  aa_order <- c(
    "R", "K", "H",
    "D", "E",
    "S", "T", "N", "Q",
    "C", "U",
    "G", "A", "V", "L", "I", "M", "F", "W", "Y", "P"
  )

  present   <- intersect(aa_order, unique(composition$amino_acid))
  extra     <- setdiff(unique(composition$amino_acid), aa_order)
  ordered_aa <- c(present, extra)

  composition <- composition |>
    dplyr::mutate(
      amino_acid = factor(amino_acid, levels = rev(ordered_aa))
    )

  ggplot2::ggplot(composition,
                  ggplot2::aes(x = position, y = amino_acid, fill = frequency)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradient(
      low    = "#F7FBFF",
      high   = "#08519C",
      name   = "Frequency",
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = scales::percent
    ) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_width(1)) +
    ggplot2::labs(title = title, x = "Position", y = "Amino Acid") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid      = ggplot2::element_blank(),
      plot.title      = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}
