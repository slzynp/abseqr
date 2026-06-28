#' Plot CDR3 length distribution
#'
#' Produces a bar chart showing the frequency of each CDR3 length in the
#' dataset. Sequences where CDR3 could not be extracted are excluded with
#' a message. Use \code{extract_cdr3()} to add CDR3 data first.
#'
#' @param data A tibble containing a \code{cdr3_length} column, as returned
#'   by \code{extract_cdr3()}.
#' @param fill Character. Bar fill colour. Default \code{"#4C8CB5"}.
#' @param title Character. Plot title. Default \code{"CDR3 Length Distribution"}.
#' @param subtitle Character or \code{NULL}. Optional subtitle.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @export
#'
#' @examples
#' fasta_file <- system.file("extdata", "example_nanobodies.fasta",
#'                            package = "abseqr")
#' seqs <- read_fasta(fasta_file) |> extract_cdr3()
#' cdr3_length_dist(seqs)
cdr3_length_dist <- function(data,
                              fill     = "#4C8CB5",
                              title    = "CDR3 Length Distribution",
                              subtitle = NULL) {

  if (!"cdr3_length" %in% names(data)) {
    rlang::abort("Input must contain a 'cdr3_length' column. Run extract_cdr3() first.")
  }

  n_missing <- sum(is.na(data$cdr3_length))
  if (n_missing > 0) {
    message(n_missing, " sequence(s) with no CDR3 detected were excluded.")
  }

  plot_data <- data[!is.na(data$cdr3_length), ]

  ggplot2::ggplot(plot_data, ggplot2::aes(x = cdr3_length)) +
    ggplot2::geom_bar(fill = fill, colour = "white", width = 0.8) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_width(1)) +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      x        = "CDR3 Length (residues)",
      y        = "Count"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold")
    )
}
