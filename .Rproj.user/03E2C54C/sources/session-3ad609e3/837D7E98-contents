## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 7,
  fig.height = 4,
  dpi = 150
)

## ----load---------------------------------------------------------------------
library(abseqr)

fasta_file <- system.file("extdata", "example_nanobodies.fasta",
                           package = "abseqr")
seqs <- read_fasta(fasta_file)
seqs

## ----cdr3---------------------------------------------------------------------
seqs <- extract_cdr3(seqs)
seqs |> dplyr::select(header, cdr3, cdr3_length)

## ----cdr3-summary-------------------------------------------------------------
n_total <- nrow(seqs)
n_found <- sum(!is.na(seqs$cdr3))
cat(sprintf("CDR3 extracted: %d / %d sequences (%.0f%%)\n",
            n_found, n_total, 100 * n_found / n_total))

