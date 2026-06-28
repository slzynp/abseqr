#' abseqr: Antibody and Nanobody Sequence Analysis Tools
#'
#' A toolkit for parsing, summarising, and visualising antibody and nanobody
#' sequences. Provides functions for reading FASTA files into tidy data frames,
#' extracting and analysing CDR regions both from sequences and from
#' Chothia-renumbered PDB structures, computing amino acid composition,
#' and producing publication-ready plots.
#'
#' @section FASTA sequence analysis:
#' - \code{read_fasta()}: Parse a FASTA file into a tidy tibble
#' - \code{extract_cdr3()}: Extract CDR3 regions from sequence strings
#' - \code{cdr3_length_dist()}: Plot CDR3 length distribution
#' - \code{aa_composition()}: Compute per-position amino acid composition
#' - \code{aa_composition_plot()}: Heatmap of amino acid composition
#'
#' @section Chothia PDB structure analysis:
#' - \code{read_chothia_pdb()}: Parse a SAbDab Chothia-renumbered PDB file
#' - \code{extract_cdr_chothia()}: Extract CDR1/2/3 loops by Chothia residue ranges
#' - \code{summarise_cdrs()}: Collapse CDR atoms to one sequence per loop
#'
#' @docType package
#' @name abseqr-package
"_PACKAGE"

## quiets R CMD CHECK notes for dplyr column names
utils::globalVariables(c(
  "header", "sequence", "length", "cdr3", "cdr3_length",
  "position", "amino_acid", "frequency", "n",
  "chain_id", "chain_type", "cdr", "res_num", "ins_code", "aa1", "atom_name"
))
