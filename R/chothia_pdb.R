# Chothia CDR residue ranges (inclusive) per chain type
# Heavy chain (H): CDR1 26-32, CDR2 52-56, CDR3 95-102
# Light chain (L): CDR1 24-34, CDR2 50-56, CDR3 89-97
# Insertion codes (e.g. 100A-100K) are treated as part of the CDR3 loop
.chothia_cdr_ranges <- list(
  H = list(
    CDR1 = c(26L, 32L),
    CDR2 = c(52L, 56L),
    CDR3 = c(95L, 102L)
  ),
  L = list(
    CDR1 = c(24L, 34L),
    CDR2 = c(50L, 56L),
    CDR3 = c(89L, 97L)
  )
)

# Three-letter to one-letter amino acid conversion
.aa3to1 <- c(
  ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E",
  GLY="G", HIS="H", ILE="I", LEU="L", LYS="K", MET="M", PHE="F",
  PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V",
  MSE="M", SEP="S", TPO="T", CSO="C"   # common modified residues
)


#' Read a Chothia-renumbered PDB file into a tidy tibble
#'
#' Parses ATOM records from a SAbDab Chothia-renumbered PDB file into a tidy
#' tibble. Each row represents one atom. Residue numbers with insertion codes
#' (e.g. \code{100A}, \code{82B}) are preserved as character strings, which is
#' essential for correct Chothia CDR boundary handling.
#'
#' Chain type (heavy/light/other) is inferred from the SAbDab \code{REMARK 5}
#' header lines when present.
#'
#' @param path Path to a Chothia-renumbered PDB file.
#' @param atoms Character vector of atom names to retain. Default
#'   \code{"CA"} (alpha carbons only). Pass \code{NULL} to keep all atoms.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{record}{Character. Record type: \code{"ATOM"} or \code{"HETATM"}.}
#'     \item{atom_serial}{Integer. Atom serial number.}
#'     \item{atom_name}{Character. Atom name (e.g. \code{"CA"}).}
#'     \item{res_name}{Character. Three-letter residue name.}
#'     \item{chain_id}{Character. Chain identifier.}
#'     \item{res_seq}{Character. Residue sequence number, including any
#'       insertion code (e.g. \code{"100A"}).}
#'     \item{res_num}{Integer. Numeric part of residue sequence number.}
#'     \item{ins_code}{Character. Insertion code, or \code{""} if none.}
#'     \item{x, y, z}{Numeric. Cartesian coordinates in Angstroms.}
#'     \item{aa1}{Character. One-letter amino acid code, or \code{NA} for
#'       non-standard residues.}
#'     \item{chain_type}{Character. \code{"H"} (heavy), \code{"L"} (light),
#'       or \code{"other"} — inferred from REMARK 5 lines.}
#'   }
#'
#' @export
#'
#' @examples
#' pdb_file <- system.file("extdata", "1f2x.pdb", package = "abseqr")
#' atoms <- read_chothia_pdb(pdb_file)
#' atoms
read_chothia_pdb <- function(path, atoms = "CA") {
  if (!file.exists(path)) {
    rlang::abort(paste0("File not found: ", path))
  }

  lines <- readLines(path, warn = FALSE)

  # --- Parse REMARK 5 chain-type annotations ---
  remark_lines <- lines[grepl("^REMARK   5", lines)]
  chain_type_map <- .parse_sabdab_remarks(remark_lines)

  # --- Parse ATOM / HETATM records (PDB fixed-width format) ---
  atom_lines <- lines[grepl("^(ATOM|HETATM)", lines)]
  if (length(atom_lines) == 0) {
    rlang::abort("No ATOM or HETATM records found in file.")
  }

  record     <- trimws(substr(atom_lines,  1,  6))
  atom_ser   <- suppressWarnings(as.integer(trimws(substr(atom_lines,  7, 11))))
  atom_name  <- trimws(substr(atom_lines, 13, 16))
  res_name   <- trimws(substr(atom_lines, 18, 20))
  chain_id   <- trimws(substr(atom_lines, 22, 22))
  res_seq_raw <- trimws(substr(atom_lines, 23, 27))  # includes insertion code

  x <- suppressWarnings(as.numeric(trimws(substr(atom_lines, 31, 38))))
  y <- suppressWarnings(as.numeric(trimws(substr(atom_lines, 39, 46))))
  z <- suppressWarnings(as.numeric(trimws(substr(atom_lines, 47, 54))))

  # Split residue sequence into numeric part and insertion code
  res_num  <- suppressWarnings(as.integer(gsub("[A-Za-z]", "", res_seq_raw)))
  ins_code <- gsub("[^A-Za-z]", "", res_seq_raw)

  aa1 <- .aa3to1[res_name]
  names(aa1) <- NULL

  chain_type <- chain_type_map[chain_id]
  chain_type[is.na(chain_type)] <- "other"
  names(chain_type) <- NULL

  result <- tibble::tibble(
    record      = record,
    atom_serial = atom_ser,
    atom_name   = atom_name,
    res_name    = res_name,
    chain_id    = chain_id,
    res_seq     = res_seq_raw,
    res_num     = res_num,
    ins_code    = ins_code,
    x           = x,
    y           = y,
    z           = z,
    aa1         = aa1,
    chain_type  = chain_type
  )

  # Filter to requested atom types
  if (!is.null(atoms)) {
    result <- result[result$atom_name %in% atoms, ]
  }

  result
}


#' Extract CDR loops from a Chothia-renumbered PDB tibble
#'
#' Uses standard Chothia residue number ranges to extract CDR1, CDR2, and CDR3
#' residues from the output of \code{read_chothia_pdb()}. Insertion codes
#' (e.g. residues \code{100A}–\code{100K} in the H-CDR3 loop) are
#' automatically included within the CDR3 range.
#'
#' @param pdb_atoms A tibble as returned by \code{read_chothia_pdb()}.
#' @param chain_types Character vector of chain types to extract CDRs for.
#'   Default \code{c("H", "L")}. Use \code{"H"} for heavy chain only or
#'   \code{"L"} for light chain only.
#'
#' @return A tibble with all columns from \code{pdb_atoms}, plus:
#'   \describe{
#'     \item{cdr}{Character. CDR label (\code{"CDR1"}, \code{"CDR2"},
#'       \code{"CDR3"}), or \code{NA} if the residue is not in a CDR.}
#'   }
#'   Only rows that fall within a CDR are returned.
#'
#' @export
#'
#' @examples
#' pdb_file <- system.file("extdata", "1f2x.pdb", package = "abseqr")
#' atoms <- read_chothia_pdb(pdb_file)
#' cdrs  <- extract_cdr_chothia(atoms)
#' cdrs
extract_cdr_chothia <- function(pdb_atoms,
                                chain_types = c("H", "L")) {

  if (!inherits(pdb_atoms, "data.frame")) {
    rlang::abort("`pdb_atoms` must be a tibble from read_chothia_pdb().")
  }

  required <- c("chain_type", "res_num", "ins_code")
  missing  <- setdiff(required, names(pdb_atoms))
  if (length(missing) > 0) {
    rlang::abort(paste0("Missing columns: ", paste(missing, collapse = ", "),
                        ". Run read_chothia_pdb() first."))
  }

  results <- lapply(chain_types, function(ct) {
    ranges <- .chothia_cdr_ranges[[ct]]
    if (is.null(ranges)) {
      rlang::warn(paste0("No Chothia ranges defined for chain type '", ct,
                         "'. Skipping."))
      return(NULL)
    }

    chain_atoms <- pdb_atoms[pdb_atoms$chain_type == ct, ]
    if (nrow(chain_atoms) == 0) return(NULL)

    cdr_label <- rep(NA_character_, nrow(chain_atoms))

    for (cdr_name in names(ranges)) {
      lo <- ranges[[cdr_name]][1]
      hi <- ranges[[cdr_name]][2]

      # A residue is in the CDR if its numeric part is within [lo, hi],
      # OR if it has an insertion code AND its numeric part equals hi
      # (handles e.g. 100A-100K insertions at the tip of H-CDR3)
      in_cdr <- (chain_atoms$res_num >= lo & chain_atoms$res_num <= hi) |
                (nzchar(chain_atoms$ins_code) &
                   chain_atoms$res_num == hi)

      cdr_label[in_cdr] <- cdr_name
    }

    chain_atoms$cdr <- cdr_label
    chain_atoms[!is.na(cdr_label), ]
  })

  dplyr::bind_rows(results)
}


#' Summarise CDR sequences from a Chothia PDB tibble
#'
#' Collapses per-atom CDR data (as returned by \code{extract_cdr_chothia()})
#' into one row per CDR loop, reporting the amino acid sequence, length, and
#' chain information.
#'
#' @param cdr_atoms A tibble as returned by \code{extract_cdr_chothia()}.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{chain_id}{Character. Chain identifier.}
#'     \item{chain_type}{Character. \code{"H"} or \code{"L"}.}
#'     \item{cdr}{Character. \code{"CDR1"}, \code{"CDR2"}, or \code{"CDR3"}.}
#'     \item{sequence}{Character. One-letter amino acid sequence of the loop.}
#'     \item{length}{Integer. Number of residues in the loop.}
#'   }
#'
#' @export
#'
#' @examples
#' pdb_file <- system.file("extdata", "1f2x.pdb", package = "abseqr")
#' atoms <- read_chothia_pdb(pdb_file)
#' cdrs  <- extract_cdr_chothia(atoms)
#' summarise_cdrs(cdrs)
summarise_cdrs <- function(cdr_atoms) {

  required <- c("chain_id", "chain_type", "cdr", "res_num",
                "ins_code", "aa1", "atom_name")
  missing  <- setdiff(required, names(cdr_atoms))
  if (length(missing) > 0) {
    rlang::abort(paste0("Missing columns: ", paste(missing, collapse = ", "),
                        ". Run extract_cdr_chothia() first."))
  }

  # Work with CA atoms only to get one row per residue
  ca <- cdr_atoms[cdr_atoms$atom_name == "CA", ]

  if (nrow(ca) == 0) {
    rlang::abort("No CA atoms found. Re-run read_chothia_pdb() with atoms = 'CA' (default).")
  }

  # Deduplicate to one row per unique residue position
  ca_unique <- ca[!duplicated(
    paste(ca$chain_id, ca$res_num, ca$ins_code)
  ), ]

  # Build sequences per chain × CDR
  groups <- split(ca_unique,
                  list(ca_unique$chain_id,
                       ca_unique$chain_type,
                       ca_unique$cdr),
                  drop = TRUE)

  rows <- lapply(groups, function(g) {
    # Sort by res_num then ins_code for correct ordering
    g <- g[order(g$res_num, g$ins_code), ]
    seq_str <- paste(g$aa1, collapse = "")
    tibble::tibble(
      chain_id   = g$chain_id[1],
      chain_type = g$chain_type[1],
      cdr        = g$cdr[1],
      sequence   = seq_str,
      length     = nchar(seq_str)
    )
  })

  result <- dplyr::bind_rows(rows)
  result[order(result$chain_id, result$cdr), ]
}


# Internal: parse SAbDab REMARK 5 lines into a named vector
# mapping chain_id -> chain_type ("H" or "L")
.parse_sabdab_remarks <- function(remark_lines) {
  chain_map <- character(0)

  for (line in remark_lines) {
    # Match e.g. "HCHAIN=K" or "LCHAIN=K"
    h_match <- regmatches(line, regexpr("HCHAIN=([A-Z])", line))
    l_match <- regmatches(line, regexpr("LCHAIN=([A-Z])", line))

    if (length(h_match) > 0) {
      chain_id <- sub("HCHAIN=", "", h_match)
      chain_map[chain_id] <- "H"
    }
    if (length(l_match) > 0) {
      chain_id <- sub("LCHAIN=", "", l_match)
      chain_map[chain_id] <- "L"
    }
  }

  chain_map
}
