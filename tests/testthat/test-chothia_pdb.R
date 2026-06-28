pdb_file <- system.file("extdata", "1f2x.pdb", package = "abseqr")

test_that("read_chothia_pdb returns correct columns", {
  atoms <- read_chothia_pdb(pdb_file)
  expected <- c("record", "atom_serial", "atom_name", "res_name",
                "chain_id", "res_seq", "res_num", "ins_code",
                "x", "y", "z", "aa1", "chain_type")
  expect_named(atoms, expected)
})

test_that("read_chothia_pdb returns only CA atoms by default", {
  atoms <- read_chothia_pdb(pdb_file)
  expect_true(all(atoms$atom_name == "CA"))
})

test_that("read_chothia_pdb parses insertion codes correctly", {
  atoms <- read_chothia_pdb(pdb_file, atoms = NULL)
  insertion_rows <- atoms[nzchar(atoms$ins_code), ]
  expect_gt(nrow(insertion_rows), 0)
  # Insertion codes should be single letters
  expect_true(all(nchar(insertion_rows$ins_code) == 1))
})

test_that("read_chothia_pdb assigns chain types from REMARK 5", {
  atoms <- read_chothia_pdb(pdb_file)
  expect_true("H" %in% atoms$chain_type)
  expect_true("L" %in% atoms$chain_type)
})

test_that("read_chothia_pdb errors on missing file", {
  expect_error(read_chothia_pdb("nonexistent.pdb"))
})

test_that("extract_cdr_chothia returns only CDR residues", {
  atoms <- read_chothia_pdb(pdb_file)
  cdrs  <- extract_cdr_chothia(atoms)
  expect_true(all(cdrs$cdr %in% c("CDR1", "CDR2", "CDR3")))
  expect_false(any(is.na(cdrs$cdr)))
})

test_that("extract_cdr_chothia finds all three CDRs for both chains", {
  atoms <- read_chothia_pdb(pdb_file)
  cdrs  <- extract_cdr_chothia(atoms)
  found_cdrs <- unique(cdrs$cdr)
  expect_setequal(found_cdrs, c("CDR1", "CDR2", "CDR3"))
})

test_that("extract_cdr_chothia includes insertion code residues in CDR3", {
  atoms <- read_chothia_pdb(pdb_file, atoms = NULL)
  cdrs  <- extract_cdr_chothia(atoms)
  cdr3_ins <- cdrs[cdrs$cdr == "CDR3" & nzchar(cdrs$ins_code), ]
  expect_gt(nrow(cdr3_ins), 0)
})

test_that("summarise_cdrs returns one row per chain x CDR", {
  atoms   <- read_chothia_pdb(pdb_file)
  cdrs    <- extract_cdr_chothia(atoms)
  summary <- summarise_cdrs(cdrs)

  expect_named(summary, c("chain_id", "chain_type", "cdr", "sequence", "length"))
  # Should have 3 CDRs x 2 chains = 6 rows
  expect_equal(nrow(summary), 6L)
})

test_that("summarise_cdrs length matches nchar(sequence)", {
  atoms   <- read_chothia_pdb(pdb_file)
  cdrs    <- extract_cdr_chothia(atoms)
  summary <- summarise_cdrs(cdrs)
  expect_equal(summary$length, nchar(summary$sequence))
})
