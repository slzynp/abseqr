test_that("extract_cdr3 adds cdr3 and cdr3_length columns", {
  fasta_file <- system.file("extdata", "example_nanobodies.fasta",
                             package = "abseqr")
  seqs   <- read_fasta(fasta_file)
  result <- extract_cdr3(seqs)

  expect_true("cdr3" %in% names(result))
  expect_true("cdr3_length" %in% names(result))
})

test_that("extract_cdr3 extracts non-NA CDR3 from known nanobody sequences", {
  fasta_file <- system.file("extdata", "example_nanobodies.fasta",
                             package = "abseqr")
  seqs   <- read_fasta(fasta_file)
  result <- extract_cdr3(seqs)

  n_found <- sum(!is.na(result$cdr3))
  expect_gt(n_found, 0L)
})

test_that("extract_cdr3 cdr3_length matches nchar(cdr3)", {
  fasta_file <- system.file("extdata", "example_nanobodies.fasta",
                             package = "abseqr")
  result <- read_fasta(fasta_file) |> extract_cdr3()
  found  <- result[!is.na(result$cdr3), ]

  expect_equal(found$cdr3_length, nchar(found$cdr3))
})

test_that("extract_cdr3 works on plain character vector", {
  seq <- "QVQLVESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGDYYCSSTSCPIWYDYWGQGTQVTVSS"
  result <- extract_cdr3(seq)
  expect_s3_class(result, "tbl_df")
  expect_true("cdr3" %in% names(result))
})

test_that("extract_cdr3 respects min_length and max_length", {
  fasta_file <- system.file("extdata", "example_nanobodies.fasta",
                             package = "abseqr")
  result <- read_fasta(fasta_file) |>
    extract_cdr3(min_length = 10L, max_length = 10L)

  found <- result[!is.na(result$cdr3_length), ]
  if (nrow(found) > 0) {
    expect_true(all(found$cdr3_length == 10L))
  }
})

test_that("extract_cdr3 errors if sequence column missing from tibble", {
  bad_tbl <- tibble::tibble(seq = "ACDEFG")
  expect_error(extract_cdr3(bad_tbl))
})
