test_that("read_fasta returns correct structure", {
  fasta_file <- system.file("extdata", "example_nanobodies.fasta",
                             package = "abseqr")
  result <- read_fasta(fasta_file)

  expect_s3_class(result, "tbl_df")
  expect_named(result, c("header", "sequence", "length"))
  expect_equal(nrow(result), 20L)
})

test_that("read_fasta strips > from headers", {
  fasta_file <- system.file("extdata", "example_nanobodies.fasta",
                             package = "abseqr")
  result <- read_fasta(fasta_file)
  expect_false(any(grepl("^>", result$header)))
})

test_that("read_fasta length column matches sequence length", {
  fasta_file <- system.file("extdata", "example_nanobodies.fasta",
                             package = "abseqr")
  result <- read_fasta(fasta_file)
  expect_equal(result$length, nchar(result$sequence))
})

test_that("read_fasta sequences are uppercase", {
  fasta_file <- system.file("extdata", "example_nanobodies.fasta",
                             package = "abseqr")
  result <- read_fasta(fasta_file)
  expect_equal(result$sequence, toupper(result$sequence))
})

test_that("read_fasta errors on missing file", {
  expect_error(read_fasta("nonexistent_file.fasta"))
})

test_that("read_fasta handles multiline sequences", {
  tmp <- tempfile(fileext = ".fasta")
  writeLines(c(">seq1", "ACDEF", "GHIKL"), tmp)
  result <- read_fasta(tmp)
  expect_equal(result$sequence, "ACDEFGHIKL")
  expect_equal(result$length, 10L)
  unlink(tmp)
})
