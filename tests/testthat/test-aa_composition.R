test_that("aa_composition returns correct columns", {
  fasta_file <- system.file("extdata", "example_nanobodies.fasta",
                             package = "abseqr")
  result <- read_fasta(fasta_file) |>
    extract_cdr3() |>
    aa_composition(length_filter = 12L)

  expect_named(result, c("position", "amino_acid", "n", "frequency"))
})

test_that("aa_composition frequencies sum to 1 per position", {
  fasta_file <- system.file("extdata", "example_nanobodies.fasta",
                             package = "abseqr")
  result <- read_fasta(fasta_file) |>
    extract_cdr3() |>
    aa_composition(length_filter = 12L)

  sums <- tapply(result$frequency, result$position, sum)
  expect_true(all(abs(sums - 1) < 1e-9))
})

test_that("aa_composition errors on mixed-length sequences without filter", {
  seqs <- c("ACDE", "ACDEF")
  expect_error(aa_composition(seqs))
})

test_that("aa_composition errors if length_filter yields no sequences", {
  fasta_file <- system.file("extdata", "example_nanobodies.fasta",
                             package = "abseqr")
  seqs <- read_fasta(fasta_file) |> extract_cdr3()
  expect_error(aa_composition(seqs, length_filter = 999L))
})

test_that("aa_composition_plot returns ggplot object", {
  fasta_file <- system.file("extdata", "example_nanobodies.fasta",
                             package = "abseqr")
  comp <- read_fasta(fasta_file) |>
    extract_cdr3() |>
    aa_composition(length_filter = 12L)

  p <- aa_composition_plot(comp)
  expect_s3_class(p, "gg")
})

test_that("cdr3_length_dist returns ggplot object", {
  fasta_file <- system.file("extdata", "example_nanobodies.fasta",
                             package = "abseqr")
  seqs <- read_fasta(fasta_file) |> extract_cdr3()
  p <- cdr3_length_dist(seqs)
  expect_s3_class(p, "gg")
})
