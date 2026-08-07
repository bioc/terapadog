library(testthat)
library(terapadog)

# Test: esetm is a matrix
test_that("esetm must be a matrix", {
  # load data
  rna_file <- system.file("extdata", "rna_counts.tsv", package = "terapadog")
  ribo_file <- system.file("extdata", "ribo_counts.tsv", package = "terapadog")
  sample_file <- system.file("extdata", "sample_info.tsv", package = "terapadog")
 # Use the paths to load the files.
  prepared_data <- prepareTerapadogData(rna_file, ribo_file, sample_file, "1", "2")
  # Unpacks the expression.data and exp_de from the output
  esetm <- "not_a_matrix"
  exp_de <- prepared_data$exp_de
  # create fake esetm

  expect_error(terapadog(esetm = esetm, exp_de = exp_de), "Error: esetm must be a matrix")
})

# Test: esetm has dimensions greater than 4
test_that("esetm must have dimensions greater than 4", {
  rna_file <- system.file("extdata", "rna_counts.tsv", package = "terapadog")
  ribo_file <- system.file("extdata", "ribo_counts.tsv", package = "terapadog")
  sample_file <- system.file("extdata", "sample_info.tsv", package = "terapadog")
  # Use the paths to load the files.
  prepared_data <- prepareTerapadogData(rna_file, ribo_file, sample_file, "1", "2")
  # Unpacks the expression.data and exp_de from the output
  esetm <- matrix(1:16, nrow=4, ncol=4)
  exp_de <- prepared_data$exp_de
  expect_error(terapadog(esetm = esetm, exp_de = exp_de),
               "Error: esetm must have dimensions greater than 4")
})

# Test: NI must be numeric
test_that("NI must be numeric", {
  rna_file <- system.file("extdata", "rna_counts.tsv", package = "terapadog")
  ribo_file <- system.file("extdata", "ribo_counts.tsv",package = "terapadog")
  sample_file <- system.file("extdata", "sample_info.tsv", package = "terapadog")
  # Use the paths to load the files.
  prepared_data <- prepareTerapadogData(rna_file, ribo_file, sample_file, "1", "2")
  # Unpacks the expression.data and exp_de from the output
  esetm <- prepared_data$expression.data
  exp_de <- prepared_data$exp_de
  NI <- "ten"
  expect_error(terapadog(esetm = esetm, exp_de = exp_de, NI = NI),
               "Error: NI must be numeric")
})

# Test: NI must be greater than 5
test_that("NI must be numeric", {
  rna_file <- system.file("extdata", "rna_counts.tsv", package = "terapadog")
  ribo_file <- system.file("extdata", "ribo_counts.tsv",package = "terapadog")
  sample_file <- system.file("extdata", "sample_info.tsv", package = "terapadog")
  # Use the paths to load the files.
  prepared_data <- prepareTerapadogData(rna_file, ribo_file, sample_file, "1", "2")
  # Unpacks the expression.data and exp_de from the output
  esetm <- prepared_data$expression.data
  exp_de <- prepared_data$exp_de
  NI <- 5
  expect_error(terapadog(esetm = esetm, exp_de = exp_de, NI = NI),
               "Error: NI must be greater than 5")
})

# Test: nno dupicate genes are in esetm
test_that("esetm must not have duplicate row names", {
  rna_file <- system.file("extdata", "rna_counts.tsv", package = "terapadog")
  ribo_file <- system.file("extdata", "ribo_counts.tsv",package = "terapadog")
  sample_file <- system.file("extdata", "sample_info.tsv", package = "terapadog")
  # Use the paths to load the files.
  prepared_data <- prepareTerapadogData(rna_file, ribo_file, sample_file, "1", "2")
  # Unpacks the expression.data and exp_de from the output
  esetm <- prepared_data$expression.data
  exp_de <- prepared_data$exp_de
  # Adding duplicate genes to esetm
  esetm <- rbind(esetm, esetm[c(2, 4), , drop = FALSE])
  expect_error(terapadog(esetm = esetm, exp_de = exp_de), "Error: Duplicate row names found in esetm")
})

# Test: terapadog runs without raising errors
test_that("terapadog runs without errors", {
  rna_file <- system.file("extdata", "rna_counts.tsv", package = "terapadog")
  ribo_file <- system.file("extdata", "ribo_counts.tsv",package = "terapadog")
  sample_file <- system.file("extdata", "sample_info.tsv", package = "terapadog")
  # Use the paths to load the files.
  prepared_data <- prepareTerapadogData(rna_file, ribo_file, sample_file, "1", "2")
  # Unpacks the expression.data and exp_de from the output
  esetm <- prepared_data$expression.data
  # Slices the data to reduce the computational costs of testing terapadog.
  exp_de <- prepared_data$exp_de
  # converts ids
  esetm <- id_converter(esetm, "ensembl_gene_id")
  # Assemble a mock list of pathways - to reduce access to KEGG servers and network when testing.
  genes <- rownames(esetm)
  mock_gslist <- list(
    "00001" = genes[1000:2250], # Big-sized fake-pathway
    "00002" = genes[3000:3200], # Middle-sized fake-pathway
    "00003" = genes[4000:4047], # Small-sized fake-pathway
    "00004" = genes[7000:7025] # Small-sized fake-pathway
  )
  mock_gs.names <- c("00001" = "Mock Set 1", "00002" = "Mock Set 2", "00003" = "Mock Set 3", "00004" = "Mock Set 4")
  # runs terapadog with reduced genes, iterations, and thresholds
  # To evade bioconductor's timeout issues.
  expect_silent(terapadog(esetm = esetm, exp_de = exp_de, gslist = mock_gslist, gs.names = mock_gs.names, NI = 7, Nmin = 2, verbose = FALSE))
})
