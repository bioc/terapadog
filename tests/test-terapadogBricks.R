library(testthat)
library(terapadog) # Adjust this path

##### Tests for prepareGeneSets() #####

set.seed(42)
mock_esetm <- matrix(runif(100), nrow = 20, dimnames = list(paste0("Gene", 1:20), NULL))

mock_gslist <- list(
  Pathway1 = paste0("Gene", 1:10),
  Pathway2 = paste0("Gene", 5:15),
  Pathway3 = paste0("Gene", 10:20)
)
mock_gs.names <- c("Pathway 1", "Pathway 2", "Pathway 3")

# Test: Organism KEGG ID must be three characters
test_that("Organism ID must be three characters", {
  expect_error(terapadog:::prepareGeneSets(mock_esetm, organism = "hs"),
               "Error: 'organism' must be a three-letter string.")
})

# Test: gslist must be a list or compatible structure
test_that("gslist must be a list", {
  expect_error(terapadog:::prepareGeneSets(mock_esetm, gslist = "not_a_list"),
               "Error: 'gslist' must be a list or compatible structure.")
})

# Test: Custom gslist requires gs.names
test_that("Custom gslist requires gs.names", {
  expect_error(terapadog:::prepareGeneSets(mock_esetm, gslist = mock_gslist, gs.names = NULL),
               "For a custom gslist, you must provide gs.names!")
})

# Test: gs.names length must match gslist
test_that("gs.names length must match gslist", {
  expect_error(terapadog:::prepareGeneSets(mock_esetm, gslist = mock_gslist,
                               gs.names = c("Only one name")),
               "Error: Length of 'gslist' and 'gs.names' must be the same.")
})

# Test: all genes have their weight
test_that("Gene frequency weighting factor calculation is correct", {
  result <- terapadog:::prepareGeneSets(mock_esetm, gslist = mock_gslist,
                            gs.names = mock_gs.names)
  expect_true(all(names(result$gf) %in% rownames(mock_esetm)))
})

# Test: At least 10 genes from esetm must be in gslist
test_that("At least 10 genes from esetm must be in gslist and no duplicates in rownames", {
  expect_error(terapadog:::prepareGeneSets(mock_esetm[1:5, , drop = FALSE],
                               gslist = mock_gslist, gs.names = mock_gs.names),
               "Error: Less than 10 genes in 'esetm' match the gene set list.")})

# Test: Gene sets smaller than Nmin are removed
test_that("Gene sets smaller than Nmin are removed", {

  valid_gslist <- list(
    Pathway1 = paste0("Gene", 1:5),
    Pathway2 = paste0("Gene", 6:7),
    Pathway3 = paste0("Gene", 9:12),
    Pathway4 = paste0("Gene", 13:20)
  )
  result <- terapadog:::prepareGeneSets(mock_esetm, gslist = valid_gslist, gs.names = c("P1", "P2", "P3", "P4"), Nmin = 3)

  # Check that small gene sets were removed
  expect_false("Pathway2" %in% names(result$gslist))
  expect_true(length(result$gslist) >= 3)
})


##### Tests for generate_permutation_matrix() ######

# Test: the matrix of permutations is not empty
test_that("combidx is not empty", {
  G <- factor(rep(c("A", "B"), each = 3))
  block <- factor(rep(1:3, 2))
  res <- terapadog:::generate_permutation_matrix(G, NI = 10, paired = FALSE, block = block, verbose = FALSE)

  combidx <- res$combidx
  expect_true(!is.null(combidx) && nrow(combidx) > 0 && ncol(combidx) > 0)
})

# Test: no duplicated permutations are kept
test_that("No duplicated combinations exist in combidx", {
  G <- factor(rep(c("A", "B"), each = 3))
  block <- factor(rep(1:3, 2))
  res <- terapadog:::generate_permutation_matrix(G, NI = 10, paired = FALSE, block = block, verbose = FALSE)

  combidx <- res$combidx
  expect_equal(ncol(combidx), length(unique(apply(combidx, 2, paste, collapse = ","))))
})

# Test: the code handles an incorrect number of permutations
test_that("Invalid number of permutations handled correctly", {
  G <- factor(rep(c("A", "B"), each = 3))
  block <- factor(rep(1:3, 2))
  res <- terapadog:::generate_permutation_matrix(G, NI = 5, paired = FALSE, block = block, verbose = FALSE)

  combidx <- res$combidx
  expect_true(ncol(combidx) == 5) # Ensure `NI` permutations are returned
})

# CTest: permutations respect "block" values, in paired experimental designs
test_that("Permutations respect block assignments", {
  G <- factor(rep(c("Wt", "Mut"), each = 3))
  block <- factor(rep(c("Pat1", "Pat2", "Pat3", "Pat1", "Pat2", "Pat3")))
  res <- terapadog:::generate_permutation_matrix(G, NI = 10, paired = TRUE, block = block, verbose = FALSE)

  combidx <- res$combidx

  # Define forbidden index pairs (which would happen if "block" not respected)
  forbidden_pairs <- list(c(1, 4), c(2, 5), c(3, 6))

  # Iterate through each permutation (column in combidx)
  for (i in seq_len(ncol(combidx))) {
    permuted_indices <- combidx[, i]

    # Check if any forbidden pair appears in this permutation
    for (pair in forbidden_pairs) {
      if (all(pair %in% permuted_indices)) {
        fail(paste("Permutation", i, "contains forbidden pair", pair[1], "and", pair[2]))
      }
    }
  }

  # If no forbidden pairs were found, the test passes
  expect_true(TRUE)
})

##### Tests for gsScoreFun() ######

# This internal function has internal checks (if ()...{stop()}).
