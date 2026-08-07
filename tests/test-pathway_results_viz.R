# Unit tests for get_FCs
library(testthat)
library(terapadog)

# Tests for get_gene_from_kegg_pathway()
test_that("get_gene_from_kegg_pathway errors when 'Identifier' column is missing", {
  df <- data.frame(NotIdentifier = c(1, 2, 3))
  expect_error(
    terapadog::get_gene_from_kegg_pathway(df, pathway_id = "4350"),
    "must contain the column 'Identifier'",
    fixed = TRUE
  )
})

test_that("get_gene_from_kegg_pathway errors when 'Identifier' is not numeric or character", {
  df <- data.frame(Identifier = c(TRUE, FALSE))
  expect_error(
    terapadog::get_gene_from_kegg_pathway(df, pathway_id = "4350"),
    "must be either numeric or a character string",
    fixed = TRUE
  )
})

test_that("get_gene_from_kegg_pathway errors when organism_prefix is not 3 letters", {
  df <- data.frame(Identifier = c(1, 2, 3))
  expect_error(
    terapadog::get_gene_from_kegg_pathway(df, pathway_id = "4350",
                                          organism_prefix = "human"),
    "organism prefix is a 3-letter acronym",
    fixed = TRUE
  )
})

test_that("get_gene_from_kegg_pathway errors when pathway_id is not a valid number", {
  df <- data.frame(Identifier = c(1, 2, 3))
  expect_error(
    terapadog::get_gene_from_kegg_pathway(df, pathway_id = "abc"),
    "Provide a pathway ID as reported in the terapadog results",
    fixed = TRUE
  )
})

test_that("get_gene_from_kegg_pathway errors when pathway_id has more than 5 digits", {
  df <- data.frame(Identifier = c(1, 2, 3))
  expect_error(
    terapadog::get_gene_from_kegg_pathway(df, pathway_id = "123456"),
    "Provide a pathway ID as reported in the terapadog results",
    fixed = TRUE
  )
})

# Tests for plot_pathway_changes
test_that("plot_pathway_changes gives an error when one required column is missing", {
  # All required columns except TE_padj.
  df <- data.frame(
    Identifier = c(1, 2, 3),
    RNA_FC     = c(0, 0, 0),
    RNA_padj   = c(0.5, 0.5, 0.5),
    RIBO_FC    = c(0, 0, 0),
    RIBO_padj  = c(0.5, 0.5, 0.5),
    TE_FC      = c(0, 0, 0)
  )
  expect_error(
    terapadog::plot_pathway_changes(df, pathway_id = "4350"),
    "The output for get_FCs() is required as input dataframe.",
    fixed = TRUE
  )
})

test_that("plot_pathway_changes gives an error when 'Identifier' is not numeric or character", {
    df <- data.frame(
    Identifier = c(TRUE, FALSE, TRUE),
    RNA_FC     = c(0, 0, 0),
    RNA_padj   = c(0.5, 0.5, 0.5),
    RIBO_FC    = c(0, 0, 0),
    RIBO_padj  = c(0.5, 0.5, 0.5),
    TE_FC      = c(0, 0, 0),
    TE_padj = c(0.5, 0.5, 0.5)
  )
  expect_error(
    terapadog::plot_pathway_changes(df, pathway_id = "4350"),
    "must be either numeric or a character string",
    fixed = TRUE
  )
})

test_that("plot_pathway_changes gives an error when organism_prefix is not 3 letters", {
	df <- data.frame(
    Identifier = c(111, 112, 113),
    RNA_FC     = c(0, 0, 0),
    RNA_padj   = c(0.5, 0.5, 0.5),
    RIBO_FC    = c(0, 0, 0),
    RIBO_padj  = c(0.5, 0.5, 0.5),
    TE_FC      = c(0, 0, 0),
    TE_padj = c(0.5, 0.5, 0.5)
  )
  expect_error(
    terapadog::plot_pathway_changes(df, pathway_id = "4350",
                                          organism_prefix = "human"),
    "organism prefix is a 3-letter acronym",
    fixed = TRUE
  )
})

test_that("plot_pathway_changes gives an error when pathway_id is not a valid number", {
	df <- data.frame(
    Identifier = c(111, 112, 113),
    RNA_FC     = c(0, 0, 0),
    RNA_padj   = c(0.5, 0.5, 0.5),
    RIBO_FC    = c(0, 0, 0),
    RIBO_padj  = c(0.5, 0.5, 0.5),
    TE_FC      = c(0, 0, 0),
    TE_padj = c(0.5, 0.5, 0.5)
  )
  expect_error(
    terapadog::plot_pathway_changes(df, pathway_id = "abc"),
    "Provide a pathway ID as reported in the terapadog results",
    fixed = TRUE
  )
})

test_that("plot_pathway_changes gives an error when pathway_id has more than 5 digits", {
  	df <- data.frame(
    Identifier = c(111, 112, 113),
    RNA_FC     = c(0, 0, 0),
    RNA_padj   = c(0.5, 0.5, 0.5),
    RIBO_FC    = c(0, 0, 0),
    RIBO_padj  = c(0.5, 0.5, 0.5),
    TE_FC      = c(0, 0, 0),
    TE_padj = c(0.5, 0.5, 0.5)
  )
  expect_error(
    terapadog::plot_pathway_changes(df, pathway_id = "123456"),
    "Provide a pathway ID as reported in the terapadog results",
    fixed = TRUE
  )
})
