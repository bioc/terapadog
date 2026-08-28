# R/terapadog.R

#' Performs the main Gene Set Enrichment Analysis, by applying a modified
#' version of the PADOG algorithm to genes undergoing changes in TE.
#' @importFrom KEGGREST keggLink keggList
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom utils combn
#' @importFrom methods is
#' @importFrom stats var quantile sd na.omit
#' @param esetm A matrix containing the counts from RNA and RIBO samples.
#' Rownames must be ensembl GENEIDs, while column names must be sample names.
#' Refer to prepareTerapadogData.R to prepare input data.
#' @param exp_de A dataframe containing information regarding the samples.
#' It has number of rows equal to the columns of esetm.
#' It has a formatted vocabulary, but can be obtained by running
#' prepareTerapadogData.R.
#' @param paired Logical. Specify is the study has a paired design or not.
#' If it does, be sure that the pairs are specified in the "Block" column of the
#' exp_de dataframe.
#' @param gslist A list of named character vectors. Each vector is named after a
#' KEGG pathway ID and each element within the vector is an ENSEMBL gene ID for
#' a gene part of said pathway.
#' @param organism A three letter string giving the name of the organism
#' supported by the "KEGGREST" package.
#' @param gs.names Character vector with the names of the gene sets.
#' If specified, must have the same length as gslist.
#' @param NI Number of iterations allowed to determine the gene set score
#' significance p-values.
#' @param Nmin The minimum size of gene sets to be included in the analysis.
#' @param verbose Logical. If true, shows number of iterations done.
#' @param ranks Logical. If true, returns also rank of each gene set.
#' @return A dataframe with the PADOG score for each pathway in exam, alongside 
#' rank value - how they rank in enrichment across all gene sets present.
#' @export
#'
terapadog <- function (esetm = NULL, exp_de = NULL, paired = FALSE,
                       gslist = "KEGGRESTpathway", organism = "hsa" ,
                       gs.names = NULL, NI = 1000, Nmin = 3, verbose = TRUE,
                       ranks = FALSE) {

  ####### UNIT TESTING - checks on inputs #####
  if (!is.matrix(esetm)) stop("Error: esetm must be a matrix")
  if (!all(dim(esetm) > 4)) stop("Error: esetm must have dimensions greater than 4")
  if (!is.numeric(NI)) stop("Error: NI must be numeric")
  if (!(NI > 5)) stop("Error: NI must be greater than 5")
  if (any(duplicated(rownames(esetm)))) stop("Error: Duplicate row names found in esetm")


  # First function call - retrieving gene set info and calculating weights
  geneSets <- terapadog:::prepareGeneSets(esetm, gslist, organism, gs.names,
                                          Nmin, verbose)
  # Extracting results
  gslist <- geneSets$gslist
  gs.names <- geneSets$gs.names
  gf <- geneSets$gf

  # This section groups matching RNA and RIBO samples under a single index (to keep the pair together during permutations).

  # Orders matching RIBO and RNA samples by the SampleName value (which is provided by the user when submitting data)
  exp_de_ordered <- exp_de[do.call(order, exp_de["SampleName"]),]
  # Retrieves the indexes (from exp_de) for each RNA and matching RIBO Sample.
  ordered_indexes <- rownames(exp_de_ordered)
  grouped_indexes <- as.data.frame(matrix(ordered_indexes, ncol = 2, byrow = TRUE))
  # Adds to the data frame info on the Group for each couple of indices
  grouped_indexes$Group <- exp_de_ordered[exp_de_ordered$SeqType != "RIBO", ]$Group
  # Sets block to NULL
  block <- NULL
  # If the study has a paired experimental design, reassigns block
  if (paired) {
    # Retrieves the info for paired samples from the dataframe
    grouped_indexes$Block <- exp_de_ordered[exp_de_ordered$SeqType != "RIBO", ]$Block
    block <- factor(grouped_indexes$Block)
  }


  # Extracts the group vector from the dataframe, to use it for permutations as PADOG does
  group <- grouped_indexes$Group
  G <- factor(group)

  # Creates permutations matrix, retrieves least and most abundant levels for G, and recomputes iterations (NI) needed
  permutations_info <- terapadog:::generate_permutation_matrix(G, NI, paired,
                                                               block, verbose)
  # Extracts results
  minG <- permutations_info$minG
  bigG <- permutations_info$bigG
  combidx <- permutations_info$combidx
  NI <- permutations_info$NI


  # Filter. Contains genes that are both in the gene sets, as well in esetm
  deINgs <- intersect(rownames(esetm), unlist(gslist))
  # For each gene in gslist, check ID against deINGs.
  # Results in a list where each gene set contains indexes corresponding to genes existing both in gslist and esetm.
  gslistINesetm <- lapply(gslist, match, table = deINgs, nomatch = 0)
  # Preparing variables to store results of iterations
  MSabsT <- MSTop <- matrix(NA, length(gslistINesetm), NI + 1)

  # Calculates scores
  for (ite in 1:(NI + 1)) {
    Sres <- terapadog:::gsScoreFun(G, block, ite, exp_de, esetm, paired,
                                   grouped_indexes, minG, bigG, gf, combidx,
                                   deINgs, gslistINesetm)
    # Extracts raw p-values
    MSabsT[, ite] <- Sres["MeanAdjP", ]
    # Extracts weighted p-values (gene set scores weighted down)
    MSTop[, ite] <- Sres["WeightedAdjP", ]

    if (verbose && (ite%%10 == 0)) {
      message(ite, "/", NI)
    }
  }

  # Extracts scores for real experiment
  meanAbsT0 <- MSabsT[, 1]
  padog0 <- MSTop[, 1]
  MSabsT <- scale(MSabsT)
  MSTop <- scale(MSTop) # Here second standardisation

  # mff is what checks the real score against the iterative ones.
  # Originally, PADOG compares t-scores (the bigger, the more significant).
  # Since we are using p-value, the smaller the more significant, so the comparison sign was changed.

  mff <- function(x) {
    mean(x[-1] < x[1], na.rm = TRUE) # Original mean(x[-1] > x[1], na.rm = TRUE)
  }
  PSabsT <- apply(MSabsT, 1, mff)
  PSTop <- apply(MSTop, 1, mff)
  PSabsT[PSabsT == 0] <- 1/(NI +1) #sets lowest p-value from 0 to 1/(NI+1) 
  PSTop[PSTop == 0] <- 1/(NI +1)

  # Compiles the final result as a dataframe
  if (!is.null(gs.names)) {
    myn <- gs.names
  }
  else {
    myn <- names(gslist)
  }
  SIZE <- unlist(lapply(gslist, function(x) {
    length(intersect(rownames(esetm), x))
  }))
  
  res <- data.frame(Name = myn, ID = names(gslist), Size = SIZE,
                    meanAbsT0, padog0, PmeanAbsT = PSabsT, Ppadog = PSTop,
                    stringsAsFactors = FALSE)
  
  if (ranks == TRUE) {
    res$S0star <- MSTop[, 1]
    res$rank <- rank(res$S0star, ties.method = "average")
    res$S0star <- NULL
  }

  ord <- order(res$Ppadog, -res$padog0)
  res <- res[ord, ]
  return(res)
}
