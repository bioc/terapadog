# R/terapadogBricks.R

#' This function retrieves KEGG-based gene sets, if 'gslist' is set to "KEGGRESTpathway".
#' oOtherwise, it uses the user-supplied list of gene sets, after
#' verifying that enough genes overlap with the provided expression data.
#' It also computes a weighting factor 'gf' to down-weight genes that appear
#' in many sets.
#'
#' @importFrom KEGGREST keggLink keggList
#' @importFrom stats var quantile sd
#' @param esetm A matrix of expression data, whose rownames are gene IDs.
#'   Required for overlap checks with the sets.
#' @param gslist Either a user-supplied list of pathways,
#'   or the default string "KEGGRESTpathway", indicating that KEGG gene sets
#'   should be retrieved via the KEGGREST package.
#' @param organism A three-letter string giving the organism code for KEGG.
#' e.g. "hsa" for human. Defaults is "hsa".
#' @param gs.names Character vector with the names of the gene sets.
#' If specified, must have the same length as gslist.
#' @param Nmin  The minimum size of gene sets to be included in the analysis.
#' @param verbose Logical. If true, shows number of gene sets being worked upon.
#'
#' @return A list with three elements: gslist, gs.names, gf.
#' gslist is the list of pathways/sets.
#' gs.names is the names of each gene pathway/set.
#' gf is the vector of gene weights, computed from how frequents genes are across sets.
#'
#' @examples
#'   # Suppose 'esetm' is a matrix of counts, rownames are gene IDs.
#'   esetm <- matrix(
#'   data = c(10, 3, 25, 12, 8, 14, 7, 1, 5, 7, 2, 10, 11, 2, 27, 11, 8, 10, 6, 2, 4, 2),
#'   nrow =11,
#'   ncol = 2,
#'   dimnames = list (
#'   c("GeneA","GeneB","GeneC", "GeneD","GeneE", "GeneF", "GeneV", "GeneX","GeneY","GeneZ", "GeneW"),
#'   c("Sample1", "Sample2"))
#'   )
#'   # If we have our own sets in a list:
#'   mySets <- list(
#'     Path1 = c("GeneA","GeneB","GeneC"),
#'     Path2 = c("GeneX","GeneY","GeneZ", "GeneW"),
#'     Path3 = c("GeneA", "GeneZ", "GeneD", "GeneV"),
#'     Path4 = c("GeneA","GeneB","GeneY", "GeneE", "GeneF"),
#'     Path5 = c("GeneC", "GeneV")
#'   )
#'   gs.names <- c(Path1 = "Pathway_1",
#'   Path2 = "Pathway_2",
#'   Path3 = "Pathway_3",
#'   Path4 = "Pathway_4",
#'   Path5 = "Anne_PHathway")
#'
#'   geneSets <- terapadog:::prepareGeneSets(esetm, gslist = mySets, gs.names = gs.names)
#'
#' @keywords internal

prepareGeneSets <- function(esetm,
                            gslist = "KEGGRESTpathway",
                            organism = "hsa",
                            gs.names = NULL,
                            Nmin = 3,
                            verbose = FALSE) {

  # Retrieve KEGG sets, if user specified "KEGGRESTpathways"
  if (length(gslist) == 1 && gslist == "KEGGRESTpathway") {

    # Checks organism is given in the KEGG format (3 letters)
    if (nchar(organism) != 3) {
      stop("Error: 'organism' must be a three-letter string.")
    }

    # Retrieve pathways for the specified organism
    res.KEGG <- KEGGREST::keggLink("pathway", organism)

    # Formats output as dataframe
    a <- data.frame(path = gsub(paste0("path:", organism),"", res.KEGG),
                    gns = gsub(paste0(organism, ":"),"", names(res.KEGG))
                    )

    # maps pathways to vector of genes
    gslist <- tapply(a$gns, a$path, function(x) {
      as.character(x)
    })

    # Creates optional names for each pathways from KEGG
    gs.names <- KEGGREST::keggList("pathway", organism)[paste0(organism, names(gslist))]
    names(gs.names) <- names(gslist)

    # Checks at least 3 pathwyas are present
    if (length(gslist) < 3) {
      stop("Error: 'gslist' must contain at least 3 gene sets.")
    }

    # Housekeeping, removes variables
    rm(res.KEGG, a)
  }

  # Checks if the user provided gs.list as a list
  if (mode(gslist) != "list") {
    stop("Error: 'gslist' must be a list or compatible structure.")
  }
  # Checks if the user-provided gs.list has at least 3 pathways/sets
  if (length(gslist) < 3) {
    stop("Error: 'gslist' must contain at least 3 gene sets.")
  }
  # Checks if the user has also provided gs.names
  if (!identical(gslist, "KEGGRESTpathway") && is.null(gs.names)) {
    stop("For a custom gslist, you must provide gs.names!")
  }
  # Checks if the user-provided gs.names matches gslist
  if (!is.null(gs.names) && length(gslist) != length(gs.names)) {
    stop("Error: Length of 'gslist' and 'gs.names' must be the same.")
  }
  # Check esetm matches the user-given list
  if (sum(rownames(esetm) %in% as.character(unlist(gslist))) <= 10) {
    stop("Error: Less than 10 genes in 'esetm' match the gene set list.")
  }

  #Build the gf weighting factor (down-weight genes that appear in many sets)

  # Create table of frequencies for genes across all sets in gslist
  gf <- table(unlist(gslist))

  # If all genes DO NOT appear in the same exact number of pathways
  if (!(stats::var(gf) == 0)) {

    # Cap extreme frequencies above the 99th percentile
    if (stats::quantile(gf, 0.99) > mean(gf) + 3 * stats::sd(gf)) {
      gf[gf > stats::quantile(gf, 0.99)] <- stats::quantile(gf, 0.99)
    }

    # Creates the weight
    gff <- function(x) {
      1 + ((max(x) - x)/(max(x) - min(x)))^0.5
    }

    # Rewrites gf, so that now holds the new weights per genes
    gf <- gff(gf)

  }

  else {
    # If everything is uniform (var == 0), assign weight 1
    uniqueGenes <- unique(unlist(gslist)) # Extracts all unique genes
    gf <- rep(1, length(uniqueGenes)) # Creates a vector with same length, but just 1s
    names(gf) <- uniqueGenes # applies gene names to vector
  }

  # I think this manages genes that are not within pathways... but I want an external opinion
  allGallP <- unique(unlist(gslist))

  # Checks if some genes in esetm are missing from gf
  missedGenes <- setdiff(rownames(esetm), names(gf))
  # if missing genes are found, they are appended to gf, with a weight of 1 each.
  extraGenes <- rep(1, length(missedGenes))
  names(extraGenes) <- missedGenes
  gf <- c(gf, extraGenes)

  # Checks sets remaining match
  if (sum(rownames(esetm) %in% allGallP) <= 10) {
    stop("Error: Less than 10 genes in 'esetm' match the 'allGallP' gene set.")
  }

  # Filtering out gene sets that are too small
  gslist <- gslist[unlist(lapply(gslist, function(x) {
    length(intersect(rownames(esetm), x)) >= Nmin
  }))]
  # Retrieve gs.names
  gs.names <- gs.names[names(gslist)]

  # Sanity check on the remaining gene sets after the filtering
  if (length(gslist) < 3) {
    stop("Error: After filtering, 'gslist' must still contain at least 3 gene sets.")
  }

  if (verbose) {
    message("Analyzing ", length(gslist), " gene sets with ", Nmin, " or more genes!")
  }

  # Return the outputs: gslist, gs.names, gf
  list(
    gslist = gslist,
    gs.names = gs.names,
    gf = gf)

}

#' Generate a Permutation Matrix for Group Assignments
#'
#' This function generates a permutation matrix based on given sample groups,
#' handling paired and unpaired designs. It also returns the group in the experiment
#' with the least and most samples. Finally, it recomputes the necessary iterations,
#' based on the permutations done.
#' @importFrom utils combn
#' @param G A factor vector indicating the group assignment for each sample.
#' @param NI Integer. The maximum number of permutations allowed.
#' @param paired Logical. If TRUE, the function accounts for paired samples.
#' @param block A factor vector specifying the blocking variable for paired samples.
#' @param verbose Logical. If TRUE, the function prints the number of permutations used.
#'
#' @return A list containing:
#'   \item{combidx}{A matrix where each column is a permutation of sample indices.}
#'   \item{NI}{The number of permutations used.}
#'   \item{bigG}{A vector of group labels excluding the smallest group.}
#'   \item{minG}{The group label with the smallest number of samples.}
#'
#' @examples
#' G <- factor(rep(c("Wt", "Mut"), each = 3))
#' block <- factor(rep(c("Pat1", "Pat2", "Pat3", "Pat1", "Pat2", "Pat3")))
#' result <- terapadog:::generate_permutation_matrix(G, NI = 1000, paired = TRUE, block = block, verbose = FALSE)
#'
#' @keywords internal


generate_permutation_matrix <- function(G, NI, paired, block, verbose) {

  Glen <- length(G) # number of samples (since each factor is tied to a sample)
  tab <- table(G) # counts frequency for each factor
  idx <- which.min(tab) # finds factors with least number of samples assigned
  minG <- names(tab)[idx] # retrieve its name
  minGSZ <- tab[idx] # retrieves number of samples assigned to it
  bigG <- rep(setdiff(levels(G), minG), length(G))

  # Takes as input gr.indices, the indices of G retrieved from seq_along(G)
  combFun <- function(gr.indices, countn = TRUE) {
    # Counts instances for each level of "group" (needs to be recalculated)
    g <- G[gr.indices]
    tab <- table(g)

    if (countn) {
      # Determines size smallest level
      minsz <- min(tab)
      # If smallest group is >10, returns -1 (too many permutations to compute)
      # Otherwise, calculates all ways tio choose minsz samples out of g
      ifelse(minsz > 10, -1, choose(length(g), minsz))
    }
    else {
      # retrieves indices in g who originally belonged to minG (level with least members)
      dup <- which(g == minG)
      # generates all possible selections of minG samples out of samples.
      cms <- utils::combn(length(g), tab[minG])
      # Flags duplicated combinations and removes them (like the original)
      del <- apply(cms, 2, setequal, dup)

      if (paired) {
        cms <- cms[, order(del, decreasing = TRUE), drop = FALSE]
        cms[] <- gr.indices[c(cms)]
        cms
      }
      else {
        cms[, !del, drop = FALSE]
      }
    }
  }

  if (paired) {
    # generates the possible permutations, awared of the paired design (block)
    bct <- tapply(seq_along(G), block, combFun, simplify = TRUE)
    nperm <- ifelse(any(bct < 0), -1, prod(bct))
    # Checks if permutations is less than zero or more than the maximum allowed (NI)
    if (nperm < 0 || nperm > NI) {
      # Alternative way to generate permutations.
      # btab is a list, each element are indices of samples belonging to same block
      btab <- tapply(seq_along(G), block, `[`, simplify = FALSE)
      # Generates permutations within each block, and randomly samples instead of
      # generating all permutations
      bSamp <- function(gr.indices) {
        g <- G[gr.indices]
        tab <- table(g)
        bsz <- length(g)
        minsz <- tab[minG]
        cms <- do.call(cbind, replicate(NI, sample.int(bsz, minsz), simplify = FALSE))
        cms[] <- gr.indices[c(cms)]
        cms
      }
      # Final result is a matrix of (sampled) possible combinations.
      combidx <- do.call(rbind, lapply(btab, bSamp))
    }
    # If number of permutations is valid
    else {
      # generates valid permutations within each block.
      bcomb <- tapply(seq_along(G), block, combFun, countn = FALSE, simplify = FALSE)
      # Generates all possible combinations of different block permutations.
      colb <- expand.grid(lapply(bcomb, function(x) 1:ncol(x)))[-1, , drop = FALSE]
      # Iterates over block permutations and all their combinations. Generating combidx
      combidx <- mapply(function(x, y) x[, y, drop = FALSE], bcomb, colb, SIMPLIFY = FALSE)
      # Merges into a single matrix
      combidx <- do.call(rbind, combidx)
    }
  } else { # Handles the cares for unpaired samples
    nperm <- combFun(seq_along(G))
    #  n. permutations is not valid
    if (nperm < 0 || nperm > NI) {
      combidx <- do.call(cbind, replicate(NI, sample.int(Glen, minGSZ), simplify = FALSE))
    } else { # generates matrix of permutations
      combidx <- combFun(seq_along(G), countn = FALSE)
    }
  }

  NI <- ncol(combidx)
  if (verbose) {
    message("# of permutations used:", NI)
  }

  return(list(combidx = combidx, NI = NI, bigG = bigG, minG = minG))
}


#' Compute Gene Set Scores after Translational Efficiency Analysis
#'
#' This internal function performs gene set scoring by applying a modified
#' version of the PADOG algorithm to genes undergoing changes in translational efficiency.
#' As reported by the results of the DeltaTE package
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#'
#' @param G A factor vector representing the group for each sample.
#' @param block A factor indicating paired samples.
#' @param ite Integer, indicating the current iteration number.
#' @param exp_de A dataframe containing metadata for each sample, including grouping information.
#' @param esetm A matrix containing RNA and RIBO count data, where rows correspond to genes and columns to samples.
#' @param paired Logical, indicating whether the study design is paired.
#' @param grouped_indexes A dataframe mapping RNA and RIBO samples to their corresponding indices.
#' @param minG A character value representing the smallest group.
#' @param bigG A character vector representing the larger group.
#' @param gf A vector containing gene weighting factors.
#' @param combidx A matrix storing all possible permutations for group shuffling.
#' @param deINgs A vector with genes that are both in the gene sets, as well in esetm.
#' @param gslistINesetm A list of indices mapping genes existing both in gslist and esetm.
#'
#' @details
#' The function performs differential translational analysis using DESeq2 and calculates gene set scores
#' based on adjusted p-values from the differential analysis.
#'
#' @return A named matrix containing two rows:
#' \describe{
#'   \item{MeanAdjP}{Mean adjusted p-value for each gene set.}
#'   \item{WeightedAdjP}{Weighted adjusted p-value incorporating gene weighting factors.}
#' }
#'
#' @keywords internal

gsScoreFun <- function(G, block, ite, exp_de, esetm, paired, grouped_indexes,
                       minG, bigG, gf, combidx, deINgs, gslistINesetm) {

  force(G) # Forces evaluation of G
  force(block) # Forces evaluation of block

  # This bit re-assigns the "group" label to the samples according to the combidx matrix, if past the first iteration

  if (ite > 1) {
    G <- bigG # Need to check the original PADOG for bigG
    G[combidx[, ite - 1]] <- minG
    G <- factor(G)

    # Unpacks grouped_indexes, so to have "group" labels for each sample (RNA count and RIBO count) in a vector
    for (i in seq_along(G)) {

      # Get the new label from G
      value <- as.character(G[i])

      # Updates exp_de according to the retrieved indexes (from grouped_indexes) with the new "group" label
      exp_de[grouped_indexes[i,]$V1, ]$Group <- value
      exp_de[grouped_indexes[i,]$V2, ]$Group <- value
    }
  }

  if (paired) {
    design_TE <- ~ Block + Group + SeqType+ Group:SeqType # Paired designs do not have batch effect correction!
  }
  else { # Add if statement, if there is batch, then do this, else do not
    design_TE <- ~ Group + SeqType+ Group:SeqType
  }

  # Setup the ddsMat object
  ddsMat <- DESeq2::DESeqDataSetFromMatrix(
    countData = esetm,
    colData = exp_de,
    design = design_TE
  )

  # Calculate results (without printing messages on console)
  ddsMat <- suppressMessages(DESeq2::DESeq(ddsMat))

  # Extract specific comparison of interest
  res.TE <- DESeq2::results(ddsMat, name = "Groupd.SeqTypeRIBO",
                            cooksCutoff = FALSE,  # New
                            independentFiltering = FALSE # New
                            )

  # Checks the result of the Differential Translational Analysis contains padjs
  if (!"padj" %in% colnames(res.TE)) {
    stop("Error: DESeq2 results did not contain 'padj' column.")
  }

  # originally, moderated t-values (abs value) were retrieved and stored in a df.
  # Now it just extracts adjusted p-values (and the gene ID, obviously)
  de <- res.TE$padj
  names(de) <- rownames(res.TE)

  # The scaling happens. Matrix is generated with original padj (de) and scaled one (de*gf)
  degf <- scale(cbind(de, de / gf[names(de)]))   # degf <- scale(cbind(de, de * gf[names(de)]))
  rownames(degf) <- names(de)
  # Drops genes that are not in deINgs
  degf <- degf[deINgs, , drop = FALSE]

  # Checks if degf is empty after filtering
  if (nrow(degf) == 0) {
    stop("Error: No genes remain in 'degf' after filtering. Check genes in 'esetm' and 'gslist'.")
  }


  # Aggregates gene-level info into gene set scores. Iterates over each set in gslistINesetm
  # Each set has row indices for degf
  res <- sapply(gslistINesetm, function(z) {
    # Extracts relevant data from degf
    X <- stats::na.omit(degf[z, , drop = FALSE])
    # Computes the score
    colMeans(X, na.rm = TRUE) * sqrt(nrow(X))
  })

  # Ensure `res` is a matrix (important if `sapply` simplifies output)
  res <- as.matrix(res)

  # Checks that res has a row named "de" (more for future upkeepers, since "de" is key to assign new row names)
  if (!"de" %in% rownames(res)) {
    stop("Error: Expected row name 'de' not found in res. Check input data and scoring process.")
  }

  # Rename the first row from "de" to "MeanAdjP"
  rownames(res)[rownames(res) == "de"] <- "MeanAdjP"
  # Assign the second row as "WeightedAdjP" directly
  rownames(res)[rownames(res) != "MeanAdjP"] <- "WeightedAdjP"

  # Return the properly structured matrix
  return(res)
}
