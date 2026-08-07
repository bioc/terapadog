# R/pathway_results_viz.R

#' This function retrieves, for a given pathway id, the member genes
#' GeneID (ensembl) and gene name.
#' It then extracts their rows from the results produced by get_FC, to allow
#' the manual inspection of the genes within a pathway of interest.
#' @importFrom biomaRt getBM useEnsembl
#' @importFrom dplyr %>% inner_join
#' @importFrom KEGGREST keggLink
#' @param fc_results A dataframe containing the results of the get_FC() function.
#' @param pathway_id The kegg ID for the pathway of interest, as reported in the
#' terapadog() results "ID" field.
#' @param organism_prefix Character. The three-letter code for the organism
#' of interest. default is H. sapiens ("hsa").
#' @param ensembl_dataset Character. Defines the ensembl dataset to use. Default
#' is "hsapiens_gene_ensembl".
#' @return A dataframe with the results of a Differential Translation Analysis
#' for only the genes belonging to the specified pathway.
#' The Ensembl GeneID and the gene Symbol are also provided.
#' @examples
#' \dontrun{
#' # DO NOT RUN: 1 - Queries KEGG and Ensembl - requires network access.
#' # 2 - KEGG does not allow to redistribute their data. Only research purposes.
#' # As such, this example cannot be run to build documentation.
#' fc_results_path <- system.file("extdata", "fc_results_example.csv",
#' package = "terapadog")
#' fc_results <- read.table(fc_results_path, sep = ",", header = TRUE)
#' mapped <- get_gene_from_kegg_pathway(fc_results, pathway_id = "4350")
#' }
#' @export


get_gene_from_kegg_pathway <- function(fc_results,
                                       pathway_id,
                                       organism_prefix = "hsa",
                                       ensembl_dataset = "hsapiens_gene_ensembl") {


  ####### UNIT TESTING - checks on inputs #####
  # Check formatting input dataframe
  if (!"Identifier" %in% names(fc_results)) {
    stop("Error: The dataframe must contain the column 'Identifier' as given in output by get_FC() function")
  }
  # Check columns is formatted properly
  if (!is.numeric(fc_results[["Identifier"]]) &&
      !is.character(fc_results[["Identifier"]])) {
    stop("Error: The 'Identifier' column must be either numeric or a character string.")
  }
  # Check organism prefix is correct
  if (!nchar(organism_prefix)==3) stop("Error: the organism prefix is a 3-letter acronym, like hsa")
  # Check ID for pathway is correct
  if (!grepl("^[0-9]{1,5}$", pathway_id)) stop("Error: Provide a pathway ID as reported in the terapadog results. It's a number of 5 or less digits.")


  ####### CODE ######
  # Assemble pathway ID
  pathway_id <- formatC(as.numeric(pathway_id), width = 5, flag = "0")
  pathway_id <- paste0("path:",organism_prefix, pathway_id)

  # Retrieve gene IDs for pathway and format them
  message(paste0("Retrieving information for KEGG pathway: ", pathway_id))
  genes <- keggLink(organism_prefix, pathway_id)
  kegg_genes <- unname(genes)
  entrez_ids <- sub(paste0("^", organism_prefix,":"), "", kegg_genes)

  # Check entrez IDs have been retrieved
  if (length(entrez_ids) == 0) {
    stop("No genes found for pathway ", pathway_id,
         "- worth checking inputs again.")
  }

  # Connect to the human Ensembl BioMart
  mart <- useEnsembl(biomart = "genes", dataset = ensembl_dataset)

  # Query to map entrez gene id, ensembl gene Id and gene symbol
  message("Using BioMart to retrieve ID conversion table")
  mapping_df <- tryCatch(
    {
    getBM(
    attributes = c("entrezgene_id", "ensembl_gene_id", "external_gene_name"),
    filters = "entrezgene_id",
    values = entrez_ids,
    mart = mart)
    }, error = function(msg) {stop("Ensembl (BioMart) is currently down. Try again later.")}
  )

  # Ensures compatibility when merging later
  mapping_df$entrezgene_id <- as.integer(mapping_df$entrezgene_id)
  if (is.character(fc_results[["Identifier"]])) {
    fc_results$Identifier <- as.integer(as.character(fc_results$Identifier))
    }


  # Subsect the results, keep ONLY the data for the pathway of interest. Discard the rest.
  merged_df <- fc_results %>%
    inner_join(mapping_df, by = c("Identifier" = "entrezgene_id"),
               relationship = "many-to-many")
  message("Extraction and Conversion completed! :)")

  return(merged_df)
}



#' This function plots, for a given pathway id, the member genes RNA, RIBO, and
#' TE fold changes on the relative pathway (3 plots).
#' Only significant genes are plotted.
#' Plots are saved as .png files.
#' @importFrom pathview pathview
#' @param fc_results A dataframe containing the results of the get_FC() function.
#' @param pathway_id The kegg ID for the pathway of interest, as reported in the
#' terapadog() results "ID" field.
#' @param alpha_threshold float. Defines the threshold for alpha. Values above
#' are considered non-significant.
#' @param organism_prefix Character. The three-letter code for the organism
#' of interest. default is H. sapiens ("hsa").
#' @return Saves three pictures in the working directory. Behaviour inherited
#' by pathview. Changing the workdir changes the output directory.
#' @examples
#' \dontrun{
#' # DO NOT RUN: 1 - Queries KEGG - requires network access.
#' # 2 - KEGG does not allow to redistribute their data in vignettes or docs.
#' # Only research purposes.
#' fc_results_path <- system.file("extdata", "fc_results_example.csv",
#' package = "terapadog")
#' fc_results <- read.table(fc_results_path, sep = ",", header = TRUE)
#' plot_pathway_changes(fc_results, "4350")
#' }
#' @export


plot_pathway_changes <- function(fc_results,
                                 pathway_id,
                                 alpha_threshold = 0.05,
                                 organism_prefix = "hsa") {

  # Input validation
  required_columns <- c("Identifier", "RNA_FC", "RNA_padj", "RIBO_FC",
                         "RIBO_padj", "TE_FC",  "TE_padj")
  missing_columns <- setdiff(required_columns, colnames(fc_results))
  if (length(missing_columns) > 0) stop("The output for get_FCs() is required as input dataframe.")
    # Check columns is formatted properly
  if (!is.numeric(fc_results[["Identifier"]]) &&
      !is.character(fc_results[["Identifier"]])) {
    stop("Error: The 'Identifier' column must be either numeric or a character string.")
  }
  # Check ID for pathway is correct
  if (!grepl("^[0-9]{1,5}$", pathway_id)) stop("Error: Provide a pathway ID as reported in the terapadog results. It's a number of 5 or less digits.")
  # Check organism prefix is correct
  if (!nchar(organism_prefix)==3) stop("Error: the organism prefix is a 3-letter acronym, like hsa")

  pathway_id <- formatC(as.numeric(pathway_id), width = 5, flag = "0")

  # Create mask for n.s. fold changes
  rna_not_sig <- fc_results$RNA_padj > alpha_threshold | is.na(fc_results$RNA_padj)
  ribo_not_sig  <- fc_results$RIBO_padj > alpha_threshold  | is.na(fc_results$RIBO_padj)
  te_not_sig <- fc_results$TE_padj > alpha_threshold  | is.na(fc_results$TE_padj)

  # Extract values from results
  rna_fc <- fc_results$RNA_FC
  ribo_fc <- fc_results$RIBO_FC
  te_fc<- fc_results$TE_FC
  # Apply names (Entrez gene ID)
  names(rna_fc) <- fc_results$Identifier
  names(ribo_fc) <- fc_results$Identifier
  names(te_fc) <- fc_results$Identifier

  rna_fc[rna_not_sig] <- 0
  ribo_fc[ribo_not_sig] <- 0
  te_fc[te_not_sig] <- 0

  # Generate picture for RNA
  pv_rna <- pathview::pathview(
    gene.data   = rna_fc,
    pathway.id  = pathway_id,
    species     = organism_prefix,
    gene.idtype = "ENTREZ",
    multi.state = FALSE,
    out.suffix  = "terapadog_RNA_significant_Fold_Changes",
    low  = list(gene = "gold", cpd = "blue"),
    mid  = list(gene = "gray", cpd = "gray"),
    high = list(gene = "blue", cpd = "yellow")
  )

  # Generate Picture for RIBO
  pv_ribo <- pathview::pathview(
    gene.data   = ribo_fc,
    pathway.id  = pathway_id,
    species     = organism_prefix,
    gene.idtype = "ENTREZ",
    multi.state = FALSE,
    out.suffix  = "terapadog_RIBO_significant_Fold_Changes",
    low  = list(gene = "gold", cpd = "blue"),
    mid  = list(gene = "gray", cpd = "gray"),
    high = list(gene = "blue", cpd = "yellow")
  )

  # Generate Picture for TE
  pv_te <- pathview::pathview(
    gene.data   = te_fc,
    pathway.id  = pathway_id,
    species     = organism_prefix,
    gene.idtype = "ENTREZ",
    multi.state = FALSE,
    out.suffix  = "terapadog_TE_significant_Fold_Changes",
    low  = list(gene = "gold", cpd = "blue"),
    mid  = list(gene = "gray", cpd = "gray"),
    high = list(gene = "blue", cpd = "yellow")
  )

}

