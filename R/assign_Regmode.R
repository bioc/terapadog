# R/assign_Regmode.R

#'
#' This function reads the results of the function getFCs and assigns a
#' regulatory mode to each gene based on the Fold Change (FC) for RNA-Seq Counts,
#' Ribo-Seq Counts, or TE. It is intended for internal use and its output is
#' the output of get_FCs.R
#'
#' @param res_df A dataframe, output of the function getFCs
#' @return A dataframe, with two extra columns with info on the Regulatory Mode
#' @examples
#' # Internal function, code cannot be run from here.
#' \dontrun{
#' mockdata <- data.frame(
#'   Identifier = c("ENSG00000248713", "ENSG00000125780"),
#'   log2FoldChange = c(-0.69, 2),
#'   padj = c(0.16, 0.001),
#'   RIBO_FC = c(0.27, 3),
#'   RIBO_padj = c(0.45, 0.002),
#'   RNA_FC = c(1, 0.56),
#'   RNA_padj = c(0.001, 0.65)
#' )
#' result <- assign_Regmode(mockdata)
#' # Only the head of the result file will be returned
#' print(head(result))
#' }
#' @keywords internal

assign_Regmode <- function(res_df) {
  # Function checks padj (the padj value for the TE change), the RIBO_padj (same but for the FC from Ribo-Seq counts)
  # and the RNA_padj (same as above, but for RNA FC).

  # ---- Input Checks ---- #
  # Check presence of required columns
  required_columns <- c("Identifier", "padj", "RIBO_padj", "RNA_padj" , "log2FoldChange", "RIBO_FC", "RNA_FC")
  missing_columns <- setdiff(required_columns, colnames(res_df))
  if (length(missing_columns) > 0) {
    stop("Input is missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  # Empty dataframe - no rows case
  if (nrow(res_df) == 0) {
    stop("Input dataframe has no rows.")
  }

  # ---- ---- #

  res_df$RegMode <- NA
  res_df$RegModeExplicit <- NA


  #Applying the "Forwarded" RegMode (genes regulated by mRNA abundance, with no significant changes in Translational Efficiency)
  forwarded <- res_df$padj >= 0.05 & res_df$RIBO_padj < 0.05 & res_df$RNA_padj < 0.05
  res_df$RegMode[forwarded] <- "Forwarded"
  # Checks directionality of FC and assigns a more "explicit" RegMode value
  up_forwarded <- res_df$RegMode == "Forwarded" & res_df$RNA_FC > 0
  res_df$RegModeExplicit[up_forwarded] <- "(Up)regulated, driven by mRNA transcription only"
  down_forwarded <- res_df$RegMode == "Forwarded" & res_df$RNA_FC < 0
  res_df$RegModeExplicit[down_forwarded] <- "(Down)regulated, driven by mRNA transcription only"


  # Applying "Exclusive" RegMode (genes regulated by translation efficiency, with no change in mRNA abundance.
  exclusive <- res_df$padj < 0.05 & res_df$RIBO_padj < 0.05 & res_df$RNA_padj >= 0.05
  res_df$RegMode[exclusive] <- "Exclusive"
  # Checks directionality of FC and assigns a more "explicit" RegMode value
  up_exclusive <- res_df$RegMode == "Exclusive" & res_df$log2FoldChange > 0
  res_df$RegModeExplicit[up_exclusive] <- "(Up)regulated, driven by mRNA translation only"
  down_exclusive <- res_df$RegMode == "Exclusive" & res_df$log2FoldChange < 0
  res_df$RegModeExplicit[down_exclusive] <- "(Down)regulated, driven by mRNA translation only"


  # Applying "Buffered" RegMode. This requires a more complex condition. All padjs must be significant, but
  # The directionality of the Fold change between TE and RNA must be opposite.
  buffered <- (res_df$padj < 0.05 & res_df$RIBO_padj < 0.05 & res_df$RNA_padj < 0.05) &
    res_df$log2FoldChange * res_df$RNA_FC < 0
  res_df$RegMode[buffered] <- "Buffered"
  # Buffered (special case, when the effect of TE and RNA cancels out the change in RIBO)
  buffered_special <- res_df$padj < 0.05 & res_df$RIBO_padj >= 0.05 & res_df$RNA_padj < 0.05
  res_df$RegMode[buffered_special] <- "Buffered"
  # Checks directionality of FC and assigns a more "explicit" RegMode value
  buffered_mRNA_down <- res_df$RegMode == "Buffered" & res_df$RNA_FC < 0
  res_df$RegModeExplicit[buffered_mRNA_down] <- "Buffered, decrease in mRNA levels counteracted by increase in translation"
  buffered_mRNA_up <- res_df$RegMode == "Buffered" & res_df$RNA_FC > 0
  res_df$RegModeExplicit[buffered_mRNA_up] <- "Buffered, increase in mRNA levels counteracted by decrease in translation"


  # Applying "Intensified" RegMode. This requires a more complex condition. All padjs must be significant, but
  # The directionality of the Fold change between TE and RNA must be the same.
  intensified <- (res_df$padj < 0.05 & res_df$RIBO_padj < 0.05 & res_df$RNA_padj < 0.05) &
    res_df$log2FoldChange * res_df$RNA_FC > 0
  res_df$RegMode[intensified] <- "Intensified"
  # Checks directionality of FC and assigns a more "explicit" RegMode value
  intensified_down <- res_df$RegMode == "Intensified" & res_df$RNA_FC < 0
  res_df$RegModeExplicit[intensified_down] <- "Synergic decrease in both transcription and translation"
  insensified_up <- res_df$RegMode == "Intensified" & res_df$RNA_FC > 0
  res_df$RegModeExplicit[insensified_up] <- "Synergic increase in both transcription and translation"


  # No change case
  no_change <- res_df$padj >= 0.05 & res_df$RIBO_padj >= 0.05 & res_df$RNA_padj >= 0.05
  res_df$RegMode[no_change] <- "No Change"
  res_df$RegModeExplicit[no_change] <- "No significant change detected in transcription or translation"


  # Case in which any of the padj values is NA -> this is a case where a padj was not calculated because the gene
  # Did not pass a quality filter enforced by DESeq2.
  no_data <- is.na(res_df$padj) | is.na(res_df$RIBO_padj) | is.na(res_df$RNA_padj)
  res_df$RegMode[no_data] <- "Undeterminable"
  res_df$RegModeExplicit[no_data] <- "One or more adjusted p-values are missing (NA)"


  # Combination of padjs not categorised by DeltaTE (Defined as "Undetermined" by Chotani in paper 2019)
  undetermined <- (res_df$padj >= 0.05 & res_df$RIBO_padj >= 0.05 & res_df$RNA_padj < 0.05)  |
    (res_df$padj >= 0.05 & res_df$RIBO_padj < 0.05 & res_df$RNA_padj >= 0.05) |
    (res_df$padj < 0.05 & res_df$RIBO_padj >= 0.05 & res_df$RNA_padj >= 0.05)
  res_df$RegMode[undetermined] <- "Undetermined"
  res_df$RegModeExplicit[undetermined] <- "Cannot be assigned to any Regulatory Mode"

  return(res_df)
}
