# R/prepareTerapadogData.R

#' Prepare Data by Loading and Validating RNA, RIBO Counts, and Metadata.
#' This function reads RNA and RIBO count files, checks input data validity and
#' merges them into a single numerical matrix (expression.data).
#' It also prepares the metatadata needed by padog (exp_de).
#' @importFrom utils read.table
#' @param path_to_RNA_counts A string representing the file path
#' to the RNA counts data file (.csv or .tsv).
#' @param path_to_RIBO_counts A string representing the file path
#' to the RIBO counts data file (.csv or .tsv).
#' @param path_to_metadata The file path to the metadata file (.csv or .tsv).
#' @param analysis.group.1 A string specifying the baseline group in
#' the experiment (e.g., WT, control, etc.).
#' @param analysis.group.2 A string specifying the target group to compare
#' against the baseline (e.g., mutant, disease, treatment, etc.).
#' @return A list containing two data frames: expression.data and exp_de.
#' @examples
#' # Data is also available in the "extdata" folder of this package.
#' # The path will be automatically generated for the purpose of this example
#' rna_file <- system.file("extdata", "rna_counts.tsv",
#' package = "terapadog")
#' ribo_file <- system.file("extdata", "ribo_counts.tsv",
#' package = "terapadog")
#' sample_file <- system.file("extdata", "sample_info.tsv",
#'  package = "terapadog")
#'  # Use the paths to load the files.
#' prepared_data <- prepareTerapadogData(rna_file, ribo_file,
#' sample_file, "1", "2")
#' # Unpacks the expression.data and exp_de from the output
#' expression.data <- prepared_data$expression.data
#' exp_de <- prepared_data$exp_de
#' # For sake of brevity, only the data frame's head will be printed out
#' print(head(expression.data))
#' print(head(exp_de))
#' @export
#'
prepareTerapadogData <- function (path_to_RNA_counts, path_to_RIBO_counts,
                                  path_to_metadata,
                                  analysis.group.1, analysis.group.2) {

  # ---- Input Checks ---- #
  # Check existence of files given as input
  if (!file.exists(path_to_RNA_counts)) stop("RNA counts file does not exist: ", path_to_RNA_counts)
  if (!file.exists(path_to_RIBO_counts)) stop("RIBO counts file does not exist: ", path_to_RIBO_counts)
  if (!file.exists(path_to_metadata)) stop("Metadata file does not exist: ", path_to_metadata)

  # ---- ---- #


  # Read the file using the detected separator
  RNA_counts <- utils::read.table(path_to_RNA_counts, sep = detect_separator(
    path_to_RNA_counts), header = TRUE)
  RIBO_counts <- utils::read.table(path_to_RIBO_counts, sep = detect_separator(
    path_to_RIBO_counts), header = TRUE)

  # Read the metadata
  metadata <- utils::read.table(path_to_metadata, sep = detect_separator(
    path_to_metadata), header = TRUE)

  # ---- Medatata Formatting Checks ---- #
  # Checks required columns are present in metadata file
  required_columns <- c("SampleName", "Condition")
  missing_cols <- setdiff(required_columns, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop("Metadata file is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  # Check that analysis.group.1 and analysis.group.2 are actually in metadata
  unique_conditions <- unique(metadata$Condition)
  if (!(analysis.group.1 %in% unique_conditions) | !(analysis.group.2 %in% unique_conditions)) {
    stop("analysis.group.1 or analysis.group.2 do not match Condition values available in Sample_info file.")
  }

  # ---- ---- #

  # Check if dataframes are empty (no rows, no columns, or NULL)
  check_input_df(RNA_counts)
  check_input_df(RIBO_counts)

  # Check colnames are matching
  check_matching_colnames(RNA_counts, RIBO_counts)

  # Check range of data, if normalised between o and 1, stops execution.
  check_value_range(RNA_counts)
  check_value_range(RIBO_counts)

  # Check if data is integer. If floats are present, rounds them to int.
  RNA_counts <- check_integer_values(RNA_counts)
  RIBO_counts <- check_integer_values(RIBO_counts)

  # Get original colnames and concatenate (will be "SampleName" in exp_de)
  RNA_SampleNames <- colnames(RNA_counts)
  RIBO_SampleNames <- colnames(RIBO_counts)
  SampleName <- c(RNA_SampleNames, RIBO_SampleNames)

  # Assigns a unique name to each column (will be SampleID in exp_de)
  colnames(RNA_counts) <- paste0(colnames(RNA_counts), "_RNA")
  colnames(RIBO_counts) <- paste0(colnames(RIBO_counts), "_RIBO")
  SampleID <- c(colnames(RNA_counts), colnames(RIBO_counts))

  # Creates vector of SeqType
  rna_vector <- rep("RNA", ncol(RNA_counts))
  ribo_vector <- rep("RIBO", ncol(RIBO_counts))
  SeqType <- c(rna_vector, ribo_vector) #NB. check if rows names are considered a column or not!

  # Merge the dataframes and convert to matrix
  expression.data <- cbind(RNA_counts, RIBO_counts)
  expression.data <- as.matrix(expression.data)

  # Create exp_de by putting together SampleID, SampleName, SeqType
  exp_de <- data.frame(SampleID, SampleName, SeqType)
  # Assign Group/Condition to exp_de based on value of SampleName in metadata
  exp_de$Condition <- metadata$Condition[match(exp_de$SampleName, metadata$SampleName)]
  # If there is info on paired design (Block column) add it to the dataframe
  if ("Block" %in% colnames(metadata)) {
    exp_de$Block <- metadata$Block[match(exp_de$SampleName, metadata$SampleName)]
  }

  # Reformat baseline and target to "c" and "d"
  baseline <- ifelse(exp_de$Condition == analysis.group.1, 1, # was grepl
                     ifelse(exp_de$Condition == analysis.group.2, 0, NA))
  target <- ifelse(exp_de$Condition == analysis.group.2, 1, # was grepl
                   ifelse(exp_de$Condition == analysis.group.1,0, NA))

  design <- data.frame(baseline, target)

  # Defines Group from selected analysis groups
  exp_de$Group <- NA
  exp_de$Group[design$baseline == 1] <- "c"
  exp_de$Group[design$target == 1] <- "d"

  # Ensures important columns are converted to factors and order of comparison.
  exp_de$Group <- factor(exp_de$Group, levels = c("c", "d"))
  exp_de$SeqType <- factor(exp_de$SeqType, levels = c("RNA", "RIBO"))

  # Return output
  return(list(expression.data = expression.data, exp_de = exp_de))
}
