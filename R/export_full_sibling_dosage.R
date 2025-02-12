#' Export Full Sibling Dosage Data to CSV
#'
#' This function exports the dosage data for full sibling groups into separate CSV files.
#' The function takes a dosage matrix and a list of full sibling groups, and generates a CSV file for each group.
#'
#' @param dosage_matrix A matrix where rows represent markers (SNPs) and columns represent individuals. The values in the matrix are dosages (e.g., allele counts).
#' @param full_sibling_groups A named list where each element is a vector of individual IDs representing a full sibling group.
#' @param file_path A string specifying the directory where the CSV files should be saved.
#'
#' @return No return value. Writes CSV files to the specified file path.
#'
#' @examples
#' # Example dosage matrix
#' dosage_matrix <- matrix(runif(100, 0, 4), nrow = 10, ncol = 10)
#' rownames(dosage_matrix) <- paste0("SNP", 1:10)
#' colnames(dosage_matrix) <- paste0("Ind", 1:10)
#'
#' # Example full sibling groups
#' full_sibs <- list("Group1" = c("Ind1", "Ind2", "Ind3"), "Group2" = c("Ind4", "Ind5"))
#'
#' # Export the dosage data
#' export_full_sibling_dosage(dosage_matrix, full_sibs, file_path = tempdir())
#'
#' @export
export_full_sibling_dosage <- function(dosage_matrix, full_sibling_groups, file_path = NULL) {

  # Ensure file_path is specified
  if (is.null(file_path)) {
    stop("Please specify a valid file path.")
  }

  founders <- full_sibling_groups[["Founders"]]
  full_sibling_groups <- full_sibling_groups[-which(names(full_sibling_groups)=="Founders")]

  # Loop through the full sibling groups
  for (group_name in names(full_sibling_groups)) {
    # Extract individuals belonging to this sibling group
    individuals <- full_sibling_groups[[group_name]]


    # Filter dosage_matrix for the individuals in this sibling group
    if (all(individuals %in% colnames(dosage_matrix))) {
      sib_group_dosage <- dosage_matrix[, individuals, drop = FALSE]
    } else {
      warning(paste("Some individuals from group", group_name, "are not found in the dosage matrix. Skipping group."))
      next
    }

    # Generate random values for other columns (these are just placeholders; adapt to your data)
    n_markers <- nrow(sib_group_dosage)
    chrom <- rep("Chr_1", n_markers)   # Placeholder chromosome data
    genome_pos <- seq(1, n_markers)    # Placeholder genome positions
    ref <- sample(c("A", "T", "G", "C"), n_markers, replace = TRUE)   # Placeholder reference alleles
    alt <- sample(c("A", "T", "G", "C"), n_markers, replace = TRUE)   # Placeholder alternate alleles

    # Create the data frame for export
    export_df <- data.frame(
      snp_id = rownames(sib_group_dosage),
      dosage_matrix[, strsplit(group_name, "_x_")[[1]], drop = FALSE],
      chrom = chrom,
      genome_pos = genome_pos,
      ref = ref,
      alt = alt,
      sib_group_dosage  # Add the filtered dosage data for the sibling group
    )

    # Construct the filename using the group name
    file_name <- paste0(group_name, "_dosage.csv")

    # Full file path
    full_file_path <- file.path(file_path, file_name)

    # Write the CSV file
    write.csv(export_df, file = full_file_path, row.names = FALSE)

    # Inform the user
    cat("Exported CSV for group:", group_name, "to", full_file_path, "\n")
  }
}
