#' Extract Full Sibling Groups
#'
#' This function extracts full sibling groups from a given pedigree matrix. Full siblings are defined as individuals
#' sharing the same parents.
#'
#' @param pedigree A data frame or matrix where columns represent Parent 1, Parent 2, and Offspring. If founders have no parents, they should be represented with NA values.
#'
#' @return A named list where each entry represents a full sibling group, named after the parent pair (or 'Founders' if both parents are NA).
#'         The elements of the list are vectors containing the IDs of the offspring that belong to each sibling group.
#'
#' @examples
#' # Define a simple pedigree matrix
#' pedigree <- data.frame(
#'   Parent1 = c(NA, NA, "F1", "F1", "O1"),
#'   Parent2 = c(NA, NA, "F2", "F2", "F3"),
#'   Offspring = c("F1", "F2", "O1", "O2", "O3")
#' )
#'
#' # Extract full sibling groups
#' extract_full_sibs(pedigree)
#'
#' @export
extract_full_sibs <- function(pedigree) {
  pedigree <- as.data.frame(pedigree)

  # Combine Parent1 and Parent2 to form a unique identifier for parent pairs
  # Ensure that Parent1 and Parent2 are in consistent order
  parent_pairs <- apply(pedigree[, c("Parent1", "Parent2")], 1, function(parents) {
    paste(sort(parents), collapse = "_x_")
  })

  # Add the parent_pairs to the pedigree
  pedigree$ParentPair <- parent_pairs

  # If no parents are found (founders), assign a 'Founders' label
  pedigree$ParentPair[pedigree$ParentPair == "NA_x_NA"] <- "Founders"
  pedigree$ParentPair[pedigree$ParentPair == ""] <- "Founders"

  # Find unique parent pairs
  unique_pairs <- unique(pedigree$ParentPair)

  # Initialize a list to store full sibling groups
  full_sib_groups <- list()

  # Loop through each unique parent pair to extract full siblings
  for (pair in unique_pairs) {
    sibs <- pedigree$Offspring[pedigree$ParentPair == pair]
    if (length(sibs) > 1) {
      full_sib_groups[[pair]] <- sibs
    }
  }

  return(full_sib_groups)
}
