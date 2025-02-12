#' Parse Pedigree Matrix
#'
#' This function parses a pedigree matrix and identifies founders and non-founders
#' (crosses) in the dataset.
#'
#' @param pedigree Matrix: Pedigree matrix with three columns: Parent 1, Parent 2, Offspring.
#'
#' @return A list containing:
#'   \describe{
#'     \item{founders}{The row indices of founders in the pedigree.}
#'     \item{crosses}{The row indices of non-founders (crosses) in the pedigree.}
#'   }
#' @export
parse_pedigree <- function(pedigree) {
  # Identify founders (rows where both Parent 1 and Parent 2 are NA)
  founders <- which(is.na(pedigree[, 1]) & is.na(pedigree[, 2]))

  # Identify non-founders (crosses) where both parents are specified
  crosses <- which(!is.na(pedigree[, 1]) & !is.na(pedigree[, 2]))

  return(list(founders = founders, crosses = crosses))
}
