#' Simulate a Pedigree for Multiparental Crosses
#'
#' This function creates a pedigree data frame that includes both founders and offspring.
#' Founders are extracted from the unique parent IDs in the parental cross matrix (with their
#' ploidy provided in \code{ploidy.vec}). Offspring are generated for each cross defined in
#' \code{parents.mat} using the corresponding number of offspring provided in \code{n_offspring}.
#' The resulting pedigree contains columns for \code{parent_1}, \code{parent_2}, \code{individual},
#' and \code{cross} (a concatenation of \code{parent_1} and \code{parent_2}).
#'
#' @param ploidy.vec A numeric vector specifying the ploidy for each parent. The vector should be named
#'   with the parent IDs (e.g., \code{c(P1 = 4, P2 = 2, P3 = 4, P4 = 2, P5 = 4, P6 = 4)}).
#' @param parents.mat A two-column matrix where each row represents a parental cross.
#'   The first column corresponds to \code{parent_1} and the second to \code{parent_2}.
#' @param n_offspring A numeric vector specifying the number of offspring for each cross.
#'   Its length must equal the number of rows in \code{parents.mat}.
#'
#' @return A data frame with the following columns:
#'   \describe{
#'     \item{parent_1}{Character; Parent 1 for the cross (or \code{NA} for founders).}
#'     \item{parent_2}{Character; Parent 2 for the cross (or \code{NA} for founders).}
#'     \item{individual}{Character; The individual ID. For founders, this is the parent's ID.
#'           For offspring, the ID is a concatenation of the parental IDs and an index (e.g., \code{"P1xP2_1"}).}
#'     \item{cross}{Character; A concatenation of \code{parent_1} and \code{parent_2} separated by "x".
#'           Note that founders will have \code{"NAxNA"} as their \code{cross} value.}
#'   }
#'
#' @examples
#' \dontrun{
#'   ploidy.vec <- c(4, 2, 4, 2, 4, 4)
#'   names(ploidy.vec) <- c("P1", "P2", "P3", "P4", "P5", "P6")
#'   parents.mat <- matrix(c("P1", "P2",
#'                           "P1", "P3",
#'                           "P2", "P2",
#'                           "P3", "P4",
#'                           "P4", "P5",
#'                           "P5", "P6",
#'                           "P5", "P2"),
#'                         ncol = 2, byrow = TRUE)
#'   n_offspring <- c(200, 30, 200, 4, 5, 3, 2)
#'   pedigree <- simulate_pedigree(ploidy.vec, parents.mat, n_offspring)
#'   head(pedigree)
#' }
#'
#' @export
simulate_pedigree <- function(ploidy.vec, parents.mat, n_offspring) {

  # For founders, both parents are NA.
  founders <- unique(as.vector(parents.mat))
  founders_df <- data.frame(
    parent_1   = rep(NA, length(founders)),
    parent_2   = rep(NA, length(founders)),
    individual = founders,
    stringsAsFactors = FALSE
  )

  # Create rows for each offspring.
  offspring_list <- mapply(function(p1, p2, n_off) {
    data.frame(
      parent_1   = rep(p1, n_off),
      parent_2   = rep(p2, n_off),
      individual = paste0(p1, "x", p2, "_", seq_len(n_off)),
      stringsAsFactors = FALSE
    )
  }, parents.mat[, 1], parents.mat[, 2], n_offspring, SIMPLIFY = FALSE)

  offspring_df <- do.call(rbind, offspring_list)

  # Combine founders and offspring into a single pedigree.
  pedigree <- rbind(founders_df, offspring_df)

  # Create the "cross" column by concatenating parent_1 and parent_2.
  pedigree$cross <- apply(pedigree[, 1:2], 1, paste0, collapse = "x")

  return(pedigree)
}
