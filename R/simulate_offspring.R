#' Simulate Offspring from Crosses in Pedigree
#'
#' This function simulates offspring based on the cross information in a pedigree matrix,
#' using the homology groups of founders or previous generations.
#'
#' @param pedigree Matrix: Pedigree matrix with three columns: Parent 1, Parent 2, Offspring.
#' @param founders List: Homology groups of founders.
#' @param n.mrk Integer: Number of markers.
#' @param map.length Numeric: Genetic map length.
#'
#' @return A list of homology groups for each offspring.
#' @export
#' @importFrom gtools permutations
simulate_offspring <- function(pedigree, founders, n.mrk, map.length) {
  parsed_data <- parse_pedigree(pedigree)
  crosses <- parsed_data$crosses
  homology_groups <-vector("list", nrow(pedigree))
  names(homology_groups) <- pedigree[,3]
  homology_groups[names(founders)] <- founders
  for (i in crosses) {
    offspring_name <- pedigree[i, 3]
    parent1_name <- pedigree[i, 1]
    parent2_name <- pedigree[i, 2]

    parent1_homology <- homology_groups[[parent1_name]]
    parent2_homology <- homology_groups[[parent2_name]]

    cross_result <- simulate_cross(n.ind = 1,
                                   h1 = parent1_homology,
                                   h2 = parent2_homology,
                                   cm.map = seq(1, map.length, length.out = n.mrk))

    homology_groups[[offspring_name]] <- cross_result$offspring[, , 1]
  }

  return(homology_groups)
}
