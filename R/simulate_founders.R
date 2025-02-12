#' Simulate Founders Based on Pedigree Information
#'
#' This function simulates the homology groups of founder individuals in a pedigree,
#' using specified ploidy levels, marker counts, and allele distributions.
#'
#' @param pedigree Matrix: Pedigree matrix with three columns: Parent 1, Parent 2, Offspring.
#' @param ploidy.vec Vector: A vector specifying the ploidy levels for each founder.
#' @param n.mrk Integer: Number of markers.
#' @param alleles List: A list of allele sets for each founder.
#' @param lambda Numeric: Poisson parameter for allele distribution.
#' @param shuffle.homolog Logical: Whether to shuffle homologs.
#'
#' @return A list of homology groups for each founder.
#' @export
simulate_founders <- function(pedigree, ploidy.vec, n.mrk, alleles, lambda = 1, shuffle.homolog = TRUE) {
  parsed_data <- parse_pedigree(pedigree)
  founders <- parsed_data$founders

  homology_groups <- list()

  for (i in founders) {
    founder_name <- pedigree[i, 3]
    ploidy <- ploidy.vec[founder_name]

    homology_group <- simulate_multiallelic_homology_group(ploidy = ploidy,
                                                           n.mrk = n.mrk,
                                                           alleles = alleles[[founder_name]],
                                                           lambda = lambda,
                                                           shuffle.homolog = shuffle.homolog)

    homology_groups[[founder_name]] <- homology_group
  }

  return(homology_groups)
}
