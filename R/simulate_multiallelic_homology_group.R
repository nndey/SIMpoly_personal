#' Simulate Multiallelic Homology Groups
#'
#' This function simulates homology groups for a given ploidy level and number of markers.
#' It assigns alleles to homologous chromosomes and optionally shuffles the homologs.
#'
#' @param ploidy Integer: The ploidy level of the individual.
#' @param n.mrk Integer: The number of markers to simulate.
#' @param alleles Vector: A set of alleles to assign to homologs (default is 0 to ploidy-1).
#' @param lambda Numeric: Poisson parameter controlling allele distribution.
#' @param shuffle.homolog Logical: Whether to shuffle homologs (default is FALSE).
#' @param seed Numeric: Optional seed for reproducibility.
#'
#' @return A matrix with rows as markers and columns as homologs, showing the allele assignments.
#' @export
simulate_multiallelic_homology_group <- function(ploidy, n.mrk, alleles = 0:(ploidy-1),
                                                 lambda = median(alleles), shuffle.homolog = FALSE,
                                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  P <- matrix(NA, n.mrk, ploidy)

  # Simulate allele distribution across markers
  for (i in 1:n.mrk) {
    P[i, ] <- sort(sample(alleles, size = ploidy, prob = dpois(x = alleles, lambda = lambda),
                          replace = TRUE))
  }

  # Optionally shuffle homologs
  if (shuffle.homolog) {
    P <- t(apply(P, 1, sample))
  } else {
    P <- t(apply(P, 1, sort))
  }

  # Set row and column names for the matrix
  dimnames(P) <- list(paste0("M", 1:nrow(P)), paste0("H", 1:ncol(P)))
  return(P)
}
