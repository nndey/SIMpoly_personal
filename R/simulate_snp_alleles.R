#' Simulate Reference and Alternate Alleles for SNPs
#'
#' Generates two character vectors of size `n`, representing reference and alternate alleles
#' for a set of SNPs. Optionally introduces transition/transversion bias for realistic simulation.
#'
#' @param n Integer. Number of SNPs to simulate.
#' @param biased Logical. If TRUE, uses realistic base frequencies and transition/transversion bias (default: TRUE).
#' @param transition_prob Numeric. Probability of choosing a transition over a transversion when `biased = TRUE` (default: 0.7).
#' @param seed Integer or NULL. Random seed for reproducibility. If NULL, randomness is not fixed.
#'
#' @return A data frame with two columns: `ref` (reference allele) and `alt` (alternate allele),
#' each of length `n`. Each element is one of "A", "T", "C", or "G", with `alt` always different from `ref`.
#'
#' @examples
#' # Simulate 10 SNPs with realistic bias
#' simulate_snp_alleles(10)
#'
#' # Simulate 5 SNPs with purely random alleles
#' simulate_snp_alleles(5, biased = FALSE)
#'
#' @export
simulate_snp_alleles <- function(n, biased = TRUE, transition_prob = 0.7, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  bases <- c("A", "T", "C", "G")

  # Choose reference alleles
  if (biased) {
    ref <- sample(bases, n, replace = TRUE, prob = c(0.3, 0.3, 0.2, 0.2))
  } else {
    ref <- sample(bases, n, replace = TRUE)
  }

  # Define transitions
  transitions <- list("A" = "G", "G" = "A", "C" = "T", "T" = "C")

  # Choose alternate alleles
  alt <- character(n)
  for (i in seq_len(n)) {
    if (biased && runif(1) < transition_prob) {
      alt[i] <- transitions[[ref[i]]]
    } else {
      alt[i] <- sample(setdiff(bases, c(ref[i], transitions[[ref[i]]])), 1)
    }
  }

  return(data.frame(ref = ref, alt = alt, stringsAsFactors = FALSE))
}
