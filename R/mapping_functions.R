#' Inverse Haldane Mapping Function (returns cM)
#'
#' @description
#' Converts recombination frequency (\eqn{r}) into genetic distance in centiMorgans (cM)
#' according to the Haldane mapping function, assuming no crossover interference.
#'
#' @param r Numeric vector. Recombination frequencies (\eqn{0 \leq r \leq 0.5}).
#'
#' @return Numeric vector of genetic distances in centiMorgans (cM).
#'
#' @details
#' The Haldane mapping function assumes no crossover interference and relates
#' recombination fraction (\eqn{r}) to map distance (\eqn{d}) by:
#' \deqn{d = -\frac{1}{2} \log(1 - 2r)}{d = -0.5 * log(1 - 2r)}
#' in Morgans. Multiplying by 100 gives centiMorgans.
#'
#' @examples
#' imf_haldane_cM(0.1)   # 11.06 cM
#' imf_haldane_cM(c(0.01, 0.1, 0.2))  # multiple recombination frequencies
#'
#' @export
imf_haldane_cM <- function(r) {
  if (any(r < 0 | r > 0.5, na.rm = TRUE))
    stop("Input recombination fractions must be between 0 and 0.5")
  -50 * log(1 - 2 * r)
}
