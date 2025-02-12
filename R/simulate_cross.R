#' Simulate Biparental Crosses
#'
#' This function simulates offspring from a biparental cross by combining gametes
#' from two parents, each represented by their homology groups.
#'
#' @param n.ind Integer: The number of offspring to generate.
#' @param h1 Matrix: Homology group for parent 1.
#' @param h2 Matrix: Homology group for parent 2.
#' @param cm.map Vector: Genetic map positions in centiMorgans.
#' @param prob1 Vector: Probability distribution for parent 1 bivalent configuration.
#' @param prob2 Vector: Probability distribution for parent 2 bivalent configuration.
#' @param seed Numeric: Optional seed for reproducibility.
#'
#' @return A list containing:
#'   \describe{
#'     \item{offspring}{An array of offspring formed by combining gametes from both parents.}
#'     \item{rf.calc}{Recombination fractions calculated for the cross.}
#'     \item{ph1}{The homology group of parent 1.}
#'     \item{ph2}{The homology group of parent 2.}
#'   }
#' @export
simulate_cross <- function(n.ind, h1, h2, cm.map, prob1 = NULL, prob2 = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Simulate gametes for both parents
  data.P1 <- simulate_gamete(h1, n.ind, cm.map, prob1)
  data.P2 <- simulate_gamete(h2, n.ind, cm.map, prob2)

  # Calculate recombination fractions
  rf.calc <- (data.P1$c.o.count + data.P2$c.o.count) / (n.ind * (ncol(h1) + ncol(h2)) / 2)

  # Combine gametes to form offspring
  offspring <- array(NA, dim = c(nrow(h1), (ncol(h1) + ncol(h2)) / 2, n.ind))
  for (i in 1:n.ind) {
    offspring[, , i] <- cbind(data.P1$gamete[, , i], data.P2$gamete[, , i])
  }

  return(list(offspring = offspring, rf.calc = rf.calc, ph1 = data.P1$homology.group, ph2 = data.P2$homology.group))
}
