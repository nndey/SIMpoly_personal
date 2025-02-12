#' Simulate Gametes from Homology Groups
#'
#' This function simulates the formation of gametes by modeling recombination events
#' during meiosis, using a provided homology group for a given individual.
#'
#' @param homology.group Matrix: The homology group of an individual.
#' @param n.ind Integer: The number of individuals (gametes) to simulate.
#' @param cm.map Vector: A vector of genetic map positions in centiMorgans.
#' @param prob Vector: Probability distribution for bivalent chromosome configurations (default NULL).
#' @param seed Numeric: Optional seed for reproducibility.
#'
#' @return A list containing:
#'   \describe{
#'     \item{gamete}{An array of simulated gametes.}
#'     \item{c.o.count}{A vector of crossover counts.}
#'     \item{homology.group}{The original homology group used in the simulation.}
#'   }
#' @export
simulate_gamete <- function(homology.group, n.ind, cm.map, prob = NULL, seed = NULL) {
  ploidy <- ncol(homology.group)
  n.mrk <- nrow(homology.group)
  if (!is.null(seed)) set.seed(seed)

  # Calculate recombination fractions using Haldane's mapping function
  dist.vec <- diff(cm.map)
  rf.vec <- 0.5 * (1 - exp(-dist.vec / 50))
  if (length(rf.vec) == 1) rf.vec <- rep(rf.vec, n.mrk - 1)

  res <- array(NA, c(n.mrk, ploidy / 2, n.ind))
  rf.res <- numeric(n.mrk - 1)

  # Generate all possible bivalent configurations
  a <- permutations(ploidy, ploidy, 1:ploidy)
  bv.conf <- vector("list", nrow(a))
  for (i in 1:nrow(a)) {
    temp <- apply(matrix(a[i, ], 2, ploidy / 2), 2, sort)
    bv.conf[[i]] <- temp[, order(temp[1, ]), drop = FALSE]
  }
  bv.conf <- unique(bv.conf)

  # Default probabilities if not provided
  if (is.null(prob)) prob <- rep(1 / length(bv.conf), length(bv.conf))

  # Simulate gametes for each individual
  for (k in 1:n.ind) {
    gen.1 <- matrix(1:ploidy, ploidy, n.mrk)
    choosed_biv <- sample(bv.conf, 1, prob = prob)[[1]]
    for (i in 1:ncol(choosed_biv)) {
      choosed_biv[, i] <- sample(choosed_biv[, i])
    }
    pole.1 <- choosed_biv[1, , drop = FALSE]
    pole.2 <- choosed_biv[2, , drop = FALSE]
    set.1 <- gen.1[pole.1, , drop = FALSE]
    set.2 <- gen.1[pole.2, , drop = FALSE]

    for (i in 1:(ploidy / 2)) {
      a <- set.1[i, ]
      b <- set.2[i, ]
      for (j in 1:(n.mrk - 1)) {
        if (runif(1) < rf.vec[j]) {
          which.swap <- c((j + 1):n.mrk)
          temp <- a[which.swap]
          a[which.swap] <- b[which.swap]
          b[which.swap] <- temp
        }
      }
      set.1[i, ] <- a
      set.2[i, ] <- b
    }

    # Sample a gamete product from the two poles
    gam <- if (sample(0:1, 1)) set.1 else set.2
    for (i in 1:n.mrk) {
      res[i, , k] <- as.numeric(homology.group[i, gam[, i]])
    }
  }
  return(list(gamete = res, c.o.count = rf.res, homology.group = homology.group))
}
