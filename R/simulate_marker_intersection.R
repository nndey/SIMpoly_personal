#' @title Generate Correlated Set Membership Data
#' @description Simulates set memberships using a correlated Bernoulli approach,
#' ensuring that every element belongs to at least one set. Returns both the Euler-compatible
#' counts and the logical membership matrix.
#' @param n Integer. The total number of elements to distribute across sets.
#' @param x Integer. The number of sets.
#' @param p Numeric vector of length `x`. The marginal probability of each element being in a given set.
#' @param rho Numeric. The correlation between sets (higher means more overlap).
#' @param multi_prob Numeric vector. Probability distribution for assigning "None" elements to 1, 2, or 3 sets.
#' @return A list containing:
#' \item{counts}{A named numeric vector of set intersections, formatted for eulerr.}
#' \item{membership_matrix}{A logical matrix (n by x), indicating set membership.}
#' @import mvtnorm
#' @import eulerr
#' @importFrom stats qnorm
#' @export
generate_correlated_sets <- function(n = 1000,
                                     x = 5,
                                     p = rep(0.3, x),
                                     rho = 0.5,
                                     multi_prob = c(0.5, 0.3, 0.2)) {
  # Create correlation matrix
  Sigma <- matrix(rho, nrow = x, ncol = x)
  diag(Sigma) <- 1  # Set variances to 1

  # Sample from correlated multivariate normal
  X <- rmvnorm(n, mean = rep(0, x), sigma = Sigma)

  # Convert to Bernoulli (TRUE/FALSE)
  membership_matrix <- X < qnorm(p)

  # Ensure every element belongs to at least one set
  none_indices <- which(rowSums(membership_matrix) == 0)  # Find "None" elements

  for (i in none_indices) {
    num_groups <- sample(1:3, 1, prob = multi_prob)  # Assign to 1, 2, or 3 sets
    random_sets <- sample(1:x, num_groups)  # Choose random sets
    membership_matrix[i, random_sets] <- TRUE
  }

  # Convert to data frame
  membership_df <- as.data.frame(membership_matrix)
  colnames(membership_df) <- paste0("Set", 1:x)

  # Generate all possible set combinations
  all_combinations <- expand.grid(rep(list(c(FALSE, TRUE)), x))

  # Count occurrences of each combination
  counts <- sapply(1:nrow(all_combinations), function(i) {
    sum(apply(membership_df, 1, function(row) all(row == unlist(all_combinations[i,]))))
  })

  # Create labels for each intersection
  region_labels <- apply(all_combinations, 1, function(row) {
    included_sets <- colnames(membership_df)[as.logical(row)]
    if (length(included_sets) == 0) return("None")  # Avoid empty names
    paste(included_sets, collapse = "&")
  })

  # Assign names properly
  names(counts) <- region_labels
  counts <- counts[counts > 0]  # Remove empty regions

  return(list(counts = counts, membership_matrix = membership_matrix))
}

#' @title Plot Euler Diagram for Correlated Sets
#' @description Plots an Euler diagram from the simulated correlated set membership data.
#' @param results A list returned by `generate_correlated_sets()`, containing `counts` and `membership_matrix`.
#' @return An Euler diagram plot.
#' @import eulerr
#' @export
plot_correlated_sets <- function(results) {

  # Extract counts from the results
  counts <- results$counts

  # Generate Euler diagram
  fit <- euler(counts)

  # Plot with numbers inside regions
  plot(fit,
       fills = list(alpha = 0.7),
       labels = list(cex = 1.2),  # Adjust text size
       quantities = TRUE)  # Show numbers in regions
}
