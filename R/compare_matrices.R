#' Compare Matrices for Maximum Similarity
#'
#' This function compares two matrices by finding the column permutation of the second matrix
#' that maximizes the similarity score with the first matrix. The similarity is measured as
#' the proportion of identical elements between the two matrices.
#'
#' @param matrix1 A matrix of reference values.
#' @param matrix2 A matrix to be compared, with the same row names as `matrix1`.
#'
#' @return A list containing:
#' \item{max_similarity}{The highest similarity score found.}
#' \item{best_permutation}{The column permutation of `matrix2` that maximizes similarity.}
#' @export
#' @importFrom gtools permutations
compare_matrices <- function(matrix1, matrix2) {
  # Get common rows by row names
  row.id <- intersect(rownames(matrix1), rownames(matrix2))
  matrix1 <- matrix1[row.id,]
  matrix2 <- matrix2[row.id,]

  # Ensure matrices have the same dimensions
  if (!all(dim(matrix1) == dim(matrix2))) {
    stop("Matrices must have the same dimensions")
  }

  n <- ncol(matrix2)
  # Generate all permutations of columns
  permutations <- gtools::permutations(n, n)

  max_similarity <- -Inf  # Track the maximum similarity score
  best_permutation <- NULL

  # Iterate over each permutation of columns
  for (i in 1:nrow(permutations)) {
    perm <- permutations[i, ]
    permuted_matrix2 <- matrix2[, perm]

    # Calculate similarity as the proportion of identical elements
    similarity <- sum(matrix1 == permuted_matrix2) / length(matrix1)

    # Update max similarity and best permutation if highest similarity so far
    if (similarity > max_similarity) {
      max_similarity <- similarity
      best_permutation <- perm
    }
  }

  # Return results: max similarity score and best permutation of columns for matrix2
  list(max_similarity = max_similarity, best_permutation = best_permutation)
}

#' Plot Difference Matrix Between Two Matrices
#'
#' This function visualizes differences between two matrices by creating a difference matrix plot.
#' The plot shows concordant and discordant elements after aligning the columns of `matrix2`
#' with the best matching permutation found by `compare_matrices`.
#'
#' @param matrix1 A reference matrix.
#' @param matrix2 A matrix to be compared with `matrix1`.
#'
#' @return A ggplot object displaying the difference matrix.
#' @export
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual labs theme_minimal theme element_text
plot_difference_matrix <- function(matrix1, matrix2) {
  # Get common rows by row names
  row.id <- intersect(rownames(matrix1), rownames(matrix2))
  matrix1 <- matrix1[row.id,]
  matrix2 <- matrix2[row.id,]

  # Use compare_matrices to align columns of matrix2
  id <- compare_matrices(matrix1, matrix2)

  # Create a difference matrix
  diff_matrix <- matrix("Concordant", nrow = nrow(matrix1), ncol = ncol(matrix1))

  # Assign values based on concordance or type of mismatch
  diff_matrix[matrix1 == 0 & matrix2[, id$best_permutation] == 1] <- "Mismatch 0->1"
  diff_matrix[matrix1 == 1 & matrix2[, id$best_permutation] == 0] <- "Mismatch 1->0"

  # Melt the difference matrix for ggplot
  df_diff <- reshape2::melt(diff_matrix)
  colnames(df_diff) <- c("Row", "Column", "Status")

  # Plot the difference matrix
  ggplot2::ggplot(df_diff, ggplot2::aes(x = Row, y = Column, fill = Status)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_manual(values = c("Concordant" = "gray", "Mismatch 0->1" = "blue", "Mismatch 1->0" = "red")) +
    ggplot2::labs(title = paste("Difference Matrix - Similarity:", round(id$max_similarity, 2)), fill = "Status") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 90, hjust = 1))
}
