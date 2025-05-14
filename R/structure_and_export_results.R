#' Structure and Export Simulation Results
#'
#' This function structures the simulation results into a data frame where each row
#' corresponds to a marker for an individual. Optionally, the results can be exported as CSV files.
#'
#' @param homology_groups List: Homology groups of all individuals.
#' @param pedigree Matrix: Pedigree matrix with three columns: Parent 1, Parent 2, Offspring.
#' @param map.length Numeric: Genetic map length.
#' @param file_path String: Optional file path for exporting the results as CSV files.
#'
#' @return A data frame containing the structured simulation results.
#' @export
structure_and_export_results <- function(homology_groups, pedigree, map.length, file_path = NULL) {
  results_list <- list()

  positions <- seq(1, map.length, length.out = nrow(homology_groups[[1]]))

  for (i in seq_along(homology_groups)) {
    individual <- names(homology_groups)[i]
    homology_group <- homology_groups[[individual]]

    homolog_cols <- paste0("Homolog_", 1:ncol(homology_group))

    individual_df <- data.frame(
      individual = rep(individual, nrow(homology_group)),
      parent_1 = pedigree[i,1],
      parent_2 = pedigree[i,2],
      marker = paste0("M", 1:nrow(homology_group)),
      position = positions,
      homology_group
    )

    colnames(individual_df)[6:ncol(individual_df)] <- homolog_cols

    results_list[[i]] <- individual_df
  }
  results_list <- pad_homologs(results_list, fixed_cols = 5)
  results_df <- do.call(rbind, results_list)
  if (!is.null(file_path)) {
    write.csv(results_df, file = file.path(file_path, "homology_groups.csv"), row.names = FALSE)
    write.csv(pedigree, file = file.path(file_path, "pedigree.csv"), row.names = FALSE)
  }
  return(results_df)
}
