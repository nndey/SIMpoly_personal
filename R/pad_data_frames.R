#' Pad Data Frames with Missing Homolog Columns
#'
#' This function takes a list of data frames that share 4 fixed columns
#' (assumed to be the first 4 columns) and additional homolog columns
#' (named "homolog_1", "homolog_2", ...). It pads each data frame with
#' extra homolog columns filled with NA so that every data frame ends up
#' having the same total number of columns.
#'
#' @param df_list A list of data frames. Each data frame must have at least 4 columns.
#' @param fixed_cols The number of fixed columns at the beginning of each data frame.
#'        Defaults to 4.
#' @param homolog_prefix The prefix for the homolog columns. Defaults to "homolog_".
#'
#' @return A list of data frames, each padded with NA columns so that they all have the same number of columns.
pad_homologs <- function(df_list, fixed_cols = 5, homolog_prefix = "Homolog_") {
  # Determine the maximum number of homolog columns present in any data frame.
  max_homolog_count <- max(sapply(df_list, function(df) {
    if (ncol(df) < fixed_cols) {
      stop("One of the data frames has fewer than the required fixed columns.")
    }
    ncol(df) - fixed_cols
  }))

  # For each data frame, add extra homolog columns (filled with NA) if needed.
  padded_list <- lapply(df_list, function(df) {
    current_homolog_count <- ncol(df) - fixed_cols
    if (current_homolog_count < max_homolog_count) {
      new_cols_count <- max_homolog_count - current_homolog_count
      # Create a data frame with the new columns, filled with NA.
      new_cols <- as.data.frame(matrix(NA, nrow = nrow(df), ncol = new_cols_count))
      # Set the column names for the new homolog columns.
      new_names <- paste0(homolog_prefix, (current_homolog_count + 1):max_homolog_count)
      names(new_cols) <- new_names
      # Append the new columns to the original data frame.
      df <- cbind(df, new_cols)
    }
    return(df)
  })

  return(padded_list)
}
