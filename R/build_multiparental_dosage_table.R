#' Build a Multiparental Dosage Data Frame from Cross Results and Pedigree
#'
#' This function processes cross results and pedigree information to create a wide-format dosage data frame.
#' In the output, rows correspond to markers (with their map positions), and columns contain dosage values
#' for each individual. The function computes a "cross" column for biparental crosses, calculates a dosage column
#' (summing across homolog columns), splits the data by group (founders or crosses), and then pivots the data
#' to a wide format. Finally, it reorders the columns so that parental columns (individuals with no recorded
#' parents in the pedigree) appear first.
#'
#' @param cross.results A data frame containing cross results. It must include at least the following columns:
#'   `parent_1`, `parent_2`, `individual`, `marker`, `map_position`, and homolog columns starting with `"Homolog_"`.
#' @param pedigree A data frame containing pedigree information with columns `parent_1`, `parent_2`, and `individual`.
#'   Individuals with `NA` for both `parent_1` and `parent_2` are considered parents.
#'
#' @return A wide-format data frame with markers as rows and dosage values as columns.
#'   The first columns are `marker` and `map_position`, followed by the parental dosage columns, and then the remaining individuals.
#'   Missing dosage values are represented as `NA`.
#'
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' \dontrun{
#'   # Assuming cross.results and pedigree are defined:
#'   wide_df <- build_multiparental_dosage(cross.results, pedigree)
#' }
#'
#' @export
build_multiparental_dosage <- function(cross.results, pedigree) {

  # Add a "cross" column (for biparental crosses) to cross.results.
  cross.results <- cross.results %>%
    mutate(cross = if_else(is.na(parent_1) | is.na(parent_2),
                           NA_character_,
                           paste0(parent_1, "x", parent_2)))

  # Create a split list by grouping individuals.
  list.of.parents.and.biaparentals <- cross.results %>%
    # Add a dosage column that sums all homolog columns
    mutate(dosage = rowSums(select(., starts_with("Homolog_")), na.rm = TRUE),
           # Create a grouping key: if cross is NA, use individual; otherwise, use cross.
           group_key = if_else(is.na(cross), as.character(individual), cross)) %>%
    group_by(group_key) %>%
    group_split()

  # Optionally, name the list elements by the group_key.
  list.of.parents.and.biaparentals <- setNames(
    list.of.parents.and.biaparentals,
    sapply(list.of.parents.and.biaparentals, function(x) unique(x$group_key))
  )

  ### Part 1: Build the Wide Dosage Data Frame

  # 1. Create a master (union) list of markers and their map positions from all tibbles.
  markers_union <- bind_rows(lapply(list.of.parents.and.biaparentals, function(df) {
    df %>% select(marker, map_position)
  })) %>% distinct(marker, map_position)

  # 2. Extract dosage information (marker, map_position, individual, dosage) from every tibble.
  all_dosage <- bind_rows(lapply(list.of.parents.and.biaparentals, function(df) {
    df %>% select(marker, map_position, individual, dosage)
  }))

  # 3. Join the master marker list with the dosage data to ensure every marker is represented.
  combined <- markers_union %>%
    full_join(all_dosage, by = c("marker", "map_position"))

  # 4. Pivot to wide format so that each unique individual becomes a column.
  wide_dosage <- combined %>%
    pivot_wider(
      id_cols = c(marker, map_position),  # each row is a unique marker
      names_from = individual,            # create one column per individual
      values_from = dosage                # fill in the dosage values
    ) %>%
    arrange(map_position)

  ### Part 2: Reorder Columns to Place Parents First

  # Extract parent's IDs from the pedigree. Parents are those individuals with NA for both parent_1 and parent_2.
  parent_ids <- pedigree %>%
    filter(is.na(parent_1) & is.na(parent_2)) %>%
    pull(individual)

  # Identify all individual dosage columns present in wide_dosage (excluding the marker info columns).
  all_indivs <- setdiff(colnames(wide_dosage), c("marker", "map_position"))

  # Determine which of these columns correspond to parents.
  parent_cols <- intersect(parent_ids, all_indivs)

  # The remaining columns (non-parents) are:
  other_cols <- setdiff(all_indivs, parent_ids)

  # Reorder the columns: marker, map_position, then parent's dosage columns, then the others.
  wide_dosage_ordered <- wide_dosage %>%
    select(marker, map_position, all_of(parent_cols), all_of(other_cols))

  return(wide_dosage_ordered)
}
