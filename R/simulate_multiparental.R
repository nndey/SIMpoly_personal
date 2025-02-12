#' Simulate Multiparental Dosage Data and Process Biparental Crosses
#'
#' This function simulates multiparental dosage data across multiple chromosomes.
#' For each chromosome, it simulates founder haplotypes, generates offspring data,
#' structures the results into a wide-format data frame (with an added chromosome identifier),
#' and then assembles a list of data frames for biparental crosses enriched with variant allele
#' information.
#'
#' @param n.chr Integer specifying the number of chromosomes to simulate.
#' @param map.len Numeric vector specifying the map length for each chromosome. Its length must equal \code{n.chr}.
#' @param pedigree A data frame containing pedigree information. It must include at least the columns
#'   \code{parent_1}, \code{parent_2}, \code{individual} and (for offspring) a \code{cross} identifier.
#'   Founders are identified by having \code{NA} in both \code{parent_1} and \code{parent_2}.
#' @param ploidy.vec Numeric vector specifying the ploidy for each founder. The names of the vector
#'   should correspond to founder IDs in the pedigree.
#' @param n_mrk Integer specifying the number of markers to simulate per chromosome.
#' @param alleles A list of allele codes for the founders (e.g. \code{list(P1 = 0:1, P2 = 0:1, ...)}).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{wide_df}}{A wide-format data frame combining dosage data from all chromosomes.
#'           Rows correspond to markers (with a \code{Chrom} column) and columns contain dosage values for individuals.
#'           The marker names are modified to indicate the chromosome.}
#'     \item{\code{dat}}{A named list of data frames, one per biparental cross, each containing:
#'           \code{snp_id} (the marker name), the dosage for the two parents (\code{P1} and \code{P2}),
#'           \code{chrom}, \code{genome_pos}, two variant allele columns (\code{alt} and \code{res}),
#'           and the dosage columns for the offspring in that cross.}
#'   }
#'
#' @import dplyr
#' @import tidyr
#' @import stringr
#'
#' @examples
#' \dontrun{
#'   # Example usage (assuming required helper functions exist):
#'   result <- simulate_multiparental_data(
#'     n.chr = 3,
#'     map.len = c(100, 120, 90),
#'     pedigree = my_pedigree,
#'     ploidy.vec = c(P1 = 4, P2 = 2, P3 = 4, P4 = 2, P5 = 4, P6 = 4),
#'     n_mrk = 200,
#'     alleles = list(P1 = 0:1, P2 = 0:1, P3 = 0:1, P4 = 0:1, P5 = 0:1, P6 = 0:1)
#'   )
#'   wide_df <- result$wide_df
#'   dat <- result$dat
#' }
#'
#' @export
simulate_multiparental_data <- function(n.chr, map.len, pedigree, ploidy.vec, n_mrk, alleles) {
  wide_df <- tibble()
  for (i in 1:n.chr) {
    set.seed(i)

    # Simulate founder haplotypes
    parent.homologs <- simulate_founders(pedigree, ploidy.vec, n_mrk, alleles)
    # Optionally, plot parent phases:
    # plot_parent_phase_gg(parent.homologs)

    # Simulate offspring for chromosome i using the corresponding map length
    offspring_sim <- simulate_offspring(pedigree, parent.homologs, n_mrk, map.len[i])

    # Structure and export the simulation results.
    cross.results <- structure_and_export_results(
      homology_groups = offspring_sim,
      pedigree = pedigree,
      map.length = map.len[i]
    )

    # Build chromosome-specific dosage data and add a Chrom column.
    chr_df <- build_multiparental_dosage(cross.results, pedigree) %>%
      mutate(Chrom = paste0("Ch_", i)) %>%
      relocate(Chrom, .after = marker) %>%
      mutate(marker = paste0(marker, "_", str_replace(Chrom, "_", "")))

    # Append the chromosome-specific data to the consolidated wide_df.
    wide_df <- bind_rows(wide_df, chr_df)
  }

  # Create a matrix of variant alleles.
  # For each marker (row in wide_df), randomly sample two alleles.
  variant <- sapply(1:nrow(wide_df), function(x) sample(c("A", "T", "C", "G"), 2))

  # Identify biparental crosses.
  # For pedigree rows where neither parent is NA, mark them as biparental.
  id <- apply(pedigree[, 1:2], 1, function(x) !any(is.na(x)))
  biparental.crosses <- unique(pedigree$cross[id])

  dat <- vector("list", length(biparental.crosses))
  names(dat) <- biparental.crosses

  for (i in seq_along(biparental.crosses)) {
    # For each biparental cross, get the offspring individual names.
    ind.names <- pedigree$individual[pedigree$cross == biparental.crosses[i]]
    # Get the two parental IDs.
    P1 <- unique(pedigree$parent_1[pedigree$cross == biparental.crosses[i]])
    P2 <- unique(pedigree$parent_2[pedigree$cross == biparental.crosses[i]])

    dat[[i]] <- data.frame(
      snp_id     = wide_df$marker,
      P1         = wide_df[[P1]],
      P2         = wide_df[[P2]],
      chrom      = wide_df$Chrom,
      genome_pos = wide_df$map_position,
      alt        = variant[1,],
      res        = variant[2,],
      wide_df[, ind.names, drop = FALSE],
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  return(list(wide_df = wide_df, dat = dat, parent_homologs = parent.homologs))
}
