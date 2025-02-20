#' @title Simulate Multiparental Dosage Data and Process Biparental Crosses
#'
#' @description Simulates multiparental dosage data across multiple chromosomes.
#' It generates founder haplotypes, simulates offspring, structures the results
#' into a wide-format data frame, and extracts biparental cross data enriched
#' with variant allele information.
#'
#' @param n.chr Integer. Number of chromosomes to simulate.
#' @param map.len Numeric vector of length \code{n.chr}. Specifies the map length for each chromosome.
#' @param pedigree Data frame. Must contain at least the columns \code{parent_1}, \code{parent_2},
#'   \code{individual}, and \code{cross}. Founders are identified by \code{NA} in both parental columns.
#' @param ploidy.vec Named numeric vector. Specifies the ploidy level for each founder.
#'   Names should correspond to founder IDs in \code{pedigree}.
#' @param n_mrk Integer vector. Specifies the number of markers per chromosome.
#' @param alleles List. Contains allele codes for the founders (e.g., \code{list(P1 = 0:1, P2 = 0:1, ...)}).
#' @param missing Numeric. Proportion of missing data to introduce. Must be in the range \code{[0.0, 1.0)}.
#'   Default is \code{0.0} (no missing data).
#' @param p Numeric. Probability of an element being in each biparental cross when generating correlated sets.
#'   Default is \code{0.3}.
#' @param rho Numeric. Correlation parameter for set overlap when generating correlated biparental crosses.
#'   Higher values increase overlap. Default is \code{0.1}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{wide_df}}{Wide-format data frame containing dosage data from all chromosomes.
#'         Rows correspond to markers (with a \code{Chrom} column) and columns contain dosage
#'         values for individuals. Marker names indicate the chromosome.}
#'   \item{\code{dat}}{Named list of data frames, one per biparental cross, each containing:
#'         \code{snp_id} (marker name), dosage for the two parents (\code{P1} and \code{P2}),
#'         \code{chrom}, \code{genome_pos}, two variant allele columns (\code{alt} and \code{res}),
#'         and the dosage columns for offspring in that cross.}
#'   \item{\code{parent_homologs}}{List of founder haplotype simulations for each chromosome.}
#'   \item{\code{counts}}{Named numeric vector of set intersections (formatted for Euler diagrams).}
#' }
#'
#' @import dplyr
#' @import tidyr
#' @import stringr
#'
#' @examples
#' \dontrun{
#'   result <- simulate_multiparental_data(
#'     n.chr = 3,
#'     map.len = c(100, 120, 90),
#'     pedigree = my_pedigree,
#'     ploidy.vec = c(P1 = 4, P2 = 2, P3 = 4, P4 = 2, P5 = 4, P6 = 4),
#'     n_mrk = c(200, 250, 180),
#'     alleles = list(P1 = 0:1, P2 = 0:1, P3 = 0:1, P4 = 0:1, P5 = 0:1, P6 = 0:1),
#'     missing = 0.05,
#'     p = 0.3,
#'     rho = 0.2
#'   )
#'   wide_df <- result$wide_df
#'   dat <- result$dat
#' }
#'
#' @export
simulate_multiparental_data <- function(n.chr,
                                        map.len,
                                        pedigree,
                                        ploidy.vec,
                                        n_mrk,
                                        alleles,
                                        missing = 0.0,
                                        p = 0.3,
                                        rho = 0.1) {
  wide_df <- tibble()
  parent_homologs <- vector("list", n.chr)
  stopifnot("'missing' must be in [0.0, 1.0)." = missing >= 0.0 & missing < 1.0)
  for (i in seq_len(n.chr)) {
    set.seed(i)  # Ensure reproducibility

    # Simulate founder haplotypes
    parent_homologs[[i]] <- simulate_founders(pedigree, ploidy.vec, n_mrk[i], alleles)

    # Simulate offspring for chromosome i using the corresponding map length
    offspring_sim <- simulate_offspring(pedigree, parent_homologs[[i]], n_mrk[i], map.len[i])

    # Structure and export the simulation results
    cross_results <- structure_and_export_results(
      homology_groups = offspring_sim,
      pedigree = pedigree,
      map.length = map.len[i]
    )

    # Build chromosome-specific dosage data
    chr_df <- build_multiparental_dosage(cross_results, pedigree) %>%
      mutate(Chrom = paste0("Ch_", i)) %>%
      relocate(Chrom, .after = marker) %>%
      mutate(marker = paste0(marker, "_", str_replace(Chrom, "_", "")))

    # Append chromosome-specific data to consolidated wide_df
    wide_df <- bind_rows(wide_df, chr_df)
  }

  # Create a matrix of variant alleles
  variant <- sapply(1:nrow(wide_df), function(x) sample(c("A", "T", "C", "G"), 2))

  # Identify biparental crosses
  biparental_idx <- apply(pedigree[, 1:2], 1, function(x) !any(is.na(x)))
  pedigree_biparental <- pedigree[biparental_idx, ]
  unique_crosses <- which(!duplicated(pedigree_biparental$cross))

  biparental_crosses <- pedigree_biparental$cross[unique_crosses]
  P1s <- pedigree_biparental$parent_1[unique_crosses]
  P2s <- pedigree_biparental$parent_2[unique_crosses]

  # Generate correlated biparental cross memberships
  correlated_sets <- generate_correlated_sets(
    n = nrow(wide_df),
    x = length(biparental_crosses),
    p = rep(p, length(biparental_crosses)),
    rho = rho
  )

  dat <- vector("list", length(biparental_crosses))
  names(dat) <- biparental_crosses

  for (i in seq_along(biparental_crosses)) {
    # Get offspring individual names for each biparental cross
    offspring_names <- pedigree$individual[pedigree$cross == biparental_crosses[i]]
    marker_indices <- which(correlated_sets$membership_matrix[, i])
    geno.temp <- wide_df[marker_indices, offspring_names, drop = FALSE]
    geno <- as.matrix(geno.temp)

    if(missing > 0.0){
      id <- sample(length(geno), size = length(geno) * missing)
      geno[id] <- NA
    }
    colnames(geno) <- colnames(geno.temp)
    dat[[i]] <- data.frame(
      snp_id     = wide_df$marker[marker_indices],
      P1         = wide_df[[P1s[i]]][marker_indices],
      P2         = wide_df[[P2s[i]]][marker_indices],
      chrom      = wide_df$Chrom[marker_indices],
      genome_pos = wide_df$map_position[marker_indices],
      alt        = variant[1, ][marker_indices],
      res        = variant[2, ][marker_indices],
      geno,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  return(list(
    wide_df = wide_df,
    dat = dat,
    parent_homologs = parent_homologs,
    counts = correlated_sets
  ))
}
