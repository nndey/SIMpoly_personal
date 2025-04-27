#' Simulate multiple selfing‐derived recombinant inbred populations
#'
#' @description
#'   From two inbred founders (all-1 vs all-0), simulate an F₂ via Poisson crossovers
#'   and then advance by selfing to arbitrary generations.  Return both genotype dosages
#'   and the historical counts of crossovers between adjacent markers for each generation.
#'
#' @param n_markers        Integer ≥ 2; number of equally spaced markers per chromosome.
#' @param map_length_cm    Numeric > 0; total genetic map length in centiMorgans.
#' @param n_individuals    Integer ≥ 1; individuals per generation.
#' @param generations      Integer vector ≥ 2; generations to simulate and return (e.g. c(2,6)).
#'
#' @return A list with two elements:
#'   - `genotypes`: named list of data.frames of dosage genotypes (0/1/2) for each generation.
#'   - `rec_counts`: named list of integer vectors (length = markers - 1) of total historical crossovers
#'                  observed between each adjacent marker in that generation.
#'
#' @examples
#' res <- simulate_selfing_multi(
#'   n_markers = 50,
#'   map_length_cm = 100,
#'   n_individuals = 20000,
#'   generations = c(2,6)
#' )
#' round(res$rec_counts$F2/(2*20000),3)
#' res$genotypes$F2     # genotype DF for F2
#' res$rec_counts$F2    # crossover counts between markers in F2
#' @export
simulate_selfing_multi <- function(n_markers,
                                   map_length_cm,
                                   n_individuals,
                                   generations = c(2, 6)) {
  # Validate inputs
  stopifnot(
    is.numeric(n_markers), n_markers >= 2,
    is.numeric(map_length_cm), map_length_cm > 0,
    is.numeric(n_individuals), n_individuals >= 1,
    is.numeric(generations), all(generations >= 2),
    all(generations == as.integer(generations))
  )

  # Helper: single-crossover swap
  switch_ch <- function(G, bp) {
    rbind(
      c(G[1, 1:bp], G[2, (bp+1):ncol(G)]),
      c(G[2, 1:bp], G[1, (bp+1):ncol(G)])
    )
  }

  # Wrapper to simulate a gamete and count historical crossovers
  gamete_gen_count <- function(G, len_cm, counts) {
    nxo <- rpois(1, len_cm / 100)
    if (nxo > 0) {
      pos <- sort(runif(nxo, 0, len_cm))
      markers <- seq(0, len_cm, length.out = ncol(G))
      for (p in findInterval(pos, markers)) {
        if (p >= 1 && p < ncol(G)) counts[p] <<- counts[p] + 1L
        G <- switch_ch(G, p)
      }
    }
    G[sample.int(2, 1), , drop = TRUE]
  }

  # Function to build a diploid and update counts
  ind_gen_count <- function(G, len_cm, counts) {
    g1 <- gamete_gen_count(G, len_cm, counts)
    g2 <- gamete_gen_count(G, len_cm, counts)
    rbind(g1, g2)
  }

  # Initialize storage
  gens <- sort(unique(generations))
  max_gen <- max(gens)
  genotype_list <- setNames(vector('list', length(gens)), paste0('F', gens))
  rec_counts_list <- setNames(vector('list', length(gens)), paste0('F', gens))

  # Founder chromatids for two inbreds
  G_founders <- rbind(rep(1L, n_markers), rep(0L, n_markers))

  # Simulate F2
  G_current <- vector('list', n_individuals)
  counts <- integer(n_markers - 1)
  for (i in seq_len(n_individuals)) {
    G_current[[i]] <- ind_gen_count(G_founders, map_length_cm, counts)
  }
  if (2 %in% gens) {
    dos <- t(vapply(G_current, function(m) colSums(m), numeric(n_markers)))
    df <- setNames(as.data.frame(dos), paste0('M_', seq_len(n_markers)))
    df$Generation <- 2L
    genotype_list[['F2']] <- df[, c('Generation', paste0('M_', seq_len(n_markers)))]
    rec_counts_list[['F2']] <- counts
  }

  # Advance by selfing
  for (gen in seq(3, max_gen)) {
    counts <- integer(n_markers - 1)
    for (i in seq_len(n_individuals)) {
      G_current[[i]] <- ind_gen_count(G_current[[i]], map_length_cm, counts)
    }
    if (gen %in% gens) {
      dos <- t(vapply(G_current, function(m) colSums(m), numeric(n_markers)))
      df <- setNames(as.data.frame(dos), paste0('M_', seq_len(n_markers)))
      df$Generation <- gen
      genotype_list[[paste0('F', gen)]] <- df[, c('Generation', paste0('M_', seq_len(n_markers)))]
      rec_counts_list[[paste0('F', gen)]] <- counts
    }
  }

  list(
    genotypes  = genotype_list,
    rec_counts = rec_counts_list
  )
}
