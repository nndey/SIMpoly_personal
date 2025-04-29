#' Simulate genotype formation and crossover counts over multiple generations
#'
#' @description
#' Simulates a population of diploid individuals undergoing repeated selfing.
#' At each meiosis, crossovers are placed according to a Poisson process along
#' a linear genetic map. Returns per-generation genotype sums and crossover counts.
#'
#' @param n.mrk Integer. Number of markers (loci) per chromosome. Default: 20.
#' @param map.len Numeric. Total map length in centimorgans. Default: 100.
#' @param n.ind Integer. Number of individuals in the population. Default: 2000.
#' @param F.generations Integer ≥ 2. Number of generations to simulate (F1 → F2 … → F\{F.generations\}). Default: 2.
#' @param plot Logical. If TRUE, calls \code{print_geno()} to plot each generation’s heatmap. Default: TRUE.
#' @param sleep Numeric. Seconds to pause between plots. Default: 0.5.
#'
#' @return A list with components:
#' \describe{
#'   \item{geno}{List of data.frames, one per generation from F2 to F\{F.generations\}. Each has rows = individuals and columns = “F_gen” plus marker sums.}
#'   \item{map}{Numeric vector of marker positions (length = \code{n.mrk}), equally spaced from 0 to \code{map.len}.}
#'   \item{count.CO}{List of numeric vectors (length = \code{n.mrk}-1), giving counts of crossovers between adjacent markers for each generation.}
#' }
#'
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#'
#' @examples
#' \dontrun{
#'   res <- simulate_selfing_multi(n.mrk = 10, map.len = 50, n.ind = 100, F.generations = 4, plot = FALSE)
#'   str(res)
#' }
#' @export simulate_selfing_multi
simulate_selfing_multi <- function(n.mrk = 20,
                          map.len = 100,
                          n.ind = 2000,
                          F.generations = 2,
                          plot = TRUE,
                          sleep = 0.5) {
  # helper: single-crossover swap
  switch_ch <- function(G, bp) {
    rbind(
      c(G[1, 1:bp], G[2, (bp+1):ncol(G)]),
      c(G[2, 1:bp], G[1, (bp+1):ncol(G)])
    )
  }

  # simulate one gamete: returns one haplotype + CO counts per interval
  gamete_gen <- function(G, map.len) {
    # expected number of crossovers (in Morgans)
    lambda <- map.len / 100
    x <- rpois(1, lambda)
    co.pos <- sort(runif(n = x, min = 0, max = map.len))
    map <- seq(0, map.len, length.out = n.mrk)

    # count crossovers in each interval
    lbl <- paste0("M", seq_len(n.mrk-1), "-M", seq(2, n.mrk))
    counts <- sapply(seq_len(n.mrk-1), function(i) {
      sum(co.pos >= map[i] & co.pos < map[i+1])
    })
    names(counts) <- lbl

    # apply each crossover in order
    G.out <- G
    for (pos in co.pos) {
      co.idx <- which(diff((seq(0, map.len, length.out = ncol(G)) - pos) < 0) == -1)
      G.out <- switch_ch(G.out, co.idx)
    }

    # randomly pick one of the two chromatids
    gam <- G.out[sample(1:2, 1), ]
    list(gamete = gam, co.pos = counts)
  }

  # generate one individual by fusing two gametes
  ind_gen <- function(G, map.len) {
    g1 <- gamete_gen(G, map.len)
    g2 <- gamete_gen(G, map.len)
    list(
      ind    = rbind(g1$gamete, g2$gamete),
      co.pos = g1$co.pos + g2$co.pos
    )
  }

  # initial “parents” for F1: one homologue all 1s, the other all 0s
  G.cur <- vector("list", n.ind)
  for (i in seq_len(n.ind)) {
    G.cur[[i]] <- rbind(rep(1, n.mrk), rep(0, n.mrk))
  }

  # prepare outputs
  marker_names <- paste0("M_", seq_len(n.mrk))
  geno     <- vector("list", F.generations-1)
  count.CO <- vector("list", F.generations-1)
  names(geno)     <- names(count.CO) <- paste0("F_", seq(2, F.generations))

  # simulate from F2 to F_F.generations
  for (j in seq(2, F.generations)) {
    # accumulate crossovers for this generation
    total_co <- numeric(n.mrk-1)

    for (i in seq_len(n.ind)) {
      tmp <- ind_gen(G.cur[[i]], map.len)
      total_co <- total_co + tmp$co.pos
      G.cur[[i]] <- tmp$ind
    }

    count.CO[[j-1]] <- total_co

    # build genotype sum matrix
    dat.cur <- t(sapply(G.cur, function(x) apply(x, 2, sum)))
    dat.cur <- data.frame(F_gen = j, dat.cur, check.names = FALSE)
    rownames(dat.cur) <- paste0("Ind_", seq_len(nrow(dat.cur)))
    colnames(dat.cur) <- c("F_gen", marker_names)

    # optionally plot
    if (plot) {
      print_geno(dat.cur, j)
      Sys.sleep(sleep)
    }

    geno[[j-1]] <- dat.cur
  }

  list(
    geno     = geno,
    map      = round(seq(0, map.len, length.out = n.mrk), 1),
    count.CO = count.CO
  )
}


#' Plot genotype heatmaps by generation
#'
#' @description
#' Given a data frame of genotype sums per individual (rows) and marker (columns),
#' this function drops the first column (assumed to be generation ID), reshapes the
#' data into long form, recodes numeric genotype calls (0,1,2) into factors (“aa”,
#' “Aa”, “AA”), and draws a tile‐based heatmap with a discrete legend.
#'
#' @param dat.cur A data.frame or matrix. First column is a generation label (will be dropped);
#'   remaining columns are numeric sums (0–2) for each marker, with rownames as individual IDs.
#' @param j Integer or character. Generation identifier used in the plot title (e.g. 2 for “F_2”).
#'
#' @return Invisibly returns the ggplot object.
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column
#' @export
#'
#' @examples
#' # assume `dat.cur` has columns: F_gen, M_1, M_2, …, M_20 and rownames Ind_1, Ind_2, …
#' print_geno(dat.cur, j = 2)
print_geno <- function(dat.cur, j) {
  # 1) Drop the first column (generation ID)
  m <- dat.cur[, -1, drop = FALSE]

  # 2) Preserve original ordering
  indiv_order  <- rownames(m)
  marker_order <- colnames(m)

  # 3) Reshape and recode
  df_long <- as.data.frame(m) %>%
    tibble::rownames_to_column("Individual") %>%
    tidyr::pivot_longer(
      cols      = -Individual,
      names_to  = "Marker",
      values_to = "Genotype"
    ) %>%
    dplyr::mutate(
      Individual = factor(Individual, levels = indiv_order),
      Marker     = factor(Marker,    levels = marker_order),
      Genotype   = factor(
        Genotype,
        levels = c("0", "1", "2"),
        labels = c("aa", "Aa", "AA")
      )
    )

  # 4) Build and print the tile plot
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Marker, y = Individual, fill = Genotype)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(
      values = c(aa = "red2",
                 Aa = "forestgreen",
                 AA = "cornflowerblue"),
      name   = "Genotype"
    ) +
    ggplot2::labs(
      title = paste0("F_", j),
      x     = NULL,
      y     = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid  = ggplot2::element_blank()
    )

  print(p)
  invisible(p)
}
#' Inverse Haldane Mapping Function (returns cM)
#'
#' @description
#' Converts recombination frequency (\eqn{r}) into genetic distance in centiMorgans (cM)
#' according to the Haldane mapping function, assuming no crossover interference.
#'
#' @param r Numeric vector. Recombination frequencies (\eqn{0 \leq r \leq 0.5}).
#'
#' @return Numeric vector of genetic distances in centiMorgans (cM).
#'
#' @details
#' The Haldane mapping function assumes no crossover interference and relates
#' recombination fraction (\eqn{r}) to map distance (\eqn{d}) by:
#' \deqn{d = -\frac{1}{2} \log(1 - 2r)}{d = -0.5 * log(1 - 2r)}
#' in Morgans. Multiplying by 100 gives centiMorgans.
#'
#' @examples
#' imf_haldane_cM(0.1)   # 11.06 cM
#' imf_haldane_cM(c(0.01, 0.1, 0.2))  # multiple recombination frequencies
#'
#' @export
imf_haldane_cM <- function(r) {
  if (any(r < 0 | r > 0.5, na.rm = TRUE))
    stop("Input recombination fractions must be between 0 and 0.5")
  -50 * log(1 - 2 * r)
}
