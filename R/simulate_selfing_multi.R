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
                                   F.generations = 2){
  return(simulate_selfing_multi_rcpp(n_mrk = n.mrk,
                                     map_len = map.len,
                                     n_ind = n.ind,
                                     F_generations = F.generations))
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
#' plot_geno(dat.cur, j = 2)
plot_geno <- function(dat.cur, j) {
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

