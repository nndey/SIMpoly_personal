#' Compute possible gamete dosages from a polyploid genotype
#'
#' Given a string representation of a polyploid genotype (e.g. “AAaa” for tetraploid),
#' returns all possible counts of the reference allele (uppercase) in the gametes.
#'
#' @param genotype Character. Even-length string of uppercase (reference) and lowercase (alternate) alleles.
#' @return Integer vector of gamete allele counts.
#' @examples
#' get_dose_segregation("AAaa")        # tetraploid simplex
#' get_dose_segregation("AAAaaa")      # hexaploid duplex
#' @export
get_dose_segregation <- function(genotype) {
  alleles <- strsplit(genotype, "")[[1]]
  doses   <- as.integer(grepl("[A-Z]", alleles))
  ploidy  <- length(doses)
  if (ploidy %% 2 != 0) {
    stop("Genotype length must be even (auto-polyploid).")
  }
  half     <- ploidy / 2
  gametes  <- combn(doses, half, sum)
  return(as.integer(gametes))
}

#' Plot a heatmap of gamete dosages plus a stripe-bar of dosage proportions
#'
#' Builds a tile plot of all parental-gamete combinations (cells labeled by dosage),
#' and below it a bar showing the overall frequency of each dosage.
#'
#' @param gamete_dosage Numeric matrix; entry [i,j] is parent1_gamete_i + parent2_gamete_j.
#' @param complete_dominance Logical; if TRUE, uses two colors (zero vs non-zero), otherwise a gradient.
#' @param col Character vector of length 2; c(low_color, high_color) for the fill scale.
#' @return A combined ggplot object (heatmap over stripe bar).
#' @import ggplot2 reshape2 dplyr patchwork
#' @examples
#' # Suppose you have a hexaploid duplex × simplex cross:
#' p1 <- get_dose_segregation("AAaaaa")
#' p2 <- get_dose_segregation("AAAAAa")
#' mat <- outer(p1, p2, "+")
#' plot_genotype_heatmap_with_stripe(mat, complete_dominance=FALSE, col=c("purple","red"))
#' @export
plot_genotype_heatmap_with_stripe <- function(
    gamete_dosage,
    complete_dominance = FALSE,
    col = c("blue", "red")
) {
  # 1) Melt to long form
  df <- reshape2::melt(
    gamete_dosage,
    varnames = c("G1", "G2"),
    value.name = "Dosage"
  )
  df$G1 <- factor(df$G1)
  df$G2 <- factor(df$G2)

  # 2) Prepare fill mapping and scales
  if (complete_dominance) {
    df$fillVal      <- ifelse(df$Dosage == 0, "zero", "nonzero")
    heat_scale      <- scale_fill_manual(
      values = c(nonzero = col[1], zero = col[2]),
      name   = NULL,
      labels = c("> 0", "0")
    )
    stripe_aes      <- aes(fill = fillVal)
    stripe_scale    <- heat_scale
  } else {
    max_val         <- max(df$Dosage)
    df$fillVal      <- df$Dosage
    heat_scale      <- scale_fill_gradient(
      low    = col[2],
      high   = col[1],
      limits = c(0, max_val),
      name   = "Dosage"
    )
    stripe_aes      <- aes(fill = Dosage)
    stripe_scale    <- heat_scale
  }

  # 3) Heatmap
  p_heat <- ggplot(df, aes(x = G2, y = G1, fill = fillVal)) +
    geom_tile() +
    heat_scale +
    geom_text(aes(label = Dosage), color = "white", size = 3) +
    theme_minimal() +
    labs(x = "Gamete dose (Parent 2)", y = "Gamete dose (Parent 1)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # 4) Stripe bar of dosage proportions
  stripe_df <- df %>%
    dplyr::count(Dosage, fillVal) %>%
    dplyr::mutate(prop = n / sum(n))

  p_stripe <- ggplot(stripe_df, aes(x = factor(Dosage), y = prop)) +
    do.call(geom_col, list(mapping = stripe_aes)) +
    stripe_scale +
    theme_minimal() +
    labs(x = "Dosage", y = "Proportion") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # 5) Stack
  combined <- p_heat / p_stripe + patchwork::plot_layout(heights = c(4, 1))
  return(combined)
}

#' Simulate segregation, return counts, proportions, and optional plot
#'
#' Crosses two polyploid genotype strings, computes full progeny allele-dosage
#' distribution, shows a heatmap+stripe plot if requested.
#'
#' @param geno_p1 Character; parent 1 genotype string (even length).
#' @param geno_p2 Character; parent 2 genotype string (even length).
#' @param plot Logical; if TRUE, displays the heatmap + stripe plot.
#' @param complete_dominance Logical; passed through to the plotting function.
#' @param col Character vector length 2; colors for the fill scale.
#' @return A list with:
#'   \item{genotype_counts}{Named integer vector of counts per dosage.}
#'   \item{expected_segregation}{Numeric vector of proportions per dosage.}
#'   \item{plot}{Invisible ggplot object if \code{plot = TRUE}.}
#' @examples
#' segregation("AAAAAa", "AAAAAA", plot = TRUE)
#' @export
segregation <- function(
    geno_p1,
    geno_p2,
    plot              = TRUE,
    complete_dominance = FALSE,
    col               = c("blue", "red")
) {
  P1 <- get_dose_segregation(geno_p1)
  P2 <- get_dose_segregation(geno_p2)
  gamete_dosage <- kronecker(P1, t(P2), "+")

  # Plot if requested
  plot_obj <- plot_genotype_heatmap_with_stripe(
    gamete_dosage,
    complete_dominance = complete_dominance,
    col                = col
  )
  if (plot) print(plot_obj)

  # Counts & proportions
  counts <- table(as.vector(gamete_dosage))
  props  <- counts / sum(counts)

  return(list(
    genotype_counts      = counts,
    expected_segregation = props,
    plot                 = invisible(plot_obj)
  ))
}

# Example:
# segregation("AAAAAa", "AAAAAA")
