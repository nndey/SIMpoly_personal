#' Plot Parent Phase with ggplot2 (Flipped and Padded)
#'
#' This function takes a named list of binary matrices (each representing a parent)
#' and produces a ggplot2 tile plot where:
#'   - The matrices are padded to have four columns (extra columns are filled with NA)
#'   - The axes are flipped relative to the original matrix (markers are on the x-axis,
#'     haplotypes on the y-axis)
#'   - Each parent's data is shown in its own row (using facet_grid)
#'
#' @param parent.phase A named list of binary matrices. Each matrix should have markers
#'   in rows and haplotypes in columns.
#' @param colors A character vector of two colors for the binary values (default: c("white", "black")).
#'        The NA fill (for padded cells) is drawn with the na.value color.
#' @param na.value Color for missing cells (default: "grey90").
#'
#' @return A ggplot object is printed.
#' @import ggplot2
#' @export
plot_parent_phase_gg <- function(parent.phase, colors = c("white", "black"), na.value = "grey90") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required but not installed.")
  }

  # Helper function to pad a matrix to have exactly 4 columns.
  pad_matrix <- function(mat, ncols = 4) {
    if (ncol(mat) < ncols) {
      extra <- matrix(NA, nrow = nrow(mat), ncol = ncols - ncol(mat))
      # Create extra column names based on existing names or defaults.
      existing <- if (!is.null(colnames(mat))) colnames(mat) else paste0("H", 1:ncol(mat))
      extra_colnames <- paste0("H", (length(existing)+1):ncols)
      colnames(extra) <- extra_colnames
      mat <- cbind(mat, extra)
    }
    return(mat)
  }

  # Convert each parent's matrix to a long-format data frame.
  df_list <- lapply(names(parent.phase), function(parentName) {
    mat <- as.matrix(parent.phase[[parentName]])
    # Ensure the matrix is numeric.
    if (!is.numeric(mat)) {
      mat <- apply(mat, c(1,2), as.numeric)
    }
    # Pad the matrix to have 4 columns.
    mat <- pad_matrix(mat, ncols = 4)

    # Get marker names (rows) and haplotype names (columns).
    markers <- if (!is.null(rownames(mat))) rownames(mat) else paste0("M", 1:nrow(mat))
    haplotypes <- if (!is.null(colnames(mat))) colnames(mat) else paste0("H", 1:ncol(mat))

    # Convert the matrix to a data frame in long format.
    df <- as.data.frame(as.table(mat))
    names(df) <- c("Marker", "Haplotype", "Value")

    # Now, since we want to flip x and y relative to the original:
    #   - x-axis: Marker (in the original matrix, these were rows)
    #   - y-axis: Haplotype (originally columns)
    # We set the factor levels to preserve the natural order.
    df$Marker <- factor(df$Marker, levels = markers)        # left-to-right: M1, M2, ...
    df$Haplotype <- factor(df$Haplotype, levels = haplotypes)  # bottom-to-top: H1, H2, H3, H4
    df$Parent <- parentName
    return(df)
  })

  df_all <- do.call(rbind, df_list)

  # Create the ggplot.
  p <- ggplot2::ggplot(df_all, ggplot2::aes(x = Marker, y = Haplotype, fill = factor(Value))) +
    ggplot2::geom_tile(color = "grey") +
    ggplot2::facet_grid(Parent ~ .) +
    ggplot2::scale_fill_manual(values = colors, name = "Value", na.value = na.value) +
    ggplot2::labs(x = "Markers", y = "Haplotypes") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "lightgrey"))

  print(p)
}
