## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)
library(knitr)

## ----install-packages---------------------------------------------------------
#  library(SIMpoly)
#  library(mappoly2)  # for consensus mapping

## ----define-pedigree, eval=TRUE-----------------------------------------------
# Define ploidy for each parent.
ploidy.vec <- c(4, 2, 4, 2, 4, 4)  # for six parents
names(ploidy.vec) <- c("P1", "P2", "P3", "P4", "P5", "P6")

# Define parental crosses.
parents.mat <- matrix(c("P1", "P2",
                        "P1", "P3",
                        "P2", "P2",
                        "P3", "P4",
                        "P4", "P5",
                        "P5", "P6",
                        "P5", "P2"),
                      ncol = 2, byrow = TRUE)

# Define number of offspring per cross.
n_offspring <- c(200, 100, 200, 100, 300, 200, 100)

# Define allele set for each founder.
alleles <- list(P1 = 0:1, P2 = 0:1, P3 = 0:1, P4 = 0:1, P5 = 0:1, P6 = 0:1)

# Generate pedigree.
pedigree <- simulate_pedigree(ploidy.vec, parents.mat, n_offspring)
head(pedigree, n = 10)

## ----simulate-data, eval=TRUE-------------------------------------------------
n.chr <- 3  # Number of chromosomes
map.len <- c(100, 120, 90)  # Chromosome lengths
n_mrk <- c(2000, 1500, 2000)  # Number of markers per chromosome

# Simulate genetic data
result <- simulate_multiparental_data(n.chr, map.len, pedigree, ploidy.vec, n_mrk, alleles, p = .3 , rho = .7)

# Extract results
wide_df <- result$wide_df
all_dat <- result$dat
parent_homologs <- result$parent_homologs[[1]]
plot_parent_phase_gg(parent_homologs)
plot_correlated_sets(result$counts)

## ----build-maps---------------------------------------------------------------
#  MAPs <- vector("list", nrow(parents.mat))
#  
#  for(i in 1:nrow(parents.mat)){
#    p1 <- parents.mat[i,1]
#    p2 <- parents.mat[i,2]
#    dat <- mappoly2:::table_to_mappoly(all_dat[[i]], ploidy.vec[p1], ploidy.vec[p2], p1, p2)
#  
#    # Filtering steps
#    dat <- filter_markers_by_chisq_pval(dat, plot = FALSE)
#    dat <- filter_markers_by_missing_rate(dat, mrk.thresh = .1, plot = FALSE)
#    dat <- filter_individuals_by_missing_rate(dat, ind.thresh = .1, plot = FALSE)
#    dat <- pairwise_rf(dat, mrk.scope = "all", ncpus = 8)
#    plot(dat)
#    dat <- rf_filter(dat, probs = c(0, 1), diagnostic.plot = FALSE)
#  
#    # Group and order markers
#    g <- group(dat, expected.groups = n.chr, comp.mat = TRUE, inter = FALSE)
#    plot(g)
#    s <- make_sequence(g, ch = list(1,2,3))
#    s <- order_sequence(s, type = "genome")
#    s <- pairwise_phasing(s, type = "genome", parent = "p1p2")
#    s <- mapping(s, type = "genome", parent = "p1p2", ncpus = n.chr)
#    s <- calc_haplotypes(s, type = "genome", ncpus = n.chr)
#  
#    MAPs[[i]] <- s
#  }
#  plot_map(MAPs[[1]], lg = 1, type = "genome")

## ----consensus-map------------------------------------------------------------
#  # Visualize multiple maps
#  plot_multi_map(MAPs)
#  
#  # Integrate maps
#  w <- prepare_to_integrate(MAPs, type = "genome")
#  plot(w)
#  
#  # Estimate consensus map
#  x <- estimate_consensus_map(w, ncpus = 8)
#  plot(x)
#  plot(x, only.consensus = TRUE, col = mappoly::mp_pallet3(3))
#  
#  # Compute haplotypes
#  system.time(x <- calc_consensus_haplo(x, ncpus = 8))
#  
#  # Visualize consensus haplotypes
#  plot_consensus_haplo(x, lg = 1, ind = "P1xP2_1")
#  plot_consensus_haplo(x, lg = 1, ind = "P5xP6_1")

