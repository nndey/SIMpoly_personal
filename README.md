---

&#x20;&#x20;

# **SIMpoly - Simulation of Genetic Marker Data in Pedigreed Populations**

### *Marcelo Mollinari*

## 🚀 **Introduction**

This vignette demonstrates how to use **SIMpoly** to simulate multiple biparental populations and how to leverage **mappoly2** to construct a consensus genetic map. The workflow includes:

1. **Installing from GitHub**
2. **Defining simulation parameters and pedigree.**
3. **Simulating multiple chromosomes** using **SIMpoly**.
4. **Constructing linkage maps** using **mappoly2**.
5. **Integrating maps** to generate a consensus map.

> ⚠️ **Note:** The consensus map generation requires **mappoly2**. If you only have **SIMpoly**, the first simulation steps will work, but the final consensus map will not be generated.

## 📦 **Installation from GitHub**

You can install the development version from GitHub. Within R, install `devtools`:

```r
install.packages("devtools")
```

If you are using Windows, please install the latest recommended version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

To install SIMpoly from GitHub, use:

```r
devtools::install_github("mmollina/SIMpoly", dependencies=TRUE)
```

---

## 🛠️ **Step 1: Installing and Loading Packages**

```r
library(SIMpoly)
library(mappoly2)  # for consensus mapping
```

---

## 🔬 **Step 2: Define Simulation Parameters and Pedigree**

```r
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
```

---

## 📊 **Step 3: Simulate Multiparental Data**

```r
n.chr <- 3  # Number of chromosomes
map.len <- c(100, 120, 90)  # Chromosome lengths
n_mrk <- c(2000, 1500, 2000)  # Number of markers per chromosome

# Simulate genetic data
result <- simulate_multiparental_data(n.chr, map.len, pedigree, ploidy.vec, n_mrk, alleles, missing = 0.1, p = .3, rho = .7)

# Extract results
wide_df <- result$wide_df
all_dat <- result$dat
parent_homologs <- result$parent_homologs[[1]]
plot_parent_phase_gg(parent_homologs)
```

---

## 🔗 **Step 4: Construct Linkage Maps with mappoly2**

```r
MAPs <- vector("list", nrow(parents.mat))

for(i in 1:nrow(parents.mat)){
  p1 <- parents.mat[i,1]
  p2 <- parents.mat[i,2]
  dat <- mappoly2:::table_to_mappoly(all_dat[[i]], ploidy.vec[p1], ploidy.vec[p2], p1, p2)
  
  # Filtering steps
  dat <- filter_markers_by_chisq_pval(dat, plot = FALSE)
  dat <- filter_markers_by_missing_rate(dat, mrk.thresh = .1, plot = FALSE)
  dat <- filter_individuals_by_missing_rate(dat, ind.thresh = .1, plot = FALSE)
  dat <- pairwise_rf(dat, mrk.scope = "all", ncpus = 8)
  dat <- rf_filter(dat, probs = c(0, 1), diagnostic.plot = FALSE)
  
  # Group and order markers
  g <- group(dat, expected.groups = n.chr, comp.mat = TRUE, inter = FALSE)
  s <- make_sequence(g, ch = list(1,2,3))
  s <- order_sequence(s, type = "genome")
  s <- mapping(s, type = "genome", parent = "p1p2", ncpus = n.chr)
  
  MAPs[[i]] <- s
}
plot_map(MAPs[[1]], lg = 1, type = "genome")
```

---

## 🏆 **Acknowledgment**

This package has been developed as part of the project [AFRI-Grant: A Genetics-Based Data Analysis System for Breeders in Polyploid Breeding Programs](https://portal.nifa.usda.gov/web/crisprojectpages/1027948-a-genetics-based-data-analysis-system-for-breeders-in-polyploid-breeding-programs.html) funded by USDA NIFA.

&#x20;

NC State University promotes equal opportunity and prohibits discrimination and harassment based on various protected characteristics.

