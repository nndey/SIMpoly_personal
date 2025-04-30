// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
double imf_haldane_cM(double r) {
  if (R_IsNA(r) || r < 0 || r > 0.5) {
    stop("Input recombination fractions must be between 0 and 0.5");
  }
  return -50.0 * log(1 - 2 * r);
}

// Helper: switch chromatids at breakpoint bp
NumericMatrix switch_chrom(const NumericMatrix& G, int bp) {
  int ncol = G.ncol();
  NumericMatrix out(2, ncol);
  for(int j = 0; j <= bp; ++j) {
    out(0, j) = G(0, j);
    out(1, j) = G(1, j);
  }
  for(int j = bp + 1; j < ncol; ++j) {
    out(0, j) = G(1, j);
    out(1, j) = G(0, j);
  }
  return out;
}

// Generate one gamete: returns haplotype and crossover counts per interval
List gamete_gen(const NumericMatrix& G, const NumericVector& map) {
  int n_mrk = G.ncol();
  double map_len = map[n_mrk - 1];
  double lambda = map_len / 100.0;
  int x = R::rpois(lambda);
  NumericVector co_pos(x);
  for(int i = 0; i < x; ++i) {
    co_pos[i] = R::runif(0, map_len);
  }
  std::sort(co_pos.begin(), co_pos.end());

  int n_int = n_mrk - 1;
  IntegerVector counts(n_int);
  for(int i = 0; i < n_int; ++i) {
    int cnt = 0;
    for(int j = 0; j < x; ++j) {
      if (co_pos[j] >= map[i] && co_pos[j] < map[i+1]) ++cnt;
    }
    counts[i] = cnt;
  }

  NumericMatrix Gout = clone(G);
  for(int i = 0; i < x; ++i) {
    double pos = co_pos[i];
    int idx = std::upper_bound(map.begin(), map.end(), pos) - map.begin() - 1;
    if (idx < 0) idx = 0;
    if (idx >= n_int) idx = n_int - 1;
    Gout = switch_chrom(Gout, idx);
  }

  int which = (R::runif(0, 1) < 0.5) ? 0 : 1;
  NumericVector gam(n_mrk);
  for(int j = 0; j < n_mrk; ++j) gam[j] = Gout(which, j);

  return List::create(_["gamete"] = gam,
                      _["co.pos"] = counts);
}

// Generate one individual by fusing two gametes
List ind_gen(const NumericMatrix& G, const NumericVector& map) {
  List g1 = gamete_gen(G, map);
  List g2 = gamete_gen(G, map);
  NumericVector gam1 = g1["gamete"];
  NumericVector gam2 = g2["gamete"];
  IntegerVector co1 = g1["co.pos"];
  IntegerVector co2 = g2["co.pos"];

  NumericMatrix ind_mat(2, G.ncol());
  for(int j = 0; j < G.ncol(); ++j) {
    ind_mat(0, j) = gam1[j];
    ind_mat(1, j) = gam2[j];
  }

  int n = co1.size();
  IntegerVector counts(n);
  for(int i = 0; i < n; ++i) counts[i] = co1[i] + co2[i];

  return List::create(_["ind"] = ind_mat,
                      _["co.pos"] = counts);
}

// [[Rcpp::export]]
List simulate_selfing_multi_rcpp(int n_mrk = 20,
                                 double map_len = 100.0,
                                 int n_ind = 2000,
                                 int F_generations = 2) {
  if (F_generations < 2) stop("F_generations must be >= 2");

  // marker positions equally spaced
  NumericVector map(n_mrk);
  for(int i = 0; i < n_mrk; ++i) map[i] = map_len * i / (n_mrk - 1);

  // initial parental matrices
  std::vector<NumericMatrix> Gcur(n_ind);
  for(int i = 0; i < n_ind; ++i) {
    NumericMatrix mat(2, n_mrk);
    for(int j = 0; j < n_mrk; ++j) {
      mat(0, j) = 1.0;
      mat(1, j) = 0.0;
    }
    Gcur[i] = mat;
  }

  int gens = F_generations - 1;
  List geno(gens), countCO(gens);
  CharacterVector gen_names(gens);
  for(int j = 0; j < gens; ++j) {
    gen_names[j] = "F_" + std::to_string(j + 2);
  }
  geno.attr("names") = gen_names;
  countCO.attr("names") = gen_names;

  CharacterVector marker_names(n_mrk);
  for(int k = 0; k < n_mrk; ++k) marker_names[k] = "M_" + std::to_string(k+1);

  for(int gen = 2; gen <= F_generations; ++gen) {
    IntegerVector total_co(n_mrk - 1);
    for(int i = 0; i < n_ind; ++i) {
      List tmp = ind_gen(Gcur[i], map);
      IntegerVector co = tmp["co.pos"];
      NumericMatrix indM = tmp["ind"];
      for(int k = 0; k < n_mrk-1; ++k) total_co[k] += co[k];
      Gcur[i] = indM;
    }
    countCO[gen - 2] = total_co;

    // build genotype sum data.frame
    NumericVector Fgen(n_ind);
    for(int i = 0; i < n_ind; ++i) Fgen[i] = gen;

    int ncol = n_mrk + 1;
    List df(ncol);

    // column 0: F_gen
    df[0] = Fgen;

    // columns 1..n_mrk: marker sums
    for(int k = 0; k < n_mrk; ++k) {
      NumericVector col(n_ind);
      for(int i = 0; i < n_ind; ++i)
        col[i] = Gcur[i](0, k) + Gcur[i](1, k);
        df[k+1] = col;
    }

    // set names
    CharacterVector nms(ncol);
    nms[0] = "F_gen";
    for(int k = 0; k < n_mrk; ++k)
      nms[k+1] = marker_names[k];
    df.attr("names") = nms;

    df.attr("class") = "data.frame";
    df.attr("row.names") = IntegerVector::create(NA_INTEGER, -n_ind);

    geno[gen - 2] = df;
  }

  return List::create(_["geno"] = geno,
                      _["map"]  = map,
                      _["count.CO"] = countCO);
}
