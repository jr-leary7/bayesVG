data {
  int<lower=1> M;  // number of spots
  int<lower=1> N;  // number of gene-spot pairs in long dataframe
  int<lower=1> G;  // number of genes
  int<lower=1> k;  // number of basis functions used to approximate GP
  array[N] int<lower=1, upper=M> spot_id;  // unique ID for each spot
  array[N] int<lower=1, upper=G> gene_id;  // unique ID for each gene
  matrix[M, k] phi;  // matrix of QR-decomposed basis functions used to approximate GP
  vector[N] y;  // vector of normalized, scaled gene expression used as response variable
}

parameters {
  real mu_beta0;
  real<lower=0> sigma_beta0;
  vector[G] z_beta0;
  real mu_amplitude;
  real<lower=0> sigma_amplitude;
  vector<lower=0>[G] amplitude;
  vector[k] mu_alpha;
  vector<lower=0>[k] sigma_alpha;
  matrix[k, G] z_alpha_t;
  real<lower=0> sigma_y;
}

transformed parameters {
  vector[G] beta0 = mu_beta0 + sigma_beta0 * z_beta0;
  vector[G] amplitude_sq = square(amplitude);
}

model {
  matrix[k, G] alpha_t;
  alpha_t = rep_matrix(mu_alpha, G) + diag_pre_multiply(sigma_alpha, z_alpha_t);
  matrix[M, G] phi_alpha = phi * alpha_t;
  vector[N] w;
  for (i in 1:N) {
    w[i] = phi_alpha[spot_id[i], gene_id[i]];
  }
  to_vector(z_alpha_t) ~ std_normal();
  z_beta0 ~ std_normal();
  mu_beta0 ~ normal(0, 2);
  sigma_beta0 ~ std_normal();
  mu_amplitude ~ normal(0, 2);
  sigma_amplitude ~ std_normal();
  amplitude ~ lognormal(mu_amplitude, sigma_amplitude);
  mu_alpha ~ normal(0, 2);
  sigma_alpha ~ std_normal();
  sigma_y ~ std_normal();
  y ~ normal(beta0[gene_id] + amplitude_sq[gene_id] .* w, sigma_y);
}
