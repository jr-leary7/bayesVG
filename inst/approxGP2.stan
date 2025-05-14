data {
  int<lower=1> M;  // number of spots
  int<lower=1> N;  // number of gene-spot pairs in long dataframe
  int<lower=1> G;  // number of genes
  int<lower=1> k;  // number of basis functions used to approximate GP
  array[N] int<lower=1, upper=M> spot_id;  // unique ID for each spot
  array[N] int<lower=1, upper=G> gene_id;  // unique ID for each gene
  matrix[M, k] phi;  // matrix of basis functions used to approximate GP
  vector[G] gene_depths;  // vector of logged gene-level sequencing depths to adjust for in the model
  vector[N] y;  // vector of normalized, scaled gene expression used as response variable
}

parameters {
  real beta0;  // global intercept
  real beta1;  // coefficient for gene library size
  matrix[k, G] alpha_t;  // transposed matrix of gene-specific coefficients for each basis function
  real<lower=0> sigma_y;  // observation noise of response variable
  vector<lower=0>[G] amplitude;  // vector of gene-specific amplitudes of the approximate GP
  real mu_amplitude;  // mean for the amplitude
  real<lower=0> sigma_amplitude;  // SD for the amplitude
  vector[k] mu_alpha;  // vector of means for the basis function coefficients
  vector<lower=0>[k] sigma_alpha;  // vector of SDs for the basis function coefficients
}

model {
  matrix[M, G] phi_alpha;
  phi_alpha = phi * alpha_t;
  vector[N] w;
  for (i in 1:N) {
    w[i] = phi_alpha[spot_id[i], gene_id[i]];
  }
  vector[G] amplitude_sq = square(amplitude);
  beta0 ~ normal(0, 2);
  beta1 ~ normal(0, 2);
  mu_alpha ~ normal(0, 2);
  sigma_alpha ~ std_normal();
  mu_amplitude ~ normal(0, 2);
  sigma_amplitude ~ std_normal();
  sigma_y ~ normal(0, 2);
  for (i in 1:k) {
    alpha_t[i] ~ normal(mu_alpha[i], sigma_alpha[i]);
  }
  amplitude ~ lognormal(mu_amplitude, sigma_amplitude);
  y ~ normal(beta0 + beta1 * gene_depths[gene_id] + amplitude_sq[gene_id] .* w, sigma_y);
}

// generated quantities {
//   matrix[M, G] phi_alpha;
//   phi_alpha = phi * alpha_t;
//   vector[G] amplitude_sq = square(amplitude);
//   array[N] real log_lik;
//   for (i in 1:N) {
//     real mu_i;
//     mu_i = beta0 + beta1 * gene_depths[gene_id[i]] + amplitude_sq[gene_id[i]] * phi_alpha[spot_id[i], gene_id[i]];
//     log_lik[i] = normal_lpdf(y[i] | mu_i, sigma_y);
//   }
// }
