data {
  int<lower=1> M;  // number of spots
  int<lower=1> N;  // number of gene-spot pairs
  int<lower=1> G;  // number of genes
  int<lower=1> k;  // number of basis functions
  array[N] int<lower=1,upper=M> spot_id;  // unique ID for each spot
  array[N] int<lower=1,upper=G> gene_id;  // unique ID for each gene
  matrix[M, k] phi;  // QR–decomposed basis functions
  vector[G] gene_depths;  // vector of logged gene-level sequencing depths to adjust for in the model
  array[N] int<lower=0> y;  // raw UMI counts
}

parameters {
  real beta0;  // global intercept
  real beta1;  // coefficient for gene library size
  real<lower=0> phi_nb;  // negative-binomial overdispersion parameter 
  vector<lower=0>[G] amplitude;  // vector of gene-specific amplitudes of the approximate GP
  real mu_amplitude;  // mean for the amplitude
  real<lower=0> sigma_amplitude;  // SD for the amplitude
  vector[k] mu_alpha;  // vector of means for the basis function coefficients
  vector<lower=0>[k] sigma_alpha;  // vector of SDs for the basis function coefficients
  matrix[k, G] z_alpha_t;  // standard normal RV for basis function coefficients
}

transformed parameters {
  vector[G] amplitude_sq = square(amplitude);
}

model {
  beta0 ~ normal(0, 2);
  beta1 ~ normal(0, 2);
  matrix[k, G] alpha_t;
  alpha_t = rep_matrix(mu_alpha, G) + diag_pre_multiply(sigma_alpha, z_alpha_t);
  matrix[M, G] phi_alpha = phi * alpha_t;
  to_vector(z_alpha_t) ~ std_normal();
  mu_alpha ~ normal(0, 2);
  sigma_alpha ~ std_normal();
  mu_amplitude ~ normal(0, 2);
  sigma_amplitude ~ std_normal();
  amplitude ~ lognormal(mu_amplitude, sigma_amplitude);
  phi_nb ~ cauchy(0, 2) T[0, ];
  vector[N] w;
  for (i in 1:N) {
    w[i] = phi_alpha[spot_id[i], gene_id[i]];
  }
  y ~ neg_binomial_2_log(beta0 + beta1 * gene_depths[gene_id] + amplitude_sq[gene_id] .* w, phi_nb);
}
