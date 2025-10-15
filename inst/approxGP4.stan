data {
  int<lower=1> M;  // number of spots
  int<lower=1> N;  // number of gene-spot pairs
  int<lower=1> G;  // number of genes
  int<lower=1> k;  // number of basis functions
  array[N] int<lower=1,upper=M> spot_id;  // unique ID for each spot
  array[N] int<lower=1,upper=G> gene_id;  // unique ID for each gene
  matrix[M, k] phi;  // QR–decomposed basis functions
  array[N] int<lower=0> y;  // raw UMI counts
}

parameters {
  real mu_beta0;  // mean for the gene-specific intercept 
  real<lower=0> sigma_beta0;  // SD for the gene-specific intercept 
  vector[G] z_beta0;  // vector of standard normal RVs for the global intercepts 
  real<lower=0> phi_nb;  // negative-binomial overdispersion parameter 
  vector<lower=0>[G] amplitude;  // vector of gene-specific amplitudes of the approximate GP
  real mu_amplitude;  // mean for the amplitude
  real<lower=0> sigma_amplitude;  // SD for the amplitude
  vector[k] mu_alpha;  // vector of means for the basis function coefficients
  vector<lower=0>[k] sigma_alpha;  // vector of SDs for the basis function coefficients
  matrix[k, G] z_alpha_t;  // standard normal RV for basis function coefficients
}

transformed parameters {
  vector[G] beta0 = mu_beta0 + sigma_beta0 * z_beta0;
  vector[G] amplitude_sq = square(amplitude);
}

model {
  matrix[k, G] alpha_t;
  alpha_t = rep_matrix(mu_alpha, G) + diag_pre_multiply(sigma_alpha, z_alpha_t);
  matrix[M, G] phi_alpha = phi * alpha_t;
  to_vector(z_alpha_t) ~ std_normal();
  mu_beta0 ~ normal(0, 2);
  sigma_beta0 ~ std_normal();
  z_beta0 ~ std_normal();
  mu_amplitude ~ normal(0, 2);
  sigma_amplitude ~ std_normal();
  amplitude ~ lognormal(mu_amplitude, sigma_amplitude);
  phi_nb ~ cauchy(0, 2) T[0, ];
  mu_alpha ~ normal(0, 2);
  sigma_alpha ~ std_normal();
  vector[N] w;
  for (i in 1:N) {
    w[i] = phi_alpha[spot_id[i], gene_id[i]];
  }
  y ~ neg_binomial_2_log(beta0[gene_id] + amplitude_sq[gene_id] .* w, phi_nb);
}
