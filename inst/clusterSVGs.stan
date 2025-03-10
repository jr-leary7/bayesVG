data {
  int<lower=1> N;  // number of SVGs
  int<lower=1> D;  // number of principal components 
  int<lower=1> K;  // number of clusters
  matrix[N, D] X;  // matrix of principal components 
}

parameters {
  simplex[K] theta;  // mixing proportions
  matrix[K, D] mu;  // cluster centers
  array[K] vector<lower=0>[D] sigma;  // cluster-specific variances (diagonal covariance)
}

model {
  theta ~ dirichlet(rep_vector(1.0, K));
  for (i in 1:K) {
    mu[i] ~ normal(0, 10);
    sigma[i] ~ std_normal();
  }
  for (i in 1:N) {
    vector[K] log_component;
    for (j in 1:K) {
      matrix[D, D] cov_j = diag_matrix(square(sigma[j]));
      log_component[j] = log(theta[j]) + multi_normal_lpdf(to_vector(X[i]) | to_vector(mu[j]), cov_j);
    }
    target += log_sum_exp(log_component);
  }
}

generated quantities {
  array[N] real log_lik;
  for (i in 1:N) {
    vector[K] log_component;
    for (j in 1:K) {
      matrix[D, D] cov_j = diag_matrix(square(sigma[j]));
      log_component[j] = log(theta[j]) + multi_normal_lpdf(to_vector(X[i]) | to_vector(mu[j]), cov_j);
    }
    log_lik[i] = log_sum_exp(log_component);
  }
}
