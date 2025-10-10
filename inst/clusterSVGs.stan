data {
  int<lower=1> N;  // number of SVGs
  int<lower=1> D;  // number of principal components
  int<lower=1> K;  // number of clusters
  matrix[N, D] X;  // matrix of principal components
}

parameters {
  simplex[K] theta;  // mixing proportions
  matrix[K, D] mu;  // cluster centers
  array[K] vector<lower=0>[D] sigma;  // cluster-specific standard deviations (diagonal covariance)
}

model {
  to_vector(mu) ~ normal(0, 7);
  for (i in 1:K) {
    sigma[i] ~ normal(0, 3);
  }
  theta ~ dirichlet(rep_vector(2, K));
  vector[K] log_theta_const;
  for (i in 1:K) {
    log_theta_const[i] = log(theta[i]) - 0.5 * (D * log(2 * pi()) + 2 * sum(log(sigma[i])));
  }
  for (i in 1:N) {
    vector[K] log_component;
    for (j in 1:K) {
      vector[D] scaled = (X[i]' - mu[j]') ./ sigma[j];
      log_component[j] = log_theta_const[j] - 0.5 * dot_self(scaled);
    }
    target += log_sum_exp(log_component);
  }
}

generated quantities {
  array[N] real log_lik;
  for (i in 1:N) {
    vector[K] log_component;
    for (j in 1:K) {
      vector[D] scaled = (X[i]' - mu[j]') ./ sigma[j];
      log_component[j] = log(theta[j]) - 0.5 * (D * log(2 * pi()) + 2 * sum(log(sigma[j])) + dot_self(scaled));
    }
    log_lik[i] = log_sum_exp(log_component);
  }
}
