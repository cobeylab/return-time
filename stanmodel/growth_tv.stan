data {
  int N;
  int cases[N];
}

parameters {
  real r[N-1];
  real<lower=0> r_sd;
  real<lower=0> C0;
  real<lower=0> phi;
}

transformed parameters {
  vector<lower=0>[N] C;
  
  C[1] = C0;
  
  for (i in 2:(N)) {
    C[i] = C[i-1] * exp(r[i-1]);
  }
}

model {
  r[1] ~ normal(0, 1);
  
  for (i in 2:(N-1)) {
    r[i] ~ normal(r[i-1], r_sd);
  }
  C0 ~ cauchy(0, 1);
  
  phi ~ normal(0, 10);
  
  for (i in 1:N) {
    cases[i] ~ neg_binomial_2(C[i], phi);
  }
}
