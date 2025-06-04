data {
  int N;
  int Npred;
  int Nnpi;
  int npistart;
  int npiend;
  int npiwhich[npiend-npistart+1];
  int week[N+Npred];
  int cases[N];
  real mu;
  real pop;
  real gamma;
}

parameters {
  real<lower=0, upper=1> S0;
  real<lower=0, upper=1> I0;
  real<lower=0, upper=1> rho;
  real<lower=0> beta[52];
  real<lower=0> sigma;
  real<lower=0> tau;
  real<lower=0> phi;
  real<lower=0> omega;
  real<lower=0> npi[Nnpi];
}

transformed parameters {
  vector<lower=0>[N+Npred] S;
  vector<lower=0>[N+Npred] I;
  vector<lower=0>[N+Npred] R;
  vector<lower=0>[N+Npred] C;
  vector<lower=0>[N+Npred] replenish;
  vector<lower=0>[N+Npred] npieff;
  real delta = 1/tau;
  
  S[1] = S0 * pop;
  I[1] = I0 * pop;
  R[1] = (1-S0-I0) * pop;
  C[1] = 0;
  
  replenish[1] = 0;
  
  for (i in 1:(npistart-1)) {
    npieff[i] = 1;
  }
  
  for (i in npistart:npiend) {
    npieff[i] = npi[npiwhich[i-npistart+1]];
  }
  
   for (i in (npiend+1):(N+Npred)) {
     npieff[i] = 1;
   }
  
  for (i in 2:(N+Npred)) {
    real foi = beta[week[i]] * (I[i-1]+omega)/pop * npieff[i];
    real Sout = (1-exp(-(foi+mu))) * S[i-1];
    real StoI = foi/(foi+mu) * Sout;
    real Iout = (1-exp(-(gamma+mu))) * I[i-1];
    real ItoR = gamma/(gamma+mu) * Iout;
    real Rout = (1-exp(-(delta+mu))) * R[i-1];
    real RtoS = delta/(delta+mu) * Rout;
    
    S[i] = S[i-1] - Sout + mu * pop + RtoS;
    I[i] = I[i-1] + StoI - Iout;
    R[i] = R[i-1] + ItoR - Rout;
    C[i] = StoI * rho;
    
    replenish[i] = (mu * pop + RtoS)/S[i];
  }
}

model {
  for (i in 2:52) {
    beta[i] ~ normal(beta[i-1], sigma);
  }
  beta[1] ~ normal(beta[52], sigma);
  
  npi ~ normal(1, 0.25);
  
  S0 ~ uniform(0, 1);
  I0 ~ normal(0, 0.001);
  rho ~ normal(0, 0.02);
  phi ~ normal(0, 10);
  omega ~ normal(0, 200);
  tau ~ normal(104, 26);
  sigma ~ normal(0, 1);
  
  for (i in 2:N) {
    cases[i] ~ neg_binomial_2(C[i], phi);
  }
}
