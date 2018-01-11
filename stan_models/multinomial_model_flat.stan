# this is the same as: threshold_test/model_flat_informative_priors_tight.stan

functions {
  real beta_ccdf(real x, real a, real b) {
    return 1 - beta_cdf(x, a, b);
  }
  
  // E[X | X > x] for X ~ beta(a,b)
  real beta_conditional_mean(real x, real a, real b) {
    return beta_ccdf(x, a+1, b) / beta_ccdf(x, a, b) * a / (a+b);
  }
}

data {
  int<lower=1> N; // number of observations
  int<lower=1> R; // number of suspect races
  int<lower=1> D; // number of counties
  
  int<lower=1,upper=R> r[N]; // race of suspect
  int<lower=1,upper=D> d[N]; // county where stop occurred
  
    int<lower=0> s[N]; // # of searches
    int<lower=0> h[N]; // # of successful searches (hits)
    
  real<lower=0,upper=1> base_population_proportions[N]; // the proportions in the base population. 
}

parameters {	
  // hyperparameters
  vector<lower=0>[R] sigma_t;
  //   real mu_phi_d;
  //   real<lower=0> sigma_phi_d;
  // real mu_lambda_d;
  // real<lower=0> sigma_lambda_d;
  
  // search thresholds
  vector[R] t_r;
  vector[N] t_i_raw;
  
  // parameters for signal distribution
  vector[R] phi_r;
  // vector[D] phi_d;
  vector[D-1] phi_d_raw;
  
  vector[R] lambda_r; 
  // vector[D] lambda_d; 
  vector[D-1] lambda_d_raw;
}

transformed parameters {
  vector[D] phi_d;
  vector[D] lambda_d;
  vector[N] phi;
  vector[N] lambda;
  vector<lower=0, upper=1>[N] search_rate;
  vector<lower=0, upper=1>[N] hit_rate;
  vector<lower=0, upper=1>[N] t_i;
  vector<lower=0, upper=1>[N] theta;
  
  phi_d[1]      <- 0;
  phi_d[2:D]    <- phi_d_raw;
  lambda_d[1]   <- 0;
  lambda_d[2:D] <- lambda_d_raw;
  
  for (i in 1:N) {	
    real a;
    real b;
    
    // implies t_i[i] ~ logit-normal(t_r[r[i]], sigma_t[r[i]])
    t_i[i] <- inv_logit(t_r[r[i]] + t_i_raw[i] * sigma_t[r[i]]);
    
    // signal distribution parameters	
    phi[i]    <- inv_logit(phi_r[r[i]] + phi_d[d[i]]);
    lambda[i] <- exp(lambda_r[r[i]] + lambda_d[d[i]]);
    
    // transformed signal distribution parameters
    a <- lambda[i] * phi[i];
    b <- lambda[i] * (1 - phi[i]);
    
    // implied search and hit rates
    search_rate[i] <- beta_ccdf(t_i[i], a, b);
    hit_rate[i]    <- beta_conditional_mean(t_i[i], a, b);
    
    theta[i] <- base_population_proportions[i] * search_rate[i];
  }
   
}

model {  
  // Draw threshold hyperparameters
  sigma_t ~ normal(0, 0.25);
  t_r     ~ normal(-2, 1);  
  
  // Draw race and department specific thresholds
  t_i_raw ~ normal(0, 1);           // implies t_i ~ logit-normal(t_r, sigma_t)
  
  // Draw race parameters
  phi_r    ~ normal(-3.5, 0.25);
  lambda_r ~ normal(2, 0.5);
  
  // Draw department parameters (for un-pinned departments)
  phi_d_raw    ~ normal(0, 0.25);    
  lambda_d_raw ~ normal(0, 0.5);     
  
  // Draw search and hit observations
  //multinomial probability
  for(i in 1:D){
    s[(R * (i - 1) + 1):(R * i)] ~ multinomial(theta[(R * (i - 1) + 1):(R * i)] / sum(theta[(R * (i - 1) + 1):(R * i)]));
  }

  //hit rate
  h ~ binomial(s, hit_rate);
}



generated quantities {
  // Stop-weighted per-race parameters
}
