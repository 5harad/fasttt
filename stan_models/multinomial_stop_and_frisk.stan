data {
  int<lower=1> N; // number of observations
  int<lower=1> R; // number of suspect races
  int<lower=1> D; // number of counties
  
  int<lower=1,upper=R> r[N]; // race of suspect
  int<lower=1,upper=D> d[N]; // county where stop occurred
  
  real<lower=0,upper=1> base_population_proportions[N];
    int<lower=0> s[N]; // # of searches
    int<lower=0> h[N]; // # of successful searches (hits)
}

parameters {	
  // hyperparameters
  real<lower=0> sigma_t; #standard deviation for the normal the thresholds are drawn from. 
  
  // search thresholds
  vector[N] t_i;
  
  // parameters for signal distribution
  vector[R] phi_r;
  vector[D-1] phi_d_raw;
  real mu_phi;
  
  vector[R] delta_r; 
  vector[D-1] delta_d_raw;
  real mu_delta;
}

transformed parameters {
  vector[D] phi_d;
  vector[D] delta_d;
  vector[N] phi;
  vector[N] delta;
  vector<lower=0, upper=1>[N] search_rate;
  vector<lower=0, upper=1>[N] hit_rate;
  real innocent_p_above_threshold;
  real guilty_p_above_threshold;
  vector<lower=0, upper=1>[N] theta;

  phi_d[1]      <- 0;
  phi_d[2:D]    <- phi_d_raw;
  delta_d[1]   <- 0;
  delta_d[2:D] <- delta_d_raw;
  
  for (i in 1:N) {	
    // phi is the fraction of people of race r, d who are guilty (ie, carrying contraband)
    phi[i]    <- inv_logit(phi_r[r[i]] + phi_d[d[i]]);
    
    // mu is the center of the guilty distribution. 
    delta[i] <- exp(delta_r[r[i]] + delta_d[d[i]]);
    
    guilty_p_above_threshold <- phi[i] * (1 - normal_cdf(t_i[i], delta[i], 1));
    innocent_p_above_threshold <- (1 - phi[i]) * (1 - normal_cdf(t_i[i], 0, 1));
    search_rate[i] = (guilty_p_above_threshold + innocent_p_above_threshold);
    hit_rate[i] = guilty_p_above_threshold / (guilty_p_above_threshold + innocent_p_above_threshold);
    theta[i] = base_population_proportions[i] * search_rate[i];
   }
}

model {  
  // Draw threshold hyperparameters
  sigma_t ~ normal(0, 1);
  
  // Draw race parameters. Each is centered at a mu, and we allow for inter-race heterogeneity. 
  mu_phi ~ normal(0, 1);#normal(-3, 1);
  mu_delta ~ normal(0, 1);
  
  phi_r    ~ normal(mu_phi, 1);
  delta_r ~ normal(mu_delta, 1);
  
  // Draw department parameters (for un-pinned departments)
  phi_d_raw    ~ normal(0, 0.1);    
  delta_d_raw ~ normal(0, 0.1); 
  
  //thresholds
  t_i ~ normal(0, sigma_t);
  
  //multinomial probability
  for(i in 1:D){
    s[(R * (i - 1) + 1):(R * i)] ~ multinomial(theta[(R * (i - 1) + 1):(R * i)] / sum(theta[(R * (i - 1) + 1):(R * i)]));
  }
  
  //hit rate
  h ~ binomial(s, hit_rate);
  
}



generated quantities {
}
