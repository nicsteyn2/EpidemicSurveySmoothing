functions {
  #include gaussian_process.stan
  #include abbott_functions.stan
}

data {
  
  int<lower = 0> t; // number of time points to model
  array[t] int<lower = 0> nSamples; // number of samples taken

  int<lower = 0> pbt; // maximum detection time

  int<lower = 0> gtmax; // maximum number of days to consider for the generation time
  array[2] real<lower = 0> gtm; // mean and standard deviation (sd) of the mean generation time
  array[2] real<lower = 0> gtsd; // mean and sd of the sd of the generation time

  int<lower = 1> M; // approximate gp dimensions
  real L; // approximate gp boundary
  int<lower = 0, upper = 2> diff_order; // Order of differencing to use for the GP (0 = cases -> mean, 1 = growth_t -> 0, 2 = growth_t -> growth_{t-1})
  
  // Fix many of the parameters
  real<lower = 0> ell; // length scale of gp (used to be called rho, but we have changed it for consistency with our methods)
  real<lower = 0> alpha; // scale of gp
  real<lower = 0> rho; // dispersion parameter for beta-binomial
  real init_inc; // Initial infections
  vector[diff_order] init_growth;
  vector<lower = 0, upper = 1>[pbt] prob_detect; // probability of detection as a function of time since infection
  
}

transformed data {
  
  matrix[t - diff_order, M] PHI = setup_gp(M, L, t - diff_order);
  
}


generated quantities {
  
  vector[M] eta; // eta of gp
  vector[t] gp; // value of gp at time t + initialisation
  vector<lower = 0, upper = 1>[t] incidence; // incident infections at time t
  vector<lower = 0, upper = 1>[t] prevalence; // detectable cases at time t
  vector<lower = 0>[t] al;
  vector<lower = 0>[t] be;
  vector<lower = 0>[t] R;
  vector[t] r;
  vector[t] rprev;
  array[t] int predictive_positives;
  real<lower = 0> gtm_sample;
  real<lower = 0> gtsd_sample;
  
  // sample gp eta
  for (j in 1:M) {
    eta[j] = normal_rng(0, 1);
  }

  // update gaussian process
  gp[(1 + diff_order):t] = update_gp(PHI, M, L, alpha, ell, eta, 0);
  
  // setup differencing of the GP
  if (diff_order) {
    gp[1:diff_order] = init_growth;
    for (j in 1:diff_order) {
      gp = cumulative_sum(gp);
    }
  }
  
  // Calculate incidence and prevalence
  incidence = inv_logit(init_inc + gp);
  prevalence = fmin(convolve(incidence, prob_detect), 0.999999);
  
  // Binomial observation model
  al = prevalence .* ((1 / rho) - 1);
  be = (1 - prevalence) .* ((1 / rho) - 1);

  {
    // sample generation time
    gtm_sample = normal_rng(gtm[1], gtm[2]);
    gtsd_sample = normal_rng(gtsd[1], gtsd[2]);

    // calculate Rt using infections and generation time
    R[1] = positive_infinity();
    R[2:t] = calculate_Rt(incidence, gtm_sample, gtsd_sample, gtmax, 0);

    // calculate growth rates
    r[1] = positive_infinity();
    r[2:t] = calculate_growth(incidence);
    
    rprev[1] = positive_infinity();
    rprev[2:t] = calculate_growth(prevalence);
    
    // generate predictive samples
    predictive_positives = beta_binomial_rng(nSamples, al, be);
    
  }

}

