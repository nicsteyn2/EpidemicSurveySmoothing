// A collection of useful functions from EpiNow2 and inc2prev

// General functions
vector convolve(vector vec1, vector vec2) {
  int lvec1 = num_elements(vec1);
  int lvec2 = num_elements(vec2);
  vector[lvec1] conv = rep_vector(1e-8, lvec1);
  for (s in 1:lvec1) {
    conv[s] += dot_product(vec1[max(1, (s - lvec2 + 1)):s], 
                           tail(vec2, min(lvec2, s)));
  }
  return(conv);
}

real normal_lub_rng(real mu, real sigma, real lb, real ub) {
  real p_lb = normal_cdf(lb | mu, sigma);
  real p_ub = normal_cdf(ub | mu, sigma);
  real u = uniform_rng(p_lb, p_ub);
  real y = mu + sigma * inv_Phi(u);
  return y;
}


vector discretised_gamma_pmf(array[] int y, real mu, real sigma, int max_val) {
  int n = num_elements(y);
  vector[n] pmf;
  real trunc_pmf;
  // calculate alpha and beta for gamma distribution
  real small = 1e-5;
  real large = 1e8;
  real c_sigma = sigma < small ? small : sigma;
  real c_mu = mu < small ? small : mu;
  real alpha = ((c_mu) / c_sigma)^2;
  real beta = (c_mu) / (c_sigma^2);
  // account for numerical issues
  alpha = alpha < small ? small : alpha;
  alpha = alpha > large ? large : alpha;
  beta = beta < small ? small : beta;
  beta = beta > large ? large : beta;
  // calculate pmf
  trunc_pmf = gamma_cdf(max_val + 1 | alpha, beta) - gamma_cdf(1 | alpha, beta);
  for (i in 1:n){
    pmf[i] = (gamma_cdf(y[i] + 1 | alpha, beta) - gamma_cdf(y[i] | alpha, beta)) /
    trunc_pmf;
  }
  return(pmf);
}


// Epi-specific functions

real update_infectiousness(vector infections, vector gt_pmf,
                           int max_gt, int index){
  int inf_start = max(1, (index - max_gt));
  int inf_end = (index - 1);
  int pmf_accessed = min(max_gt, index - 1);
  real new_inf = dot_product(infections[inf_start:inf_end], tail(gt_pmf, pmf_accessed));
  return(new_inf);
}


vector calculate_growth(vector infections) {
  int t = num_elements(infections);
  vector[t-1] growth = log(infections[2:t]) - log(infections[1:(t-1)]);
  return(growth);
}


vector calculate_Rt(vector infections,
                    real gt_mean, real gt_sd, int max_gt,
                    int smooth) {
  vector[max_gt] gt_pmf;
  array[max_gt] int gt_indexes;
  int t = num_elements(infections);
  vector[t-1] R;
  vector[t-1] sR;
  vector[t-1] infectiousness = rep_vector(1e-5, t);
  
  // calculate PMF of the generation time
  for (i in 1:(max_gt)) {
    gt_indexes[i] = max_gt - i + 1;
  }
  gt_pmf = discretised_gamma_pmf(gt_indexes, gt_mean, gt_sd, max_gt);
  
  // calculate Rt using the renewal model
  for (s in 1:(t-1)) {
    infectiousness[s] += update_infectiousness(infections, gt_pmf, max_gt, s);
    R[s] = infections[s+1] / infectiousness[s];
  }
  
  // This is unecessary given we are passing in smoothed data. A relic from old code? Disabling by default but leaving in to test.
  if (smooth) {
    for (s in 1:(t-1)) {
      real window = 0;
      sR[s] = 0;
      for (i in max(1, s - smooth):min(t-1, s + smooth)) {
        sR[s] += R[i];
        window += 1;
      }
      sR[s] = sR[s] / window;
    }
  } else{
    sR = R;
  }
  return(sR);
}
