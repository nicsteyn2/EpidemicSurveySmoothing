// Copyright 2021 Oliver Eales
// Copyright 2017 Milad Kharratzadeh
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

data {
  int n;             // number of data points
  array[n] int Y;
  array[n] int N;
}

parameters {
  vector[n] a_raw;
  real<lower=0> sigma;
  real<lower=0, upper=1> rho;
}


transformed parameters {
  vector[n] a;
  a[1] = a_raw[1];
  a[2] = a_raw[2];
  for (i in 3:n)
    a[i] = 2*a[i-1] - a[i-2] + a_raw[i]*sigma;

  real<lower=0> effective_n = (1/rho) - 1;
  vector[n] P = inv_logit(a);
  vector[n] alpha = P * effective_n;
  vector[n] beta = (1 - P) * effective_n;
  
}

model {
  // Priors
  sigma ~ inv_gamma(0.0001, 0.0001);
  a_raw[3:n] ~ normal(0,1);

  //Likelihood
  Y ~ beta_binomial(N, alpha, beta);
  
}

generated quantities {
  
  array[n] int predictive_positives;
  predictive_positives = beta_binomial_rng(N, alpha, beta);
  
}
