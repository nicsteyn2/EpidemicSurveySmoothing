
library(tidyverse)
library(cmdstanr)
library(splines)
library(extraDistr)


# # For testing
# rm(list=ls())
# df = read.csv("data/example.csv")
# 
# method = "overdispersed"
# spline_degree = 3
# target_dist_between_knots = 5
# iter_sampling = 300
# iter_warmup = 200
# max_treedepth = 15
# chains = 3
# calculateRt = TRUE
# tau_max = 14
# gentime_a = 2.29
# gentime_b = 0.36
# init_fn = function() { list(sigma=0.05, rho=0.0002) }
# 
# df = read.csv("data/reactdata.csv") %>%
#   filter(round >= 17 & round <= 19) %>%
#   mutate(date=as.Date(date))


runEalesMethod = function(df,
                          method = "overdispersed",
                          spline_degree=3, target_dist_between_knots=5,
                          iter_warmup=200, iter_sampling=300, chains=3,
                          max_treedepth = 15,
                          calculateRt = FALSE, tau_max = 14, gentime_a = 2.29, gentime_b = 0.36,
                          init_fn = "random") {
  
  # Process data
  Y = df$nPos
  N = df$nSamples
  if ("date" %in% colnames(df)) {
    X = as.numeric(df$date - min(df$date) + 1)
  } else if ("t" %in% colnames(df)) {
    X = df$t
  } else {
    stop("No date or t column found in the input dataframe.")
  }
  
  # Setup parallel processing
  options(mc.cores = chains)
  
  # Specify helper functions
  inv_logit <- function(Num) { 1/(1+exp(-Num)) }
  
  gammaDist <- function(x, a, b) { (b^a) * (x^(a-1)) * exp(-b*x) / gamma(a) }
  
  # Pre-process data
  min_X = min(X)
  max_X = max(X)
  num_knots <- ceiling((max_X - min_X)/target_dist_between_knots) + 7
  days_per_knot <- (max_X - min_X)/(num_knots - 7)
  print(paste0("Realised days per knot = ", days_per_knot))
  
  # Prepare data for stan
  data_list = list(
    num_data = length(X),
    num_knots = num_knots,
    knots = unname(seq(min_X - 3*days_per_knot, max_X + 3*days_per_knot, length.out = num_knots)),
    X = X,
    Y = Y,
    N = N,
    spline_degree = spline_degree
  )
  
  # Load stan model
  if (method == "original" | method=="binomial") {
    stan_model = cmdstanr::cmdstan_model("OtherMethods/Eales/eales.stan",
                                         cpp_options = list(stan_threads=FALSE),
                                         stanc_options = list("O1"))
  } else if (method == "overdispersed" | method=="betabinom") {
    stan_model = cmdstanr::cmdstan_model("OtherMethods/Eales/eales_betabinom.stan",
                                         cpp_options = list(stan_threads=FALSE),
                                         stanc_options = list("O1"))
  } else {
    stop("Method not recognized")
  }
  
  # Run stan model
  fit = stan_model$sample(data=data_list,
                          chains=chains, parallel_chains=chains,
                          iter_warmup=iter_warmup, iter_sampling=iter_sampling,
                          adapt_delta = 0.9, max_treedepth=max_treedepth,
                          init = init_fn)
  
  # Extract MCMC summary
  summary = fit$summary(NULL,
                        posterior::default_summary_measures()[1:3],
                        quantile = ~quantile(.x, probs = c(0.025, 0.975), na.rm=TRUE),
                        posterior::default_convergence_measures()[1:2]) %>%
    rename(lower = `2.5%`, upper = `97.5%`)
  
  # Extract samples
  samples = fit$draws(format="df")
  
  # Extract logit Pt samples
  a = as.matrix(samples[,grep("^a\\[", colnames(samples))])
  if (method=="overdispersed" | method == "betabinom") { rho = as.matrix(samples[,"rho"]) }
  nsamples = nrow(a)
  
  # Setup for posterior estimation
  num_basis = num_knots + spline_degree - 1
  Xbasis = seq(min_X - 3*days_per_knot, max_X + 3*days_per_knot, 0.1)
  Xest = seq(min_X, max_X, by=1)
  Nest = rep(0, length(Xest))
  Nest[Xest %in% X] = N
  
  # Calculate basis function matrix
  B = splines::bs(Xbasis, df=num_basis, degree=spline_degree, intercept=TRUE)
  B = predict(B, Xest)
  
  # Extract prevalence (on logit scale) and predictive positives 
  Y_array = matrix(NA, nrow=nsamples, ncol=length(Xest))
  Y_pred = matrix(NA, nrow=nsamples, ncol=length(Xest))
  for (i in 1:nsamples) {
    a_i = as.numeric(a[i,])
    Y_array[i,] = as.vector(B %*% a_i)
    if (method == "original" | method=="binomial") {
      Y_pred[i,] = rbinom(length(Xest), size=Nest, prob=inv_logit(Y_array[i,]))
    } else if (method == "overdispersed" | method == "betabinom") {
      P = inv_logit(Y_array[i,])
      al = P * ((1/rho[i]) - 1)
      be = (1 - P) * ((1/rho[i]) - 1)
      Y_pred[i,] = rbbinom(length(Xest), size=Nest, alpha=al, beta=be)
    }
  }
  
  # Back-calculate growth-rates
  Y_growth = matrix(NA, nrow=nsamples, ncol=length(Xest))
  for (i in 1:nsamples) {
    first_d = diff(Y_array[i,])/diff(Xest)
    exp_term = exp(Y_array[i,])
    first_d = first_d[1:length(first_d)]
    exp_term = exp_term[2:length(exp_term)]
    Y_growth[i,] <- c(NA, (first_d/(exp_term+1)))
  }
  
  df_prev = data.frame(t=Xest,
                       mean = apply(inv_logit(Y_array), 2, mean),
                       lower = apply(inv_logit(Y_array), 2, function(x) quantile(x, 0.025)),
                       upper = apply(inv_logit(Y_array), 2, function(x) quantile(x, 0.975)),
                       variable = "Pt")
  
  df_growth = data.frame(t=Xest,
                         mean = apply(Y_growth, 2, mean),
                         lower = apply(Y_growth, 2, function(x) quantile(x, 0.025, na.rm=TRUE)),
                         upper = apply(Y_growth, 2, function(x) quantile(x, 0.975, na.rm=TRUE)),
                         variable = "rt")
  
  # Back-calculate reproduction numbers
  if (calculateRt) {
    
    Y_prev_tmp = inv_logit(Y_array)
    Y_Rt = matrix(NA, nrow=nsamples, ncol=length(Xest))
    gamma_vals = gammaDist(seq_len(tau_max), gentime_a, gentime_b)
    for (i in 1:nsamples) {
      for (tt in (tau_max+1):length(Xest)) {
        integral <- sum(Y_prev_tmp[i,(tt-tau_max):(tt-1)] * gamma_vals)
        Y_Rt[i,tt] <- Y_prev_tmp[i,tt] / (integral / sum(gamma_vals))
      }
    }
    
    df_repro = data.frame(t=Xest,
                          mean = apply(Y_Rt, 2, mean),
                          lower = apply(Y_Rt, 2, function(x) quantile(x, 0.025, na.rm=TRUE)),
                          upper = apply(Y_Rt, 2, function(x) quantile(x, 0.975, na.rm=TRUE)),
                          variable = "Rt")
  } else {
    df_repro = data.frame()
  }
  
  nPosForPred = numeric(length(Xest)) * NA
  nPosForPred[X] = df$nPos
  df_pred = data.frame(t=Xest,
                       mean = apply(Y_pred, 2, mean),
                       lower = apply(Y_pred, 2, function(x) quantile(x, 0.025)),
                       upper = apply(Y_pred, 2, function(x) quantile(x, 0.975)),
                       variable = "nPos")
  
  df_out = rbind(df_prev, df_growth, df_repro, df_pred)
  
  # Append date if available
  if ("date" %in% colnames(df)) {
    df_out$date = df_out$t + min(df$date) - 1
  }
  
  return(list(df_out, summary))
  
}





