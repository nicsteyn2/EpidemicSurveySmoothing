
library(tidyverse)
library(cmdstanr)
library(splines)


# # # For testing
# df = read.csv("data/reactdata.csv") %>%
#   filter(round >= 17 & round <= 19) %>%
#   select(round, date, nSamples, nPos) %>%
#   mutate(date=as.Date(date)) %>%
#   arrange(date) %>%
#   complete(date = seq.Date(min(date)-3, max(date)+3, by="day"), fill=list(nSamples=0, nPos=0))
# 
# iter_warmup = 200
# iter_sampling = 300
# chains = 3
# max_treedepth = 15
# init_fn = function() { list(sigma=0.05, rho=0.0002) }

runEalesMethodNoSpline = function(df,
                          iter_warmup=200, iter_sampling=300, chains=3,
                          max_treedepth = 15,
                          init_fn = "random") {
  
  warning("This function requires input df to be \"completed\". There should not be any temporal gaps.")
  
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
  
  # Prepare data for stan
  data_list = list(n=length(Y), Y=Y, N=N)
  
  # Load stan model
  stan_model = cmdstanr::cmdstan_model("OtherMethods/Eales/eales_nospline.stan",
                                       cpp_options = list(stan_threads=FALSE),
                                       stanc_options = list("O1"))
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
    rename(lower = `2.5%`, upper = `97.5%`) %>%
    mutate(t = as.numeric(str_extract(variable, "(?<=\\[)\\d+(?=\\])")),
           variable = sub("\\[.*", "", variable))
  
  return(summary)
  
}





