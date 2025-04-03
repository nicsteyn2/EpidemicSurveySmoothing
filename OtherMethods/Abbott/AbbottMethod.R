
library(tidyverse)
library(cmdstanr)

# # For testing
# windin = 50
# windout = 21
# start_round = 17
# end_round = 19
# 
# iter_warmup = 200
# iter_sampling = 300
# chains = 3
# gp_m = 0.3
# alpha = 0.00001
# beta = 0.00001
# differencing = 1
# gentime_mean = c(mean=3.64, mean_sd=0.71)
# gentime_sd = c(sd=3.08, sd_sd=0.77)
# gentime_max = c(max=15)
# init_inc_mean = 0.01
# inc2prevmethod = "Hellewell"
# init_fn = function() { list(ell=10, alpha=0.07, rho=0.0002) } # Prevent issues with warmup
# 
# df = read.csv("data/reactdata.csv") %>%
#   rename(nSamples=n) %>%
#   select(date, round, nSamples, nPos) %>%
#   filter(round >= start_round & round <= end_round) %>%
#   mutate(date=as.Date(date))
# 
# df = read.csv("paper/outputs/simulations_steyn.csv") %>%
#   mutate(simulation="Steyn") %>%
#   filter(iter==1) %>% 
#   pivot_wider(names_from="variable", values_from="value")


# df is a sorted dataframe containing the data. It must contain columns nPos, nSamples, and either t or date
runAbbottMethod = function(df,
                           windin = 50, windout = 21,
                           iter_warmup = 200, iter_sampling = 300, chains = 3,
                           gp_m = 0.3,
                           alpha = 0.0001, # Lengthscale prior parameter, for default informative use 6.81
                           beta = 0.0001, # Lengthscale prior parameter, for default informative use 200.12
                           differencing = 1,
                           gentime_mean = c(mean=3.64, mean_sd=0.71), # Prior on the mean of the generation time distribution
                           gentime_sd = c(sd=3.08, sd_sd=0.77), # Prior on the standard deviation of the generation time distribution
                           gentime_max = c(max=15), # Assumed maximum generation times
                           init_inc_mean = 0.01,
                           init_fn = NULL,
                           inc2prevmethod = "Hellewell"
) {
  
  # Incorporate windin and windout periods into the data
  if ("date" %in% colnames(df)) {
    df_extended = df %>%
      complete(date=seq(min(date)-windin, max(date)+windout, by="day"), fill=list(nSamples=0, nPos=0)) %>%
      arrange(date)
  } else {
    df_extended = df %>%
      complete(t=seq(min(t)-windin, max(t)+windout), fill=list(nSamples=0, nPos=0)) %>%
      arrange(t)
  }
  
  # Load probability of detection curves
  if (inc2prevmethod == "Hellewell") {
    df_probdetect = read.csv("data/PCR-curves/hellewell.csv") %>%
      pivot_longer(-sample, names_to = "time", values_to = "p") %>%
      mutate(time = as.numeric(sub("X", "", time))) %>%
      group_by(time) %>%
      summarise(mean=signif(mean(p)), median=signif(median(p)), sd=signif(sd(p)))
  } else if (inc2prevmethod == "Binny") {
    df_probdetect = read.csv("data/PCR-curves/binny.csv") %>%
      rename(time=days_since_infection, lower=lower_95, upper=upper_95) %>%
      filter(time %in% seq(0, floor(max(time)))) %>%
      mutate(mean = (upper+lower)/2, sd = (upper-lower)/(2*1.96))
  } else {
    stop("Inc2prev method not recognised.")
  }
  
  # Specify data for Stan
  dat = list(
    t = nrow(df_extended), nPos = df_extended$nPos, nSamples = df_extended$nSamples,
    prob_detect_mean = rev(df_probdetect$mean), prob_detect_sd = rev(df_probdetect$sd),
    pbt = max(df_probdetect$time) + 1,
    gtm = gentime_mean, gtsd = gentime_sd,
    gtmax = gentime_max,
    lengthscale_alpha = alpha, lengthscale_beta = beta,
    M = ceiling(nrow(df_extended) * gp_m),
    L = 2,
    diff_order = differencing,
    init_inc_mean = log(init_inc_mean / (1 - init_inc_mean))
  )

  
  # Load the model
  mod = cmdstanr::cmdstan_model("OtherMethods/Abbott/abbott_model.stan",
                                include_path = "OtherMethods/Abbott/",
                                cpp_options = list(stan_threads=FALSE),
                                stanc_options = list("O1"))
  
  # Fit the model
  fit = mod$sample(data=dat,
                   chains=chains, parallel_chains=chains,
                   iter_warmup=iter_warmup, iter_sampling=iter_sampling,
                   adapt_delta = 0.9, max_treedepth=15,
                   init = init_fn)
  
  # Extract fit summary and tidy up
  df_out = fit$summary(NULL,
                       posterior::default_summary_measures()[1:3],
                       quantile = ~quantile(.x, probs = c(0.025, 0.975), na.rm=TRUE),
                       posterior::default_convergence_measures()[1:2]) %>%
    rename(lower = `2.5%`, upper = `97.5%`) %>%
    # Extract time-indices from variable names and clean-up variable names
    mutate(t = as.numeric(str_extract(variable, "(?<=\\[)\\d+(?=\\])")),
           variable = sub("\\[.*", "", variable)) %>%
    # Separate out different kinds of time
    mutate(days = ifelse(variable=="prob_detect", t, NA),
           t_approximation = ifelse(variable=="eta", t, NA),
           t = ifelse(variable %in% c("prob_detect", "eta", "init_growth"), NA, t)) %>%
    # Remove windin and windout periods and renormalise time
    mutate(t = t - windin) %>%
    filter((t >= 1 & t <= max(t, na.rm=TRUE) - windout) | is.na(t)) %>%
    # Append nice variable names
    mutate(variable = case_when(
      variable == "predictive_positives" ~ "nPos",
      variable == "incidence" ~ "It",
      variable == "prevalence" ~ "Pt",
      variable == "rprev" ~ "rt",
      variable == "r" ~ "rtinc",
      variable == "R" ~ "Rt",
      TRUE ~ variable
    ))
  
  # Append dates if available
  if ("date" %in% colnames(df)) {
    df_out = df_out %>% mutate(date = min(df$date) + t - 1)
  }


  # Return output
  return(df_out)
  
}
