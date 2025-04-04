
rm(list=ls())

library(tidyverse)
library(cmdstanr)
library(posterior) # Used to extract prior predictive samples in a convenient format

logit = function(x) { return(log(x/(1-x)))}
invlogit = function(x) { return(1/(1 + exp(-x))) }

# Set options
T = 100
nSamples = 5000

# Setup windin and windout periods
windin = 50
windout = 21

# Load probability of detection curves
df_probdetect = read.csv("data/PCR-curves/hellewell.csv") %>%
  pivot_longer(-sample, names_to = "time", values_to = "p") %>%
  mutate(time = as.numeric(sub("X", "", time))) %>%
  group_by(time) %>%
  summarise(mean=signif(mean(p)), median=signif(median(p)), sd=signif(sd(p)))

# Specify data for Stan
dat = list(
  
  # Specify length of data and assumed daily samples
  t = T+windin+windout,
  nSamples = rep(nSamples, T+windin+windout),
  
  # Gaussian process parameters (similar to those estimated from REACT-1 data)
  ell = 10,
  alpha = 0.08,
  rho = 0.0002,
  init_inc = logit(0.001),
  init_growth = 0,
  
  # Fix the probability of detection to the pointwise mean
  pbt = max(df_probdetect$time) + 1,
  prob_detect = rev(df_probdetect$mean),
  
  # Remaining hyperparameters
  gtm = c(mean=3.64, mean_sd=0.71), gtsd = c(sd=3.08, sd_sd=0.77),
  gtmax = c(max=15),
  M = ceiling((T+windin+windout) * 0.3), # The 0.3 is gp_m from the full model
  L = 2,
  diff_order = 1
  
)

# Load the model
mod = cmdstanr::cmdstan_model("OtherMethods/Abbott/abbott_simulator.stan",
                              include_path = "OtherMethods/Abbott/",
                              cpp_options = list(stan_threads=FALSE),
                              stanc_options = list("O1"))

# Fit the model
fit = mod$sample(data=dat,
                 chains=1, parallel_chains=1,
                 fixed_param=TRUE,
                 iter_sampling = 10000, iter_warmup = 0) # Sample a large number to guarantee at least 100 good iterations


# Extract output in a tidy format
df_samples = as_draws_df(fit$draws()) %>%
  select(-.draw, -.chain) %>%
  pivot_longer(cols=-c(.iteration), names_to="variable", values_to="value") %>%
  rename(iter=.iteration) %>%
  mutate(t = as.numeric(str_extract(variable, "(?<=\\[)\\d+(?=\\])")),
         variable = sub("\\[.*", "", variable)) %>%
  mutate(variable = case_when(
    variable == "predictive_positives" ~ "nPos",
    variable == "incidence" ~ "It",
    variable == "prevalence" ~ "Pt",
    variable == "rprev" ~ "rt",
    variable == "r" ~ "rtinc",
    variable == "R" ~ "Rt",
    TRUE ~ variable)
  )

# Find valid simulations
df_valid = df_samples %>%
  filter(variable %in% c("rt", "Pt"), t>=windin, t<=max(t, na.rm=TRUE)-windout) %>%
  group_by(iter, variable) %>%
  summarise(minval = min(value), maxval = max(value)) %>%
  mutate(isValid = ifelse( (variable=="rt" & minval >= -0.3 & maxval <= 0.3) |
                             (variable=="Pt" & minval >= 0.001 & maxval <= 0.3),
                           1, 0)) %>%
  group_by(iter) %>%
  summarise(isValid = all(isValid==1)) %>%
  filter(isValid==1)

# Randomly select 100 valid simulations
df_selected = df_valid %>% filter(iter %in% sample(unique(df_valid$iter), 100, replace=FALSE)) %>%  select(-isValid) %>% mutate(newiter=1:100)

# Extract selected valid simulations
df_out = df_samples %>%
  filter(iter %in% df_selected$iter) %>% left_join(df_selected, by="iter") %>% select(-iter) %>% rename(iter=newiter) %>%
  mutate(t = t - windin + 1) %>%
  filter(t >= 1, t <= T)

# Append nSamples
df_nSamples = df_out %>% filter(variable=="rt") %>% mutate(variable="nSamples", value=nSamples)
df_out = rbind(df_out, df_nSamples)

# Plot output for checking
plt = ggplot(df_out %>% filter(variable %in% c("rt", "Pt"), iter<=10)) +
  geom_line(aes(x=t, y=value, color=factor(iter))) +
  facet_wrap(~variable, scales="free_y", ncol=1)
plt

# Save output
df_out_clean = df_out %>% filter(variable %in% c("It", "Pt", "Rt", "rtinc", "rt", "nPos", "nSamples"))
write.csv(df_out_clean, "paper/outputs/2-simcompare/simulations_abbott.csv", row.names=FALSE)

