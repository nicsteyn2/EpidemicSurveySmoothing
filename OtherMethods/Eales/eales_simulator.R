
rm(list=ls())

library(tidyverse)
library(cmdstanr)
library(splines)

nIterations = 10000

# Set options
T = 50
nSamples = 5000

spline_degree=3
target_dist_between_knots=5

X = seq(1, T)
N = rep(nSamples, T)
num_knots <- ceiling((max(X) - min(X))/target_dist_between_knots) + 7
days_per_knot <- (max(X) - min(X))/(num_knots - 7)

# Specify helper functions
inv_logit <- function(Num) { 1/(1+exp(-Num)) }
gammaDist <- function(x, a, b) { (b^a) * (x^(a-1)) * exp(-b*x) / gamma(a) }

# Prepare data for stan
dat = list(
  num_data = length(X),
  num_knots = num_knots,
  knots = unname(seq(min(X) - 3*days_per_knot, max(X) + 3*days_per_knot, length.out = num_knots)),
  X = X,
  N = N,
  spline_degree = spline_degree,
  sigma = 0.1, rho = 0.0002
)

# Load the model
mdl = cmdstanr::cmdstan_model("OtherMethods/Eales/eales_simulator.stan",
                                     cpp_options = list(stan_threads=FALSE),
                                     stanc_options = list("O1"))

# Fit the model
fit = mdl$sample(data=dat,
                 chains=1, parallel_chains=1,
                 fixed_param=TRUE,
                 iter_sampling = nIterations, iter_warmup = 0) # Sample a large number to guarantee at least 100 good iterations


# Extract samples
samples = fit$draws(format="df")

# Extract logit Pt samples
a = as.matrix(samples[,grep("^a\\[", colnames(samples))])
nsamples = nrow(a)

# Setup for posterior estimation
num_basis = num_knots + spline_degree - 1
Xbasis = seq(min(X) - 3*days_per_knot, max(X) + 3*days_per_knot, 0.1)
Xest = seq(min(X), max(X), by=1)
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
  Y_pred[i,] = rbinom(length(Xest), size=Nest, prob=inv_logit(Y_array[i,]))
}

# Back-calculate growth-rates
Y_growth = matrix(NA, nrow=nsamples, ncol=length(Xest))
for (i in 1:nsamples) {
  first_d = diff(Y_array[i,])/diff(Xest)
  exp_term = exp(Y_array[i,])
  exp_term = exp_term[2:length(exp_term)]
  Y_growth[i,] <- c(NA, (first_d/(exp_term+1)))
}

# Convert into long-form dataframes
df_prev = Y_array %>%
  data.frame() %>%
  mutate(iter=seq(1,nIterations)) %>%
  pivot_longer(cols=-iter, names_to="V", values_to="value") %>%
  mutate(t = as.integer(sub("^X", "", V)), variable="Pt", value=inv_logit(value)) %>%
  select(-V)

df_growth = Y_growth %>%
  data.frame() %>%
  mutate(iter=seq(1,nIterations)) %>%
  pivot_longer(cols=-iter, names_to="V", values_to="value") %>%
  mutate(t = as.integer(sub("^X", "", V)), variable="rt") %>%
  select(-V)

df_nPos = Y_pred %>%
  data.frame() %>%
  mutate(iter=seq(1,nIterations)) %>%
  pivot_longer(cols=-iter, names_to="V", values_to="value") %>%
  mutate(t = as.integer(sub("^X", "", V)), variable="nPos") %>%
  select(-V)

df_nSamples = df_nPos %>% mutate(value=nSamples, variable="nSamples")

df_all = rbind(df_prev, df_growth, df_nPos, df_nSamples)

na_iters = df_all %>% filter(!(t==1 & variable=="rt")) %>% group_by(iter) %>% summarise(hasNA = any(is.na(value)))

# Find valid iterations
df_valid = df_all %>%
  filter(variable %in% c("Pt", "rt"), !(iter %in% na_iters$iter[na_iters$hasNA])) %>%
  group_by(iter, variable) %>%
  summarise(minval = min(value, na.rm=TRUE), maxval = max(value, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(isValid = ifelse( (variable=="rt" & minval >= -0.3 & maxval <= 0.3) |
                             (variable=="Pt" & minval >= 0.001 & maxval <= 0.3),
                           1, 0)) %>%
  group_by(iter) %>%
  summarise(isValid = all(isValid==1)) %>%
  filter(isValid==1)

# Randomly select 100 valid iterations
df_selected = df_valid %>% filter(iter %in% sample(unique(df_valid$iter), min(nrow(df_valid), 100), replace=FALSE)) %>%  select(-isValid) %>% mutate(newiter=1:min(nrow(df_valid), 100))

# Extract selected valid iterations
df_out = df_all %>% filter(iter %in% df_selected$iter) %>% left_join(df_selected, by="iter") %>% select(-iter) %>% rename(iter=newiter)

# Plot output for checking
plt = ggplot(df_out %>% filter(variable %in% c("Pt", "rt"))) +
  geom_line(aes(x=t, y=value, color=factor(iter))) +
  facet_wrap(~variable, scales="free_y", ncol=1)
plt

# Save output
write.csv(df_out, sprintf("paper/outputs/simulations_eales_T%d_nSamples%d.csv", T, nSamples), row.names=FALSE)






