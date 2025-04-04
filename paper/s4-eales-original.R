
# Code to fit original (take modified results from main)

rm(list=ls())

source("OtherMethods/Eales/EalesMethod.R")

start_round = 17
end_round = 19

# Set options
method = "binomial"
target_knot_spacing = 5
max_treedepth = 10
iter_warmup = 1000
iter_sampling = 19000
init_fn = function() { list(sigma=0.05) }


# Fit original model, target_dist = 5
df = read.csv("data/reactdata.csv") %>%
  filter(round >= start_round & round <= end_round) %>%
  select(round, date, nSamples, nPos) %>%
  mutate(date=as.Date(date)) %>%
  arrange(date)

start = proc.time()
out = runEalesMethod(df, method="original", target_dist_between_knots=5,
                     max_treedepth=max_treedepth, iter_warmup=iter_warmup, iter_sampling=iter_sampling,
                     init_fn=init_fn)
end = proc.time()
time_taken = end - start
time_taken[[3]]

results = out[[1]]
summary = out[[2]]

summary$timetofit = time_taken[[3]]

# Save results
write.csv(results, sprintf("paper/outputs/s4-eales-original/original_r%dto%d_results.csv", start_round, end_round), row.names=FALSE)
write.csv(summary, sprintf("paper/outputs/s4-eales-original/original_r%dto%d_summary.csv", start_round, end_round), row.names=FALSE)



# Code to fit modified
rm(list=ls())

source("OtherMethods/Eales/EalesMethod.R")

start_round = 17
end_round = 19

# Set options
method = "betabinom"
target_knot_spacing = 5
max_treedepth = 15
iter_warmup = 200
iter_sampling = 300
init_fn = function() { list(sigma=0.05, rho=2e-4) }


# Fit original model, target_dist = 5
df = read.csv("data/reactdata.csv") %>%
  filter(round >= start_round & round <= end_round) %>%
  select(round, date, nSamples, nPos) %>%
  mutate(date=as.Date(date)) %>%
  arrange(date)

start = proc.time()
out = runEalesMethod(df, method=method, target_dist_between_knots=5,
                     max_treedepth=max_treedepth, iter_warmup=iter_warmup, iter_sampling=iter_sampling,
                     init_fn=init_fn)
end = proc.time()
time_taken = end - start
time_taken[[3]]

results = out[[1]]
summary = out[[2]]

summary$timetofit = time_taken[[3]]

# Save results
write.csv(results, sprintf("paper/outputs/s4-eales-original/modified_r%dto%d_results.csv", start_round, end_round), row.names=FALSE)
write.csv(summary, sprintf("paper/outputs/s4-eales-original/modified_r%dto%d_summary.csv", start_round, end_round), row.names=FALSE)

