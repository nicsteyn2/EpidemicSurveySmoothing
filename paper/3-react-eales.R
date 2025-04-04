
rm(list=ls())

library(tidyverse)

source("OtherMethods/Eales/EalesMethod.R")

# Set options
start_round = 17
end_round = 19
iter_warmup = 200
iter_sampling = 300
max_treedepth = 15 # Set to 16 for round 1-to-19, else set to 15

# Load data
df = read.csv("data/reactdata.csv") %>%
  filter(round >= start_round & round <= end_round) %>%
  select(round, date, nSamples, nPos) %>%
  mutate(date=as.Date(date)) %>%
  arrange(date)

init_fn = function() { list(sigma=0.1, rho=0.0002) } # Prevent issues with warmup

start = proc.time()
results = runEalesMethod(df,
                         method="overdispersed",
                         iter_warmup=iter_warmup, iter_sampling=iter_sampling, max_treedepth = max_treedepth,
                         calculateRt=TRUE,
                         init_fn = init_fn)
end = proc.time() 
time_taken = end - start
time_taken[[3]]

states = results[[1]]
summary = results[[2]]

max(summary$rhat)
min(summary$ess_bulk)

# Save outputs
write.csv(states, paste0("paper/outputs/3-react/eales_r", start_round, "to", end_round, "_states.csv"))
write.csv(summary, paste0("paper/outputs/3-react/eales_r", start_round, "to", end_round, "_summary.csv"))

# Process coverage and widths
tmp = states %>% filter(variable=="nPos") %>% left_join(df, by="date") %>% filter(nSamples>0) %>% mutate(inPredInterval = (nPos >= lower) & (nPos <= upper))
coverage = mean(tmp$inPredInterval)

tmp = states %>% filter(variable=="Pt") %>% left_join(df,by="date") %>% filter(nSamples>0) %>% mutate(width = upper - lower)
meanwidth = 100*mean(tmp$width)

# Produce summary text
s = summary
outfile = file(paste0("paper/outputs/3-react/eales_r", start_round, "to", end_round, "_summary.txt"), "w")
cat("Summary of results...\n\n", file=outfile)
cat(paste0("Time taken = ", round(time_taken[3], digits=2), "s\n\n"), file=outfile)
cat(paste0("Sampling iterations = ", iter_sampling, "\n"), file=outfile)
cat(paste0("Warmup iterations = ", iter_warmup, "\n\n"), file=outfile)
cat(paste0("Max Rhat = ", round(max(s$rhat), digits=2), "\n"), file=outfile)
cat(paste0("Min ESS = ", round(min(s$ess_bulk), digits=2), "\n\n"), file=outfile)
cat(paste0("sigma = ", signif(s[s$variable=="sigma", "mean"], digits=2), " (", signif(s[s$variable=="sigma", "lower"], digits=2), ", ", signif(s[s$variable=="sigma", "upper"], digits=2), ")\n"), file=outfile)
cat(paste0("rho = ", 1e4*signif(s[s$variable=="rho", "mean"], digits=2), " (", 1e4*signif(s[s$variable=="rho", "lower"], digits=2), ", ", 1e4*signif(s[s$variable=="rho", "upper"], digits=2), ")\n\n"), file=outfile)
cat(paste0("Coverage of nPos = ", round(100*coverage, digits=2), "%\n"), file=outfile)
cat(paste0("Mean width of credible intervals on Pt = ", round(meanwidth, digits=2), "%\n"), file=outfile)
close(outfile)
