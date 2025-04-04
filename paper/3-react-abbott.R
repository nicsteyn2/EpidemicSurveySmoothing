
rm(list=ls())

source("OtherMethods/Abbott/AbbottMethod.R")

# Set options
start_round = 17
end_round = 19
iter_warmup = 200
iter_sampling = 300

# Load data and append Abbott method-specific windin and windout periods
df = read.csv("data/reactdata.csv") %>%
  filter(round >= start_round & round <= end_round) %>%
  select(date, round, nSamples, nPos) %>%
  mutate(date=as.Date(date)) %>%
  arrange(date)

# Fit the model
start = proc.time()
results = runAbbottMethod(df, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                          gentime_mean = c(mean=6.3611, mean_sd=0.001), gentime_sd = c(sd=4.2035, sd_sd=0.001)) # Set sd of gentime params to be near zero to treat this as fixed
end = proc.time()
time_taken = end - start

# Save results
write.csv(results, paste0("paper/outputs/3-react/abbott_r", start_round, "to", end_round, "_results.csv"), row.names=FALSE)

# Process coverage and widths
tmp = results %>% filter(variable=="nPos") %>% left_join(df, by="date") %>% filter(nSamples>0) %>% mutate(inPredInterval = (nPos >= lower) & (nPos <= upper))
coverage = mean(tmp$inPredInterval)

tmp = results %>% filter(variable=="Pt") %>% left_join(df, by="date") %>% filter(nSamples>0) %>% mutate(width = upper - lower)
meanwidth = 100*mean(tmp$width)


# Produce summary text
outfile = file(paste0("paper/outputs/3-react/abbott_r", start_round, "to", end_round, "_summary.txt"), "w")
cat("Summary of results...\n\n", file=outfile)
cat(paste0("Time taken = ", round(time_taken[3], digits=1), "s\n\n"), file=outfile)
cat(paste0("Sampling iterations = ", iter_sampling, "\n"), file=outfile)
cat(paste0("Warmup iterations = ", iter_warmup, "\n\n"), file=outfile)
cat(paste0("Max Rhat = ", round(max(results$rhat, na.rm=TRUE), digits=2), "\n"), file=outfile) #TODO: Fix error that requires na.rm
cat(paste0("Min ESS = ", round(min(results$ess_bulk, na.rm=TRUE), digits=2), "\n\n"), file=outfile)#TODO: Fix error that requires na.rm

sig_mean = signif(results %>% filter(variable=="alpha") %>% pull(mean), digits=2)
sig_lower = signif(results %>% filter(variable=="alpha") %>% pull(lower), digits=2)
sig_upper = signif(results %>% filter(variable=="alpha") %>% pull(upper), digits=2)
cat(paste0("sig2gp = ", sig_mean, " (", sig_lower, ", ", sig_upper, ")\n"), file=outfile)

ell_mean = signif(results %>% filter(variable=="ell") %>% pull(mean), digits=2)
ell_lower = signif(results %>% filter(variable=="ell") %>% pull(lower), digits=2)
ell_upper = signif(results %>% filter(variable=="ell") %>% pull(upper), digits=2)
cat(paste0("ell = ", ell_mean, " (", ell_lower, ", ", ell_upper, ")\n"), file=outfile)

rho_mean = 1e4*signif(results %>% filter(variable=="rho") %>% pull(mean), digits=2)
rho_lower = 1e4*signif(results %>% filter(variable=="rho") %>% pull(lower), digits=2)
rho_upper = 1e4*signif(results %>% filter(variable=="rho") %>% pull(upper), digits=2)
cat(paste0("rho = ", rho_mean, " (", rho_lower, ", ", rho_upper, ")\n\n"), file=outfile)

cat(paste0("Coverage of nPos = ", round(100*coverage, digits=2), "%\n"), file=outfile)
cat(paste0("Mean width of credible intervals on Pt = ", round(meanwidth, digits=2), "\n"), file=outfile)

close(outfile)

