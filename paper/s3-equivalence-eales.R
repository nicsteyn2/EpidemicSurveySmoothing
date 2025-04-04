
rm(list=ls())

source("OtherMethods/Eales/EalesMethod.R")

# ---------------------------- Fit splines-based models -----------------------------

realised_days_per_knot = 4.95121951219512 # Just define this for saving .txt file
# realised_days_per_knot = 5

# Set options (edit these)
method = "betabinom"
target_knot_spacing = 5

# Fit original model, target_dist = 5
init_fn = function() { list(sigma=0.05, rho=0.0002) }

df = read.csv("data/reactdata.csv") %>%
  filter(round >= 14 & round <= 19) %>%
  select(round, date, nSamples, nPos) %>%
  mutate(date=as.Date(date)) %>%
  arrange(date)

start = proc.time()
out = runEalesMethod(df, method=method, target_dist_between_knots=target_knot_spacing, init_fn=init_fn)
end = proc.time()
time_taken = end - start
time_taken[[3]]

results = out[[1]]
summary = out[[2]]

write.csv(results, sprintf("paper/outputs/s3-equivalence/eales-%s-targetdist%d-results.csv", method, target_knot_spacing), row.names=FALSE)
write.csv(summary, sprintf("paper/outputs/s3-equivalence/eales-%s-targetdist%d-summary.csv", method, target_knot_spacing), row.names=FALSE)

# Give duration in m and s
write(sprintf("\n\nEales method, %s, target_dist = %d...\n", method, target_knot_spacing), file="paper/outputs/s3-equivalence/s3-equivalence-eales-results.txt", append=TRUE)
write(sprintf("Fitting took %dm %ds", floor(time_taken[[3]]/60), round(time_taken[[3]] %% 60)), file="paper/outputs/s3-equivalence/s3-equivalence-eales-results.txt", append=TRUE)
write(sprintf("Max Rhat: %f", max(summary$rhat)), file="paper/outputs/s3-equivalence/s3-equivalence-eales-results.txt", append=TRUE)
write(sprintf("Min ESS: %f", min(summary$ess_bulk)), file="paper/outputs/s3-equivalence/s3-equivalence-eales-results.txt", append=TRUE)
write("\n", file="paper/outputs/s3-equivalence/s3-equivalence-eales-results.txt", append=TRUE)
write(capture.output(summary %>% filter(variable %in% c("sigma", "rho"))), file="paper/outputs/s3-equivalence/s3-equivalence-eales-results.txt", append=TRUE)
write("\n", file="paper/outputs/s3-equivalence/s3-equivalence-eales-results.txt", append=TRUE)
write(capture.output(summary %>% filter(variable == "sigma") %>% mutate(mean=mean/sqrt(realised_days_per_knot),
                                                                        lower=lower/sqrt(realised_days_per_knot),
                                                                        upper=upper/sqrt(realised_days_per_knot))),
      file="paper/outputs/s3-equivalence/s3-equivalence-eales-results.txt", append=TRUE)
write("\n", file="paper/outputs/s3-equivalence/s3-equivalence-eales-results.txt", append=TRUE)




# ---------------------------- Fit no-splines model -----------------------------

rm(list=ls())

init_fn = function() { list(sigma=0.05, rho=0.0002) }

source("OtherMethods/Eales/EalesMethodNoSpline.R")

df = read.csv("data/reactdata.csv") %>%
  filter(round >= 14 & round <= 19) %>%
  select(round, date, nSamples, nPos) %>%
  mutate(date=as.Date(date)) %>%
  arrange(date) %>%
  complete(date = seq.Date(min(date)-3, max(date)+3, by="day"), fill=list(nSamples=0, nPos=0))

start = proc.time()
results = runEalesMethodNoSpline(df, init_fn=init_fn)
end = proc.time()
time_taken = end - start
time_taken[[3]]

write.csv(results, "paper/outputs/s3-equivalence/eales-nosplines-results.csv", row.names=FALSE)

# Give duration in m and s
write("\n\nEales method (no splines)...\n", file="paper/outputs/s3-equivalence/s3-equivalence-eales-results.txt", append=TRUE)
write(sprintf("Fitting took %dm %ds", floor(time_taken[[3]]/60), round(time_taken[[3]] %% 60)), file="paper/outputs/s3-equivalence/s3-equivalence-eales-results.txt", append=TRUE)
write(sprintf("Max Rhat: %f", max(results$rhat, na.rm=TRUE)), file="paper/outputs/s3-equivalence/s3-equivalence-eales-results.txt", append=TRUE)
write(sprintf("Min ESS: %f", min(results$ess_bulk, na.rm=TRUE)), file="paper/outputs/s3-equivalence/s3-equivalence-eales-results.txt", append=TRUE)
write("\n", file="paper/outputs/s3-equivalence/s3-equivalence-eales-results.txt", append=TRUE)
write(capture.output(results %>% filter(variable %in% c("sigma", "rho"))), file="paper/outputs/s3-equivalence/s3-equivalence-eales-results.txt", append=TRUE)
write("\n", file="paper/outputs/s3-equivalence/s3-equivalence-eales-results.txt", append=TRUE)
