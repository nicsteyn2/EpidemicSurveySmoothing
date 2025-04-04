
# This script fits the Abbott method to 10 simulations from each method
# Progress can be tracked by checking the tmp.txt file in the root directory
# Code to produce the Eales simulations can be found in /OtherMethods/Eales/eales_simulator.R

rm(list=ls())

source("OtherMethods/Eales/EalesMethod.R")

# Load data
df = rbind(
  read.csv("paper/outputs/2-simcompare/simulations_steyn.csv") %>% mutate(simulation="Steyn"),
  read.csv("paper/outputs/2-simcompare/simulations_abbott.csv") %>% mutate(simulation="Abbott"),
  read.csv("paper/outputs/2-simcompare/simulations_eales.csv") %>% mutate(simulation="Eales")
)

combinations = expand.grid(
  iter = seq(1,10),
  simulation = c("Steyn", "Abbott", "Eales")
)

results_all = data.frame()
summary_all = data.frame()

init_fn = function() { list(sigma=0.1, rho=0.0002) }

st = proc.time()

for (ii in seq(1, nrow(combinations))) {
  
  print(paste0("Simulation: ", combinations$simulation[ii], ", Iteration: ", combinations$iter[ii]))
  
  # Process data for fitting
  df_fit = df %>%
    filter(simulation==combinations$simulation[ii], iter==combinations$iter[ii], variable %in% c("rt", "Pt", "nPos", "nSamples")) %>%
    pivot_wider(names_from=variable, values_from=value)
  
  # Fit model
  counter = 1
  write(sprintf("ii=%d, simulation=%s, iteration=%d, counter=%d", ii, combinations$simulation[ii], combinations$iter[ii], counter), file="tmp.txt", append=TRUE)
  out = runEalesMethod(df_fit, init_fn = init_fn)
  results = out[[1]]
  summary = out[[2]]
  check = (max(summary$rhat, na.rm=TRUE) < 1.05)
  
  # Fit again until it works, attempt 20 times
  while (!check & counter < 20) {
    
    counter = counter + 1
    write(sprintf("ii=%d, simulation=%s, iteration=%d, counter=%d", ii, combinations$simulation[ii], combinations$iter[ii], counter), file="tmp.txt", append=TRUE)
    
    out = runEalesMethod(df_fit, init_fn = init_fn)
    results = out[[1]]
    summary = out[[2]]
    check = (max(summary$rhat, na.rm=TRUE) < 1.05)
    
  }
  
  # Tidy results
  df_true = df %>%
    filter(simulation==combinations$simulation[ii], iter==combinations$iter[ii]) %>%
    select(t, value, variable) %>%
    rename(TrueValue = value)
  
  results_tidy = results %>%
    mutate(iter=combinations$iter[ii], simulation=combinations$simulation[ii], method="Eales", fitAttempts = counter) %>%
    left_join(df_true, by=c("t", "variable"))
  
  summary_tidy = summary %>%
    mutate(iter=combinations$iter[ii], simulation=combinations$simulation[ii], method="Eales", fitAttempts = counter)
  
  results_all = rbind(results_all, results_tidy)
  summary_all = rbind(summary_all, summary_tidy)
  
  write.csv(results_all, paste0("paper/outputs/2-simcompare/simcompare_eales_results.csv"), row.names=FALSE)
  write.csv(summary_all, paste0("paper/outputs/2-simcompare/simcompare_eales_summary.csv"), row.names=FALSE)
  
}

en = proc.time()
time_taken = en - st
time_taken[[3]]