
# This script fits the Abbott method to 10 simulations from each method
# Progress can be tracked by checking the tmp.txt file in the root directory
# Code to produce the Abbott simulations can be found in /OtherMethods/Abbott/abbott_simulator.R

rm(list=ls())

source("OtherMethods/Abbott/AbbottMethod.R")

# Load simulations
df = rbind(
  read.csv("paper/outputs/2-simcompare/simulations_steyn.csv") %>% mutate(simulation="Steyn"),
  read.csv("paper/outputs/2-simcompare/simulations_abbott.csv") %>% mutate(simulation="Abbott"),
  read.csv("paper/outputs/2-simcompare/simulations_eales.csv") %>% mutate(simulation="Eales")
)

# Specify fittinng combinations (10 iterations of each method)
combinations = expand.grid(
  iter = seq(1,10),
  simulation = c("Steyn", "Abbott", "Eales")
)

# Define empty dataframe to store results
results_all = data.frame()

# Reduce MCMC warmup issues by specifying an function of initial values
init_fn = function() { list(ell=10, alpha=0.07, rho=0.0002) }

st = Sys.time()
for (ii in seq(1, nrow(combinations))) {
  
  print(paste0("Simulation: ", combinations$simulation[ii], ", Iteration: ", combinations$iter[ii]))

  # Process data for fitting
  df_fit = df %>%
    filter(simulation==combinations$simulation[ii], iter==combinations$iter[ii], variable %in% c("rt", "Pt", "nPos", "nSamples")) %>%
    pivot_wider(names_from=variable, values_from=value)

  # Fit model
  counter = 1
  write(sprintf("ii=%d, simulation=%s, iteration=%d, counter=%d", ii, combinations$simulation[ii], combinations$iter[ii], counter), file="tmp.txt", append=TRUE)
  results = runAbbottMethod(df_fit, init_fn=init_fn)
  check = (max(results$rhat, na.rm=TRUE) < 1.05)
  
  # Fit again until it works, attempt 20 times
  while (!check & counter < 20) {
    
    counter = counter + 1
    write(sprintf("ii=%d, simulation=%s, iteration=%d, counter=%d", ii, combinations$simulation[ii], combinations$iter[ii], counter), file="tmp.txt", append=TRUE)
    
    results = runAbbottMethod(df_fit, init_fn=init_fn)
    
    check = (max(results$rhat, na.rm=TRUE) < 1.05)
    
  }
  
  # Tidy results
  df_true = df %>%
    filter(simulation==combinations$simulation[ii], iter==combinations$iter[ii]) %>%
    select(t, value, variable) %>%
    rename(TrueValue = value)
  
  results_tidy = results %>%
    mutate(iter=combinations$iter[ii], simulation=combinations$simulation[ii], method="Abbott", fitAttempts = counter) %>%
    filter(variable %in% c("ell", "alpha", "rho", "init_inc", "init_growth", "It", "Pt", "Rt", "rtinc", "rt", "nPos")) %>%
    left_join(df_true, by=c("t", "variable"))
  
  results_all = rbind(results_all, results_tidy)
  
  write.csv(results_all, "paper/outputs/2-simcompare/simcompare_abbott.csv", row.names=FALSE)
  
}
et = Sys.time()
time_taken = et - st
time_taken
