
rm(list=ls())

source("OtherMethods/Abbott/AbbottMethod.R")

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

init_fn = function() { list(ell=10, alpha=0.07, rho=0.0002) } # Prevent issues with warmup

for (ii in seq(1, nrow(combinations))) {
  
  print(paste0("Simulation: ", combinations$simulation[ii], ", Iteration: ", combinations$iter[ii]))
  
  # Process data for fitting
  df_fit = df %>%
    filter(simulation==combinations$simulation[ii], iter==combinations$iter[ii], variable %in% c("rt", "Pt", "nPos", "nSamples")) %>%
    pivot_wider(names_from=variable, values_from=value)
  
  # Fit model
  counter = 1
  write(sprintf("ii=%d, simulation=%s, iteration=%d, counter=%d", ii, combinations$simulation[ii], combinations$iter[ii], counter), file="tmp.txt", append=TRUE)
  results = runAbbottMethod(df_fit, inc2prevmethod="Binny", init_fn=init_fn)
  check = (max(results$rhat, na.rm=TRUE) < 1.05)
  
  # Fit again until it works, attempt 20 times
  while (!check & counter < 20) {
    
    counter = counter + 1
    write(sprintf("ii=%d, simulation=%s, iteration=%d, counter=%d", ii, combinations$simulation[ii], combinations$iter[ii], counter), file="tmp.txt", append=TRUE)
    
    results = runAbbottMethod(df_fit, inc2prevmethod="Binny", init_fn=init_fn)
    
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
  
  write.csv(results_all, "paper/outputs/s5-pcrcurves/simulations_abbott_binnypcr.csv", row.names=FALSE)
  
}
