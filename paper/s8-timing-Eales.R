
rm(list=ls())

source("OtherMethods/Eales/EalesMethod.R")

# Set global options
iter_warmup = 200
iter_sampling = 300
max_treedepth = 16
init_fn = function() { list(sigma=0.1, rho=0.0002) } # Prevent issues with warmup

# Set options
start_round = 1
end_round = 19
n = 3

fname = paste0("paper/outputs/sx-timing/eales_R",start_round, "to", end_round, "_warmup", iter_warmup, "_sampling", iter_sampling, "_treedepth", max_treedepth, ".csv")


# Function to load data
load_data = function(start_round, end_round) {
  df = read.csv("data/reactdata.csv") %>%
    filter(round >= start_round & round <= end_round) %>%
    select(round, date, nSamples, nPos) %>%
    mutate(date=as.Date(date)) %>%
    arrange(date)
  return(df)
}


# Function to run Eales method
run_eales_method = function(df) {
  
  start = proc.time()
  results = runEalesMethod(df,
                           method="overdispersed",
                           iter_warmup=iter_warmup, iter_sampling=iter_sampling, max_treedepth = max_treedepth,
                           calculateRt=TRUE,
                           init_fn = init_fn)
  end = proc.time() 
  time_taken = end - start
  
  states = results[[1]]
  summary = results[[2]]
  
  return(list(states, summary, time_taken[[3]]))
  
}


# Function to run multiple methods
run_eales_method_ntimes = function(df, n, fname) {
  
  # Dataframe to store the results
  df_out = data.frame()
  
  # Start a counter
  n_attempts = 1
  n_successful = 0
  
  # Loop over the number of times to run the method, up to a maximum of 50 attempts
  while ((n_successful < n) & (n_attempts <= 50)) {
    
    # Run the method
    results = run_eales_method(df)
    states = results[[1]]
    summary = results[[2]]
    time_taken = results[[3]]
    
    # Check if successful
    is_success = (max(summary$rhat) <= 1.05) & (min(summary$ess_bulk) >= 100)
    
    # Store key results
    df_temp = data.frame(
      start_round = min(df$round),
      end_round = max(df$round),
      iter = n_attempts,
      success = is_success,
      maxrhat = max(summary$rhat),
      miness = min(summary$ess_bulk),
      time_taken = time_taken,
      num_obs = nrow(df),
      num_days = as.numeric(max(df$date) - min(df$date)) + 1,
      iter_sampling = iter_sampling,
      iter_warmup = iter_warmup,
      max_treedepth = max_treedepth
    )
    
    # Append onto the output dataframe and save
    df_out = rbind(df_out, df_temp)
    write.csv(df_out, fname, row.names=FALSE)
    
    # Increment the counter
    n_attempts = n_attempts + 1
    n_successful = n_successful + is_success
    
  }
  
  return(df_out)
  
}


df = load_data(start_round, end_round)
test = run_eales_method_ntimes(df, n, fname)
