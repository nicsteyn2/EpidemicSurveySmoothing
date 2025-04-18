
rm(list=ls())

source("OtherMethods/Abbott/AbbottMethod.R")

# Set global options
iter_warmup = 200
iter_sampling = 300

# Set options
start_round = 1
end_round = 2
n = 10

fname = paste0("paper/outputs/s8-timing/abbott_R",start_round, "to", end_round, "_warmup", iter_warmup, "_sampling", iter_sampling, ".csv")


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
run_abbott_method = function(df) {
  
  start = proc.time()
  results = runAbbottMethod(df, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                            gentime_mean = c(mean=6.3611, mean_sd=0.001), gentime_sd = c(sd=4.2035, sd_sd=0.001)) %>%
    filter(variable != "nPos") # Removing predictive results as they can have NA Rhat values
  end = proc.time() 
  time_taken = end - start
  
  return(list(results, time_taken[[3]]))
  
}


# Function to run multiple methods
run_abbott_method_ntimes = function(df, n, fname) {
  
  # Dataframe to store the results
  df_out = data.frame()
  
  # Start a counter
  n_attempts = 1
  n_successful = 0
  
  # Loop over the number of times to run the method, up to a maximum of 50 attempts
  while ((n_successful < n) & (n_attempts <= 50)) {
    
    # Run the method
    results = run_abbott_method(df)
    df_results = results[[1]]
    time_taken = results[[2]]
    
    # Check if successful
    maxrhat = max(df_results$rhat)
    miness = min(df_results$ess_bulk)
    
    is_success = (maxrhat <= 1.05) & (miness >= 100)
    # Check if NA
    if (is.na(maxrhat) | is.na(miness)) {
      is_success = FALSE
    }
    
    # Store key results
    df_temp = data.frame(
      start_round = min(df$round),
      end_round = max(df$round),
      iter = n_attempts,
      success = is_success,
      maxrhat = maxrhat,
      miness = miness,
      time_taken = time_taken,
      num_obs = nrow(df),
      num_days = as.numeric(max(df$date) - min(df$date)) + 1,
      iter_sampling = iter_sampling,
      iter_warmup = iter_warmup
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
run_abbott_method_ntimes(df, n, fname)
