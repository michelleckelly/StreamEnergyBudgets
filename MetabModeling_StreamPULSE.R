# Modeling metabolism for sites with model records not yet
# compiled on StreamPULSE
# NOTE: You will have to roll back your R version to 3.6.3 or earlier before 
# proceeding, as streamMetabolizer isn't yet built for R 4.0.0

# Load libraries
library(tidyverse)
library(StreamPULSE)
library(streamMetabolizer)

# Set parallelization
options(mc.cores = parallel::detectCores())

# Choose model type and framework
model_type <- "bayes"
model_name <- "streamMetabolizer"

# Request variables
vars <- c('DO_mgL','DOsat_pct','satDO_mgL','WaterPres_kPa',
          'Depth_m','WaterTemp_C','Light_PAR','AirPres_kPa',
          'Discharge_m3s', "ChlorophyllA_ugL", "Nitrate_mgL")

# Define list of sites for metabolism modeling
siteList <- c("AZ_WB")
# Define dates, YYYY-MM-DD enter NULL for both if wish is to pull period of record
startDate <- "2019-05-01"
endDate <- "2019-06-01"

# This loop will iterate through the list of sites, pulling data from 
# StreamPULSE and modeling metabolism for the period of record. Metabolism 
# models for each site will be saved to RData files within /model_outputs/
for (i in 1:length(siteList)){
  # Define site
  site <- siteList[i]
  print(site)
  flush.console()
  # Start timer for run
  startTime <- proc.time()
  
  ######################## Pull data from StreamPULSE ########################
  data <- request_data(sitecode = site, variables = vars, 
                       startdate = startDate, enddate = endDate)
  # If error in pulling data, throw warning and move on to next site in list
  if (!exists("data")){
    warning(paste0("Data could not be pulled for site ", site))
    flush.console()
    next
  }
  # Display success message
  message(paste0("Data pulled for site ", site, ". \nElapsed time: ", 
                 (proc.time()[3]-startTime[3])/60, " minutes"))
  flush.console()
  
  ######################## Format data for modeling ########################
  data_prep <- prep_metabolism(d = data, type = model_type, model = model_name,
                               rm_flagged = list("Bad Data", "Questionable"),
                               interval = "15 min",
                               estimate_areal_depth = TRUE,
                               retrieve_air_pres = TRUE)
  # Display success message
  message(paste0("Data formatted for site ", site, ". \nElapsed time: ", 
                 (proc.time()[3]-startTime[3])/60, " minutes"))
  flush.console()
  
  ############################ Model metabolism ############################
  # If no discharge data, set poolK_600='none'
  if (grepl("Discharge_m3s", data_prep$specs$missing_variables)){
    message(paste0("Discharge data missing for site ", site, 
                   ". Will use fit_metabolism with poolK_600='none'."))
    flush.console()
    model_fit <- fit_metabolism(data_prep, pool_K600 = "none")
  } else {
    # Otherwise, execute as normal
    model_fit <- fit_metabolism(data_prep)
  }
  # Display success message
  message(paste0("Metabolism modeled for site ", site, ". \nElapsed time: ", 
                 (proc.time()[3]-startTime[3])/60, " minutes"))
  flush.console()
}

# Visualize results
library(ggplot2)

ggplot(data = model_fit$predictions) +
  geom_point(aes(x = date, y = GPP)) +
  geom_ribbon(aes(ymin = GPP.lower, ymax = GPP.upper, x = date))
ggplot(data = model_fit$predictions) +
  geom_point(aes(x = date, y = ER)) +
  geom_ribbon(aes(ymin = ER.lower, ymax = ER.upper, x = date))

ggplot(data = model_fit[["fit"]]@fit[["daily"]]) +
  geom_line(aes(x = date, y = K600_daily_mean)) +
  geom_ribbon(aes(ymin = K600_daily_25pct, ymax = K600_daily_75pct, x = date))