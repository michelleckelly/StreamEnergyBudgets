# Function name: request_NEON
# Function purpose: pull data products necessary for two-station metabolism 
# modeling at NEON sites
#
# Arguments:
# NEONsites     character string specifying 4-letter NEON site code to request 
#               data from (ex. `"HOPB"`). Can be more than one site (ex. 
#               `c("HOPB", "BLDE")` but be warned data pull will take longer)
# startdate     YYYY-MM character string defining start month and year for 
#               data request
# enddate       YYY-MM character string defining end month and year for 
#               data request

# For troubleshooting:
# NEONsites <- c("HOPB")
# startdate <- "2018-01"
# enddate <- "2018-12"

request_NEON <- function(NEONsites){
  
  #### Load libraries ########################################################
  # NOTE: remove and add as importFrom when reformatting function into package
  library(tidyverse)
  library(neonUtilities) # Built under R version 3.6.2
  # Since neonUtilities is an old build, we need to tell R to ignore warnings
  # thrown when downloading the NEON-water-quality localPressureDO package, 
  # otherwise the install will fail
  #Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
  #devtools::install_github(repo = "NEONScience/NEON-water-quality/localPressureDO")
  library(localPressureDO)
  library(streamMetabolizer)
  #devtools::install_github(repo = "michelleckelly/NEON-reaeration/reaRate")
  library(reaRate)
  
  #### Input parameters ######################################################
  # Define parameters of interest necessary for metabolism modeling
  params <- c("DP1.20288.001", "DP1.20053.001", "DP1.00024.001", 
              "DP1.20033.001", "DP4.00130.001")
  names(params) <- c("WaterQual", "Temp", "PAR",
                     "NO3", "Discharge")
  
  #### Pull NEON data from api ##################################################
  for (i in seq_along(params)){
    # Set dpID to i-th parameter
    dpID <- params[i]
    
    # Pull NEON data for parameter, saving the data to the variable name from
    # names(params)
    # Pull basic only to save time
    assign(names(dpID), 
           value = neonUtilities::loadByProduct(dpID = dpID, site = NEONsites, 
                                                startdate = startdate, enddate = enddate, 
                                                package = "basic", check.size = F))
  }
  
  # Pull reaeration data from full period of record
  Reaeration <- neonUtilities::loadByProduct(dpID = "DP1.20190.001", 
                                             site = NEONsites, 
                                             package = "expanded", 
                                             check.size = F)
  FieldDischarge <- neonUtilities::loadByProduct(dpID = "DP1.20048.001", 
                                                 site = NEONsites, 
                                                 package = "expanded", 
                                                 check.size = F)
  
  #### Correct dissolved oxygen percent saturation ##############################
  # Pull DO-related products for sites of interest
  DOdata <- localPressureDO::getAndFormatData(siteName = NEONsites, 
                                              startDate = startdate, 
                                              endDate = enddate)
  # Calculate DO percent saturation at local using the Benson-Krause equation,
  # Which is the same method as used by the USGS since 2011
  DOcalcd <- localPressureDO::calcBK_eq(DOdata)
  # Trim DOcalcd dataset
  DOcalcd <- 
    DOcalcd %>%
    select(horizontalPosition, startDateTime, siteID, 
           dissolvedOxygenSatCorrected)
  # Remove DO percent saturation from WaterQual
  # Merge corrected DO percent saturation with WaterQual dataset
  WaterQual$waq_instantaneous <- 
    right_join(WaterQual$waq_instantaneous, DOcalcd, 
               by = c("siteID", "startDateTime", "horizontalPosition"))
  
  #### Format and merge dataframes ##############################################
  # S1 = upstream, S2 = downstream
  ### Nitrate ###
  NO3_data <- 
    NO3$NSW_15_minute %>%
    select(siteID, startDateTime, surfWaterNitrateMean) %>%
    mutate(startDateTime = lubridate::with_tz(startDateTime, tz = "UTC")) %>%
    rename(DateTime_UTC = startDateTime, Nitrate_uMolL = surfWaterNitrateMean)
  
  ### Water quality ###
  WQ_data <- 
    WaterQual$waq_instantaneous %>%
    select(siteID, horizontalPosition, startDateTime, specificConductance, 
           dissolvedOxygen, dissolvedOxygenSatCorrected, pH, chlorophyll, 
           turbidity) %>%
    mutate(horizontalPosition = recode(horizontalPosition, "101" = "S1", 
                                       "102" = "S2"),
           startDateTime = lubridate::with_tz(startDateTime, tz = "UTC")) %>%
    rename(DateTime_UTC = startDateTime, DO_mgL = dissolvedOxygen, 
           DOsat_pct = dissolvedOxygenSatCorrected, 
           chlorophyll_ugL = chlorophyll, turbidity_NTU = turbidity,
           specificConductance_uScm = specificConductance) %>%
    filter(lubridate::minute(DateTime_UTC) %in% c(0, 15, 30, 45))
  
  ### Temp ###
  Temp_data <-
    Temp$TSW_5min %>%
    select(siteID, horizontalPosition, startDateTime, surfWaterTempMean) %>%
    mutate(horizontalPosition = recode(horizontalPosition, "101" = "S1", 
                                       "102" = "S2"),
           startDateTime = lubridate::with_tz(startDateTime, tz = "UTC")) %>%
    rename(DateTime_UTC = startDateTime, WaterTemp_C = surfWaterTempMean) %>%
    filter(lubridate::minute(DateTime_UTC) %in% c(0, 15, 30, 45))
  
  ### Air pressure ###
  AirPres_data <- 
    BP_1min %>%
    select(siteID, startDateTime, staPresMean) %>%
    mutate(startDateTime = lubridate::with_tz(startDateTime, tz = "UTC")) %>%
    rename(DateTime_UTC = startDateTime, AirPres_kPa = staPresMean) %>%
    filter(lubridate::minute(DateTime_UTC) %in% c(0, 15, 30, 45))
  
  ### Discharge & Depth ###
  Discharge_data <-
    Discharge$csd_continuousDischarge %>%
    select(siteID, endDate, calibratedPressure, equivalentStage, maxpostDischarge) %>%
    mutate(endDate = lubridate::with_tz(endDate, tz = "UTC"),
           maxpostDischarge = maxpostDischarge / 1000) %>% # convert discharge, measured as L/s, to m3/s
    rename(DateTime_UTC = endDate, WaterPres_kPa = calibratedPressure,
           Depth_m = equivalentStage, Discharge_m3s = maxpostDischarge) %>%
    filter(lubridate::minute(DateTime_UTC) %in% c(0, 15, 30, 45))
  
  ### Light ###
  PAR_data <- 
    PAR$PARPAR_1min %>%
    select(siteID, startDateTime, PARMean) %>%
    mutate(startDateTime = lubridate::with_tz(startDateTime, tz = "UTC")) %>% 
    rename(DateTime_UTC = startDateTime, Light_PAR = PARMean) %>%
    filter(lubridate::minute(DateTime_UTC) %in% c(0, 15, 30, 45))
  
  ### Merge all into data dataframe ###
  data <- full_join(WQ_data, NO3_data, 
                    by = c("siteID","DateTime_UTC"))
  data <- left_join(data, Temp_data, 
                    by = c("siteID","horizontalPosition","DateTime_UTC"))
  data <- left_join(data, AirPres_data,
                    by = c("siteID","DateTime_UTC"))
  data <- left_join(data, Discharge_data,
                    by = c("siteID","DateTime_UTC"))
  data <- left_join(data, PAR_data,
                    by = c("siteID","DateTime_UTC"))
  
  ### Convert datetime to chron solartime ###
  # Grab longitude of currently in use sensor stations
  sensorPos <- 
    WaterQual$sensor_positions_20288[end == ""] %>%
    select(siteID, HOR.VER, referenceLongitude) %>%
    mutate(HOR.VER = str_extract(HOR.VER, regex("^\\d{3}")),
           HOR.VER = recode(HOR.VER, "101" = "S1", "102" = "S2")) %>%
    rename(horizontalPosition = HOR.VER)
  
  # Add longitude for use in solartime conversion
  data <- right_join(data, sensorPos, by = c("siteID","horizontalPosition"))
  
  # Convert from UTC to solar time 
  data$solarTime <- streamMetabolizer::convert_UTC_to_solartime(data$DateTime_UTC, 
                                                                longitude = data$referenceLongitude)
  
  #################### Reaeration Rate (K) Calculations #########################
  # Format reaeration data product
  Reaeration_data <- 
    def.format.reaeration(rea_backgroundFieldCondData = Reaeration$rea_backgroundFieldCondData,
                          rea_backgroundFieldSaltData = Reaeration$rea_backgroundFieldSaltData,
                          rea_fieldData = Reaeration$rea_fieldData,
                          rea_plateauMeasurementFieldData = Reaeration$rea_plateauMeasurementFieldData,
                          rea_externalLabDataSalt = Reaeration$rea_externalLabDataSalt,
                          rea_externalLabDataGas = Reaeration$rea_externalLabDataGas,
                          rea_widthFieldData = Reaeration$rea_widthFieldData,
                          dsc_fieldData = FieldDischarge$dsc_fieldData,
                          dsc_individualFieldData = FieldDischarge$dsc_individualFieldData)
  
  k600_expanded <- def.calc.reaeration(inputFile = Reaeration_data,
                                       loggerData = Reaeration$rea_conductivityFieldData,
                                       namedLocation = "namedLocation",
                                       injectionTypeName = "injectionType",
                                       eventID = "eventID",
                                       stationToInjectionDistance = "stationToInjectionDistance",
                                       plateauGasConc = "plateauGasConc",
                                       corrPlatSaltConc = "corrPlatSaltConc",
                                       hoboSampleID = "hoboSampleID",
                                       discharge = "fieldDischarge",
                                       waterTemp = "waterTemp",
                                       wettedWidth = "wettedWidth",
                                       plot = TRUE,
                                       savePlotPath = NULL,
                                       processingInfo = NULL)
  
  k600 <- k600_expanded$outputDF
  k600_clean <- k600[k600$k600 > 0 & k600$travelTime > 0,]
  
  # Linear model of Q vs K600
  lmk600 <- lm(k600 ~ meanQ, data = k600_clean)
  
  ################### Output data to user #######################################
  output <- list(data = data,
                 k600_clean = k600_clean,
                 k600_fit = lmk600,
                 k600_expanded = k600_expanded)
  
  return(output)
}
