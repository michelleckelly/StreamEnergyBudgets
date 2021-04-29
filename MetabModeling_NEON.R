###################### Load Data from NEON ####################################
# Load libraries
library(neonUtilities) # Built under R version 3.6.2
# Since neonUtilities is an old build, we need to tell R to ignore warnings
# thrown when downloading the NEON-water-quality localPressureDO package, 
# otherwise the install will fail
#Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
#devtools::install_github(repo = "NEONScience/NEON-water-quality/localPressureDO")
library(localPressureDO)
library(streamMetabolizer)

## See if you can pull in NEON discharge data with NEON utilities package
# Data to be pulled from NEON:
# Depth_m, WaterTemp_C, Nitrate_mgL, ChlorophyllA_ugL, DO_mgL, DOsat_pct, 
# satDO_mgL, WaterPres_kPa, Light_PAR, AirPres_kPa, Discharge_m3s 

# Streampulse variable    NEON product ID     NEON name
# Depth_m                 DP1.20016.001       Elevation of surface water
# WaterTemp_C             DP1.20053.001       Temperature (PRT) in surface water
# DO_mgL                  DP1.20288.001       Water quality
# DOsat_pct
# satDO_mgL
# WaterPres_kPa
# Light_PAR               DP1.00024.001       Photosynthetically active radiation (PAR)
# AirPres_kPa
# Discharge_m3s
# Nitrate_mgL
# ChlorophyllA_ugL      DP1.20288.001       Water quality

# Ahh. NEON collects field-based surveys of discharge, intermittant

## Pull raw NEON dataset -----------------------------------------------------
# Define sites of interest
NEONsites <- c("ARIK")
#NEONsites <- c("ARIK", "REDB", "PRIN", "MCDI", "KING", "HOPB", "GUIL", "CUPE", 
#               "SYCA", "CARI", "BLDE", "MCRA", "MART")

# Define parameters of interest
params <- c("DP1.20288.001", "DP1.20053.001", #"DP1.20267.001", 
            "DP1.20033.001")
names(params) <- c("WaterQual", "Temp", #"GaugeHt", 
                   "NO3")

# Define start and end dates
startdate <- "2019-01"
enddate <- "2019-12"

# Pull data from NEON
for (i in seq_along(params)){
  # Set dpID to i-th parameter
  dpID <- params[i]
  
  # Pull NEON data for parameter, saving the data to the variable name from
  # names(params)
  assign(names(dpID), 
         value = loadByProduct(dpID = dpID, site = NEONsites, 
                               startdate = startdate, enddate = enddate, 
                               package = "basic"))
}

## Correct dissolved oxygen percent saturation -------------------------------
# Pull DO-related products for sites of interest
DOdata <- getAndFormatData(siteName = NEONsites, startDate = startdate, 
                           endDate = enddate)
# Calculate DO percent saturation at local using the Benson-Krause equation,
# Which is the same method as used by the USGS since 2011
DOcalcd <- calcBK_eq(DOdata)
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

# Save raw NEON datafiles, with corrected DO sat, to external harddrive -------
#saveRDS(WaterQual, 
#        file = "/Volumes/Ext001_MCK/Datasets/CAREER_NEON/WaterQual.RData")
#saveRDS(Temp, file = "/Volumes/Ext001_MCK/Datasets/CAREER_NEON/Temp.RData")
#saveRDS(NO3, file = "/Volumes/Ext001_MCK/Datasets/CAREER_NEON/NO3.RData")

# Load raw NEON datafiles from external harddrive
#WaterQual <- 
#  readRDS(file = "/Volumes/Ext001_MCK/Datasets/CAREER_NEON/WaterQual.RData")
#Temp <- readRDS(file = "/Volumes/Ext001_MCK/Datasets/CAREER_NEON/Temp.RData")
#NO3 <- readRDS(file = "/Volumes/Ext001_MCK/Datasets/CAREER_NEON/NO3.RData")

# Create metadata tibble for site information
# May be unneeded - consult the CAREER site list metadata file
WQ_meta <- 
  WaterQual$sensor_positions_20288 %>%
  select(siteID, referenceDescription, referenceStart, referenceEnd,
         referenceLatitude, referenceLongitude, referenceElevation) %>%
  # Filter to active stations
  filter(referenceEnd == "")
# Add horizontal position column
WQ_meta$horizontalPosition <- c("S1", "S2", "S2")

# S1 = upstream, S2 = downstream
WQ_data <- 
  WaterQual$waq_instantaneous %>%
  select(siteID, horizontalPosition, startDateTime, specificConductance, 
         dissolvedOxygen, dissolvedOxygenSatCorrected, pH, chlorophyll, 
         turbidity) %>%
  mutate(horizontalPosition = recode(horizontalPosition, "101" = "S1", 
                                     "102" = "S2"))
# Add longitude to WQ_data for use in solartime conversion
WQ_data <- 
  right_join(WQ_data, WQ_meta %>% select(horizontalPosition, referenceLongitude), 
             by = "horizontalPosition")

# Convert datetime to chron solartime -----------------------------------------
# Convert to UTC timezone
WQ_data$startDateTime <- with_tz(WQ_data$startDateTime, tz = "UTC")
# Convert from UTC to solar time using streamMetabolizer
WQ_data$solarTime <- 
  convert_UTC_to_solartime(WQ_data$startDateTime, 
                           longitude = WQ_data$referenceLongitude)

# Save trimmed WQ_data to local folder ----------------------------------------
saveRDS(WQ_data, file = "./data/NEON_formatted/NEON_WaterQual.RData")
saveRDS(WQ_meta, file = "./data/NEON_formatted/NEON_WaterQualMeta.RData")

# Stack to format for 2-station metabolism, then save formatted stacked file

#################### Reaeration Rate (K) Calculations #########################
# Load libraries 
#devtools::install_github(repo = "NEONScience/NEON-reaeration/reaRate")
library(reaRate)

# Download NEON reaeration data product
rea_data <- def.format.reaeration(site = NEONsites)
