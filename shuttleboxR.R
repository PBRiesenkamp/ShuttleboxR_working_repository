# Script for R package shuttleboxR
#
# 25/026/2024
#
# Shaun Killen and Emil Christensen
#


# TO DO:  
  # Edit inspect function for non shuttlesoft data to check it contains the right data
  # Add acclimation exclusion to each function
  # Remove dependencies as much as possible
    # Check if segmented regression can be done without segmented package

  # Make overall diagnostic function
    # Check diagnostic functions
    # Add option to highlight specific individuals in histograms

rm(list=ls())

# renv:: restore checks whether the packages you are using are the same versions as in the script
# if collaborators update a library and then edit the script, renv::restore will pick up on that
renv::restore()

# Confirm that R is looking in the right place
getwd()


#### File definitions ####
# proj_data <- read.csv("project_database.csv") # This should be a summary data file with data from numerous fish (e.g. all fish in a study)

file <- "Sys1_20240613_parcatus_rr.txt"
# data<-read.csv("data.csv")
# data$date <- as.Date(data$date, format = "%d/%m/%Y")

metadata<-read.csv("additional_data.csv", na.strings = c("NaN", ""))


#### Read shuttlesoft files function ####

# Collect existing information into one data file:
# date, time, zone, INCR_T, DECR_T, coordinates/distance, pixel ratio, trial_start, mass, a_value, b_value, acclimation temp

read.shuttlesoft <- function (file, metadata, multidat = F){
  # Load .txt data file that shuttlesoft produces into R
  # Don't load headers, because the additional info at the top of the .txt file will make a confusing dataframe
  # Rename data columns
  data<-read.delim(file,
                   header = F,
                   col.names = (c("time", "zone", "core_T", "Tpref_loligo", "INCR_T",
                                  "DECR_T", "x_pos", "y_pos", "velocity", "distance", "time_in_INCR", "time_in_DECR",
                                  "delta_T", "dyn_hysteresis", "stat_T_INCR", "stat_hyst_INCR", "stat_T_DECR",
                                  "stat_hyst_DECR", "k", "max_T", "min_T", "change_rate", "avoidance_upper",
                                  "avoidance_upper_core", "avoidance_lower", "avoidance_lower_core")),

                   na.strings = c("NaN", ""))
  
  # Extract notes, the file created info and the pixel ratio from the dataframe
  data<-data[!is.na(data$time), ]
  notes <- data$zone[data$time == "Notes"]
  filecreated <- as.POSIXct(data$zone[data$time == "File created"], format = "%d/%m/%Y; %H:%M")
  pixel_ratio <- data$zone[data$time == "Pixel Ratio [cm/pix]"]
  
   # Remove the additional info and reset rownames
  data<-data[!is.na(data$INCR_T), ]
  numeric_rows <- !is.na(suppressWarnings(as.numeric(data$INCR_T)))
  data<-data[numeric_rows, ]
  rownames(data)<-NULL 
  
  # Add notes and pixel ratio
  data$notes <- notes
  data$pixel_ratio <- pixel_ratio
  
  # Add the date of the trial by extracting it from the filecreated object
  data$date <- as.Date(filecreated)
  fileID <- gsub(paste0("^", getwd(), "/?"), "", file)
  data$fileID <- fileID
  
  trialstart <-  function(metadata){
    
    if (missing(metadata)) {
      trial_start <- data$time[1]
      warning("Metadata not provided; cannot define 'trial_start'")
    } else {
      
      # If the metadata is present, ensure necessary columns exist in the metadata
      if (!all(c("file_name", "trial_start", "mass", "a_value", "b_value") %in% colnames(metadata))) {
        warning("The metadata does not contain one or more necessary columns: 'file_name', 'trial_start', 'mass', 'a_value', 'b_value'.")
      }
      
      # Find rows with NA in specified columns
      columns_to_check <- c("file_name", "trial_start", "mass", "a_value", "b_value")
      rows_with_na <- apply(metadata[columns_to_check], 1, function(row) any(is.na(row)))
      
      # If there are any rows with NA, throw a warning mentioning file names
      if (any(rows_with_na)) {
        names_with_na <- metadata$file_name[rows_with_na]
        warning (paste("The metadata contains NA values in one or more of these columns: 'file_name', 'trial_start', 'mass', 'a_value', 'b_value'
",paste(names_with_na, collapse = "\n ")))
      }
      
      #Ensure necessary columns are in the right format
      if (all(grepl("^\\d{2}:\\d{2}:\\d{2}$", metadata$trial_start))==F){
        message("Warning: One or more entries in trial_start are not in an hh:mm:ss format")
      }
      
      if(!("trial_start" %in% names(metadata))|is.na(metadata$trial_start[metadata$file_name==fileID])){
        trial_start <- data$time[1]
      }else {
        trial_start<-metadata$trial_start[metadata$file_name==fileID]
      }
    }
    return(trial_start)}
  
  if (multidat) {
    data$trial_start <- suppressWarnings(trialstart (metadata = metadata))} else {
      data$trial_start <- trialstart(metadata = metadata)}
 
  data$a_value <- NA
  data$b_value <- NA
  data$mass <- NA
  data$initial_T <- NA
  data$a_value <- metadata$a_value[metadata$file_name == fileID]
  data$b_value <- metadata$b_value[metadata$file_name == fileID]
  data$mass <- metadata$mass[metadata$file_name == fileID]
  data$initial_T <- metadata$initial_T[metadata$file_name == fileID] 
  
  # # Remove unnecessary columns
  # data = subset(data, select = -c(k, Tpref_loligo, avoidance_upper, avoidance_upper_core, avoidance_lower, avoidance_lower_core))
  
  # # Reorder columns
  # col_order <- c("date", "time", "zone", "INCR_T", "DECR_T", "core_T",  "x_pos", "y_pos", "velocity", "distance",
  #                "time_in_INCR", "time_in_DECR", "max_T", "min_T", "change_rate", "delta_T", "dyn_hysteresis", "stat_T_INCR", "stat_hyst_INCR", "stat_T_DECR",
  #                "stat_hyst_DECR", "trial_start", "a_value", "b_value", "mass", "initial_T", "pixel_ratio", "notes")
  # data <- data[, col_order]
  
  return(data)
}

data <- read.shuttlesoft(file, metadata)

#### Inspect data function ####
inspect <- function(data){
  
  if(!all(c("date", "time", "zone", "INCR_T", "DECR_T") %in% colnames(data))){
    stop("Dataset is missing one or more of the following necessary columns:  'date', 'time', 'zone', 'INCR_T', 'DECR_T'")
  }
  else if(!all(c("x_pos", "y_pos", "velocity", "distance", "time_in_INCR", "time_in_DECR",
            "delta_T", "dyn_hysteresis", "stat_T_INCR", "stat_hyst_INCR", "stat_T_DECR",
            "stat_hyst_DECR", "max_T", "min_T", "change_rate")%in% colnames(data))) {
    warning("Dataset is missing one or more of the following columns:'x_pos', 'y_pos', 'velocity', 'distance', 'time_in_INCR', 'time_in_DECR', 'delta_T', 'dyn_hysteresis', 'stat_T_INCR', 'stat_hyst_INCR', 'stat_T_DECR', 'stat_hyst_DECR', 'max_T', 'min_T', 'change_rate'")
  } else {praise()}
   
  
}
inspect(data)




#### Function to prepare data file for further processing and detect shuttles ####

file_prepare <- function(data) {
  
  # Add time in seconds and hours assuming a sampling frequency of 1 Hz
  data$time_sec <- seq(0, by = 1, length.out = nrow(data))
  data$time_h <- data$time_sec / 3600
  
  if (all(c("x_pos", "y_pos") %in% colnames(data))) {
  
  data$x_pos[data$x_pos == "No object"] <- NA
  data$y_pos[data$y_pos == "No object"] <- NA
  
  data$x_pos <- as.numeric(data$x_pos)
  data$y_pos <- as.numeric(data$y_pos)
  }
  
  # Make sure the date is correct in case the experiment goes on past midnight
  # Extract the second when the experiment passed midnight
  # Make all observations past that date the next day
  if(length(unique(data$date)) == 1){
  
  midnight_observation <- data$time_sec[data$time == "00:00:00"]
  if(length(midnight_observation)>0){
    data$date<-ifelse(data$time_sec>=midnight_observation, 
                      data$date+1, 
                      data$date)
  }
  }
  
  # Create datetime object
  data$datetime<-as.POSIXct(paste(data$date, data$time), format = "%Y-%m-%d %H:%M:%S")
  
  # Add a column that tells you whether the system is in static or dynamic
  if ("delta_T" %in% colnames(data))  data$dyn_stat <- ifelse(is.na(data$delta_T), "static", "dynamic")
  
  # Add the phase of the experiment (e.g. acclimation and trial)
  trial_start_second<-data$time_sec[data$time == data$trial_start]
  data$trial_phase<-ifelse(data$time_sec>=trial_start_second, "trial", "acclimation")
  
  # Detect a shuttle (chamber side change) in new column
  data$shuttle <- 0
  for (i in 2:nrow(data)) {
    if (data$zone[i] != data$zone[i-1]) {
      data$shuttle[i] <- 1
    } else {
      data$shuttle[i] <- 0
    }
  }
  
  # Ensure the shuttle column is numeric
  data$shuttle <- as.numeric(data$shuttle)
  
  # Check if necessary columns exist in the data. If not, return NA for those columns
  
  if (!any(c("x_pos", "y_pos", "distance") %in% colnames(data))) {
    warning("The dataset does not contain a 'distance' column nor 'x_pos' and 'y_pos' to calculate 'distance', adding empty column")
    data$distance <- NA
  }
  
  if (!all(c("max_T", "min_T") %in% colnames(data))) {
    warning("The dataset does not contain 'min_T' and/or 'max_T' columns, adding empty columns")
    data$min_T <- NA
    data$max_T <- NA
  }
  
  data<-type.convert(data, as.is = T)
  data$datetime<-as.POSIXct(data$datetime)
  return(data)
}

data <- file_prepare(data)

#### Function to calculate body core temperature if not already done in Shuttlesoft ####

calc_coreT <- function(data) {
  
  # Initialize the core_T column
  if ("initial_T" %in% colnames(data) && any(!is.na(data$initial_T))) {
    data$core_T <- NA
    data$core_T[1] <- unique(data$initial_T)} else {
    if(!"core_T" %in% colnames(data) || all(is.na(data[["core_T"]])))
    stop("Initial temperature not defined, cannot calculate core body temperature")
    }
  
  # Ensure necessary columns exist
  if (!any(c("a_value", "b_value", "mass", "shuttle") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'a_value', 'b_value', 'mass','shuttle'.")
  }
  
  #Ensure necessary columns are in the right format
  if (!all(is.numeric(c(data$a_value, data$b_value, data$mass, data$core_T, data$shuttle)))){
    stop("One or more necessary variables are not in the right format")
  }
  
  if (any(c(length(unique(data$a_value)),length(unique(data$b_value)), length(unique(data$mass)))>1)){
    stop("One or more of the following parameters have more than 1 value: 'a_value', 'b_value', 'mass' ")
  }

  a_value <- unique(data$a_value)
  b_value <- unique(data$b_value)
  mass <- unique(data$mass)
  
  k <- a_value * mass^b_value
  
  data$ambient_T<-ifelse(data$zone == "INCR", data$INCR_T, data$DECR_T)
  
  # Calculate body core temperature from second to second
  for (i in 2:nrow(data)) {
    if (data$shuttle[i] == 1) {
      data$core_T[i] <- data$core_T[i-1]  # Carry forward the body core temperature at the moment of shuttling
    } else {
      data$core_T[i] <- data$ambient_T[i] + 
        (data$core_T[i-1] - data$ambient_T[i]) * exp(-k * 1 / 60)
    }
  }
  
  return(data)
}

data <- calc_coreT(data)

#### Function to calculate Tpref using either the mean or median of all body core temperature values ####

calc_Tpref <- function(data, method = c("mean", "median", "mode"), exclude_start_minutes = 0, exclude_end_minutes = 0, print_results = T) {
  method <- match.arg(method)
  
  # Convert the time to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Convert the core_T column to numeric
  data$core_T <- as.numeric(data$core_T)
  
  # Calculate Tpref based on the specified method
  if (method == "mean") {
    Tpref <- mean(data$core_T, na.rm = TRUE)
  } else if (method == "median") {
    Tpref <- median(data$core_T, na.rm = TRUE)
  } else if (method == "mode") {
    Tpref <- as.numeric(names(sort(table(data$core_T), decreasing = TRUE))[1])
  }
  
  if (print_results) {
  print(paste("Tpref:", Tpref))
  }
  return(Tpref)

}

# Example usage
# 
Tpref <- calc_Tpref(data, method = "median", exclude_start_minutes = 0, exclude_end_minutes = 0)
# 
# print(Tpref)

#### Function to calculate upper and lower avoidance temperatures ####

calc_Tavoid <- function(data, percentiles = c(0.05, 0.95), exclude_start_minutes = 0, exclude_end_minutes = 0, print_results = T) {
  # Ensure the percentiles are valid
  if (length(percentiles) != 2 || any(percentiles < 0) || any(percentiles > 1)) {
    stop("Tavoid percentiles should be a vector of two values between 0 and 1")
  }
  
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Calculate the percentiles
  Tavoid_lower <- quantile(data$core_T, percentiles[1], na.rm = TRUE)
  Tavoid_upper <- quantile(data$core_T, percentiles[2], na.rm = TRUE)
  
  if (print_results) {
  # Print values
  print(paste("Tavoid Lower:", Tavoid_lower))
  print(paste("Tavoid Upper:", Tavoid_upper)) 
  }
  
  return(c(Tavoid_lower, Tavoid_upper))

}

# Example usage
# Calculate Tavoid values

Tavoid_values <-calc_Tavoid(data, percentiles = c(0.05, 0.95), exclude_start_minutes = 0, exclude_end_minutes = 0, print_results = F)

#### Function to calculate total distance moved ####

# calc_distance1 <- function(data){
#   data$x_pos[data$x_pos == "No object"] <- NA
#   data$y_pos[data$y_pos == "No object"] <- NA
#   
#   data$x_pos <- as.numeric(data$x_pos)*data$pixel_ratio
#   data$y_pos <- as.numeric(data$y_pos)*data$pixel_ratio
#   
#   data$delta_distance <- NA
#   data$delta_distance [1] <- 0
#   
#   for (i in 2:nrow(data)) {
#     data$delta_distance[i] <- sqrt((data$x_pos[i] - data$x_pos[i-1])^2 + (data$y_pos[i]-data$y_pos[i-1])^2)
#   }
#   data$delta_distance[is.na(data$delta_distance)]<-0
#   
#   data$cum_dist <- 0
#   for (i in 2:nrow(data)){
#     data$cum_dist[i]<- data$cum_dist[i-1]+data$delta_distance[i-1]
#   }
#   
#   data = subset(data, select = -c(delta_distance))
#   
#   
#   return(data)
#   
# }

# data<-calc_distance1(data)

calc_distance <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0, print_results = T) {
  # Ensure the Distance moved column exists
  if (!"distance" %in% colnames(data)) {
    stop("The dataset does not contain a 'distance' column")
  }
  
  # Convert exclude minutes to seconds
  # Exclude initial and final data if necessary
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - exclude_end_minutes * 60
  
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  initial_distance <- data$distance[1]
  
  data$distance <- data$distance - initial_distance

  # Calculate net distance moved
  total_distance <- max(data$distance, na.rm = TRUE)
  
  # Print the result
  if (print_results) {
  print(paste("Distance Moved:", total_distance, "cm"))
  }
  
  return(total_distance)
}

# Example usage
dist <- calc_distance(data)

#### Function to calculate shuttling frequency ####

calc_shuttling_frequency <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0, print_results = T) {
  # Ensure the 'shuttle' and 'time_sec' columns exist
  if (!all(c("shuttle", "time_sec") %in% colnames(data))) {
    stop("The dataset does not contain 'shuttle' and/or 'time_sec' columns. Run file_prepare on your data.")
  }
  
  # Convert time_sec to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Calculate shuttling frequency
  shuttles <- sum(data$shuttle, na.rm = TRUE)

  if (print_results) {
  # Print the result
  print(paste("Shuttling frequency:", shuttles))
  }
  
  return(shuttles)
}

# Example usage
shut_freq <- calc_shuttling_frequency(data)

#### Function to calculate occupancy time in each chamber ####

calc_occupancy_time <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0, print_results = T) {
  # Ensure the 'zone' and 'time_sec' columns exist
  if (!all(c("zone", "time_sec") %in% colnames(data))) {
    stop("The dataset does not contain 'zone' and/or 'time_sec' columns.")
  }
  
  # Convert time_sec to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Calculate occupancy time in each chamber
  time_in_chamberDECR <- sum(data$zone == "DECR", na.rm = TRUE)
  time_in_chamberINCR <- sum(data$zone == "INCR", na.rm = TRUE)
  
  if (print_results) {
  # Print the results
  print(paste("Time in DECR Chamber:", time_in_chamberDECR))
  print(paste("Time in INCR Chamber:", time_in_chamberINCR))
  }
  
  return(c(time_in_chamberDECR, time_in_chamberINCR))
}

# Example usage
occ_time <- calc_occupancy_time(data)

#### Function to plot a histogram of time spent at different core body temperatures ####

library(ggplot2)

plot_core_T_histogram <- function(data, bin_size, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure the body core temperature column exists
  if (!"core_T" %in% colnames(data)) {
    stop("The dataset does not contain a 'core_T' column. Please run file_prepare and calc_coreT")
  }
  
  # Convert time_sec to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Calculate the duration of each bin in minutes
  data$Time_min <- data$time_sec / 60
  total_time <- max(data$Time_min, na.rm = TRUE)
  
  # Create the histogram data
  hist_data <- hist(data$core_T, breaks = seq(floor(min(data$core_T, na.rm = TRUE)),
                                                      ceiling(max(data$core_T, na.rm = TRUE)),
                                                      by = bin_size), plot = FALSE)
  
  # Find the peak value
  Tpref_peak <- hist_data$mids[which.max(hist_data$counts)]
  
  # Create the histogram plot
  hist_plot <- ggplot(data, aes(x = core_T)) +
    geom_histogram(binwidth = bin_size, aes(y = (..count.. / sum(..count..)) * 100), fill = "#F79518", color = "black") +
    geom_vline(xintercept = Tpref_peak, color = "#9EF1A4", linetype = "dashed", linewidth = 2) +
    labs(title = "Histogram of Percent of Time Spent at Different Body Core Temperatures",
         subtitle = paste("Peak Tpref:", round(Tpref_peak, 2), "°C"),
         x = "Body Core Temperature (°C)",
         y = "Percent of Total Time (%)") +
    theme_minimal()
  
  print(hist_plot)
  print(paste("Peak Tpref:", Tpref_peak))
  
  return(Tpref_peak)
}

# # Example usage
# Tpref_peak <- plot_core_T_histogram(data, bin_size = 0.1, exclude_start_minutes = 240, exclude_end_minutes = 5)
# print(paste("Peak Tpref:", Tpref_peak))

#### Function to calculate variance in core body temperature experienced throughout the trial #####

calc_variance <- function(data, variance_type = c("std_error", "std_deviation", "coeff_variation"), exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure the body core temperature column exists
  if (!"core_T" %in% colnames(data)) {
    stop("The dataset does not contain a 'core_T' column. Please run file_prepare and calc_coreT")
  }
  
  # Convert exclude minutes to seconds
  exclude_start_seconds <- exclude_start_minutes * 60
  exclude_end_seconds <- exclude_end_minutes * 60
  
  # Exclude initial and final data if necessary
  start_time <- exclude_start_seconds
  end_time <- max(data$time_sec) - exclude_end_seconds
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Select the variance type
  variance_type <- match.arg(variance_type)
  
  # Calculate the variance measure
  if (variance_type == "std_error") {
    variance <- sd(data$core_T, na.rm = TRUE) / sqrt(length(na.omit(data$core_T)))
  } else if (variance_type == "std_deviation") {
    variance <- sd(data$core_T, na.rm = TRUE)
  } else if (variance_type == "coeff_variation") {
    variance <- sd(data$core_T, na.rm = TRUE) / mean(data$core_T, na.rm = TRUE)
  }
  
  return(variance)
}

# # Example usage
variance_se <- calc_variance(data, variance_type = "std_error", 240, 5)
# variance_sd <- calc_variance(data, variance_type = "std_deviation", 240, 5)
# variance_cv <- calc_variance(data, variance_type = "coeff_variation", 240, 5)
# 
# print(paste("Standard Error:", variance_se))
# print(paste("Standard Deviation:", variance_sd))
# print(paste("Coefficient of Variation:", variance_cv))

#### Function to plot cumulative changes in distance moved over time ####

  plot_cumulative_distance <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0) {
    if (!("time_sec" %in% colnames(data) && "distance" %in% colnames(data))) {
      stop("The dataset does not contain 'time_sec' and/or 'distance' columns.")
    }
    
    # Convert time_sec to numeric if not already
    data$time_sec <- as.numeric(data$time_sec)
    
    # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
    start_time <- exclude_start_minutes * 60
    end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
    data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
    
    # Convert time in seconds to minutes
    data$Time_min <- data$time_sec / 60
    #substract first distance value to reset cumulative distance to 0
    data$distance<-data$distance-data$distance[[1]]
    
    # Plot cumulative distance
    plot <- ggplot(data, aes(x = Time_min, y = (distance))) +
      geom_line() +
      labs(title = "Cumulative Distance Moved During Trial",
           x = "Time (minutes)",
           y = "Distance Moved (cm)") +
      theme_minimal()
    
    print(plot)
  }

# Example usage
plot_cumulative_distance(data, 240, 5)

#### Function to calculate minimum gravitation time ####

calc_min_gravitation <- function(start_T, target_T, rate_of_change, print_results = T) {
  # Calculate the change in temperature
  delta_T <- target_T - start_T
  
  # Calculate the minimum gravitation time
  gravitation_time <- delta_T / (rate_of_change/60)
  
  if (print_results) {
  # Print the result
  print(paste("Minimum Gravitation Time:", gravitation_time, "minutes"))
  }
  
  return(gravitation_time)
}

# # Example usage
start_T <- 27  # Starting mean temperature in °C
target_T <- 35  # Target temperature in °C, e.g. approximation of anticipated Tpref
rate_of_change <- 5  # Rate of temperature change in °C/h

min_grav_time <- calc_min_gravitation(start_T, target_T, rate_of_change)

#### Function to calculate time spent at extreme temperatures, near set limits ####

calc_extremes <- function(data, threshold = 0.2*(max(data$max_T)-max(data$min_T)), exclude_start_minutes = 0, exclude_end_minutes = 0, print_results = T) {
  # Ensure necessary columns exist
  if (!("time_sec" %in% colnames(data) && "core_T" %in% colnames(data))) {
    stop("The dataset does not contain 'time_sec' and/or 'core_T' columns.")
  }
  
  # Convert time_sec to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Convert time in secxonds to minutes
  data$Time_min <- data$time_sec / 60
  
  
  upper_limit<-max(data$max_T)
  lower_limit<-max(data$min_T)
  
  # Determine time spent near upper and lower extreme temperatures with the threshold
  upper_threshold <- upper_limit - threshold
  lower_threshold <- lower_limit + threshold
  data$Upper_Extreme_T <- ifelse(data$core_T > upper_threshold, 1, 0)
  data$Lower_Extreme_T <- ifelse(data$core_T < lower_threshold, 1, 0)
  
  time_near_upper_extreme <- sum(data$Upper_Extreme_T) / nrow(data) * 100
  time_near_lower_extreme <- sum(data$Lower_Extreme_T) / nrow(data) * 100
  
  if (print_results) {
  # Print the results
  print(paste("Time spent near upper extreme temperatures (%):", time_near_upper_extreme))
  print(paste("Time spent near lower extreme temperatures (%):", time_near_lower_extreme))
  }
  
  return(c(time_near_lower_extreme, time_near_upper_extreme))
}

# # Example usage
extreme_time_percent <- calc_extremes(data)
calc_extremes(data, exclude_start_minutes = 0, exclude_end_minutes = 0)[2]


#### Function to plot histogram of movement speeds ####

library(ggplot2)

plot_speed_histogram <- function(data, binwidth = 0.1, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary column exists
  if (!"velocity" %in% colnames(data)) {
    stop("The dataset does not contain a 'velocity' column.")
  }
  
  # Convert time_sec to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec, na.rm = TRUE) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Create the histogram plot
  plot <- ggplot(data, aes(x = velocity)) +
    geom_histogram(binwidth = binwidth, fill = "#F79518", color = "black") +
    labs(title = "Histogram of Movement Speeds",
         x = "Movement Speed (cm/s)",
         y = "Frequency") +
    # scale_y_log10()+
    theme_classic()
  
  print(plot)
}

# Example usage
# plot_speed_histogram(data, binwidth = 1, exclude_start_minutes = 240, exclude_end_minutes = 5)

#### Function to plot movement speed against body core temperature ####

plot_speed_vs_core_T <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!all(c("core_T", "distance", "time_sec") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'core_T', 'distance', 'time_sec'. Please run file_prepare function.")
  }
  
  # Convert time_sec to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Calculate movement speed as the difference in distance moved over time
  data$Movement_Speed <- c(NA, diff(data$distance) / diff(data$time_sec))
  
  # Remove rows with NA values
  data <- na.omit(data[, c("core_T", "Movement_Speed")])
  
  # Create the scatter plot
  plot <- ggplot(data, aes(x = core_T, y = Movement_Speed)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "#F79518", se = FALSE) +
    labs(title = "Movement Speed vs. Core Body Temperature",
         x = "Core Body Temperature (°C)",
         y = "Movement Speed (cm/s)") +
    theme_classic()
  
  return(plot)
}

# Example usage
# Create the plot
# speed_temp_plot <- plot_speed_vs_core_T(data, exclude_start_minutes = 240)
# speed_temp_plot

#### Function to plot temperatures in each side of the shuttlebox over time ####

library(ggplot2)

plot_temperature_gradient <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Convert time_sec to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Convert time in seconds to hours
  data$time_h <- data$time_sec / 3600
  
  # Plot the temperature gradient
  plot <- ggplot(data, aes(x = time_h)) +
    geom_line(aes(y = `DECR_T`, color = "DECR side temp")) +
    geom_line(aes(y = `INCR_T`, color = "INCR side temp")) +
    scale_color_manual(values = c("DECR side temp" = "dodgerblue", "INCR side temp" = "brown1")) +
    labs(title = "Temperature Gradient Over Time",
         x = "Time (hours)",
         y = "Temperature (°C)",
         color = "Temperature Side") +
    theme_classic()
  
  print(plot)
}

# Example usage
# Assuming 'data' contains temperature columns for both chambers
# plot_temperature_gradient(data, 0, 0)

#### Function to plot core body temperature, Tpref, avoidance temperatures, segmented regression of changes in body core temperature during the trial #####

library(ggplot2)
library(segmented)

plot_T_segmented <- function(data, Tpref, Tavoid_lower, Tavoid_upper, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!all(c("time_h", "core_T") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'time_h', 'core_T'.")
  }
  
  # Convert the time to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Perform segmented regression
  fit <- lm(core_T ~ time_h, data = data)
  seg_fit <- segmented(fit, seg.Z = ~time_h, npsi = 1)
  
  # Get summary and breakpoint
  seg_summary <- summary(seg_fit)
  breakpoints <- seg_fit$psi[, "Est."]
  
    # Plot data with segmented regression
  plot <- ggplot(data, aes(x = time_h, y = core_T)) +
    geom_point(alpha = 0.5, size = 0.4) +
    geom_line(aes(y = predict(seg_fit)), color = "#F79518", size = 1) +
    geom_hline(yintercept = Tpref, linetype = "dashed", color = "#9EF1A4", size = 1.75) +
    geom_hline(yintercept = Tavoid_lower, linetype = "dashed", color = "#9ECFF1", size = 1.75) +
    geom_hline(yintercept = Tavoid_upper, linetype = "dashed", color = "#F1AD9E", size = 1.75) +
    labs(title = "Body Core Temperature vs Time with Segmented Regression",
         x = "Time (h)",
         y = "Body Core Temperature (°C)") +
    theme_classic()
  
  print(plot)
  
  return(list(summary = seg_summary, breakpoints = breakpoints))
  
  }

# Example usage
plot_T_segmented(data, Tpref = Tpref, Tavoid_values[1], Tavoid_values[2])

#### Function to plot changes in core body temperature over time, as well as Tpref, and upper and lower avoidance temperatures, but WITHOUT segmented regression lines #####

plot_T <- function(data, Tpref, Tavoid_lower, Tavoid_upper, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!all(c("time_h", "core_T") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'time_h', 'core_T'.")
  }
  
  # Convert the time to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Plot data
  plot <- ggplot(data, aes(x = time_h, y = core_T)) +
    geom_point(alpha = 0.5, size = 0.4) +
    geom_hline(yintercept = Tpref, linetype = "dashed", color = "#9EF1A4", size = 1.75) +
    geom_hline(yintercept = Tavoid_lower, linetype = "dashed", color = "#9ECFF1", size = 1.75) +
    geom_hline(yintercept = Tavoid_upper, linetype = "dashed", color = "#F1AD9E", size = 1.75) +
    labs(title = "Body Core Temperature vs Time with Segmented Regression",
         x = "Time (h)",
         y = "Body Core Temperature (°C)") +
    theme_classic()
  
  print(plot)
}

# Example usage
plot_T(data, Tpref, Tavoid_values[1], Tavoid_values[2])

##### Function to calculate actual gravitation time ####

library(segmented)

  calc_act_gravitation <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0, print_results = T) {
  
  # Convert the time to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Perform segmented regression
  fit <- lm(core_T ~ time_h, data = data)
  seg_fit <- segmented(fit, seg.Z = ~time_h, npsi = 1)
  
  # Get summary and breakpoint
  seg_summary <- summary(seg_fit)
  breakpoints <- seg_fit$psi[, "Est."]
  
  # Extract breakpoints
  breakpoints <- seg_fit$psi[, "Est."]
  
  if (length(breakpoints) == 0) {
    stop("No breakpoints found in the segmented model.")
  }
  
  # Assume the first breakpoint as the gravitation time
  gravitation_time <- breakpoints[1]
  
  if (print_results) {
  print(paste("Gravitation Time (hours):", gravitation_time))
  }
  
  return(gravitation_time)
}

# Calculate gravitation time
act_grav_time <- calc_act_gravitation(data)

##### Function to calculate and plot interval means for shuttling rate over time ####

library(ggplot2)
library(dplyr)

shuttling_aggregated <- function(data, interval_minutes = 10, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!all(c("time_sec", "shuttle") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'time_sec', 'shuttle'.")
  }
  
  # Convert time_sec to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec, na.rm = TRUE) - (exclude_end_minutes * 60)
  data <- data %>% filter(time_sec > start_time & time_sec < end_time)
  
  # Calculate time intervals
  data <- data %>% 
    mutate(Time_interval = floor((time_sec - start_time) / (interval_minutes * 60)) * interval_minutes)
  
  # Calculate mean shuttling rate for each interval
  shuttling_rate <- data %>% 
    group_by(Time_interval) %>% 
    summarise(mean_shuttle_rate = mean(shuttle, na.rm = TRUE))
  
  # Convert Time_interval to numeric for plotting
  shuttling_rate$Time_interval <- as.numeric(shuttling_rate$Time_interval)
  
  # Plot the mean shuttling rate
  plot <- ggplot(shuttling_rate, aes(x = Time_interval/60, y = mean_shuttle_rate*60)) +
    geom_line(color = "black") +
    geom_point(color = "#F79518") +
    labs(title = "Mean Shuttling Rate Over the Trial",
         x = "Time (h)",
         y = "Mean Shuttling Rate (shuttles/hour)") +
    theme_classic()
  
  print(plot)
  
  return(shuttling_rate)
}

# Example usage
shuttling_rate <- shuttling_aggregated(data, interval_minutes = 60)

##### Function to calculate and plot interval means for speed over time ####

library(ggplot2)
library(dplyr)

speed_aggregated <- function(data, interval_minutes = 10, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!all(c("time_sec", "velocity") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'time_sec', 'velocity'.")
  }
  
  # Convert time_sec to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec, na.rm = TRUE) - (exclude_end_minutes * 60)
  data <- data %>% filter(time_sec > start_time & time_sec < end_time)
  
  # Calculate time intervals
  data <- data %>% 
    mutate(Time_interval = floor((time_sec - start_time) / (interval_minutes * 60)) * interval_minutes)
  
  # Calculate mean movement speed for each interval
  movement_speed <- data %>% 
    group_by(Time_interval) %>% 
    summarise(mean_movement_speed = mean(velocity, na.rm = TRUE))
  
  # Convert Time_interval to numeric for plotting
  movement_speed$Time_interval <- as.numeric(movement_speed$Time_interval)
  
  # Plot the mean movement speed
  plot <- ggplot(movement_speed, aes(x = Time_interval/60, y = mean_movement_speed)) +
    geom_line(color = "black") +
    geom_point(color = "#F79518") +
    labs(title = "Mean Movement Speed Over the Trial",
         x = "Time (h)",
         y = "Mean Movement Speed (cm/s)") +
    theme_classic()
  
  print(plot)
  
  return(movement_speed)
}

# Example usage
movement_speed <- speed_aggregated(data, interval_minutes = 60, exclude_end_minutes = 5)

##### Function to calculate and plot interval means for body core temperature over time ####

library(ggplot2)
library(dplyr)

Tcore_aggregated <- function(data, interval_minutes = 10, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!all(c("time_sec", "core_T") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'time_sec', 'core_T'.")
  }
  
  # Convert time_sec to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec, na.rm = TRUE) - (exclude_end_minutes * 60)
  data <- data %>% filter(time_sec > start_time & time_sec < end_time)
  
  # Calculate time intervals
  data <- data %>% 
    mutate(Time_interval = floor((time_sec - start_time) / (interval_minutes * 60)) * interval_minutes)
  
  # Calculate mean preferred temperature for each interval
  Tcore_aggregated <- data %>% 
    group_by(Time_interval) %>% 
    summarise(mean_Tcore = mean(core_T, na.rm = TRUE))
  
  # Convert Time_interval to numeric for plotting
  Tcore_aggregated$Time_interval <- as.numeric(Tcore_aggregated$Time_interval)
  
  # Plot the mean preferred temperature
  plot <- ggplot(Tcore_aggregated, aes(x = Time_interval/60, y = mean_Tcore)) +
    geom_line(color = "black") +
    geom_point(color = "#F79518") +
    labs(title = "Mean Body Core Temperature Over the Trial",
         x = "Time (h)",
         y = "Mean Core Temperature (°C)") +
    theme_classic()
  
  print(plot)
  
  return(Tcore_aggregated)
}

# Example usage
Tcore_agg <- Tcore_aggregated(data, interval_minutes = 60, exclude_start_minutes = 10, exclude_end_minutes = 5)

##### Function to calculate interval means for speed, shuttling rate, and Tcore, then produces a correlation matrix ####

library(ggplot2)
library(dplyr)

mean_correlations_scatterplots <- function(data, interval_minutes = 10, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!all(c("time_sec", "velocity", "core_T", "shuttle") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'time_sec', 'velocity', 'core_T', 'shuttle'.")
  }
  
  # Convert time_sec to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec, na.rm = TRUE) - (exclude_end_minutes * 60)
  data <- data %>% filter(time_sec > start_time & time_sec < end_time)
  
  # Calculate time intervals
  data <- data %>% 
    mutate(Time_interval = floor((time_sec - start_time) / (interval_minutes * 60)) * interval_minutes)
  
  # Calculate mean values for velocity, Tpref (core_T), and shuttling rate for each interval
  summary_data <- data %>% 
    group_by(Time_interval) %>% 
    summarise(mean_velocity = mean(velocity, na.rm = TRUE),
              mean_Tpref = mean(core_T, na.rm = TRUE),
              mean_shuttle_rate = mean(shuttle, na.rm = TRUE))
  
  # Calculate correlation matrix
  correlation_matrix <- cor(summary_data[, -1], use = "complete.obs")
  
  # Print the correlation matrix
  print(correlation_matrix)
  
  # Create scatterplots for pairwise comparisons
  plot1 <- ggplot(summary_data, aes(x = mean_velocity, y = mean_Tpref)) +
    geom_point(color = "#F79518") +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    labs(title = "Mean velocity vs Mean Tpref",
         x = "Mean velocity (cm/s)",
         y = "Mean Tpref (°C)") +
    theme_classic()
  
  plot2 <- ggplot(summary_data, aes(x = mean_velocity, y = mean_shuttle_rate)) +
    geom_point(color = "#F79518") +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    labs(title = "Mean velocity vs Mean Shuttling Rate",
         x = "Mean velocity (cm/s)",
         y = "Mean Shuttling Rate") +
    theme_classic()
  
  plot3 <- ggplot(summary_data, aes(x = mean_Tpref, y = mean_shuttle_rate)) +
    geom_point(color = "#F79518") +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    labs(title = "Mean Tpref vs Mean Shuttling Rate",
         x = "Mean Tpref (°C)",
         y = "Mean Shuttling Rate") +
    theme_classic()
  
  # Print the plots
  print(plot1)
  print(plot2)
  print(plot3)
  
  return(list(summary_data = summary_data, correlation_matrix = correlation_matrix, plots = list(plot1, plot2, plot3)))
}

# Example usage
results <- mean_correlations_scatterplots(data, interval_minutes = 10, exclude_start_minutes = 10, exclude_end_minutes = 5)
summary_data <- results$summary_data
correlation_matrix <- results$correlation_matrix
plots <- results$plots

correlation_matrix


# results
##### Function to plot a 3D animation of fish movements at selected time intervals, with time as a dimension ####

library(rgl)

animate_movements <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!all(c("x_pos", "y_pos", "time_sec", "core_T") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'x_pos', 'y_pos', 'time_sec', 'core_T'.")
  }
  
  # Convert time_sec to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec, na.rm = TRUE) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Set up 3D plot
  plot3d(data$x_pos, data$y_pos, data$time_sec, type = "l", col = "black", alpha = 0.15,
         xlab = "X Position", ylab = "Y Position", zlab = "Time (seconds)")
  
  # Determine color scaling based on Tpref
  colors <- colorRampPalette(c("royalblue3", "red2"))(length(unique(data$core_T)))
  data$color <- colors[as.numeric(cut(data$core_T, breaks = length(colors)))]
  
  # Animate the plot
  for (i in seq_len(nrow(data))) {
    points3d(data$x_pos[i], data$y_pos[i], data$time_sec[i], col = data$color[i], size = 5)
    Sys.sleep(0.01)
  }
}

# Example usage
# animate_movements(data, exclude_end_minutes = 5)

##### Function to produce a heat map of fish locations over the trial ####

library(ggplot2)

plot_heatmap <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!all(c("x_pos", "y_pos", "time_sec") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'x_pos', 'y_pos', 'time_sec'.")
  }
  
  # Convert time_sec to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec, na.rm = TRUE) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Plot heatmap
  heatmap_plot <- ggplot(data, aes(x = x_pos, y = y_pos)) +
    stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_fill_viridis_c() +
    labs(title = "Heat Map of Fish Locations within the Shuttlebox",
         x = "X Position",
         y = "Y Position") +
    theme_classic() +
    coord_fixed()
  
  return(heatmap_plot)
}

# Example usage
heatmap_plot <- plot_heatmap(data)
print(heatmap_plot)

##### Function to visualise links between shuttling and activity across individuals in the entire project dataset ####

library(ggplot2)

plot_distance_vs_shuttles <- function(proj_data) {
  # Ensure necessary columns exist
  if (!all(c("ID", "distance", "shuttles") %in% colnames(proj_data))) {
    stop("The dataset does not contain one or more necessary columns: 'ID', 'distance', 'shuttles'.")
  }
  
  # Create the plot
  plot <- ggplot(proj_data, aes(x = distance, y = shuttles, label = ID)) +
    geom_point(color = "#F79518") +
    geom_text(vjust = -1, hjust = 1.5) +
    labs(title = "Relationship between Distance Moved and Shuttles across Individuals",
         x = "Distance Moved (cm)",
         y = "Number of Shuttles") +
    theme_classic()
  
  print(plot)
}

# Example usage
plot_distance_vs_shuttles(proj_data)

##### Function to visualise links between time near system temperature limits and activity across individuals in the entire project dataset ####

library(ggplot2)

plot_limits_vs_distance <- function(proj_data) {
  # Ensure necessary columns exist
  if (!all(c("ID", "distance", "time_near_limits") %in% colnames(proj_data))) {
    stop("The dataset does not contain one or more necessary columns: 'ID', 'distance', 'time_near_limits'.")
  }
  
  # Create the plot
  plot <- ggplot(proj_data, aes(x = time_near_limits, y = distance, label = ID)) +
    geom_point(color = "#F79518") +
    geom_text(vjust = -1, hjust = 1.5) +
    labs(title = "Relationship between Time Near Limits and Distance Moved",
         x = "Time Near Limits (minutes)",
         y = "Distance Moved (cm)") +
    theme_classic()
  
  print(plot)
}

# Example usage
# plot_limits_vs_distance(proj_data)


##### Function to visualise links between time near system temperature limits and shuttles across individuals in the entire project dataset ####

library(ggplot2)

plot_limits_vs_shuttles <- function(proj_data) {
  # Ensure necessary columns exist
  if (!all(c("ID", "time_near_limits", "shuttles") %in% colnames(proj_data))) {
    stop("The dataset does not contain one or more necessary columns: 'ID', time_near_limits', 'shuttles'.")
  }
  
  # Create the plot
  plot <- ggplot(proj_data, aes(x = time_near_limits, y = shuttles, label = ID)) +
    geom_point(color = "#F79518") +
    geom_text(vjust = -1, hjust = 1.5) +
    labs(title = "Relationship between Time Near Limits and Number of shuttles",
         x = "Time Near Limits (minutes)",
         y = "Shuttles") +
    theme_classic()
  
  print(plot)
}

# Example usage
# plot_limits_vs_shuttles(proj_data)

##### Function to perform PCA with shuttles, distance moved, and time near limits, then plot PC scores vs Tpref ####

library(ggplot2)
library(dplyr)
library(FactoMineR)
library(factoextra)

pca_and_plot <- function(proj_data) {
  # Ensure necessary columns exist
  if (!all(c("ID", "distance", "shuttles", "time_near_limits", "Tpref") %in% colnames(proj_data))) {
    stop("The dataset does not contain one or more necessary columns: 'ID', 'distance', 'shuttles', 'time_near_limits', 'Tpref'.")
  }
  
  # Select the columns for PCA
  pca_data <- proj_data %>%
    dplyr::select(distance, shuttles, time_near_limits) %>%
    na.omit()  # Remove any rows with missing values
  
  # Perform PCA with scaling
  pca_result <- PCA(pca_data, scale.unit = TRUE, graph = FALSE)
  
  # Extract the scores of the first principal component
  scores <- data.frame(ID = proj_data$ID, PC1 = pca_result$ind$coord[, 1], Tpref = proj_data$Tpref)
  
  # Extract the loadings
  loadings <- data.frame(Variable = rownames(pca_result$var$coord), pca_result$var$coord)
  
  # Display the loadings
  print("PCA Loadings:")
  print(loadings)
  
  # Create the plot with 95% CI
  plot <- ggplot(scores, aes(x = PC1, y = Tpref, label = ID)) +
    geom_point(color = "#F79518") +
    geom_smooth(method = "lm", color = "black") +  # Adding 95% CI shading
    geom_text(vjust = -1, hjust = 1.5) +
    labs(title = "PCA Scores from the First Component vs Tpref with 95% CI",
         x = "PCA First Component Score",
         y = "Tpref (°C)") +
    theme_classic()
  
  print(plot)
  
  return(list(pca_result = pca_result, scores = scores, loadings = loadings))
}

# # Example usage
results <- pca_and_plot(proj_data)
pca_result <- results$pca_result
scores <- results$scores
loadings <- results$loadings
# 
loadings
##### Function to create frequency distributions for key measures across individuals in the data set ####

library(ggplot2)
library(gridExtra)

plot_histograms <- function(proj_data, bin_size_Tpref = 1, bin_size_Tavoid_upper = 1, bin_size_Tavoid_lower = 1, bin_size_distance = 2000, bin_size_shuttles = 20, bin_size_time_near_limits = 5) {
  # Ensure necessary columns exist
  required_columns <- c("Tpref", "Tavoid_upper", "Tavoid_lower", "distance", "shuttles", "time_near_limits")
  if (!all(required_columns %in% colnames(proj_data))) {
    stop("The dataset does not contain one or more necessary columns: 'Tpref', 'Tavoid_upper', 'Tavoid_lower', 'distance', 'shuttles', 'time_near_limits'.")
  }
  
  # Create individual plots
  p1 <- ggplot(proj_data, aes(x = Tpref)) +
    geom_histogram(binwidth = bin_size_Tpref, fill = "#F79518", color = "black") +
    labs(title = "Distribution of Tpref", x = "Tpref (°C)", y = "Frequency") +
    theme_classic()
  
  p2 <- ggplot(proj_data, aes(x = Tavoid_upper)) +
    geom_histogram(binwidth = bin_size_Tavoid_upper, fill = "#E41A1C", color = "black") +
    labs(title = "Distribution of Tavoid_upper", x = "Tavoid_upper (°C)", y = "Frequency") +
    theme_classic()
  
  p3 <- ggplot(proj_data, aes(x = Tavoid_lower)) +
    geom_histogram(binwidth = bin_size_Tavoid_lower, fill = "#377EB8", color = "black") +
    labs(title = "Distribution of Tavoid_lower", x = "Tavoid_lower (°C)", y = "Frequency") +
    theme_classic()
  
  p4 <- ggplot(proj_data, aes(x = distance)) +
    geom_histogram(binwidth = bin_size_distance, fill = "#4DAF4A", color = "black") +
    labs(title = "Distribution of Distance", x = "Distance", y = "Frequency") +
    theme_classic()
  
  p5 <- ggplot(proj_data, aes(x = shuttles)) +
    geom_histogram(binwidth = bin_size_shuttles, fill = "#984EA3", color = "black") +
    labs(title = "Distribution of Shuttles", x = "Shuttles", y = "Frequency") +
    theme_classic()
  
  p6 <- ggplot(proj_data, aes(x = time_near_limits)) +
    geom_histogram(binwidth = bin_size_time_near_limits, fill = "#EEDD82", color = "black") +
    labs(title = "Distribution of Time Near Limits", x = "Time Near Extremes", y = "Frequency") +
    theme_classic()
  
  # Combine the plots into a multipanel plot
  multipanel_plot <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)
  
  return(multipanel_plot)
}

# Example usage
multipanel_plot <- plot_histograms(proj_data)



#### COMPILE PROJ_DATA FILE ####



# TO DO
# fix Tpref
# fix fileID

read.shuttlesoft.project <- function (metadata, directory = getwd()){
  
  # Check if metadata is specified and give warnings for what happens if data is missing
  if (missing(metadata)) {
    warning("Metadata not provided; cannot define one or more of the following variables: 'trial_start', 'a-value', 'b-value', 'initial_temp.'")
  } else {
    
    # If the metadata is present, ensure necessary columns exist in the metadata
    if (!all(c("file_name", "trial_start", "mass", "a_value", "b_value") %in% colnames(metadata))) {
      warning("The metadata does not contain one or more necessary columns: 'file_name', 'trial_start', 'mass', 'a_value', 'b_value'.")
    }
    
    # Find rows with NA in specified columns
    columns_to_check <- c("file_name", "trial_start", "mass", "a_value", "b_value")
    rows_with_na <- apply(metadata[columns_to_check], 1, function(row) any(is.na(row)))
    
    # If there are any rows with NA, throw warning mentioning file names
    if (any(rows_with_na)) {
      names_with_na <- metadata$file_name[rows_with_na]
      warning (paste("The metadata contains NA values in one or more of these columns: 'file_name', 'trial_start', 'mass', 'a_value', 'b_value'
",paste(names_with_na, collapse = "\n ")))
    }
    
    #Ensure necessary columns are in the right format
    if (all(grepl("^\\d{2}:\\d{2}:\\d{2}$", metadata$trial_start))==F){
      message("Warning: One or more entries in trial_start are not in an hh:mm:ss format")
    }
  }
  
  # Search for all the .txt files in the specified directory
  txt_files <- list.files(path = directory, pattern = "\\.txt$", full.names = TRUE)
  
  # Read those text files and compile them into a list
  data_read <- lapply(txt_files, read.shuttlesoft, metadata = metadata, multidat = T)
  
  # Name each list element as the directory of the text files
  names(data_read)<-txt_files
  
  return(data_read)
}

data_read<-read.shuttlesoft.project(metadata)

compile_project_data <- function(data_read,
                                 exclude_acclimation = F,
                                 exclude_start_minutes = 0, 
                                 exclude_end_minutes = 0,
                                 Tpref_method = "median",
                                 Tavoid_percintiles = c(0.05, 0.95),
                                 textremes_threshold = expression(0.2 * (max(df$max_T) - max(df$min_T))),
                                 core_T_variance_type = "std_error"){
  
  # Prepare the data using file_prepare. Multidat is TRUE because this will suppress warnings from preparing each data file separately
  data_list <- lapply(data_read, file_prepare)

  
  # Create a function that pulls all the information from a datafile in one go, using the functions specified earlier 
  apply_functions <- function(df, df_name,
                              exc_accl = exclude_acclimation,
                              exc_start = exclude_start_minutes, 
                              exc_end = exclude_end_minutes,
                              Tpref_met = Tpref_method,
                              Tavoid_perc = Tavoid_percintiles,
                              textremes_th = eval(textremes_threshold),
                              Tcore_variance = core_T_variance_type) {
    
    # Define the fileID the function is currently considering
    fileID <- unique(df$fileID)
    
    # Calculate core body temperature
    df<-calc_coreT(df)

    # Exclude acclimation period if requested from the start. Acclimation period is defined in file_prepare
    ifelse(exclude_acclimation == T, df<-subset(df, df$trial_phase != "acclimation"), df)
    
    # Calculate all variables using default settings (these can be edited in the main function)
    Tpref <- calc_Tpref(df, method = Tpref_met, exclude_start_minutes = exc_start, exclude_end_minutes = exc_end, print_results = F)
    grav_time <- calc_act_gravitation(df, exclude_start_minutes = exc_start, exclude_end_minutes = exc_end, print_results = F)
    distance <- calc_distance(df, exclude_start_minutes = exc_start, exclude_end_minutes = exc_end, print_results = F)
    Tavoid <- calc_Tavoid(df, percentiles = Tavoid_perc, exclude_start_minutes = exc_start, exclude_end_minutes = exc_end, print_results = F)
    Tpref_range <- Tavoid[2]-Tavoid[1]
    textremes <- calc_extremes(df, threshold = textremes_th, exclude_start_minutes = exc_start, exclude_end_minutes = exc_end, print_results = F)
    core_T_variance <- calc_variance(df, variance_type = Tcore_variance, exclude_start_minutes = exc_start, exclude_end_minutes = exc_end)
    nr_shuttles <- calc_shuttling_frequency(df, exclude_start_minutes = exc_start, exclude_end_minutes = exc_end, print_results = F)
    chamber_seconds <- calc_occupancy_time(df, exclude_start_minutes = exc_start, exclude_end_minutes = exc_end, print_results = F)
    
    # Return a list of the variables you want
    return(list(fileID = fileID,
                Tpref = Tpref,
                Tpref_range = Tpref_range,
                grav_time = grav_time, 
                distance = distance,
                Tavoid = Tavoid,
                textremes = textremes,
                core_T_variance = core_T_variance,
                nr_shuttles = nr_shuttles,
                chamber_seconds = chamber_seconds))
  }
  
  # Apply functions to each dataframe in the list and collect results
  results <- do.call(rbind, lapply(names(data_list), function(name) {
    
    # Apply the functions
    func_results <- apply_functions(df = data_list[[name]], df_name = name)
    
    # Combine results into a data frame
    data.frame(
      fileID = func_results$fileID,
      Tpref = func_results$Tpref,
      Tpref_range = func_results$Tpref_range,
      grav_time  = func_results$grav_time,
      distance  = func_results$distance,
      Tavoid_low  = func_results$Tavoid[1],
      Tavoid_high = func_results$Tavoid[2],
      t_near_low_extreme  = func_results$textremes[1],
      t_near_high_extreme  = func_results$textremes[2],
      core_T_variance  = func_results$core_T_variance,
      nr_shuttles  = func_results$nr_shuttles,
      seconds_in_DECR = func_results$chamber_seconds[1],
      seconds_in_INCR = func_results$chamber_seconds[2]
      
    )
    
  }))
  
  # Reset rownames as they have acquired names from the list
  rownames(results) <- NULL
  
  return(results)
}

results<-compile_project_data(data_read, exclude_acclimation = F)

calc_Tpref(data, method = "median")
calc_act_gravitation(data)
calc_distance(data)
Tavoid<-calc_Tavoid(data)
Tavoid[2] - Tavoid[1]
calc_extremes(data)
calc_variance(data)
calc_shuttling_frequency(data)
calc_occupancy_time(data)


write.csv(data, "data.csv")
