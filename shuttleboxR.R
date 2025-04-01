# Script for R package shuttleboxR
#
# 25/026/2024
#
# Shaun Killen and Emil Christensen
#


# TO DO:  
  # Edit inspect() to check type of data (e.g. time/numeric/whatever)
  # Write inspect.project function?
  # Add acclimation exclusion to each function
  # Remove dependencies as much as possible
    # Check if segmented regression can be done without segmented package

  # Make overall diagnostic function
    # Check diagnostic functions
    # Add option to highlight specific individuals in histograms

rm(list=ls())

# renv:: restore checks whether the packages you are using are the same versions as in the script
# if collaborators update a library and then edit the script, renv::restore will pick up on that
renv::status()
renv::restore()

# Confirm that R is looking in the right place
getwd()



library(segmented)
library(ggplot2)
library(rgl)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(dbscan)
library(gridExtra)


##### (PREPARE) Read shuttlesoft files function ####

# Collect existing information into one data file:
# date, time, zone, INCR_T, DECR_T, coordinates/distance, pixel ratio, trial_start, mass, a_value, b_value, acclimation temp
read_shuttlesoft <- function (file, metadata, multidat = F){
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
  if (length(data$zone[data$time == "Notes"]) == 0 || is.na(data$zone[data$time == "Notes"])) {
    notes <- NA
    }else{
      notes <- data$zone[data$time == "Notes"]
    }
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
  fileID <- basename(file)
  data$fileID <- fileID
  
  trialstart <-  function(metadata){
    
    if (missing(metadata)) {
      trial_start <- data$time[1]
      warning("Metadata not provided; cannot define 'trial_start'")
    } else {
      
      # If the metadata is present, ensure necessary columns exist in the metadata
      if (!all(c("file_name", "trial_start", "mass", "a_value", "b_value") %in% colnames(metadata))) {
        warning("The metadata does not contain one or more valuable columns: 'file_name', 'trial_start', 'mass', 'a_value', 'b_value'.")
      }
      
      # Find rows with NA in specified columns
      columns_to_check <- intersect(c("file_name", "trial_start", "mass", "a_value", "b_value"), colnames(metadata))
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
    
      if(!("trial_start" %in% colnames(metadata))){
        trial_start <- data$time[1]
      }else if("trial_start" %in% colnames(metadata) && is.na(metadata$trial_start[metadata$file_name==fileID])){
        trial_start <- data$time[1]
      }else {
        trial_start<-metadata$trial_start[metadata$file_name==fileID]
      }
    }
    return(trial_start)}
  
  if (multidat) {
    data$trial_start <- suppressWarnings(trialstart (metadata = metadata))} else {
      data$trial_start <- trialstart(metadata = metadata)}
 
  data$trial_start <- as.character(data$trial_start)
  
  data$a_value <- NA
  data$b_value <- NA
  data$mass <- NA
  data$initial_T <- NA
  data$a_value <- metadata$a_value[metadata$file_name == fileID]
  data$b_value <- metadata$b_value[metadata$file_name == fileID]
  data$mass <- metadata$mass[metadata$file_name == fileID]
  data$initial_T <- metadata$initial_T[metadata$file_name == fileID] 
  
   return(data)
}


##### (PREPARE) inspect data function ####
inspect <- function(data){
  
   if(!all(c("date", "time", "zone", "INCR_T", "DECR_T") %in% colnames(data))){
    stop("Dataset is missing one or more of the following necessary columns:  'date', 'time', 'zone', 'INCR_T', 'DECR_T'")
  }
  else if (!all(c("x_pos", "y_pos", "distance") %in% colnames(data))) {
    warning("Dataset is missing 'x_pos', 'y_pos', 'distance', no way to determine distance found")
  }
  else if(!all(c( "velocity", "time_in_INCR", "time_in_DECR",
            "delta_T", "dyn_hysteresis", "stat_T_INCR", "stat_hyst_INCR", "stat_T_DECR",
            "stat_hyst_DECR", "max_T", "min_T", "change_rate")%in% colnames(data))) {
    warning("Dataset is missing one or more of the following columns:'velocity', 'time_in_INCR', 'time_in_DECR', 'delta_T', 'dyn_hysteresis', 'stat_T_INCR', 'stat_hyst_INCR', 'stat_T_DECR', 'stat_hyst_DECR', 'max_T', 'min_T', 'change_rate'")
  } else {print("No errors were found in the dataset")}
}


##### (PREPARE) prepare data file for further processing and detect shuttles ####

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
  if (length(midnight_observation) > 0) {
    data$date[data$time_sec >= midnight_observation] <- data$date[data$time_sec >= midnight_observation] + 1
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

##### (CALCULATE) coreT if not already done in Shuttlesoft ####

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

##### (CALCULATE) Temperature preference (using mean or median of coreT values) ####

calc_Tpref <- function(data, method = c("median", "mean", "mode"), exclude_start_minutes = 0, exclude_end_minutes = 0, exclude_acclimation = F, print_results = T) {
  method <- match.arg(method)
  
  # Convert the time to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  if (exclude_acclimation) {
    data <- data[data$trial_phase != "acclimation", ]
  }
  
  # Convert the core_T column to numeric
  data$core_T <- as.numeric(data$core_T)
  
  # Calculate Tpref based on the specified method
  if (method == "median") {
    Tpref <- median(data$core_T, na.rm = TRUE)
  } else if (method == "mean") {
    Tpref <- mean(data$core_T, na.rm = TRUE)
  } else if (method == "mode") {
    Tpref <- as.numeric(names(sort(table(data$core_T), decreasing = TRUE))[1])
  }
  
  if (print_results) {
  print(paste("Tpref:", Tpref))
  }
  return(Tpref)

}

##### (CALCULATE) Avoidance temperature (upper and lower) ####

calc_Tavoid <- function(data, percentiles = c(0.05, 0.95), exclude_start_minutes = 0, exclude_end_minutes = 0, exclude_acclimation = F, print_results = T) {
  # Ensure the percentiles are valid
  if (length(percentiles) != 2 || any(percentiles < 0) || any(percentiles > 1)) {
    stop("Tavoid percentiles should be a vector of two values between 0 and 1")
  }
  
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  if (exclude_acclimation) {
    data <- data[data$trial_phase != "acclimation", ]
  }
  # Calculate the percentiles
  Tavoid_lower <- quantile(data$core_T, percentiles[1], na.rm = TRUE, names = F)
  Tavoid_upper <- quantile(data$core_T, percentiles[2], na.rm = TRUE, names = F)
  
  if (print_results) {
  # Print values
  print(paste("Tavoid Lower:", Tavoid_lower))
  print(paste("Tavoid Upper:", Tavoid_upper)) 
  }
  
  return(c(Tavoid_lower, Tavoid_upper))

}

##### (CALCULATE) Distance moved ####

calc_distance <- function(data, pixel_to_cm = T){
  
  # If coordinates are missing, make them the last known coordinate
  for (i in 2:nrow(data)) {
    if(is.na(data$x_pos[i])) data$x_pos[i]<-data$x_pos[i-1]
    if(is.na(data$y_pos[i])) data$y_pos[i]<-data$y_pos[i-1]
  }
  
  #Initialise the velocity column
  data$velocity<-NA
  
  # Calculate the distance covered in each observation (velocity) using Pythagoras
  for (i in 2:nrow(data)) {
    data$velocity[i] <- sqrt((data$x_pos[i] - data$x_pos[i-1])^2 + (data$y_pos[i]-data$y_pos[i-1])^2)
  }
  
  # Change any NA velocities to zero
  data$velocity[is.na(data$velocity)]<-0
  
  # If the coordinates are specified in pixels rather than cm, convert
  if (pixel_to_cm){
    data$velocity <- data$pixel_ratio*data$velocity
  }
  
  #Initialise the distance column
  data$distance <- 0
  
  # Calculate cumulative distance
  for (i in 2:nrow(data)){
    data$distance[i]<- data$distance[i-1]+data$velocity[i-1]
  }
  
  return(data)

}

calc_tot_distance <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0, exclude_acclimation = F, print_results = T) {
  # Ensure the Distance moved column exists
  if (!"distance" %in% colnames(data)) {
    stop("The dataset does not contain a 'distance' column")
  }
  
  # Convert exclude minutes to seconds
  # Exclude initial and final data if necessary
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - exclude_end_minutes * 60
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  if (exclude_acclimation) {
    data <- data[data$trial_phase != "acclimation", ]
  }
  
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

##### (CALCULATE) Shuttling frequency ####

calc_shuttles <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0, exclude_acclimation = F, print_results = T) {
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
  
  if (exclude_acclimation) {
    data <- data[data$trial_phase != "acclimation", ]
  }
  
  # Calculate shuttling frequency
  shuttles <- sum(data$shuttle, na.rm = TRUE)

  if (print_results) {
  # Print the result
  print(paste("Shuttling frequency:", shuttles))
  }
  
  return(shuttles)
}

##### (CALCULATE) Occupancy time in each chamber ####

calc_occupancy <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0, exclude_acclimation = F, print_results = T) {
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
  
  if (exclude_acclimation) {
    data <- data[data$trial_phase != "acclimation", ]
  }
  
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

##### (CALCULATE) Core temperature variance  #####

calc_coreT_variance <- function(data, variance_type = c("std_error", "std_deviation", "coeff_variation"), exclude_start_minutes = 0, exclude_end_minutes = 0, exclude_acclimation = F) {
  # Ensure the body core temperature column exists
  if (!"core_T" %in% colnames(data)) {
    stop("The dataset does not contain a 'core_T' column. Please run file_prepare and calc_coreT")
  }
  
  # Convert exclude minutes to seconds
  exclude_start_seconds <- exclude_start_minutes * 60
  exclude_end_seconds <- exclude_end_minutes * 60
  
  if (exclude_acclimation) {
    data <- data[data$trial_phase != "acclimation", ]
  }
  
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

##### (CALCULATE) Actual Gravitation time  ####

calc_gravitation <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0, exclude_acclimation = F, print_results = T) {
  
  # Convert the time to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  if (exclude_acclimation) {
    data <- data[data$trial_phase != "acclimation", ]
  }
  
  # Perform segmented regression
  fit <- lm(core_T ~ time_h, data = data)
  seg_fit <- segmented::segmented(fit, seg.Z = ~time_h, npsi = 1)
  
  # Get summary and breakpoint
  seg_summary <- summary(seg_fit)
  breakpoints <- seg_fit$psi[, "Est."]
  
  # Extract breakpoints
  breakpoints <- seg_fit$psi[, "Est."]
  
  if (length(breakpoints) == 0) {
    stop("No breakpoints found in the segmented model.")
  }
  
  # Assume the first breakpoint as the gravitation time
  gravitation_time <- breakpoints[1]-min(data$time_h)
  
  if (print_results) {
  print(paste("Gravitation Time (hours):", gravitation_time))
  }
  
  return(gravitation_time)
}

##### (CALCULATE) Time at extreme temperatures (near set limits) ####

calc_extremes <- function(data, threshold = 0.2*(max(data$max_T)-max(data$min_T)), exclude_start_minutes = 0, exclude_end_minutes = 0, exclude_acclimation = F, print_results = T) {
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
  
  if (exclude_acclimation) {
    data <- data[data$trial_phase != "acclimation", ]
  }
  
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

##### (CALCULATE) Proportion of lost coordinates ####

calc_track_accuracy <- function (data, exclude_start_minutes = 0, exclude_end_minutes = 0, exclude_acclimation = F, print_results = T){
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Exclude acclimation if requested
  if (exclude_acclimation) {
    data <- data[data$trial_phase != "acclimation", ]
  }
  
  observations <- nrow(data)
  missed_coordinates <- sum(is.na(data$x_pos))
  
  prop <- 1-(missed_coordinates/observations)
  
  if (print_results) {
  print(paste0("Proportion of time where subject was tracked: ", prop*100, "%"))
  }
  return(prop)
  
}
  
####### (INSPECT) Temperatures in each side of the shuttlebox over time ####

plot_T_gradient <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Convert time_sec to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Convert time in seconds to hours
  data$time_h <- data$time_sec / 3600
  
  # Plot the temperature gradient
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = time_h)) +
    ggplot2::geom_line(ggplot2::aes(y = `DECR_T`, color = "Cold chamber")) +
    ggplot2::geom_line(ggplot2::aes(y = `INCR_T`, color = "Warm chamber")) +
    ggplot2::scale_color_manual(values = c("Warm chamber" = "brown1", "Cold chamber" = "dodgerblue")) +
    ggplot2::labs(title = "Temperature Gradient Over Time",
         x = "Time (hours)",
         y = "Temperature (°C)",
         color = "Temperature Side") +
    ggplot2::theme_light()
  
  print(plot)
}

####### (INSPECT) Core body temperature over time #####

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
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = time_h, y = core_T)) +
    ggplot2::geom_point(alpha = 0.5, size = 0.4) +
    ggplot2::geom_hline(yintercept = Tpref, linetype = "dashed", color = "#9EF1A4", size = 1.75) +
    ggplot2::annotate("text", x = 0, hjust = 0, y = Tpref + 0.5, label = paste("Tpref:", round(Tpref, 1)), color = "#9EF1A4", size = 4.5)+
    ggplot2::geom_hline(yintercept = Tavoid_lower, linetype = "dashed", color = "#9ECFF1", size = 1.75) +
    ggplot2::annotate("text", x = 0, hjust = 0, y = Tavoid_lower + 0.5, label = paste("Low Tavoid:", round(Tavoid_lower, 1)), color = "#9ECFF1", size = 4.5)+
    ggplot2::geom_hline(yintercept = Tavoid_upper, linetype = "dashed", color = "#F1AD9E", size = 1.75) +
    ggplot2::annotate("text", x = 0, hjust = 0, y = Tavoid_upper + 0.5, label = paste("High Tavoid:", round(Tavoid_upper, 1)), color = "#F1AD9E", size = 4.5)+
    ggplot2::labs(title = "Body Core Temperature vs Time",
         x = "Time (h)",
         y = "Body Core Temperature (°C)") +
    ggplot2::theme_light()
  
  print(plot)
}

####### (INSPECT) Segmented regressions of core body temperature during trial #####

plot_T_segmented <- function(data, Tpref_method = "median", Tavoid_percentiles = c(0.05, 0.95),exclude_start_minutes = 0, exclude_end_minutes = 0, exclude_acclimation = F, overlay_chamber_temp = T) {
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
  
  if (exclude_acclimation) {
  fit <- lm(core_T ~ time_h, data = data[data$trial_phase != "acclimation", ])
  seg_fit <- segmented::segmented(fit, seg.Z = ~time_h, npsi = 1)
  predicted_values <- data.frame(
    time_h = data$time_h[data$trial_phase != "acclimation"],
    predicted = predict(seg_fit))
  } else {
    fit <- lm(core_T ~ time_h, data = data)
    seg_fit <- segmented::segmented(fit, seg.Z = ~time_h, npsi = 1)
    predicted_values <- data.frame(
      time_h = data$time_h,
      predicted = predict(seg_fit))
  }
  
  # Perform segmented regression
  acclimation_time <- min(data$time_h[data$trial_phase != "acclimation"], na.rm = TRUE)
  
  Tpref <- calc_Tpref(data, method = Tpref_method, print_results = F, exclude_start_minutes = exclude_start_minutes, exclude_end_minutes = exclude_end_minutes,  exclude_acclimation = exclude_acclimation)
  Tavoid <- calc_Tavoid(data, percentiles = Tavoid_percentiles, print_results = F, exclude_start_minutes = exclude_start_minutes, exclude_end_minutes = exclude_end_minutes,  exclude_acclimation = exclude_acclimation)
  Tavoid_lower <- Tavoid[1]
  Tavoid_upper <- Tavoid[2]
  
  # Get summary and breakpoint
  seg_summary <- summary(seg_fit)
  breakpoints <- seg_fit$psi[, "Est."]- acclimation_time
  
  # Plot data with segmented regression
  plot <-
    ggplot2::ggplot(data, ggplot2::aes(x = time_h, y = core_T)) +
    ggplot2::geom_point(alpha = 1, size = 0.4) +
    ggplot2::geom_line(data = predicted_values, ggplot2::aes(x = time_h, y = predicted), color = "#F79518", size = 1)+
    ggplot2::annotate("text", x = min(predicted_values$time_h), hjust = 0, y = min(data$core_T), label = paste("Gravitation time:", round(breakpoints, 1), "h"), color = "#F79518", size = 4.5)+
    ggplot2::geom_hline(yintercept = Tpref, linetype = "dashed", color = "#03DB6A", size = 1.75) +
    ggplot2::annotate("text", x = min(data$time_h), hjust = 0, y = Tpref + 0.4, label = paste("Preference:", round(Tpref, 1)), color = "#03DB6A", size = 4.5)+
    ggplot2::geom_hline(yintercept = Tavoid_lower, linetype = "dashed", color = "#9ECFF1", size = 1.75) +
    ggplot2::annotate("text", x = min(data$time_h) + exclude_start_minutes/60, hjust = 0, y = Tavoid_lower + 0.4, label = paste("Lower avoidance:", round(Tavoid_lower, 1)), color = "#9ECFF1", size = 4.5)+
    ggplot2::geom_hline(yintercept = Tavoid_upper, linetype = "dashed", color = "#F1AD9E", size = 1.75) +
    ggplot2::annotate("text", x = min(data$time_h) + exclude_start_minutes/60, hjust = 0, y = Tavoid_upper + 0.4, label = paste("Upper avoidance:", round(Tavoid_upper, 1)), color = "#F1AD9E", size = 4.5)+
    ggplot2::labs(title = "Body Core Temperature vs Time with Segmented Regression",
         x = "Time (h)",
         y = "Body Core Temperature (°C)")+
    ggplot2::theme_light()
  
  if(overlay_chamber_temp){
    plot <- plot + 
      ggplot2::geom_line(ggplot2::aes(y = `DECR_T`, color = "DECR side temp")) +
      ggplot2::geom_line(ggplot2::aes(y = `INCR_T`, color = "INCR side temp")) +
      ggplot2::scale_color_manual(values = c("DECR side temp" = "dodgerblue", "INCR side temp" = "brown1"), guide="none")
  }
  
  print(list(summary = seg_summary, breakpoints = breakpoints))
  
  return(plot)
  
}

####### (INSPECT) Histogram of time spent at different core body temperatures ####

plot_coreT_histogram <- function(data, bin_size = 0.1, exclude_start_minutes = 0, exclude_end_minutes = 0, exclude_acclimation = F) {
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
  
  
  if (exclude_acclimation) {
    data <- data[data$trial_phase != "acclimation", ]
  }
  
  
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
  hist_plot <- ggplot2::ggplot(data, ggplot2::aes(x = core_T)) +
    ggplot2::geom_histogram(binwidth = bin_size, ggplot2::aes(y = (ggplot2::after_stat(count) / sum(ggplot2::after_stat(count))) * 100), fill = "#F79518", color = "black") +
    ggplot2::geom_vline(xintercept = Tpref_peak, color = "#9EF1A4", linetype = "dashed", linewidth = 2) +
    ggplot2::labs(title = "Histogram of Percent of Time Spent at Different Body Core Temperatures",
         subtitle = paste("Peak Tpref:", round(Tpref_peak, 2), "°C"),
         x = "Body Core Temperature (°C)",
         y = "Percent of Total Time (%)") +
    ggplot2::theme_light()
  
  print(hist_plot)
  print(paste("Peak Tpref:", Tpref_peak))
  # return(Tpref_peak)
}

####### (INSPECT) Histogram of selected variable ####

plot_histogram <- function(data, column, binwidth = 0.1, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary column exists
  if (!paste(column) %in% colnames(data)) {
    stop(paste0("The dataset does not contain a '",column,"' column"))
  }
  
  # Convert time_sec to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec, na.rm = TRUE) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Create the histogram plot
  plot <- ggplot2::ggplot(data, ggplot2::aes_string(x = column)) +
    ggplot2::geom_histogram(binwidth = binwidth, fill = "#F79518", color = "black") +
    ggplot2::labs(title = paste("Histogram of", column),
         x = column,
         y = "Frequency") +
    ggplot2::theme_light()
  
  print(plot)
}

####### (INSPECT) Cumulative distance changes over time ####

plot_distance <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0) {
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
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = Time_min, y = (distance))) +
      ggplot2::geom_line(color = 'orange') +
      ggplot2::labs(title = "Cumulative Distance Moved During Trial",
           x = "Time (minutes)",
           y = "Distance Moved (cm)") +
      ggplot2::theme_light()
    
    print(plot)
  }

####### (INSPECT) Movement speed vs body core temperature ####

plot_speed_coreT <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0) {
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
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = core_T, y = Movement_Speed)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_smooth(method = "loess", formula = 'y ~ x', color = "#F79518", se = FALSE) +
    ggplot2::labs(title = "Movement Speed vs. Core Body Temperature",
         x = "Core Body Temperature (°C)",
         y = "Movement Speed (cm/s)") +
    ggplot2::theme_light()
  
  return(plot)
}

####### (INSPECT) Interval plot of selected variable #####

plot_interval <- function(data, column, interval_minutes = 10, exclude_start_minutes = 0, exclude_end_minutes = 0) {

    # Convert time_sec to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec, na.rm = T) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Calculate time intervals
  data$t_interval <- floor((data$time_sec - start_time) / (interval_minutes * 60)) * interval_minutes
  
  # Calculate mean rate for each interval
  rate <- aggregate(data[[column]] ~ t_interval, data = data, FUN = function(x) mean(x, na.rm = TRUE))
  names(rate)[names(rate) == "data[[column]]"] <- column
    
 rate$t_interval<-rate$t_interval/60
 
 if (column == "shuttle") {rate[[column]] <- rate[[column]]*60}
 
  # Convert Time_interval to numeric for plotting
  rate$t_interval <- as.numeric(rate$t_interval)
  
  # Plot the mean shuttling rate
  plot <- ggplot2::ggplot(rate, ggplot2::aes_string(x = "t_interval", y = column)) +
    ggplot2::geom_line(color = "black") +
    ggplot2::geom_point(color = "#F79518") +
    ggplot2::labs(title = paste("Mean ", column,  "over the trial"),
         x = "Time (h)",
         y = paste("Mean ", column))+
    ggplot2::theme_light()
  
  print(rate)
  
  return(plot)
}
####### (INSPECT) Missed tracks over trial ####

plot_tracking <- function(data, interval_minutes = 60, exclude_start_minutes = 0, exclude_end_minutes = 0, exclude_acclimation = F) {
  
  # Ensure necessary columns exist
  if (!all(c("time_h", "x_pos", "y_pos") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'time_h', 'x_pos', 'y_pos'.")
  }
  
  # Convert the time to numeric if not already
  data$time_sec <- as.numeric(data$time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$time_sec) - (exclude_end_minutes * 60)
  data <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
  
  # Exclude acclimation if requested
  if (exclude_acclimation) {
    data <- data[data$trial_phase != "acclimation", ]
  }
  
  data$mistrack <- NA
  
  for (i in seq_len(nrow(data))) {
    if (is.na(data$x_pos[i])) {
      data$mistrack[i] <- 1
    } else {
      data$mistrack[i] <- 0
    }
  }
  
  # Calculate time intervals
  data$t_interval <- floor((data$time_sec - start_time) / (interval_minutes * 60)) * interval_minutes
  
  # Calculate mean rate for each interval
  rate <- aggregate(data$mistrack ~ t_interval, data = data, FUN = function(x) sum(x, na.rm = TRUE))
  names(rate)[names(rate) == "data$mistrack"] <- "mistrack"
  
  rate$t_interval<-rate$t_interval/60
  
  # Convert Time_interval to numeric for plotting
  rate$t_interval <- as.numeric(rate$t_interval)
  
  # Plot the mean shuttling rate
  plot <- ggplot2::ggplot(rate, ggplot2::aes(x = t_interval, y = mistrack)) +
    ggplot2::geom_line(color = "black") +
    ggplot2::geom_point(color = "#F79518") +
    ggplot2::labs(title = "Number of missed tracks over the trial",
         x = "Time (h)",
         y = paste("Missed tracks/h "))+
    ggplot2::theme_light()
  
  print(rate)
  
  return(plot)
  
}

####### (INSPECT) Produce a correlation matrix ####

correlation_matrix <- function(data, columns, exclude_start_minutes = 0, exclude_end_minutes = 0, interval_minutes = 10) {
    
    # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
    start_time <- exclude_start_minutes * 60
    end_time <- max(data$time_sec, na.rm = TRUE) - (exclude_end_minutes * 60)
    data_2 <- data[data$time_sec >= start_time & data$time_sec <= end_time, ]
    
    # Create time intervals
    data_2$t_interval <- floor((data_2$time_sec - start_time) / (interval_minutes * 60)) * interval_minutes
    data
    # Compute mean values for selected columns within each interval
    mean_data_2 <- aggregate(data_2[, columns, drop = FALSE], 
                             by = list(t_interval = data_2$t_interval), 
                             FUN = function(x) mean(x, na.rm = TRUE))
    
    # Function to display correlation values
    panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor = 1.2, ...) {
      r <- cor(x, y, use = "complete.obs")  # Compute correlation
      txt <- formatC(r, format = "f", digits = digits)  # Format value
      txt <- paste0(prefix, txt)  # Append prefix if needed
      
      # Place text at the center of the panel
      text(mean(range(x)), mean(range(y)), txt, cex = cex.cor)
    }
    
    # Generate correlation matrix plot
    pairs(
      mean_data[, columns, drop = FALSE],  # Use mean-aggregated data
      lower.panel = panel.smooth,  # Scatter plots with smoothing
      upper.panel = panel.cor  # Correlation values
    )
  }

####### (INSPECT) 3D animation of fish movements at selected time intervals, with time as a dimension ####


animate_movements <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0, speed = 100) {
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
  rgl::plot3d(data$x_pos, data$y_pos, data$time_sec, type = "l", col = "darkgrey", alpha = 0.15,
         xlab = "X Position", ylab = "Y Position", zlab = "Time (seconds)")
  
  # Determine color scaling based on Tpref
  colors <- colorRampPalette(c("royalblue3", "red2"))(length(unique(data$core_T)))
  data$color <- colors[as.numeric(cut(data$core_T, breaks = length(colors)))]
  
  # # Animate the plot
  # for (i in seq(1, nrow(data))) {  # Update every 10 points
  #   rgl::points3d(data$x_pos[i], data$y_pos[i], data$time_sec[i], col = data$color[i], size = 5)
  # }
  # 
  # rgl::play3d(spin3d(axis = c(0, 0, 1), rpm = 100), duration = 5)
  
  run_speed = 10^speed^speed
  
  # Function to update the animation
  animate_step <- function(time) {
    index <- floor(time * run_speed) + 1  # Adjust speed
    if (index > nrow(data)) return()
    rgl::points3d(data$x_pos[index], data$y_pos[index], data$time_sec[index], col = data$color[index], size = 5)
  }
  
  # Use play3d with a time sequence
  rgl::play3d({
    for (i in seq_len(nrow(data))) {
      rgl::points3d(data$x_pos[i], data$y_pos[i], data$time_sec[i], col = data$color[i], size = 5)
      Sys.sleep(1 / (run_speed * 100))  # Adjust speed dynamically
    }
  }, duration = nrow(data) / (run_speed * 10))
}

####### (INSPECT) Heat map of fish locations over the trial ####

plot_heatmap <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0, exclude_acclimation = F) {
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
  
  if (exclude_acclimation) {
    data <- data[data$trial_phase != "acclimation", ]
  }
  
  # Get min and max x_pos for each temperature side
  cold_min <- min(data$x_pos[data$zone == "DECR"], na.rm = TRUE)
  cold_max <- max(data$x_pos[data$zone == "DECR"], na.rm = TRUE)
  warm_min <- min(data$x_pos[data$zone == "INCR"], na.rm = TRUE)
  warm_max <- max(data$x_pos[data$zone == "INCR"], na.rm = TRUE)
  
  # Determine the border between either zone
  split_position <- mean(c(min(cold_max, warm_max), max(cold_min, warm_min)))
  
  # Plot heatmap
  heatmap_plot <- ggplot2::ggplot(data, ggplot2::aes(x = x_pos, y = y_pos)) +
    ggplot2::stat_density2d(ggplot2::aes(fill = ggplot2::after_stat(density)), geom = "raster", contour = FALSE) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::geom_vline(xintercept = split_position, color = "white", linetype = "dashed", size = 1) +
    ggplot2::annotate("text", x = cold_min + (cold_max-cold_min)/4, y = min(data$y_pos, na.rm = TRUE)+20, 
             label = "Cold Chamber", color = "dodgerblue", size = 5, hjust = 0) +
    ggplot2::annotate("text", x = warm_min + (warm_max-warm_min)/4, y = min(data$y_pos, na.rm = TRUE)+20, 
             label = "Warm Chamber", color = "brown1", size = 5, hjust = 0) +
    ggplot2::labs(title = "Heat Map of Fish Locations within the Shuttlebox",
         x = "X Position",
         y = "Y Position") +
    ggplot2::theme_light() +
    ggplot2::coord_fixed()
  
  return(heatmap_plot)
}

####### (PREPARE) Read shuttlesoft project from directory ####

read_shuttlesoft_project <- function (metadata, directory = getwd()){
  
  # Check if metadata is specified and give warnings for what happens if data is missing
  if (missing(metadata)) {
    warning("Metadata not provided; cannot define one or more of the following variables: 'trial_start', 'a-value', 'b-value', 'initial_temp.'")
  } else {
    
    # If the metadata is present, ensure necessary columns exist in the metadata
    if (!all(c("file_name", "trial_start", "mass", "a_value", "b_value") %in% colnames(metadata))) {
      warning("The metadata does not contain one or more necessary columns: 'file_name', 'trial_start', 'mass', 'a_value', 'b_value'.")
    }
    
    # Find rows with NA in specified columns
    columns_to_check <- intersect(c("file_name", "trial_start", "mass", "a_value", "b_value"), colnames(metadata))
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
  
  # Check if all files are represented in the metadata
  file_names <- basename(txt_files)
  missing_files <- setdiff(file_names, metadata$file_name)
  if (length(missing_files) > 0) {
    stop(paste("The following files are not represented in metadata:", 
               paste(missing_files, collapse = "\n ")))
  }
  
  # Read those text files and compile them into a list
  data_read <- lapply(txt_files, read_shuttlesoft, metadata = metadata, multidat = T)
  
  # Name each list element as the directory of the text files
  names(data_read)<-txt_files
  
  return(data_read)
}

####### (PREPARE) Compile project data into single dataframe ####

compile_project_data<-function(data_read){
  
  # Prepare the data using file_prepare. Multidat is TRUE because this will suppress warnings from preparing each data file separately
  data_list <- lapply(data_read, file_prepare)
  processed_list <- lapply(data_list, calc_coreT)
  combined_df <- do.call(rbind, processed_list)
  rownames(combined_df)<-NULL
  return(combined_df)
}

##### (CALCULATE) Calculate variables of project read object ####

calc_project_results <- function(data_read,
                                 calculate_distance = F,
                                 pixel_to_cm = T,
                                 exclude_acclimation = F,
                                 exclude_start_minutes = 0, 
                                 exclude_end_minutes = 0,
                                 Tpref_method = "median",
                                 Tavoid_percentiles = c(0.05, 0.95),
                                 textremes_threshold = expression(0.2 * (max(df$max_T) - max(df$min_T))),
                                 core_T_variance_type = "std_error"){
  message("Initialising...")
  
  # Prepare the data using file_prepare. Multidat is TRUE because this will suppress warnings from preparing each data file separately
  data_list <- lapply(data_read, file_prepare)
  
  
  # Create a function that pulls all the information from a datafile in one go, using the functions specified earlier 
  apply_functions <- function(df, df_name) {
    
    # Define the fileID the function is currently considering
    fileID <- unique(df$fileID)
    
    # Calculate core body temperature
    df<-calc_coreT(df)
    
    # Calculate distance covered if true
    if (calculate_distance) df<-calc_distance (df, pixel_to_cm)
    
    # Exclude acclimation period if requested from the start. Acclimation period is defined in file_prepare
    if(exclude_acclimation) {
      df <- subset(df, df$trial_phase != "acclimation") 
    }
    
    # Subset custom time window
    df <- df[df$time_sec >= exclude_start_minutes * 60 & df$time_sec <= max(df$time_sec) - (exclude_end_minutes * 60), ]
    
    
    textremes_th <- eval(textremes_threshold) 
    
    # Calculate all variables using default settings (these can be edited in the main function)
    Tpref <- calc_Tpref(df, method = Tpref_method, print_results = F)
    grav_time <- calc_gravitation(df, print_results = F)
    tot_distance <- calc_tot_distance(df, print_results = F)
    Tavoid <- calc_Tavoid(df, percentiles = Tavoid_percentiles, print_results = F)
    Tpref_range <- Tavoid[2] - Tavoid[1]
    textremes <- calc_extremes(df, threshold = textremes_th, print_results = F)
    textremes_tot <- textremes[1]+textremes[2]
    track_accuracy <- calc_track_accuracy(df, print_results = F)
    core_T_variance <- calc_coreT_variance(df, variance_type = core_T_variance_type)
    nr_shuttles <- calc_shuttles(df, print_results = F)
    chamber_seconds <- calc_occupancy(df, print_results = F)
    
    # Return a list of the variables you want
    return(list(fileID = fileID,
                Tpref = Tpref,
                Tpref_range = Tpref_range,
                grav_time = grav_time, 
                tot_distance = tot_distance,
                Tavoid = Tavoid,
                textremes = textremes,
                textremes_tot = textremes_tot,
                track_accuracy = track_accuracy,
                core_T_variance = core_T_variance,
                nr_shuttles = nr_shuttles,
                chamber_seconds = chamber_seconds))
  }
  
  # Apply functions to each dataframe in the list and collect results
  results <- do.call(rbind, lapply(names(data_list), function(name, total) {
    
    current_iter <- which(names(data_list) == name)  
    message("Processing dataset ", current_iter, " of ", total, " (", round((current_iter / total) * 100, 2), "% complete)")  
    
    # Apply the functions
    func_results <- apply_functions(df = data_list[[name]], df_name = name)
    
    # Combine results into a data frame
    data.frame(
      fileID = func_results$fileID,
      Tpref = func_results$Tpref,
      Tpref_range = func_results$Tpref_range,
      grav_time  = func_results$grav_time,
      tot_distance  = func_results$tot_distance,
      Tavoid_lower  = func_results$Tavoid[1],
      Tavoid_upper = func_results$Tavoid[2],
      t_near_min  = func_results$textremes[1],
      t_near_max  = func_results$textremes[2],
      t_near_limits = func_results$textremes_tot,
      track_accuracy = func_results$track_accuracy,
      core_T_variance  = func_results$core_T_variance,
      nr_shuttles  = func_results$nr_shuttles,
      seconds_in_DECR = func_results$chamber_seconds[1],
      seconds_in_INCR = func_results$chamber_seconds[2]
      
    )
    
  }, total = length(data_list)))
  
  # Reset rownames as they have acquired names from the list
  rownames(results) <- NULL
  
  return(results)
}

####### (INSPECT) Link shuttling and activity across individuals in the entire project dataset ####


plot_distance_vs_shuttles <- function(proj_data) {
  # Ensure necessary columns exist
  if (!all(c("fileID", "distance", "nr_shuttles") %in% colnames(proj_data))) {
    stop("The dataset does not contain one or more necessary columns: 'fileID', 'distance', 'nr_shuttles'.")
  }
  
  # Create the plot
  plot <- ggplot2::ggplot(proj_data, ggplot2::aes(x = distance, y = nr_shuttles, label = fileID)) +
    ggplot2::geom_point(color = "#F79518") +
    ggplot2::geom_text(vjust = -1, hjust = 1.5) +
    ggplot2::labs(title = "Relationship between Distance Moved and Shuttles across Individuals",
         x = "Distance Moved (cm)",
         y = "Number of Shuttles") +
    ggplot2::theme_light()
  
  print(plot)
}

####### (INSPECT) Link time near system temperature limits and activity across individuals in the entire project dataset ####


plot_limits_vs_distance <- function(proj_data) {
  # Ensure necessary columns exist
  if (!all(c("fileID", "distance", "t_near_limits") %in% colnames(proj_data))) {
    stop("The dataset does not contain one or more necessary columns: 'fileID', 'distance', 't_near_limits'.")
  }
  
  # Create the plot
  plot <- ggplot2::ggplot(proj_data, ggplot2::aes(x = t_near_limits, y = distance, label = fileID)) +
    ggplot2::geom_point(color = "#F79518") +
    ggplot2::geom_text(vjust = -1, hjust = 1.5) +
    ggplot2::labs(title = "Relationship between Time Near Limits and Distance Moved",
         x = "Time Near Limits (minutes)",
         y = "Distance Moved (cm)") +
    ggplot2::theme_light()
  
  print(plot)
}

####### (INSPECT) Link time near system temperature limits and shuttles across individuals in the entire project dataset ####


plot_limits_vs_shuttles <- function(proj_data, id_col) {
  # Ensure necessary columns exist
  if (!all(c("fileID", "t_near_limits", "nr_shuttles") %in% colnames(proj_data))) {
    stop("The dataset does not contain one or more necessary columns: 'fileID', t_near_limits', 'nr_shuttles'.")
  }
  
  # Create the plot
  plot <- ggplot2::ggplot(proj_data, ggplot2::aes(x = t_near_limits, y = nr_shuttles, label = fileID)) +
    ggplot2::geom_point(color = "#F79518") +
    ggplot2::geom_text(vjust = -1, hjust = 1.5) +
    ggplot2::labs(title = "Relationship between Time Near Limits and Number of shuttles",
         x = "Time Near Limits (minutes)",
         y = "Shuttles") +
    ggplot2::theme_light()
  
  print(plot)
}

####### (INSPECT) PCA with shuttles, distance moved, and time near limits, then plot PC scores vs Tpref ####

pca <- function(data, mahalanobis_th = 0.7, dbscan_th = 1, print_labels = TRUE, id_col = "fileID", var_col = "Tpref", biplot_variables = T) {
  
  # Change the id column into row names 
  # Remove all non-numeric columns and perform a PCA
  pca_data <- data
  rownames(pca_data) <- pca_data[[id_col]]
  pca_data[[id_col]] <- NULL
  pca_data <- pca_data[, sapply(pca_data, is.numeric)]
  
  pca <- FactoMineR::PCA(pca_data, scale.unit = TRUE, graph = FALSE)
  pca_scores <- pca$ind$coord
  
  # Compute the Mahalanobis distance
  center <- colMeans(pca_scores)  # Mean of the PCA scores
  cov_matrix <- cov(pca_scores)  # Covariance matrix of PCA scores
  mahalanobis_dist <- mahalanobis(pca_scores, center, cov_matrix)
  mahalanobis_dist_df<-as.data.frame(mahalanobis_dist)
  
  threshold <- qchisq(mahalanobis_th, df = ncol(pca_scores))  # Chi-squared threshold
  outliers_mahalanobis <- rownames(mahalanobis_dist_df)[mahalanobis_dist_df$mahalanobis_dist > threshold]
  
  dbscan_result <- dbscan::dbscan(pca_scores[, 1:2], eps = dbscan_th, minPts = 5)  # Using PC1 and PC2 for DBSCAN
  outliers_dbscan <- rownames(data)[dbscan_result$cluster == 0]  # Points with cluster = 0 are outliers
  
  combined_outliers <- matrix(NA, nrow = 2, ncol = max(length(outliers_mahalanobis), length(outliers_dbscan)))
  combined_outliers[1, 1:length(outliers_mahalanobis)] <- outliers_mahalanobis
  combined_outliers[2, 1:length(outliers_dbscan)] <- outliers_dbscan
  combined_outliers <- as.data.frame(combined_outliers)
  rownames(combined_outliers) <- c("Mahalanobis Outliers", "DBSCAN Outliers")
  
  p1 <- factoextra::fviz_screeplot(pca, addlabels = print_labels, main = "Scree Plot")
  if (print_labels) {
    label_content <- "all"
  }else {label_content <- "var"
  }
  
  p2 <- factoextra::fviz_pca_biplot(pca, label = label_content,
                        repel = TRUE,  # Avoid overlapping labels
                        col.var = if (biplot_variables) "blue" else NA,  , # Color for variables
                        col.ind = "red",  # Color for individuals
                        title = "PCA Biplot")
  
  p3 <- factoextra::fviz_contrib(pca, choice = "var", axes = 1)+
    ggplot2::ggtitle("Contribution of variables to PC1")
  
  plot_dat <- data.frame(id = data[[id_col]], PC1 = pca$ind$coord[, 1], var_y = data[[var_col]])
  if(print_labels == F){
    plot_dat$id <- ""
  }
  plot_dat[!is.na(plot_dat$id), ]
  
  p4 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = PC1, y = var_y)) +
    ggplot2::geom_point(color = "#F79518") +
    ggplot2::geom_smooth(method = "lm", color = "black") +  # Adding 95% CI shading
    ggrepel::geom_text_repel(ggplot2::aes(label = id), vjust = -1, hjust = 1.5) +
    ggplot2::labs(title = paste0("PCA Scores from the First Component vs ",var_col," with 95% CI"),
         x = "PCA First Component Score",
         y = paste0(var_col)) +
    ggplot2::theme_light()
  
  pca_loadings <- pca$var$coord
  
  output<-list(pca = pca,
               pca_loadings = pca_loadings,
               pca_scores = pca_scores,
               outliers = combined_outliers, 
               plots = list (screeplot = p1, biplot = p2, pc1contributionplot = p3, varplot = p4))
  
  return(output)
  
}

####### (INSPECT) Frequency distributions for key measures across individuals in the data set ####

plot_histograms <- function(proj_data, bin_size_Tpref = 1, bin_size_Tavoid_upper = 1, bin_size_Tavoid_lower = 1, bin_size_distance = 2000, bin_size_shuttles = 20, bin_size_t_near_limits = 5, nr_sd = 2, iqr_multiplier = 1.5, id_col = "fileID") {
  # Ensure necessary columns exist
  required_columns <- c("Tpref", "Tavoid_upper", "Tavoid_lower", "tot_distance", "nr_shuttles", "t_near_limits")
  if (!all(required_columns %in% colnames(proj_data))) {
    stop("The dataset does not contain one or more necessary columns: 'Tpref', 'Tavoid_upper', 'Tavoid_lower', 'tot_distance', 'nr_shuttles', 't_near_limits'.")
  }
  
  # Create individual plots
  p1 <- ggplot2::ggplot(proj_data, ggplot2::aes(x = Tpref)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(count)),
                   binwidth = bin_size_Tpref, 
                   fill = "#F79518", 
                   color = "black") +
    ggplot2::stat_function(fun = function(x) dnorm(x, 
                                          mean = mean(proj_data$Tpref, na.rm = TRUE), 
                                          sd = sd(proj_data$Tpref, na.rm = TRUE)) * length(proj_data$Tpref) * bin_size_Tpref,
                  linewidth = 1, 
                  color = 'darkorange2') +
    ggplot2::labs(x = "Temperature preference (°C)", 
         y = "Frequency") +
    ggplot2::theme_light()
  
  # p2: Distribution of Tavoid_upper
  p2 <- ggplot2::ggplot(proj_data, ggplot2::aes(x = Tavoid_upper)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(count)),
                   binwidth = bin_size_Tavoid_upper, 
                   fill = "#E41A1C", 
                   color = "black") +
    ggplot2::stat_function(fun = function(x) dnorm(x, 
                                          mean = mean(proj_data$Tavoid_upper, na.rm = TRUE), 
                                          sd = sd(proj_data$Tavoid_upper, na.rm = TRUE)) * length(proj_data$Tavoid_upper) * bin_size_Tavoid_upper,
                  linewidth = 1, 
                  color = 'darkred') +  
    ggplot2::labs(x = "Upper avoidance temperature (°C)", y = "Frequency") +
    ggplot2::theme_light()
  
  # p3: Distribution of Tavoid_lower
  p3 <- ggplot2::ggplot(proj_data, ggplot2::aes(x = Tavoid_lower)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(count)),
                   binwidth = bin_size_Tavoid_lower, 
                   fill = "#377EB8", 
                   color = "black") +
    ggplot2::stat_function(fun = function(x) dnorm(x, 
                                          mean = mean(proj_data$Tavoid_lower, na.rm = TRUE), 
                                          sd = sd(proj_data$Tavoid_lower, na.rm = TRUE)) * length(proj_data$Tavoid_lower) * bin_size_Tavoid_lower,
                  linewidth = 1, 
                  color = 'darkblue') +  
    ggplot2::labs(x = "Lower avoidance temperature (°C)", y = "Frequency") +
    ggplot2::theme_light()
  
  # p4: Distribution of Distance
  p4 <- ggplot2::ggplot(proj_data, ggplot2::aes(x = tot_distance)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(count)),
                   binwidth = bin_size_distance, 
                   fill = "#4DAF4A", 
                   color = "black") +
    ggplot2::stat_function(fun = function(x) dnorm(x, 
                                          mean = mean(proj_data$tot_distance, na.rm = TRUE), 
                                          sd = sd(proj_data$tot_distance, na.rm = TRUE)) * length(proj_data$tot_distance) * bin_size_distance,
                  linewidth = 1, 
                  color = 'darkgreen') +  
    ggplot2::labs(x = "Distance", y = "Frequency") +
    ggplot2::theme_light()
  
  # p5: Distribution of Shuttles
  p5 <- ggplot2::ggplot(proj_data, ggplot2::aes(x = nr_shuttles)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(count)),
                   binwidth = bin_size_shuttles, 
                   fill = "#984EA3", 
                   color = "black") +
    ggplot2::stat_function(fun = function(x) dnorm(x, 
                                          mean = mean(proj_data$nr_shuttles, na.rm = TRUE), 
                                          sd = sd(proj_data$nr_shuttles, na.rm = TRUE)) * length(proj_data$nr_shuttles) * bin_size_shuttles,
                  linewidth = 1, 
                  color = 'purple') +  
    ggplot2::labs(x = "Nr shuttles", y = "Frequency") +
    ggplot2::theme_light()
  
  # p6: Distribution of Time Near Limits
  p6 <- ggplot2::ggplot(proj_data, ggplot2::aes(x = t_near_limits)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(count)),
                   binwidth = bin_size_t_near_limits, 
                   fill = "#EEDD82", 
                   color = "black") +
    ggplot2::stat_function(fun = function(x) dexp(x,
                                         rate = 1 / mean(proj_data$t_near_limits, na.rm = TRUE)) *
                    length(proj_data$t_near_limits) * bin_size_t_near_limits,
                  linewidth = 1, 
                  color = 'orange') +  
    ggplot2::labs(x = "Time near limits", y = "Frequency") +
    ggplot2::theme_light()
  
  variables <- c("Tpref", "Tavoid_upper", "Tavoid_lower", "tot_distance", "nr_shuttles", "t_near_limits")

  # Function to detect outliers based on Z-score for selected variables and return a dataframe with outlier ids
  detect_outliers <- function(data, variables) {
    variables <- c("Tpref", "Tavoid_upper", "Tavoid_lower", "tot_distance", "nr_shuttles", "t_near_limits")
    outlier_df <- data.frame(matrix(NA, nrow = nrow(data), ncol = length(variables)))  # Initialize dataframe with NA
    # Set column names of the dataframe to the variable names
    colnames(outlier_df) <- variables
    
    # Loop through the specified variables
    for (var in variables) {
      
      if (var == "t_near_limits") {  # Use IQR for t_near_limits
        Q1 <- quantile(data[[var]], 0.25, na.rm = TRUE)
        Q3 <- quantile(data[[var]], 0.75, na.rm = TRUE)
        IQR_value <- Q3 - Q1
        upper_bound <- Q3 + iqr_multiplier * IQR_value  # Ignore lower bound for exponential data
        
        outlier_indices <- which(data[[var]] > upper_bound)
        
      } else {
        # Calculate the Z-scores for the variable
        z_scores <- (data[[var]] - mean(data[[var]], na.rm = TRUE)) / sd(data[[var]], na.rm = TRUE)
        # Identify outliers
        outlier_indices <- which(abs(z_scores) > nr_sd)
      }
      
      # Place outlier ids in the respective columns, keeping others as NA
      outlier_df[outlier_indices, var] <- data[[id_col]][outlier_indices]
      
    }
    
    return(outlier_df)
  }
  outlier_df<-detect_outliers(proj_data, variables)
  
  names_Tpref <- proj_data[!is.na(proj_data$Tpref) & proj_data[[id_col]] %in% outlier_df$Tpref, c(id_col, "Tpref")]
  colnames(names_Tpref) <- c("id","Tpref")
  plot_build <- ggplot2::ggplot_build(p1)
  y_min <- plot_build$layout$panel_params[[1]]$y.range[1]
  y_max <- plot_build$layout$panel_params[[1]]$y.range[2]
  y_mid <- (y_min + y_max) / 2  
  
  p1 <- p1 + 
    ggplot2::geom_vline(data = names_Tpref, ggplot2::aes(xintercept = Tpref), color = "blue", linetype = "dashed") +
    ggplot2::geom_text(data = names_Tpref, ggplot2::aes(x = Tpref, y = y_mid, label = id), angle = 90, vjust = -1, color = "blue", position = ggplot2::position_jitter(width = 0, height = 3))
  
  names_Tavoid_upper <- proj_data[!is.na(proj_data$Tavoid_upper) & proj_data[[id_col]] %in% outlier_df$Tavoid_upper, c(id_col, "Tavoid_upper")]
  colnames(names_Tavoid_upper) <- c("id","Tavoid_upper")
  plot_build <- ggplot2::ggplot_build(p2)
  y_min <- plot_build$layout$panel_params[[1]]$y.range[1]
  y_max <- plot_build$layout$panel_params[[1]]$y.range[2]
  y_mid <- (y_min + y_max) / 2
  
  p2 <- p2 + 
    ggplot2::geom_vline(data = names_Tavoid_upper, ggplot2::aes(xintercept = Tavoid_upper), color = "blue", linetype = "dashed") +
    ggplot2::geom_text(data = names_Tavoid_upper, ggplot2::aes(x = Tavoid_upper, y = y_mid, label = id), angle = 90, vjust = -1, color = "blue", position = ggplot2::position_jitter(width = 0, height = 3))
  
  names_Tavoid_lower <- proj_data[!is.na(proj_data$Tavoid_lower) & proj_data[[id_col]] %in% outlier_df$Tavoid_lower, c(id_col, "Tavoid_lower")]
  colnames(names_Tavoid_lower) <- c("id","Tavoid_lower")
  plot_build <- ggplot2::ggplot_build(p2)
  y_min <- plot_build$layout$panel_params[[1]]$y.range[1]
  y_max <- plot_build$layout$panel_params[[1]]$y.range[2]
  y_mid <- (y_min + y_max) / 2
  
  p3 <- p3 + 
    ggplot2::geom_vline(data = names_Tavoid_lower, ggplot2::aes(xintercept = Tavoid_lower), color = "blue", linetype = "dashed") +
    ggplot2::geom_text(data = names_Tavoid_lower, ggplot2::aes(x = Tavoid_lower, y = y_mid, label = id), angle = 90, vjust = -1, color = "blue", position = ggplot2::position_jitter(width = 0, height = 3))
  
  names_distance <- proj_data[!is.na(proj_data$tot_distance) & proj_data[[id_col]] %in% outlier_df$tot_distance, c(id_col, "tot_distance")]
  colnames(names_distance) <- c("id","tot_distance")
  plot_build <- ggplot2::ggplot_build(p2)
  y_min <- plot_build$layout$panel_params[[1]]$y.range[1]
  y_max <- plot_build$layout$panel_params[[1]]$y.range[2]
  y_mid <- (y_min + y_max) / 2
  
  p4 <- p4 + 
    ggplot2::geom_vline(data = names_distance, ggplot2::aes(xintercept = tot_distance), color = "blue", linetype = "dashed") +
    ggplot2::geom_text(data = names_distance, ggplot2::aes(x = tot_distance, y = y_mid, label = id), angle = 90, vjust = -1, color = "blue", position = ggplot2::position_jitter(width = 0, height = 3))
 
  names_nr_shuttles <- proj_data[!is.na(proj_data$nr_shuttles) & proj_data[[id_col]] %in% outlier_df$nr_shuttles, c(id_col, "nr_shuttles")]
  colnames(names_nr_shuttles) <- c("id","nr_shuttles")
  plot_build <- ggplot2::ggplot_build(p2)
  y_min <- plot_build$layout$panel_params[[1]]$y.range[1]
  y_max <- plot_build$layout$panel_params[[1]]$y.range[2]
  y_mid <- (y_min + y_max) / 2
  
  p5 <- p5 + 
    ggplot2::geom_vline(data = names_nr_shuttles, ggplot2::aes(xintercept = nr_shuttles), color = "blue", linetype = "dashed") +
    ggplot2::geom_text(data = names_nr_shuttles, ggplot2::aes(x = nr_shuttles, y = y_mid, label = id), angle = 90, vjust = -1, color = "blue", position = ggplot2::position_jitter(width = 0, height = 3))
  
  names_t_near_limits <- proj_data[!is.na(proj_data$t_near_limits) & proj_data[[id_col]] %in% outlier_df$t_near_limits, c(id_col, "t_near_limits")]
  colnames(names_t_near_limits) <- c("id","t_near_limits")
  plot_build <- ggplot2::ggplot_build(p2)
  y_min <- plot_build$layout$panel_params[[1]]$y.range[1]
  y_max <- plot_build$layout$panel_params[[1]]$y.range[2]
  y_mid <- (y_min + y_max) / 2
  
  p6 <- p6 + 
    ggplot2::geom_vline(data = names_t_near_limits, ggplot2::aes(xintercept = t_near_limits), color = "blue", linetype = "dashed") +
    ggplot2::geom_text(data = names_t_near_limits, ggplot2::aes(x = t_near_limits, y = y_mid, label = id), angle = 90, vjust = -1, color = "blue", position = ggplot2::position_jitter(width = 0, height = 3))
  
  # Combine the plots into a multipanel plot
  multipanel_plot <- gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)
  
  print(invisible(multipanel_plot))
  return(outlier_df)
}
