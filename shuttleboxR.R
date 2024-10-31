# Script for R package shuttleboxR
#
# 25/026/2024
#
# Shaun Killen and Emil Christensen
#
# reset R
rm(list=ls())

# renv:: restore checks whether the packages you are using are the same versions as in the script
# if collaborators update a library and then edit the script, renv::restore will pick up on that
renv::restore()

# Confirm that R is looking in the right place
getwd()

proj_data <- read.csv("project_database.csv") # This should be a summary data file with data from numerous fish (e.g. all fish in a study)

additional_data<-read.csv("additional_data.csv")

#### Function to prepare data file for further processing and detect shuttles ####

file_prepare <- function(file, additional_data) {
  
  # Load .txt data file that shuttlesoft produces into R
  # Don't load headers, because the additional info at the top of the .txt file will make a confusing dataframe
  # Rename data columns
  
  data <- read.delim(file,
                     header = F,
                     col.names = (c("Clock_time_hours", "Side_presence", "Body_core_temp", "Pref_temp_loligo", "INCR_side_temp",
                                    "DECR_side_temp", "x_pos", "y_pos", "U_swim", "Dist_moved", "Time_in_INCR", "Time_in_DECR",
                                    "Side_temp_diff", "Hyst_side_temp_diff", "INCR_side_set", "Hyst_INCR_side_set", "DECR_side_set",
                                    "Hyst_DECR_side_set", "k", "Max_temp_loligo", "Min_temp", "Temp_change_rate", "Avoid_up_mean_loligo",
                                    "Unknown_1_loligo", "Unknown_2_loligo", "Unknown_3_loligo")),
                     
                     na.strings = c("NaN", ""))
  
  # Ensure necessary columns exist in the additional_data file
  if (!all(c("file_name", "trial_start", "mass", "k", "constant") %in% colnames(additional_data))) {
    stop("The additional dataset does not contain one or more necessary columns: 'file_name', 'trial_start', 'mass', 'k', 'constant'.")
  }
  
  #Ensure necessary columns are in the right format
  if (all(grepl("^\\d{2}:\\d{2}:\\d{2}$", additional_data$trial_start))==F){
    stop("One or more entries in trial_start are not in an hh:mm:ss format")
  }
  
  # Extract notes, and the file created info from the datafrmae
  notes <- data[3,2]
  filecreated <- as.POSIXct(data[1,2], format = "%d/%m/%Y; %H:%M")
  
  # Remove the additional info and reset rownames
  data<-data[-c(1:6), ]
  rownames(data)<-NULL
  
  # Add time in seconds and hours assuming a sampling frequency of 1 Hz
  data$Time_sec <- seq(0, by = 1, length.out = nrow(data))
  data$Time_h <- data$Time_sec / 3600
  
  # Add the date off the trial by extracting it from the filecreated object
  data$date <- as.Date(filecreated)
  
  # Make sure the date is correct in case the experiment goes on past midnight
  # Extract the second when the experiment passed midnight
  # Make all observations past that date the next day
  midnight_observation <- data$Time_sec[data$Clock_time_hours == "00:00:00"]
  if(length(midnight_observation)>0){
    data$date<-ifelse(data$Time_sec>=midnight_observation, 
                      data$date+1, 
                      data$date)
  }
  
  data$datetime<-as.POSIXct(paste(data$date, data$Clock_time_hours), format = "%Y-%m-%d %H:%M:%S")
  
  # Add a column that tells you whether the system is in static or dynamic
  data$dyn_stat <- ifelse(is.na(data$Side_temp_diff), "static", "dynamic")
  
  filename <- gsub(paste0("^", getwd(), "/?"), "", file)
  
  trial_start<-additional_data$trial_start[additional_data$file_name==filename]
  
  # Add the phase of the experiment (e.g. acclimation and trial)
  trial_start_second<-data$Time_sec[data$Clock_time_hours == trial_start]
  data$trial_phase<-ifelse(data$Time_sec>=trial_start_second, "trial", "acclimation")
  
  # Detect a shuttle (chamber side change) in new column
  data$shuttle <- 0
  for (i in 2:nrow(data)) {
    if (data$Side_presence[i] != data$Side_presence[i-1]) {
      data$shuttle[i] <- 1
    } else {
      data$shuttle[i] <- 0
    }
  }
  
  # Ensure the shuttle column is numeric
  data$shuttle <- as.numeric(data$shuttle)
  
  data<-type.convert(data, as.is = T)
  
  return(data)
}

data <- file_prepare(file = "C:/GitHub/ShuttleboxR/Sys1_20240613_parcatus_rr.txt", additional_data = additional_data)

#### Function to calculate body core temperature if not already done in Shuttlesoft ####

calc_core_body_temp <- function(data, BM, constant, k) {
  # Ensure necessary columns exist
  if (!all(c("Body_core_temp", "shuttle") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'Body_core_temp', 'shuttle'.")
  }
  
  #Ensure necessary columns are in the right format
  if (!all(is.numeric(c(data$Body_core_temp, data$shuttle)))){
    stop("One or more necessary variables are not in the right format")
  }

  # Initialize the Body_core_temp column
  data$Body_core_temp <- data$Body_core_temp
  
  # Calculate body core temperature from second to second
  for (i in 2:nrow(data)) {
    if (data$shuttle[i] == 1) {
      data$Body_core_temp[i] <- data$Body_core_temp[i-1]  # Carry forward the body core temperature at the moment of shuttling
    } else {
      data$Body_core_temp[i] <- data$Body_core_temp[i] + 
        (data$Body_core_temp[i-1] - data$Body_core_temp[i]) * exp(-(constant * BM^k) * 1 / 60)
    }
  }
  
  return(data)
}

# Example usage
BM <- 30  # Example body mass
constant <- 3.69  # Example constant
k <- -0.574  # Example k value

# Calculate core body temperature
data <- calc_core_body_temp(data, BM, constant, k)


#### Function to calculate Tpref using either the mean or median of all body core temperature values ####

calc_Tpref <- function(data, method = c("mean", "median", "mode"), exclude_start_minutes = 0, exclude_end_minutes = 0) {
  method <- match.arg(method)
  
  # Convert the time to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec) - (exclude_end_minutes * 60)
  data <- data[data$Time_sec > start_time & data$Time_sec < end_time, ]
  
  # Convert the Body_core_temp column to numeric
  data$Body_core_temp <- as.numeric(data$Body_core_temp)
  
  # Calculate Tpref based on the specified method
  if (method == "mean") {
    Tpref <- mean(data$Body_core_temp, na.rm = TRUE)
  } else if (method == "median") {
    Tpref <- median(data$Body_core_temp, na.rm = TRUE)
  } else if (method == "mode") {
    Tpref <- as.numeric(names(sort(table(data$Body_core_temp), decreasing = TRUE))[1])
  }
  
  print(paste("Tpref:", Tpref))
  return(Tpref)
}

# Example usage
Tpref <- calc_Tpref(data, method = "median", exclude_start_minutes = 0, exclude_end_minutes = 0)

#### Function to calculate upper and lower avoidance temperatures ####

calc_Tavoid <- function(data, percentiles = c(0.05, 0.95), exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure the percentiles are valid
  if (length(percentiles) != 2 || any(percentiles < 0) || any(percentiles > 1)) {
    stop("Percentiles should be a vector of two values between 0 and 1")
  }
  
  # Convert exclude minutes to seconds
  exclude_start_seconds <- exclude_start_minutes * 60
  exclude_end_seconds <- exclude_end_minutes * 60
  
  # Exclude initial and final data if necessary
  start_time <- exclude_start_seconds
  end_time <- max(data$Time_sec) - exclude_end_seconds
  data <- data[data$Time_sec > start_time & data$Time_sec < end_time, ]
  
  # Calculate the percentiles
  Tavoid_lower <- quantile(data$Body_core_temp, percentiles[1], na.rm = TRUE)
  Tavoid_upper <- quantile(data$Body_core_temp, percentiles[2], na.rm = TRUE)
  
  # Print values
  print(paste("Tavoid Lower:", Tavoid_lower))
  print(paste("Tavoid Upper:", Tavoid_upper)) 
  
  return(c(Tavoid_lower, Tavoid_upper))

}

# Example usage
# Calculate Tavoid values
Tavoid_values <-calc_Tavoid(data, percentiles = c(0.05, 0.95), exclude_start_minutes = 0, exclude_end_minutes = 0)

#### Function to calculate total distance moved ####

calc_distance <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure the Distance moved column exists
  if (!"Dist_moved" %in% colnames(data)) {
    stop("The dataset does not contain a 'Dist_moved' column. Please run file_prepare function.")
  }
  
  # Convert exclude minutes to seconds
  exclude_start_seconds <- exclude_start_minutes * 60
  exclude_end_seconds <- exclude_end_minutes * 60
  
  # Exclude initial and final data if necessary
  start_time <- exclude_start_seconds
  end_time <- max(data$Time_sec) - exclude_end_seconds
  initial_distance <- 0
  
  if (exclude_start_seconds > 0) {
    initial_distance <- data$Dist_moved[which(data$Time_sec <= exclude_start_seconds)]
    if (length(initial_distance) > 0) {
      initial_distance <- max(initial_distance, na.rm = TRUE)
    } else {
      initial_distance <- 0
    }
  }
  
  data <- data[data$Time_sec > start_time & data$Time_sec < end_time, ]
  
  # Calculate net distance moved
  total_distance <- max(data$Dist_moved, na.rm = TRUE)
  net_distance_moved <- total_distance - initial_distance
  
  # Print the result
  print(paste("Distance Moved:", net_distance_moved, "cm"))
  
  return(net_distance_moved)
}

# Example usage
dist <- calc_distance(data, 240)

#### Function to calculate shuttling frequency ####

calc_shuttling_frequency <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure the 'shuttle' and 'Time_sec' columns exist
  if (!all(c("shuttle", "Time_sec") %in% colnames(data))) {
    stop("The dataset does not contain 'shuttle' and/or 'Time_sec' columns.")
  }
  
  # Convert Time_sec to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec) - (exclude_end_minutes * 60)
  data <- data[data$Time_sec > start_time & data$Time_sec < end_time, ]
  
  # Calculate shuttling frequency
  shuttles <- sum(data$shuttle, na.rm = TRUE)
  
  # Print the result
  print(paste("Shuttling frequency:", shuttles))
  
  return(shuttles)
}

# Example usage
cal_shut_freq <- calc_shuttling_frequency(data, exclude_start_minutes = 240, exclude_end_minutes = 5)

#### Function to calculate occupancy time in each chamber ####

calc_occupancy_time <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure the 'Side_presence' and 'Time_sec' columns exist
  if (!all(c("Side_presence", "Time_sec") %in% colnames(data))) {
    stop("The dataset does not contain 'Side_presence' and/or 'Time_sec' columns.")
  }
  
  # Convert Time_sec to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec) - (exclude_end_minutes * 60)
  data <- data[data$Time_sec > start_time & data$Time_sec < end_time, ]
  
  # Calculate occupancy time in each chamber
  time_in_chamberDECR <- sum(data$Side_presence == "DECR", na.rm = TRUE)
  time_in_chamberINCR <- sum(data$Side_presence == "INCR", na.rm = TRUE)
  
  # Print the results
  print(paste("Time in DECR Chamber:", time_in_chamberDECR))
  print(paste("Time in INCR Chamber:", time_in_chamberINCR))
  
  return(list(time_in_chamberDECR = time_in_chamberDECR, time_in_chamberINCR = time_in_chamberINCR))
}

# Example usage
occ_time <- calc_occupancy_time(data, exclude_start_minutes = 240, exclude_end_minutes = 5)

#### Function to plot a histogram of time spent at different core body temperatures ####

library(ggplot2)

plot_body_core_temp_histogram <- function(data, bin_size, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure the body core temperature column exists
  if (!"Body_core_temp" %in% colnames(data)) {
    stop("The dataset does not contain a 'Body_core_temp' column. Please run file_prepare function.")
  }
  
  # Convert Time_sec to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec) - (exclude_end_minutes * 60)
  data <- data[data$Time_sec > start_time & data$Time_sec < end_time, ]
  
  # Calculate the duration of each bin in minutes
  data$Time_min <- data$Time_sec / 60
  total_time <- max(data$Time_min, na.rm = TRUE)
  
  # Create the histogram data
  hist_data <- hist(data$Body_core_temp, breaks = seq(floor(min(data$Body_core_temp, na.rm = TRUE)),
                                                      ceiling(max(data$Body_core_temp, na.rm = TRUE)),
                                                      by = bin_size), plot = FALSE)
  
  # Find the peak value
  Tpref_peak <- hist_data$mids[which.max(hist_data$counts)]
  
  # Create the histogram plot
  hist_plot <- ggplot(data, aes(x = Body_core_temp)) +
    geom_histogram(binwidth = bin_size, aes(y = (..count.. / sum(..count..)) * 100), fill = "#F79518", color = "black") +
    geom_vline(xintercept = Tpref_peak, color = "#9EF1A4", linetype = "dashed", size = 2) +
    labs(title = "Histogram of Percent of Time Spent at Different Body Core Temperatures",
         subtitle = paste("Peak Tpref:", round(Tpref_peak, 2), "°C"),
         x = "Body Core Temperature (°C)",
         y = "Percent of Total Time (%)") +
    theme_minimal()
  
  print(hist_plot)
  print(paste("Peak Tpref:", Tpref_peak))
  
  return(Tpref_peak)
}

# Example usage
Tpref_peak <- plot_body_core_temp_histogram(data, bin_size = 0.1, exclude_start_minutes = 240, exclude_end_minutes = 5)
print(paste("Peak Tpref:", Tpref_peak))

#### Function to calculate variance in core body temperature experienced throughout the trial #####

calc_variance <- function(data, variance_type = c("std_error", "std_deviation", "coeff_variation"), exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure the body core temperature column exists
  if (!"Body_core_temp" %in% colnames(data)) {
    stop("The dataset does not contain a 'Body_core_temp' column. Please run file_prepare function.")
  }
  
  # Convert exclude minutes to seconds
  exclude_start_seconds <- exclude_start_minutes * 60
  exclude_end_seconds <- exclude_end_minutes * 60
  
  # Exclude initial and final data if necessary
  start_time <- exclude_start_seconds
  end_time <- max(data$Time_sec) - exclude_end_seconds
  data <- data[data$Time_sec > start_time & data$Time_sec < end_time, ]
  
  # Select the variance type
  variance_type <- match.arg(variance_type)
  
  # Calculate the variance measure
  if (variance_type == "std_error") {
    variance <- sd(data$Body_core_temp, na.rm = TRUE) / sqrt(length(na.omit(data$Body_core_temp)))
  } else if (variance_type == "std_deviation") {
    variance <- sd(data$Body_core_temp, na.rm = TRUE)
  } else if (variance_type == "coeff_variation") {
    variance <- sd(data$Body_core_temp, na.rm = TRUE) / mean(data$Body_core_temp, na.rm = TRUE)
  }
  
  return(variance)
}

# Example usage
variance_se <- calc_variance(data, variance_type = "std_error", 240, 5)
variance_sd <- calc_variance(data, variance_type = "std_deviation", 240, 5)
variance_cv <- calc_variance(data, variance_type = "coeff_variation", 240, 5)

print(paste("Standard Error:", variance_se))
print(paste("Standard Deviation:", variance_sd))
print(paste("Coefficient of Variation:", variance_cv))

#### Function to plot cumulative changes in distance moved over time ####

plot_cumulative_distance <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  if (!("Time_sec" %in% colnames(data) && "Dist_moved" %in% colnames(data))) {
    stop("The dataset does not contain 'Time_sec' and/or 'Dist_moved' columns.")
  }
  
  # Convert Time_sec to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec) - (exclude_end_minutes * 60)
  data <- data[data$Time_sec > start_time & data$Time_sec < end_time, ]
  
  # Convert time in seconds to minutes
  data$Time_min <- data$Time_sec / 60
  #substract first distance value to reset cumulative distance to 0
  data$Dist_moved<-data$Dist_moved-data$Dist_moved[[1]]
  
  # Plot cumulative distance
  plot <- ggplot(data, aes(x = Time_min, y = (Dist_moved))) +
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

calc_min_gravitation <- function(start_temp, target_temp, rate_of_change) {
  # Calculate the change in temperature
  delta_temp <- target_temp - start_temp
  
  # Calculate the minimum gravitation time
  gravitation_time <- delta_temp / (rate_of_change/60)
  
  # Print the result
  print(paste("Minimum Gravitation Time:", gravitation_time, "minutes"))
  
  return(gravitation_time)
}

# Example usage
start_temp <- 20  # Starting mean temperature in °C
target_temp <- 25  # Target temperature in °C, e.g. approximation of anticipated Tpref
rate_of_change <- 5  # Rate of temperature change in °C/h

min_grav_time <- calc_min_gravitation(start_temp, target_temp, rate_of_change)

#### Function to calculate time spent at extreme temperatures, near set limits ####

calc_extremes <- function(data, threshold = 0.2*(max(data$Max_temp_loligo)-max(data$Min_temp)), exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!("Time_sec" %in% colnames(data) && "Body_core_temp" %in% colnames(data))) {
    stop("The dataset does not contain 'Time_sec' and/or 'Body_core_temp' columns.")
  }
  
  # Convert Time_sec to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec) - (exclude_end_minutes * 60)
  data <- data[data$Time_sec > start_time & data$Time_sec < end_time, ]
  
  # Convert time in secxonds to minutes
  data$Time_min <- data$Time_sec / 60
  
  upper_limit<-max(data$Max_temp_loligo)
  lower_limit<-max(data$Min_temp)
  
  # Determine time spent near upper and lower extreme temperatures with the threshold
  upper_threshold <- upper_limit - threshold
  lower_threshold <- lower_limit + threshold
  data$Upper_Extreme_Temp <- ifelse(data$Body_core_temp > upper_threshold, 1, 0)
  data$Lower_Extreme_Temp <- ifelse(data$Body_core_temp < lower_threshold, 1, 0)
  
  time_near_upper_extreme <- sum(data$Upper_Extreme_Temp) / nrow(data) * 100
  time_near_lower_extreme <- sum(data$Lower_Extreme_Temp) / nrow(data) * 100
  
  # Print the results
  print(paste("Time spent near upper extreme temperatures (%):", time_near_upper_extreme))
  print(paste("Time spent near lower extreme temperatures (%):", time_near_lower_extreme))
  
  return(c(time_near_lower_extreme, time_near_upper_extreme))
}

# Example usage
extreme_time_percent <- calc_extremes(data, exclude_start_minutes = 240, exclude_end_minutes = 5)
calc_extremes(data, exclude_start_minutes = 0, exclude_end_minutes = 0)[2]


#### Function to plot histogram of movement speeds ####

library(ggplot2)

plot_speed_histogram <- function(data, binwidth = 0.1, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary column exists
  if (!"U_swim" %in% colnames(data)) {
    stop("The dataset does not contain a 'U_swim' column.")
  }
  
  # Convert Time_sec to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec, na.rm = TRUE) - (exclude_end_minutes * 60)
  data <- data[data$Time_sec > start_time & data$Time_sec < end_time, ]
  
  # Create the histogram plot
  plot <- ggplot(data, aes(x = U_swim)) +
    geom_histogram(binwidth = binwidth, fill = "#F79518", color = "black") +
    labs(title = "Histogram of Movement Speeds",
         x = "Movement Speed (cm/s)",
         y = "Frequency") +
    # scale_y_log10()+
    theme_classic()
  
  print(plot)
}

# Example usage
plot_speed_histogram(data, binwidth = 1, exclude_start_minutes = 240, exclude_end_minutes = 5)

#### Function to plot movement speed against body core temperature ####

plot_speed_vs_core_temp <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!all(c("Body_core_temp", "Dist_moved", "Time_sec") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'Body_core_temp', 'Dist_moved', 'Time_sec'. Please run file_prepare function.")
  }
  
  # Convert Time_sec to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec) - (exclude_end_minutes * 60)
  data <- data[data$Time_sec > start_time & data$Time_sec < end_time, ]
  
  # Calculate movement speed as the difference in distance moved over time
  data$Movement_Speed <- c(NA, diff(data$Dist_moved) / diff(data$Time_sec))
  
  # Remove rows with NA values
  data <- na.omit(data[, c("Body_core_temp", "Movement_Speed")])
  
  # Create the scatter plot
  plot <- ggplot(data, aes(x = Body_core_temp, y = Movement_Speed)) +
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
speed_temp_plot <- plot_speed_vs_core_temp(data, exclude_start_minutes = 240)
speed_temp_plot

#### Function to plot temperatures in each side of the shuttlebox over time ####

library(ggplot2)

plot_temperature_gradient <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Convert Time_sec to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec) - (exclude_end_minutes * 60)
  data <- data[data$Time_sec > start_time & data$Time_sec < end_time, ]
  
  # Convert time in seconds to hours
  data$Time_h <- data$Time_sec / 3600
  
  # Plot the temperature gradient
  plot <- ggplot(data, aes(x = Time_h)) +
    geom_line(aes(y = `DECR_side_temp`, color = "DECR side temp")) +
    geom_line(aes(y = `INCR_side_temp`, color = "INCR side temp")) +
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
plot_temperature_gradient(data, 0, 0)

#### Function to plot core body temperature, Tpref, avoidance temperatures, segmented regression of changes in body core temperature during the trial #####

library(ggplot2)
library(segmented)

plot_temp_segmented <- function(data, Tpref, Tavoid_lower, Tavoid_upper, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!all(c("Time_h", "Body_core_temp") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'Time_h', 'Body_core_temp'.")
  }
  
  # Convert the time to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec) - (exclude_end_minutes * 60)
  data <- data[data$Time_sec > start_time & data$Time_sec < end_time, ]
  
  # Perform segmented regression
  fit <- lm(Body_core_temp ~ Time_h, data = data)
  seg_fit <- segmented(fit, seg.Z = ~Time_h, npsi = 1)
  
  # Get summary and breakpoint
  seg_summary <- summary(seg_fit)
  breakpoints <- seg_fit$psi[, "Est."]
  
    # Plot data with segmented regression
  plot <- ggplot(data, aes(x = Time_h, y = Body_core_temp)) +
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
plot_temp_segmented(data, Tpref = Tpref, Tavoid_lower, Tavoid_upper, exclude_start_minutes = 240, exclude_end_minutes = 5)

#### Function to plot changes in core body temperature over time, as well as Tpref, and upper and lower avoidance temperatures, but WITHOUT segmented regression lines #####

plot_temp <- function(data, Tpref, Tavoid_low, Tavoid_up, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!all(c("Time_h", "Body_core_temp") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'Time_h', 'Body_core_temp'.")
  }
  
  # Convert the time to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec) - (exclude_end_minutes * 60)
  data <- data[data$Time_sec > start_time & data$Time_sec < end_time, ]
  
  # Plot data
  plot <- ggplot(data, aes(x = Time_h, y = Body_core_temp)) +
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
plot_temp(data, Tpref, Tavoid_low, Tavoid_up, exclude_start_minutes = 240, exclude_end_minutes = 5)

##### Function to calculate actual gravitation time ####

library(segmented)

calc_act_gravitation <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  
  # Convert the time to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec) - (exclude_end_minutes * 60)
  data <- data[data$Time_sec > start_time & data$Time_sec < end_time, ]
  
  # Perform segmented regression
  fit <- lm(Body_core_temp ~ Time_h, data = data)
  seg_fit <- segmented(fit, seg.Z = ~Time_h, npsi = 1)
  
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
  
  print(paste("Gravitation Time (hours):", gravitation_time))
  
  return(gravitation_time)
}

# Calculate gravitation time
act_grav_time <- calc_act_gravitation(data)

##### Function to calculate and plot interval means for shuttling rate over time ####

library(ggplot2)
library(dplyr)

shuttling_aggregated <- function(data, interval_minutes = 10, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!all(c("Time_sec", "shuttle") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'Time_sec', 'shuttle'.")
  }
  
  # Convert Time_sec to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec, na.rm = TRUE) - (exclude_end_minutes * 60)
  data <- data %>% filter(Time_sec > start_time & Time_sec < end_time)
  
  # Calculate time intervals
  data <- data %>% 
    mutate(Time_interval = floor((Time_sec - start_time) / (interval_minutes * 60)) * interval_minutes)
  
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
shuttling_rate <- shuttling_aggregated(data, interval_minutes = 60, exclude_start_minutes = 240, exclude_end_minutes = 5)

##### Function to calculate and plot interval means for speed over time ####

library(ggplot2)
library(dplyr)

speed_aggregated <- function(data, interval_minutes = 10, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!all(c("Time_sec", "U_swim") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'Time_sec', 'U_swim'.")
  }
  
  # Convert Time_sec to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec, na.rm = TRUE) - (exclude_end_minutes * 60)
  data <- data %>% filter(Time_sec > start_time & Time_sec < end_time)
  
  # Calculate time intervals
  data <- data %>% 
    mutate(Time_interval = floor((Time_sec - start_time) / (interval_minutes * 60)) * interval_minutes)
  
  # Calculate mean movement speed for each interval
  movement_speed <- data %>% 
    group_by(Time_interval) %>% 
    summarise(mean_movement_speed = mean(U_swim, na.rm = TRUE))
  
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
movement_speed <- speed_aggregated(data, interval_minutes = 60, exclude_start_minutes = 240, exclude_end_minutes = 5)

##### Function to calculate and plot interval means for body core temperature over time ####

library(ggplot2)
library(dplyr)

Tcore_aggregated <- function(data, interval_minutes = 10, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!all(c("Time_sec", "Body_core_temp") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'Time_sec', 'Body_core_temp'.")
  }
  
  # Convert Time_sec to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec, na.rm = TRUE) - (exclude_end_minutes * 60)
  data <- data %>% filter(Time_sec > start_time & Time_sec < end_time)
  
  # Calculate time intervals
  data <- data %>% 
    mutate(Time_interval = floor((Time_sec - start_time) / (interval_minutes * 60)) * interval_minutes)
  
  # Calculate mean preferred temperature for each interval
  Tcore_aggregated <- data %>% 
    group_by(Time_interval) %>% 
    summarise(mean_Tcore = mean(Body_core_temp, na.rm = TRUE))
  
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
  if (!all(c("Time_sec", "U_swim", "Body_core_temp", "shuttle") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'Time_sec', 'U_swim', 'Body_core_temp', 'shuttle'.")
  }
  
  # Convert Time_sec to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec, na.rm = TRUE) - (exclude_end_minutes * 60)
  data <- data %>% filter(Time_sec > start_time & Time_sec < end_time)
  
  # Calculate time intervals
  data <- data %>% 
    mutate(Time_interval = floor((Time_sec - start_time) / (interval_minutes * 60)) * interval_minutes)
  
  # Calculate mean values for U_swim, Tpref (Body_core_temp), and shuttling rate for each interval
  summary_data <- data %>% 
    group_by(Time_interval) %>% 
    summarise(mean_U_swim = mean(U_swim, na.rm = TRUE),
              mean_Tpref = mean(Body_core_temp, na.rm = TRUE),
              mean_shuttle_rate = mean(shuttle, na.rm = TRUE))
  
  # Calculate correlation matrix
  correlation_matrix <- cor(summary_data[, -1], use = "complete.obs")
  
  # Print the correlation matrix
  print(correlation_matrix)
  
  # Create scatterplots for pairwise comparisons
  plot1 <- ggplot(summary_data, aes(x = mean_U_swim, y = mean_Tpref)) +
    geom_point(color = "#F79518") +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    labs(title = "Mean U_swim vs Mean Tpref",
         x = "Mean U_swim (cm/s)",
         y = "Mean Tpref (°C)") +
    theme_classic()
  
  plot2 <- ggplot(summary_data, aes(x = mean_U_swim, y = mean_shuttle_rate)) +
    geom_point(color = "#F79518") +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    labs(title = "Mean U_swim vs Mean Shuttling Rate",
         x = "Mean U_swim (cm/s)",
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

results
##### Function to plot a 3D animation of fish movements at selected time intervals, with time as a dimension ####

library(rgl)

animate_movements <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!all(c("x_pos", "y_pos", "Time_sec", "Body_core_temp") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'x_pos', 'y_pos', 'Time_sec', 'Body_core_temp'.")
  }
  
  # Convert Time_sec to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec, na.rm = TRUE) - (exclude_end_minutes * 60)
  data <- data[data$Time_sec > start_time & data$Time_sec < end_time, ]
  
  # Set up 3D plot
  plot3d(data$x_pos, data$y_pos, data$Time_sec, type = "l", col = "black", alpha = 0.15,
         xlab = "X Position", ylab = "Y Position", zlab = "Time (seconds)")
  
  # Determine color scaling based on Tpref
  colors <- colorRampPalette(c("royalblue3", "red2"))(length(unique(data$Body_core_temp)))
  data$color <- colors[as.numeric(cut(data$Body_core_temp, breaks = length(colors)))]
  
  # Animate the plot
  for (i in seq_len(nrow(data))) {
    points3d(data$x_pos[i], data$y_pos[i], data$Time_sec[i], col = data$color[i], size = 5)
    Sys.sleep(0.01)
  }
}

# Example usage
# animate_movements(data, exclude_start_minutes = 900, exclude_end_minutes = 5)

##### Function to produce a heat map of fish locations over the trial ####

library(ggplot2)

plot_heatmap <- function(data, exclude_start_minutes = 0, exclude_end_minutes = 0) {
  # Ensure necessary columns exist
  if (!all(c("x_pos", "y_pos", "Time_sec") %in% colnames(data))) {
    stop("The dataset does not contain one or more necessary columns: 'x_pos', 'y_pos', 'Time_sec'.")
  }
  
  # Convert Time_sec to numeric if not already
  data$Time_sec <- as.numeric(data$Time_sec)
  
  # Exclude rows based on the exclude_start_minutes and exclude_end_minutes parameters
  start_time <- exclude_start_minutes * 60
  end_time <- max(data$Time_sec, na.rm = TRUE) - (exclude_end_minutes * 60)
  data <- data[data$Time_sec > start_time & data$Time_sec < end_time, ]
  
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
heatmap_plot <- plot_heatmap(data, exclude_start_minutes = 240, exclude_end_minutes = 0)
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
plot_limits_vs_distance(proj_data)


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
plot_limits_vs_shuttles(proj_data)

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

# Example usage
results <- pca_and_plot(proj_data)
pca_result <- results$pca_result
scores <- results$scores
loadings <- results$loadings

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

# CONTENTS OF EXTERNAL ADDITIONAL DATA FILE
  # trial start
  # body mass
  # k value
  # constant (for body core temp)

# AIM
  # load all .txt files simultaneously
  # logical operater to filter for trial phase only
  # calculate all variables for a text file, and compile into vector
  # combine values for each .txt file into

# NECESSARY ITEMS
  # directory with .txt files
  # functions to calculate every variable
  # overall (external) data file containing data that is not in the .txt files (mass, trial start, etc.)

# VARIABLES AND ACCOMPANYING FUNCTIONS
  # # ID - from file name 
  # total_length - external data
  # # gravitation time - calc_act_gravitation
  # mass - external data
  # # distance travelled - calc_distance
  # # number of shuttles - file_prepare (so not necessary)
  # # Tpref - calc_Tpref
  # # Tavoid - calc_Tavoid
  # pref range - ?
  # # time near limits - 
  # # body core temp - calc_core_body_temp()
 
  # 

# PROBLEMS


compile_project_data <- function(directory = getwd(), 
                                 additional_data, 
                                 print_results = F,
                                 Tpref_method = "median", 
                                 Tavoid_percintiles = c(0.05, 0.95),
                                 textremes_threshold = 0.2*(max(data$Max_temp_loligo)-max(data$Min_temp)),
                                 core_temp_variance_type = "std_error"){
  
  txt_files <- list.files(path = directory, pattern = "\\.txt$", full.names = TRUE)
  data_list <- lapply(txt_files, file_prepare, additional_data = additional_data)
  
  apply_functions <- function(df, 
                              Tpref_method = "median", 
                              Tavoid_percintiles = c(0.05, 0.95),
                              textremes_threshold = 0.2*(max(df$Max_temp_loligo)-max(df$Min_temp)),
                              core_temp_variance_type = "std_error") {
    Tpref <- calc_Tpref(df, method = "median")
    grav_time <- calc_act_gravitation(df)
    distance <- calc_distance(df)
    Tavoid <- calc_Tavoid(df, percentiles = Tavoid_percintiles)
    textremes <- calc_extremes(df, threshold = textremes_threshold)
    core_temp_variance <- calc_variance(df, variance_type = core_temp_variance_type)
    nr_shuttles <- calc_shuttling_frequency(df)
    
    return(c(Tpref, grav_time, distance, Tavoid, textremes, core_temp_variance, nr_shuttles))
  }
  
  names(data_list)<-txt_files
  
  if (print_results == F) sink(tempfile())
  
  # Apply functions to each dataframe in the list and collect results
  results <- do.call(rbind, lapply(names(data_list), function(name) {
    # Apply the functions
    func_results <- apply_functions(data_list[[name]])
    # invisible(capture.output(apply_functions(), file = "NUL")) 
    
    # Combine results into a data frame
    data.frame(
      filename = gsub(paste0("^", directory, "/?"), "", name),
      Tpref = func_results[1],
      grav_time  = func_results[2],
      distance  = func_results[3],
      Tavoid_low  = func_results[4],
      Tavoid_high = func_results[5],
      t_near_low_extreme  = func_results[6],
      t_near_high_extreme  = func_results[7],
      core_temp_variance  = func_results[8],
      nr_shuttles  = func_results[9]
      
    )
    
  }))
  
  if (print_results == F) sink()
  
  rownames(results) <- NULL
  
  return(results)
  
}

results2<-compile_project_data(print_results = T, additional_data = additional_data)

class(proj_data)

