#### load ShuttleboxR #####
install.packages("ShuttleboxR")
library(ShuttleboxR)
library(tidyverse)

# Function to split project data into chunks small enough for Google Sheets
split_pd <- function(project_data, directory = getwd()){
  
  # Define the maximum number of cells allowed in a Google Sheets
  num_rows <- nrow(project_data)
  num_cols <- ncol(project_data)
  max_cells <- 10e6-(num_cols+1)
  
  # Calculate the max number of rows per sheet
  max_rows_per_sheet <- max_cells / num_cols
  
  # Round down the value to an integer, because we can't have partial rows
  max_rows_per_sheet <- floor(max_rows_per_sheet)
  
  # Split the dataset into chunks of data that do not exceed the max row limit
  split_data <- split(project_data, ceiling(seq_along(1:num_rows) / max_rows_per_sheet))
  
  object_name <- deparse(substitute(project_data))
  
  # Save each chunk into separate CSV files
  # Optionally, save the files as RDS for later use
  for (i in seq_along(split_data)) {
    # Save to CSV
    write_xlsx(split_data[[i]], paste0(directory, "/", object_name, "_chunk_", i, ".xlsx"))
  }
  
  # Print the number of chunks
  cat("Data has been split into", length(split_data), "chunks.\n")
}

# Function to read project data  Google Sheets
fetch_pd <- function(sheet_url) {
  options(timeout = 300)
  
  # Ensure the URL is in CSV format
  sheet_url <- sub("/edit.*", "/pub?output=csv", sheet_url)
  
  # Read the Google Sheet containing links and chunk numbers
  link_data <- read.csv(sheet_url, stringsAsFactors = FALSE)
  
  # Extract links (assuming they are in the first column)
  csv_links <- link_data[[2]]
  
  # Extract chunk numbers (assuming they are in the second column)
  chunk_numbers <- link_data[[1]]
  
  # Initialize an empty list to store data frames
  df_list <- list()
  
  # Loop through each link, read the CSV, and store it in the list
  for (i in seq_along(csv_links)) {
    cat("Reading chunk:", chunk_numbers[i], "\n")  # Print chunk number instead of link
    
    # Try to read the CSV, skip if there's an error
    df <- tryCatch(read.csv(csv_links[i], stringsAsFactors = FALSE), 
                   error = function(e) {
                     message("Failed to read chunk:", chunk_numbers[i])
                     return(NULL)
                   })
    
    # Add to list if the read was successful
    if (!is.null(df)) {
      df_list[[i]] <- df
    }
  }
  
  # Combine all data frames into one
  final_df <- do.call(rbind, df_list)
  
  return(final_df)
}

#### Glasgow examples ####

library(tidyverse)

getwd()

# import Glasgow metadata
metadata_glasgow<-read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vR6iHLwGLQh-2MOd9NlPjMAaWubOSsa56RQZQ5IGlqJBIr1DmKlbeHpxgPGjanCVD1hGLNRkW9PLAKL/pub?gid=44598193&single=true&output=csv")

# import Glasfow single trial
gg_file <- "Fish_1_13_3.txt"
gg <- read_shuttlesoft(gg_file, metadata_glasgow)
gg <- file_prepare(gg)
gg <- calc_coreT(gg)
gg <- gg%>%
  mutate(trial_phase = ifelse(dyn_stat == "static", "acclimation", "trial"))

# Calculate shuttle-box metrics for single trial
calc_Tpref(gg, exclude_acclimation = T, print_results = F)
calc_Tavoid(gg, exclude_acclimation = T, print_results = F)
calc_extremes(gg, exclude_acclimation = T, print_results = F)
calc_gravitation(gg, exclude_acclimation = T, print_results = F)
calc_tot_distance(gg, exclude_acclimation = T, print_results = F)
calc_shuttles(gg, exclude_acclimation = T, print_results = F)
calc_occupancy(gg, exclude_acclimation = T, print_results = F)
calc_coreT_variance(gg, exclude_acclimation = T, print_results = F, variance_type = "std_error")

# Inspect single trial
plot_T_gradient(gg)
plot_T_segmented(gg, exclude_acclimation = T)
plot_interval(gg, "shuttle", 30)
plot_coreT_histogram(gg, exclude_acclimation = T)
plot_distance(gg)
plot_histogram(gg, column = "velocity", binwidth = 1)
plot_interval(gg, column = "velocity")
plot_heatmap(gg, exclude_acclimation = T)
plot_speed_coreT(gg)

# import Glasgow project data (output of compile_project_data)
pd_glasgow <- fetch_pd(sheet_url <- "https://docs.google.com/spreadsheets/d/e/2PACX-1vSp_wq2ExdlA_fGDquavCeFIZLqnR_uLgn8G53gvUjRqOfmif9JVg_jtx6upE1OyyAk__jv_wjOF2cw/pub?gid=0&single=true&output=csv")

library(ggokabeito)

# Inspect coreT of project dataset
pd_glasgow %>%
  filter(fileID %in% unique(fileID)[1:9])%>%
  ggplot(aes(x=time_h, y = core_T, colour = id))+
  geom_line(linewidth = 1)+
  theme_light()+
  scale_color_okabe_ito()+
  labs(title = "Body Core Temperature vs Time",
       x = "Time (h)",
       y = "Body Core Temperature (°C)")

# import Glasgow results (output of calc_project_results)
  # glasgow_read <- read_shuttlesoft_project(metadata_glasgow,"C:/GitHub/ShuttleboxR/Glasgow_2011")
  # glasgow_results <- calc_project_results(glasgow_read, exclude_acclimation = T)
  # write_csv(glasgow_results, "G:/My Drive/MSc2 Frans-Polynesië/ShuttleboxR paper/Data/glasgow_results.csv")

glasgow_results<-read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTrk1nBRqdY6JZMVBYN_YtP5Dg9wmTjpH9qcIO8m4kLPM-7QdTPFzx_4YOtON_sek_H8wg8APBEsJsm/pub?gid=1375442230&single=true&output=csv")
glasgow_results<-left_join(glasgow_results,
                           metadata_glasgow%>%
                             dplyr::select(file_name, id),
                           by = join_by(fileID == file_name))%>%
  mutate(t_near_limits = t_near_min + t_near_max)

# Perform PCA on results dataset

gg_pca_dat<-glasgow_results
gg_pca<-pca(gg_pca_dat, id_col = "id")
gg_pca$outliers
gg_pca$plots$biplot
gg_pca$plots$varplot
gg_pca$plots$screeplot
gg_pca$pca_scores


# Plot frequency distributions of project dataset
plot_histograms(glasgow_results, bin_size_distance = 50000, bin_size_shuttles = 500, nr_sd = 1.5, iqr_multiplier = 1, id_col = "id")

#### Moorea examples ####

# Metadata
library(tidyverse)

# import Moorea metadata
metadata_moorea <- read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRfr5gvIdH3Z7ZSpeHmbKMfyVhVZWwPXGl6-EW5MCT8yUGDpOYtskOzCrcQPn72AibVONUhub54FYQV/pub?gid=1247062796&single=true&output=csv")
getwd()



#### Moorea SINGLE TRIALS ####



#Import Moorea single trial example
mo_file <- "Sys1_20240620_dflavicaudus_gg.d.03.txt"

mo <- read_shuttlesoft(file = mo_file, metadata = metadata_moorea)
mo <- file_prepare(mo)
mo <- calc_coreT(mo)

# Calculate shuttle-box metrics
calc_Tpref(mo, exclude_acclimation = T, print_results = F)
calc_Tpref(mo, print_results = F, exclude_start_minutes = exclude_start_minutes, exclude_end_minutes = exclude_end_minutes,  exclude_acclimation = exclude_acclimation)
calc_Tavoid(mo, exclude_acclimation = T, print_results = F)
calc_extremes(mo, exclude_acclimation = T, print_results = F)
calc_gravitation(mo, exclude_acclimation = T, print_results = F)
calc_tot_distance(mo, exclude_acclimation = T, print_results = F)
calc_shuttles(mo, exclude_acclimation = F, print_results = F)
calc_occupancy(mo, exclude_acclimation = T, print_results = T)
calc_coreT_variance(mo, exclude_acclimation = T, variance_type = "std_error")
calc_track_accuracy(mo)


# Inspect moorea single trial
plot_T_gradient(mo)
plot_T_segmented(mo, exclude_acclimation = T, Tavoid_percentiles = c(0.1, 0.9), overlay_chamber_temp = F)
plot_interval(mo, "shuttle", 30)
plot_coreT_histogram(mo_1, exclude_acclimation = T)
plot_distance(mo_2)
plot_histogram(mo_2, column = "velocity", binwidth = 1)
plot_interval(mo, column = "velocity")
plot_heatmap(mo, exclude_acclimation = F)
plot_speed_coreT(mo)

# Import project dataset (output of compile_project_data)
pd_moorea <- fetch_pd("https://docs.google.com/spreadsheets/d/e/2PACX-1vRaz3RkVSO_wFGaZE4q_aTulTMrbukd4mv1eqeQ2Ao1hGaup16i1MA5FHvaW5vToHBZFxDBSrHome_p/pub?gid=0&single=true&output=csv")
pd_moorea <- pd_moorea%>%
  mutate(species = case_when(grepl(".d.", id, ignore.case = T) ~ "D. flavicaudus",
                             grepl(".h.", id, ignore.case = T) ~ "P. arcatus",
                             grepl(".c.", id, ignore.case = T) ~ "C. maculatus"))

pd_moorea <- pd_moorea%>%
  group_by(fileID)%>%
  mutate(int_hour = ceiling(time_h))%>%
  ungroup()%>%
  group_by(fileID, int_hour)%>%
  mutate(shuttling_frequency = sum(shuttle))

cb_palette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
  "#D55E00", "#CC79A7", "#AD7700", "#1C91D4", "#28A745", 
  "#D9CE00", "#1464A0", "#A14600", "#7E0E6D", "#F4A582", 
  "#92C5DE", "#D6604D", "#4393C3", "#88CCEE", "#44AA99", 
  "#DDCC77", "#332288", "#117733", "#737373", "#661100")



# Inspect coreT of project dataset
pd_moorea%>%
  filter(trial_phase != "acclimation")%>%
  ggplot(aes(x=time_h, y = core_T, colour = id))+
  geom_line(size = 1)+
  theme_light()+
  labs(title = "Body Core Temperature vs Time",
       x = "Time (h)",
       y = "Body Core Temperature (°C)")+
  scale_colour_manual(values = cb_palette) +
  theme(legend.position = "top")+
  facet_grid(cols = vars(species))

# Import project results (output of calc_project results)
moorea_read <- read_shuttlesoft_project(metadata_moorea, "C:/GitHub/ShuttleboxR/Moorea_2024")
# moorea_results <- calc_project_results(moorea_read, exclude_acclimation = T)
# write_csv(moorea_results, "G:/My Drive/MSc2 Frans-Polynesië/ShuttleboxR paper/Data/moorea_results.csv")
moorea_results <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vSUHzgCHYVwoy1ZfIxZSFGstaN6gLBxcDeBYPv6tbeHyCCOcBiRBpoCA473n6tdxXbSEfBTfJTxijEc/pub?gid=1387120886&single=true&output=csv")%>%
  left_join(metadata_moorea %>%
              dplyr::select(file_name, id),
            by = join_by(fileID == file_name))%>%
  mutate(species = case_when(grepl(".d.", id, ignore.case = T) ~ "Damselfish",
                             grepl(".h.", id, ignore.case = T) ~ "Hawkfish",
                             grepl(".c.", id, ignore.case = T) ~ "Croucher"))

pd_moorea%>%
  dplyr::group_by(species)%>%
  dplyr::mutate(mn_shuttles = mean(nr_shuttles))


mo_pca_dat <-moorea_results%>%
  dplyr::select(Tpref, track_accuracy, grav_time, nr_shuttles, t_near_limits, id)

library(ggrepel)
mo_pca <- pca(mo_pca_dat, mahalanobis_th = 0.7,dbscan_th = 1,id_col = "id", var_col = "t_near_limits")
  
mo_pca$outliers
mo_pca$plots$biplot
mo_pca$plots$pc1contributionplot
mo_pca$plots$varplot
  # labs(y = "Temperature preference (°C)")

plot_histograms(moorea_results,nr_sd = 1, iqr_multiplier = 1 , bin_size_distance = 50000, bin_size_shuttles = 500, id_col = "id")

moorea_results%>%
  group_by(species)%>%
  summarise(limits = mean(t_near_limits))
