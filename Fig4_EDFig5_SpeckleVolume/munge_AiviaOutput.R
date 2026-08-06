library(dplyr)
library(readr)

# Set the main directory containing all image folders
main_dir <- "C:/myPath/2026-07-08-12-49-27 IF_analysis/"   # UPDATE THIS PATH

# Specify a particular measurement file name
measurement_filename <- "Measurements_Set_Meshes.csv" 


# Get all subdirectories in the main folder
folders <- list.dirs(main_dir, recursive = FALSE, full.names = TRUE)

cat("Found", length(folders), "folders to process\n\n")


results_list <- list()

for (folder in folders) {
  folder_name <- basename(folder)
  
  # Find .aivia.tif file to derive the image name
  aivia_files <- list.files(folder, pattern = "\\.aivia\\.tif$", full.names = FALSE)
  if (length(aivia_files) == 0) next
  
  # Extract image name (strip the last two extensions, e.g. .aivia.tif)
  image_name <- sub("\\.[^.]+\\.[^.]+$", "", aivia_files[1])
  
  # Locate the Measurements folder
  measurements_dir <- file.path(folder, "Measurements")
  if (!dir.exists(measurements_dir)) next
  
  # Determine which CSV file to read
  if (is.null(measurement_filename)) {
    csv_files <- list.files(measurements_dir, pattern = "\\.csv$", full.names = TRUE)
    if (length(csv_files) == 0) next
    measurement_path <- csv_files[1]
  } else {
    measurement_path <- file.path(measurements_dir, measurement_filename)
    if (!file.exists(measurement_path)) next
  }
  
  # Read the measurement table and tag rows with the image of origin
  df <- read_csv(measurement_path, show_col_types = FALSE) %>%
    mutate(Image = image_name, .before = 1)
  
  results_list[[folder_name]] <- df
}

# COMBINE RESULTS


# Check if any data was collected
# if (length(results_list) == 0) {
#   stop("No measurement tables were successfully read")
# }


# Combine all tables with rbind
master_table <- bind_rows(results_list)


# Save the combined table
setwd("C:/myPath/")
write_csv(master_table, "combined_measurements_2026-07-03-15-22-29 IF_analysis.csv")

