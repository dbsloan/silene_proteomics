library(tidyverse)
library(readxl)
library(purrr)


# Set your working directory
folder_path <- "."

# List all .xlsx files in the directory
file_list <- list.files(folder_path, pattern = "\\.xlsx$", full.names = TRUE)

# Function to process a single Excel file
process_file <- function(file_path) {
  # Read the Excel file
  df <- read_excel(file_path)
  
  # Filter out rows with semicolons in the "Master Protein Accessions" column
  df_filtered <- df %>%
    filter(!grepl(";", `Master Protein Accessions`)) %>%
    group_by(`Master Protein Accessions`) %>%
    summarize(PSM_sum = sum(`# PSMs`, na.rm = TRUE), .groups = "drop")
  
  # Rename the PSM_sum column to the file name (without extension)
  sample_name <- tools::file_path_sans_ext(basename(file_path))
  names(df_filtered)[2] <- sample_name
  
  return(df_filtered)
}

# Process all files and combine them using full_join by gene name
psm_list <- lapply(file_list, process_file)
combined_df <- reduce(psm_list, full_join, by = "Master Protein Accessions")

# Replace NA with 0 if needed
combined_df[is.na(combined_df)] <- 0




