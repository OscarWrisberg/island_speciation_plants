# Making a pipeline script which runs all the scripts in the correct order

# Start by setting the working directory as the directory in which this file is saved
# Specify the file path
file_path <- "/home/au543206/Documents/island_speciation_plants/pipeline.r"

# Get the directory name
dir_path <- dirname(file_path)

# Set the folder for the raw data
datadir <- "/home/au543206/Documents/Coryphoideae/Figures_in_r/data"

# Set the working directory to the file's directory
setwd(dir_path)

# Start by running the Data_setup_and_functions.Rmd
source("Data_setup_and_functions.r")
