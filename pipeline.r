# Making a pipeline script which runs all the scripts in the correct order

# Start by setting the working directory as the directory in which this file is saved
# Specify the file path
file_path <- "path/to/your/example_script.Rmd"
file_path <- "pipeline.r"


# Get the directory name
dir_path <- dirname("pipeline.r")

# Set the working directory to the file's directory
setwd(dir_path)



# Print the working directory to verify
print(getwd())




# Start by running the Data_setup_and_functions.Rmd
source("Data_setup_and_functions.Rmd")
