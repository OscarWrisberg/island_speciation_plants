# Making a pipeline script which runs all the scripts in the correct order

# Loading the necessary packages
#Packages
packages <- c("devtools", "picante","data.table", "phytools","RColorBrewer", "geiger", "readr", "tidyverse", "ggpubr", "ggplot2","ggnewscale", "rnaturalearth", "rnaturalearthdata", "sf", "MetBrewer", "MonoPhy","ggrepel", "modelsummary", "MuMIn","MASS","pscl","png", "magick", "cowplot", "lwgeom", "this.path") # "ggtree",

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

#Packages
if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
library(BiocManager)

if (!require("ggtree", quietly = TRUE)){
    BiocManager::install("ggtree")
}
library(ggtree)

# Start by setting the working directory as the directory in which this file is saved
# Specify the file path
file_path <- "/home/au543206/Documents/island_speciation_plants/pipeline.r"
output_folder <- "/home/au543206/Documents/island_speciation_plants/Figures"

# Get the directory name
dir_path <- dirname(file_path)

# Set the folder for the raw data
datadir <- "/home/au543206/Documents/Coryphoideae/Figures_in_r/data"

# Set the working directory to the file's directory
setwd(dir_path)

# Then we run the Coryphoideae section.
source("Finding_speciation_events_coryphoideae.r")

# Start by running the Data_setup_and_functions.Rmd
source("Data_setup_and_functions.r")

# And the same section for the phylogenetic trees of all plants
source("Finding_speciation_events_all_plants.r")

# Then we run the scripts for the model fitting of all the other models
source("Model_fitting_other_models.r")

# Running the script which creates the Coryphoideae phylogeny
source("Coryphoideae_phylogeny_plot.r")

# Running the plotting script to create all the figues
source("plotting_script.R")