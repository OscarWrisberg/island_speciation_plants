# IN this script I will try to calculate the max height for all botanical countries.


#Packages
packages <- c("devtools", "picante", "phytools","RColorBrewer", "geiger", "readr",
              "tidyverse", "ggpubr", "ggplot2","ggnewscale", "rnaturalearth", "rnaturalearthdata", "sf", "MetBrewer",
              "MonoPhy","ggtree","ggrepel", "modelsummary", "MuMIn","MASS","pscl","png", "magick", "cowplot", "raster")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

#Setting Directory
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
getwd()


datadir <- dirname(file.path("../data"))


# Loading some data
botanical_countries <- st_read(dsn = file.path("../data","wgsrpd-master/level3"), layer = "level3")