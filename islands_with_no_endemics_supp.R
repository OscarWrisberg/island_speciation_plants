#Packages
packages <- c("devtools", "picante", "phytools","RColorBrewer", "geiger", "readr", "tidyverse", "ggpubr", "ggplot2","ggnewscale", "rnaturalearth", "rnaturalearthdata", "sf", "MetBrewer", "MonoPhy","ggtree","ggrepel", "modelsummary", "MuMIn","MASS","pscl", "cowplot")

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

datadir <- dirname(file.path("../data"))


# Trying to load WCVP distribution data to see if it is better
dist_data <- read.csv(file.path("../data","checklist_distribution.txt"), sep = "|") 


#We probably need to remove the introduced locations.
dist_data_no_introduced <- dist_data[which(dist_data$introduced != 1 & dist_data$extinct != 1 & dist_data$location_doubtful != 1),]


#Soo how many different locations are there in the new database??
areas <- unique(dist_data_no_introduced$area) # 183. Seems doable.
dist_binary <- as.data.frame(dist_data_no_introduced[,1])


# Selecting only islands which I do not have any data from from Melanie.
islands_in_scatter <- output_all_sp_congregated[,2]

# I really do not like hardcoding this...
islands_not_in_scatter <- c("bouvet I.","Chagos Archipelago","Cocos (Keeling) Is.","Crozet Is.","Føroyar",
                            "Galápagos","Greenland","Great Britain","Howland-Baker Is.","Heard-McDonald Is.",
                            "Ireland","Kazan-retto","Laccadive Is.","Mozambique Channel I","Marcus I.","Maldives",
                            "Marshall Is.","Nauru", "Prince Edward I.","Phoenix Is.","Réunion","South China Sea","South Sandwich Is.","Tuvalu")

# Subsetting the distribution data so I only have Islands which are in the scatter plot
dist_data_no_introduced_no_islands_not_scatter <- dist_data_no_introduced[which(!(dist_data_no_introduced$area %in% islands_in_scatter)),]

# Subsetting the distribution data so I only have Islands which are NOT in the scatter plot
dist_data_list_mel<- dist_data_no_introduced[which((dist_data_no_introduced$area %in% islands_not_in_scatter)),]


#Converting data to wide format and 
community_matrixes_all_mel <- dist_data_list_mel %>%
  dplyr::count(plant_name_id, area) %>%
  tidyr::pivot_wider(names_from = plant_name_id, values_from = n, values_fill = 0) %>%
  as.matrix()

# Setting Column names
names(community_matrixes_all_mel) <- community_matrixes_all_mel[1,]
#community_matrixes <- community_matrixes[-1,]

#Setting row names 
row.names(community_matrixes_all_mel) <- community_matrixes_all_mel[,1]
community_matrixes_all_mel <- community_matrixes_all_mel[,-1]

#Checking data overlap between tree and distrubution data. 
community_matrixes_all_mel <- t(community_matrixes_all_mel)

# Now I should be able to filter the other dist database so it has WAY fewer entries.
# So I remove plant names which are not found on any of the islands which are not in the scatterplot
# This means that when we check wether the species on these islands are endemic or not, were not checking species which are not found on these islands.
dist_data_no_introduced_problem_islands <- dist_data_no_introduced[which(dist_data_no_introduced$plant_name_id %in% rownames(community_matrixes_all_mel)),]

#Converting data to wide format and 
community_matrixes_species_on_problem_islands <- dist_data_no_introduced_problem_islands %>%
  dplyr::count(plant_name_id, area) %>%
  tidyr::pivot_wider(names_from = plant_name_id, values_from = n, values_fill = 0) %>%
  as.matrix()

# Setting Column names
names(community_matrixes_species_on_problem_islands) <- community_matrixes_species_on_problem_islands[1,]
#community_matrixes <- community_matrixes[-1,]

#Setting row names 
row.names(community_matrixes_species_on_problem_islands) <- community_matrixes_species_on_problem_islands[,1]
community_matrixes_species_on_problem_islands <- community_matrixes_species_on_problem_islands[,-1]

#Checking data overlap between tree and distrubution data. 
community_matrixes_species_on_problem_islands <- t(community_matrixes_species_on_problem_islands)

class(community_matrixes_species_on_problem_islands) <- "numeric"

# Adding a column which counts the number of botanical countries which in the species occur.
community_matrixes_species_on_problem_islands_total <- cbind(community_matrixes_species_on_problem_islands, rowSums(community_matrixes_species_on_problem_islands[])) 

# Giving the column a name
colnames(community_matrixes_species_on_problem_islands_total)[370] <- "Total"

#Only selecting species which are endemic
community_matrixes_species_on_problem_islands_total_endems <- community_matrixes_species_on_problem_islands_total[which(community_matrixes_species_on_problem_islands_total[,370] == 1),]

#Selecting islands which are in the list of islands not in the scatterplot
community_matrixes_species_on_problem_islands_total_endems_1 <- community_matrixes_species_on_problem_islands_total_endems[,which(colnames(community_matrixes_species_on_problem_islands_total_endems) %in% colnames(community_matrixes_all_mel))]

# Combining 
community_matrixes_species_on_problem_islands_total <- rbind(community_matrixes_species_on_problem_islands_total_endems_1, colSums(community_matrixes_species_on_problem_islands_total_endems_1[])) 

community_matrixes_species_on_problem_islands_total[1175,]

islands_not_in_scatter[which(!(islands_not_in_scatter %in% names(community_matrixes_species_on_problem_islands_total[1175,])))]

# There are zero species recorded here
dist_data_no_introduced[which(dist_data_no_introduced$area == "Mozambique Channel I"),] # there are species found on "Mozambique Channel I
dist_data_no_introduced[which(dist_data_no_introduced$area == "bouvet I.")] # Seems to not exist in World checklist, apparantly it is 97% covered by ice and there is only bryophytes and lichens on the island.
dist_data_no_introduced[which(dist_data_no_introduced$area == "Cocos (Keeling) Is."),] # should be "Cocos (Keeling) Is."
dist_data_no_introduced[which(dist_data_no_introduced$area == "Nauru"),] # should be "Nauru" or "Prince Edward I." ?
dist_data_no_introduced[which(dist_data_no_introduced$area == "Prince Edward I."),]

