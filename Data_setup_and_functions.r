#Packages
packages <- c("devtools", "picante", "phytools","RColorBrewer", "geiger", "readr", "tidyverse", "ggpubr", "ggplot2","ggnewscale", "rnaturalearth", "rnaturalearthdata", "sf", "MetBrewer", "MonoPhy","ggrepel", "modelsummary", "MuMIn","MASS","pscl","png", "magick", "cowplot", "lwgeom", "this.path") # "ggtree",

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

#Setting Directory
setwd(this.path::here())

#####################################################################################################################
############################################---- Loading data ----###################################################
#####################################################################################################################

#Loading data for plotting
all_fles <- list.files(file.path(datadir,"TACT_trees"), pattern = ".rds")
dat_list <- lapply(all_fles, function(x) readRDS(file.path(datadir,"TACT_trees",x))) # Results from Melanies Tacted trees
output_with_island <- readRDS(file.path(datadir,"output_with_island_data.rds"))
#data.table::setnames(output_with_island,'GeologicalOrigin.y','GeologicalOrigin')
island_data_origin <- read_csv(file.path(datadir,"/Islands_TDWG_AllData.txt"))
botanical_countries <- st_read(dsn = file.path(datadir,"wgsrpd-master/level3"), layer = "level3")
dist_data <- read.csv(file.path(datadir,"checklist_distribution.txt"), sep = "|") 
elev_data <- read.csv(file.path(datadir,"bot_countries_maxElev.csv"))


#####################################################################################################################
#########################################---- Massaging the data ----################################################
#####################################################################################################################
# Calculating the lengths of each list
output_all_sp_uniq <- lapply(dat_list, function(x) lapply(x, function(y) lapply(y, unique)))
output_all_sp_lengths <- lapply(output_all_sp_uniq, function(x) lapply(x, function(y) lapply(y, length)))

# Adding 0's to the places where there is no radiation nodes or edges
output_all_sp_lengths <- lapply(output_all_sp_lengths, function(x) lapply(x,function(y) if(length(y)==3) {append(y,c(0,0))}else{y}))

# Here I should be adding the island botanical countries which have neither and make sure they follow the same pattern as the previous entries
# The script called Islands_with_no_endemics_supp.R examines these Islands and makes sure that none of them actually have any endemic species.
# But because none of them had any, we will just add them to the list.

# Lets create a list of the botanical countries to be added.
# Lets start doing this by adding the countries not in the scatter plot
# I really do not like hardcoding this... but I dont think I have any other method.
islands_not_in_scatter <- c("bouvet I.","Chagos Archipelago","Crozet Is.",
                            "Greenland","Great Britain","Howland-Baker Is.","Heard-McDonald Is.",
                            "Ireland","Kazan-retto","Laccadive Is.","Marcus I.","Maldives",
                            "Marshall Is.","Nauru", "Prince Edward I.","Phoenix Is.","Réunion","South China Sea","South Sandwich Is.","Tuvalu")

# This HORRIFIC piece of code selects the short 3-letter code "area_code_l3" of the islands which are missing in the dataset.
three_code_missing_islands <- unique(dist_data[which(dist_data$area %in% islands_not_in_scatter),6])[!(unique(dist_data[which(dist_data$area %in% islands_not_in_scatter),6]) %in% names(output_all_sp_lengths[[1]]))]

# Now I need to add these missing island botanical countries to all the datasets.
for (isl in three_code_missing_islands) {
  print(isl)
  output_all_sp_lengths <- lapply(output_all_sp_lengths, function(x) {
    if (!(isl %in% names(x))) {
      x[[isl]] <- list(0)  # Add a new list with the isl iterator
    }
    x
  })
}

# Now I need to add some 0's
output_all_sp_lengths <- lapply(output_all_sp_lengths, function(x) lapply(x,function(y) if(length(y)==1) {append(y,c(0,0,0,0))}else{y}))

#Converting list of list to df and naming columns
output_as_df_all_sp <- lapply(output_all_sp_lengths, function(x) as.data.frame(t(as.data.frame(do.call(cbind, x)))))

for (x in 1:length(output_as_df_all_sp)){
  colnames(output_as_df_all_sp[[x]]) <- c("colonisation_nodes","colonization_edges", "colonization_sp", "radiating_nodes","radiating_sp")
}

#Changing rownames to be a column
output_as_df_all_sp <- lapply(output_as_df_all_sp, function(x) tibble::rownames_to_column(x,"LEVEL3_COD"))

#Loading the gift dataset
gift_data <- read.csv(file.path(datadir,"env_tdwg_GIFT.csv"), sep = ",")

#Editing Gift dataset to make Haiti and Dominican republic one botanical country.
# Adding the Areas together
gift_data[which(gift_data$LEVEL3_COD == "DOM"),8] <- gift_data[which(gift_data$LEVEL3_COD == "DOM"),8] + gift_data[which(gift_data$LEVEL3_COD == "HAI"),8]
# Selecting the lowest distance to mainland
gift_data[which(gift_data$LEVEL3_COD == "DOM"),9] <- min(gift_data[which(gift_data$LEVEL3_COD == "HAI" | gift_data$LEVEL3_COD == "DOM"),9]) 
# Calculating the mean maximum elevation for the entire island.
gift_data[which(gift_data$LEVEL3_COD == "DOM"),42] <- gift_data[which(gift_data$LEVEL3_COD == "HAI"),42]*(gift_data[which(gift_data$LEVEL3_COD == "HAI"),8]/gift_data[which(gift_data$LEVEL3_COD == "DOM"),8]) + gift_data[which(gift_data$LEVEL3_COD == "DOM"),42]*(1-gift_data[which(gift_data$LEVEL3_COD == "HAI"),8]/gift_data[which(gift_data$LEVEL3_COD == "DOM"),8])


#Merging the datasets
output_as_df_all_sp <- lapply(output_as_df_all_sp, function(x) merge(x, gift_data, by="LEVEL3_COD"))


#####################################################################################################################
#####################################---- Adding Geological Origin ----##############################################
#####################################################################################################################
#Renaming column to be the same as in the gift database
names(island_data_origin)[1] <- "LEVEL3_COD"

#selecting only the geological origin and the name
island_data_origin_merge <- island_data_origin[,c(1,6)]

#Merging geological origin with island counts data.
output_all_sp_origin <- lapply(output_as_df_all_sp, function(x) merge(x, island_data_origin_merge, by="LEVEL3_COD"))


# The venezuelan antilles are not in the Geological origin dataset
which(!(output_as_df_all_sp[[1]][,1] %in% island_data_origin_merge$LEVEL3_COD))
print(output_as_df_all_sp[[1]][119,1])

#####################################################################################################################
#####################################---- Adding Elevation data ----#################################################
#####################################################################################################################

#selecting only the geological origin and the name
island_data_elevation_merge <- elev_data[,c(3,6,7)]

#Merging elevation with island counts data for all species.
output_all_sp <- lapply(output_all_sp_origin, function(x) merge(x, island_data_elevation_merge, by="LEVEL3_COD"))

#Merging elevation with island counts data for coryphoid dataset.
output_with_island <- merge(output_with_island, island_data_elevation_merge, by="LEVEL3_COD")

# Checking if any island is missing
which(!(output_all_sp[[1]][,1] %in% island_data_origin_merge$LEVEL3_COD)) # None of them is

#####################################################################################################################
#####################################---- Fixing the data names ----#################################################
#####################################################################################################################

# There are some problems with Melanies data where some of the island botanical countries have bad names
# Faulty island names
# "F\xf8royar"            "Gal\xe1pagos"          "Juan Fern\xe1ndez Is." "R\xe9union", Mozambique Channel I", "Cocos (Keeling) Is." 

#Should be Faroe, Galapagos, Juan Fernandes Islands, Reunion
# LEVEL3_COD names = FOR GAL JNF REU
bad_names <- c("FOR","GAL","JNF","REU","MCI","CKI")
nnames <- botanical_countries[which(botanical_countries$LEVEL3_COD %in% bad_names),]

# Renaming fields with bad names 
# This function might crash from time to tome because the order of the columns shift ??
for (i in 1:100) {
  # This piece of code finds the 3-letter code of the rows where the LEVEL3_NAM is not in the botanical_countries df
  problems <- output_all_sp[[i]]$LEVEL3_COD[which(!(output_all_sp[[i]]$LEVEL3_NAM %in% botanical_countries$LEVEL3_NAM))]
  #print(problems)
  # Here I change the LEVEL3_NAM of the rows where the LEVEL3_NAM is not in botanical_countries
  output_all_sp[[i]]$LEVEL3_NAM[which(!(output_all_sp[[i]]$LEVEL3_NAM %in% botanical_countries$LEVEL3_NAM))] <- nnames[[1]][which(nnames$LEVEL3_COD %in% problems)]
}

#Now I need to create a subset of the botanical countries which are the islands.
botanical_countries_islands <- botanical_countries[which(botanical_countries$LEVEL3_NAM %in% output_all_sp[[1]]$LEVEL3_NAM),]

botanical_countries_islands$LEVEL3_COD[which(botanical_countries_islands$LEVEL3_COD == "HAI")] <- "DOM" # 122

sf::sf_use_s2(FALSE) # Does not work with sf_use_s2(TRUE)

#Combining Haiti and the Dominican Republic
botanical_countries_islands %>% 
    group_by(LEVEL3_COD) %>%
    summarise(geometry = sf::st_union(geometry)) %>%
    ungroup()


#####################################################################################################################
#####################################---- Calculating fragmentation ----#############################################
#####################################################################################################################

#Setting sf to NOT use S2 distances
sf::sf_use_s2(FALSE)

near_neighbour_dist <- data.frame()

#Calculating the average neighbour distance for all the islands in each island botanical country
for (k in 1:length(botanical_countries_islands[[1]])) { # Looping through island botanical countries
  print(c("starting on ",botanical_countries_islands[[1]][k]))
  k_area <- c()
  test_cast <- st_cast(botanical_countries_islands[k,1], "POLYGON") # Changing the shape file of a island botanical country to several polygons
  for (i in 1:length(test_cast[[1]])) { # Looping through polygons
    k_area <- append(k_area,as.numeric(st_area(test_cast[i,]))/1000000) # Calculating the area of each single polygon
  }
  
  dist_bord <- st_distance(test_cast) #Calculating the distances between the borders of each polygon for each island botanical country
  test_cast_cent <- st_centroid(test_cast) # Finding the centroid for each polygon
  dist_cent <- st_distance(test_cast_cent) # Finding the distances between all centroids
  
  near_neighbour_dist_bord <- apply(dist_bord,1,function (v) {min(v[v>0])}) # Finding min value for each row which is bigger than 0
  near_neighbour_dist_cent <- apply(dist_cent,1,function (v) {min(v[v>0])})
  
  dist_bord <- units::drop_units(dist_bord)
  dist_cent <- units::drop_units(dist_cent)
  
  #Now I need to scale these to the sizes of the islands used
  near_neighbour_dist_bord_scaled <- lapply(near_neighbour_dist_bord, function (x) {x*(k_area[which(dist_bord == near_neighbour_dist_bord[2], arr.ind = TRUE)[1]] + k_area[which(dist_bord == near_neighbour_dist_bord[2], arr.ind = TRUE)[2]])})
  
   near_neighbour_dist_cent_scaled <- lapply(near_neighbour_dist_cent, function (x) {x*(k_area[which(dist_cent == near_neighbour_dist_cent[2], arr.ind = TRUE)[1]] + k_area[which(dist_cent == near_neighbour_dist_cent[2], arr.ind = TRUE)[2]])})
  
  #print(mean(unlist(near_neighbour_dist_bord_scaled)))
  #print(mean(unlist(near_neighbour_dist_cent_scaled)))
  
  near_neighbour_dist[k,1] <- botanical_countries_islands[[1]][k]
  near_neighbour_dist[k,2] <- mean(unlist(near_neighbour_dist_bord_scaled))
  near_neighbour_dist[k,3] <- mean(unlist(near_neighbour_dist_cent_scaled))
  near_neighbour_dist[k,4] <- mean(near_neighbour_dist_bord)
  near_neighbour_dist[k,5] <- mean(near_neighbour_dist_cent)
}

names(near_neighbour_dist) <- c("LEVEL3_NAM", "nearest_neighbour_distance_border_scaled", "nearest_neighbour_distance_centroid_scaled", "neareast_neighbour_bord", "neareast_neighbour_cent") 

output_all_sp <- lapply(output_all_sp, function(x) merge(x, near_neighbour_dist, by="LEVEL3_NAM"))

# - Removing haiti
output_all_sp <- lapply(output_all_sp, function(x) x[-(which(x$LEVEL3_COD == "HAI")),])


# Checking to see if all islands are recovered
which(!(output_all_sp[[1]][,8] %in% near_neighbour_dist$LEVEL3_NAM))
length(output_all_sp[[1]][,1]) # 121

# Converting some columns in my dataframe from lists to be a numeric. 
for (i in 1:100){
  for (k in 3:7)
  output_all_sp[[i]][,k] <- as.numeric(unlist(output_all_sp[[i]][,k]))
}

#####################################################################################################################
#####################################---- Fixing the data names ----#################################################
#####################################################################################################################

# Making tempoary dataframes for the different variables.
tmp_colonisation_nodes <- as.data.frame(sapply(output_all_sp, function(x) unlist(x[["colonisation_nodes"]])))
tmp_colonisation_sp <- as.data.frame(sapply(output_all_sp, function(x) unlist(x[["colonization_sp"]])))
tmp_colonisation_edges <- as.data.frame(sapply(output_all_sp, function(x) unlist(x[["colonization_edges"]])))
tmp_radiating_nodes <- as.data.frame(sapply(output_all_sp, function(x) unlist(x[["radiating_nodes"]])))
tmp_radiating_sp <- as.data.frame(sapply(output_all_sp, function(x) unlist(x[["radiating_sp"]])))

# Creating a list of the tempoary dataframes
tmp_data <- list(tmp_colonisation_edges,tmp_colonisation_nodes,tmp_colonisation_sp,tmp_radiating_nodes,tmp_radiating_sp)

#Setting rownames
tmp_names <- output_all_sp[[1]]$LEVEL3_COD
for (i in 1:length(tmp_data)){
  row.names(tmp_data[[i]]) <- tmp_names
}

# Creating a dataframe where i can congregate the mean, min and max of all the valeiues
output_all_sp_congregated<-output_all_sp[[1]]

for (i in 1:length(tmp_data)){
  loop_df <- data.frame()
  vars <- c("colonisation_nodes","colonization_sp", "colonization_edges", "radiating_nodes", "radiating_sp")
  ranges <- c("mean","min","max")
  loop_mean <- tibble::rownames_to_column(as.data.frame(apply(tmp_data[[i]],1,FUN = mean)), "LEVEL3_COD")
  print(loop_mean)
  loop_min <- tibble::rownames_to_column(as.data.frame(apply(tmp_data[[i]],1,FUN = min)), "LEVEL3_COD")
  loop_max <- tibble::rownames_to_column(as.data.frame(apply(tmp_data[[i]],1,FUN = max)), "LEVEL3_COD")
  loop_df <- cbind(loop_mean)
  loop_df <- merge(loop_df, loop_min, by="LEVEL3_COD")
  loop_df <- merge(loop_df, loop_max, by="LEVEL3_COD")
  names(loop_df) <- c("LEVEL3_COD",paste(vars[i],ranges[1], sep = "_"),paste(vars[i],ranges[2], sep = "_"),paste(vars[i],ranges[3], sep = "_"))
    
  output_all_sp_congregated <- merge(output_all_sp_congregated, loop_df, by="LEVEL3_COD")
}

output_all_sp_test <- output_all_sp

# This is where I change NA's to 0's in for nearest_neighbour_distance_border_scaled
for (i in 1:100) {
  output_all_sp_test[[i]]$nearest_neighbour_distance_border_scaled[is.na(output_all_sp_test[[i]]$nearest_neighbour_distance_border_scaled)] <- 0
}

# Removing introduced and doubtfull locations
dist_data <- dist_data[which(dist_data$introduced != 1 & dist_data$location_doubtful != 1 ),]

#So now I need to look through the islands that I have and count the total number of species on the islands
total_island_sp <- data.frame()


for (i in 1:length(botanical_countries$LEVEL3_COD)){
  total_island_sp[i,1] <- botanical_countries$LEVEL3_COD[i]
  total_island_sp[i,2] <- length(which(dist_data$area_code_l3 == botanical_countries$LEVEL3_COD[i]))
}

names(total_island_sp) <- c("LEVEL3_COD","Total_sp")

output_all_sp_test <- lapply(output_all_sp_test, function(x) merge(x, total_island_sp, by="LEVEL3_COD"))

saveRDS(output_all_sp_test, file.path(datadir,"output_for_test_all_sp")) 

output_all_sp_test_subset_log <- list()

# Creating a subset of the values I need
for (i in 1:100){
output_all_sp_test_subset_log[[i]] <- subset(output_all_sp_test[[i]], select = c(LEVEL3_NAM,radiating_sp,radiating_nodes,colonization_sp,Total_sp,mean_mx30_grd, area, dist,nearest_neighbour_distance_border_scaled,GeologicalOrigin,max_elev_30m))
output_all_sp_test_subset_log[[i]]$mean_mx30_grd <- log10(output_all_sp_test_subset_log[[i]]$mean_mx30_grd)
output_all_sp_test_subset_log[[i]]$max30_elev <- log10(output_all_sp_test_subset_log[[i]]$max_elev_30m+1)
output_all_sp_test_subset_log[[i]]$area <- log10(output_all_sp_test_subset_log[[i]]$area)
output_all_sp_test_subset_log[[i]]$dist <- log10(output_all_sp_test_subset_log[[i]]$dist+1)
output_all_sp_test_subset_log[[i]]$nearest_neighbour_distance_border_scaled <- log10(output_all_sp_test_subset_log[[i]]$nearest_neighbour_distance_border_scaled+1)
}

saveRDS(output_all_sp_test_subset_log, file.path(datadir,"output_for_test_all_sp_subset_log"))

# Creating a subset for Coryphoideae
output_coryphoideae_test_subset_log <- list()

# Creating a subset of the values I need
output_coryphoideae_test_subset_log <- subset(output_with_island, select = c(LEVEL3_NAM,`No. edges`,`No. cladogenesis`,`No. anagenesis`,mean_mx30_grd, area, dist,nearest_neighbour_distance_border_scaled,GeologicalOrigin, max_elev_30m))

output_coryphoideae_test_subset_log$max30_elev <- log10(output_coryphoideae_test_subset_log$max_elev_30m+1)
output_coryphoideae_test_subset_log$area <- log10(output_coryphoideae_test_subset_log$area)
output_coryphoideae_test_subset_log$dist <- log10(output_coryphoideae_test_subset_log$dist+1)
output_coryphoideae_test_subset_log$nearest_neighbour_distance_border_scaled <- log10(output_coryphoideae_test_subset_log$nearest_neighbour_distance_border_scaled+1)
output_coryphoideae_test_subset_log[which(is.na(output_coryphoideae_test_subset_log$nearest_neighbour_distance_border_scaled)),8] <- 0


saveRDS(output_coryphoideae_test_subset_log, file.path(datadir,"output_for_test_coryphoideae_subset_log"))


#####################################################################################################################
#####################################---- Creating the functions ----#################################################
#####################################################################################################################

# Modifying qdredge to work with quasi poisson
# modify a glm() output object so that
# it contains a quasipoisson fit but the 
# AIC (likelihood) from the equivalent regular Poisson model
x.quasipoisson <- function(...) {
res <- quasipoisson(...)
res$aic <- poisson(...)$aic
res
}

# function to extract the overdispersion parameter
# from a quasi model
dfun <- function(object) {
with(object,sum((weights * residuals^2)[weights > 0])/df.residual)
}

# function that modifies MuMIn::dredge() 
# for use with quasi GLM
qdredge <- function(model, family=x.quasipoisson, na.action=na.fail, chat = dfun(model), rank='QAIC', ...){
  model2 <- update(model, family=family, na.action=na.action)
  (dt <- dredge(model2, rank=rank, chat=chat, ...))
}

# modify a glm() output object so that
# it contains a quasipoisson fit but the 
# AIC (likelihood) from the equivalent regular binomial model
x.quasibinomial <- function(...) {
res <- quasibinomial(...)
res$aic <- binomial(...)$aic
res
}

# function to extract the overdispersion parameter
# from a quasi model
dfun <- function(object) {
with(object,sum((weights * residuals^2)[weights > 0])/df.residual)
}

# function that modifies MuMIn::dredge() 
# for use with quasi GLM
qdredge_bin <- function(model, family=x.quasibinomial, na.action=na.fail, chat = dfun(model), rank='QAIC', ...){
  model2 <- update(model, family=family, na.action=na.action)
  (dt <- dredge(model2, rank=rank, chat=chat, ...))
}

# Function which counts which model is the best
# function that modifies MuMIn::dredge() 
# for use with quasi GLM
best_model_counts <- function(input_formula, dataset_list, family_x){
print(input_formula)
best_model <- list()
best_model_vars <- data.frame()

for (i in 1:100){
test_dat_log <- dataset_list[[i]]


glm_model_for_count <- glm(input_formula, data = test_dat_log, family = family_x, na.action = na.fail)

print(glm_model_for_count)

dredge_output <- qdredge(model = glm_model_for_count)

print(dredge_output)

best_model[[i]] <- (dredge_output[which(dredge_output$delta == 0),1:10])
best_model_vars <- rbind(best_model_vars, dredge_output[which(dredge_output$delta == 0),1:10])

}

best_model_vars[which(!is.na(best_model_vars[,1:10])),1:10]
best_model_vars_test <- is.na(best_model_vars)
best_model_vars_count <- best_model_vars_test %>% as.data.frame() %>% group_by_all() %>% count()

return(best_model_vars_count)
}

# Functions which adds data to the plots.
AddVars_mean <- function(dat_set, term, type){
  output_df <- dat_set
  nr_rows <- length(dat_set[1,])
  cl_names <- names(dat_set)
  
  # I need to check if all the combinations of terms and type are in the dataset and if they are not, we add a row of the missing combinations just filled with NA
  for (i in type){
      for ( j in term) {
        exists_comb <- subset(dat_set, term == j & type == i)
        #print(exists_comb)
        #print(subset(dat_set, term == j & type == i))
        if (nrow(exists_comb) == 0) {
          #print(paste("The combination of",j,"and",i,"does not exists in the dataframe", sep = " "))
          na_row <- c(j,NA,NA,NA,NA,NA,NA,i)
          names(na_row) <- cl_names
          output_df <- rbind(output_df, na_row)
        } else {
          #print(paste("The combination of",j,"and",i,"exists in the dataframe", sep = " "))
        }
    }
  }
  return(output_df)
}

AddVars_violin <- function(dat_set, term, type){
  output_df <- dat_set
  nr_rows <- length(dat_set[1,])
  cl_names <- names(dat_set)
  
  # I need to check if all the combinations of terms and type are in the dataset and if they are not, we add a row of the missing combinations just filled with NA
  for (i in type){
      for ( j in term) {
        exists_comb <- subset(dat_set, term == j & type == i)
        #print(exists_comb)
        #print(subset(dat_set, term == j & type == i))
        if (nrow(exists_comb) == 0) {
          #print(paste("The combination of",j,"and",i,"does not exists in the dataframe", sep = " "))
          na_row <- c(j,0,0,0,0,0,0,i)
          names(na_row) <- cl_names
          for (x in range(10)) {
            output_df <- rbind(output_df, na_row)
          }
        } else {
          #print(paste("The combination of",j,"and",i,"exists in the dataframe", sep = " "))
        }
    }
  }
  return(output_df)
}


get_legend <- function(plot, legend = NULL) {
  
  gt <- ggplotGrob(plot)
  
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  
  indices <- grep(pattern, gt$layout$name)

  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}

# List of terms for the plots and models 
terms_list <- c("Area","Isolation","Isolation Squared", "Max Elevation","Mixed Origin","Continental","Fragmentation")
type_list <- c("Between-region speciation","Number Endemics","Within-region speciation")
