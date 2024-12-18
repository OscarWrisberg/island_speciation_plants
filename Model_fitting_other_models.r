# Loading Data
output_all_sp_test_subset_log <- readRDS(file.path(datadir,"output_for_test_all_sp_subset_log"))
test_dat_log_cory <- readRDS(file.path(datadir,"output_for_test_coryphoideae_subset_log"))


## Fitting the models 

# Create a dataframe which contains the mean radiating sp, colonizing sp radiating nodes and total sp.
# Making tempoary dataframes for the different variables.
tmp_colonisation_sp <- as.data.frame(sapply(output_all_sp_test_subset_log, function(x) unlist(x[["colonization_sp"]])))
tmp_radiating_nodes <- as.data.frame(sapply(output_all_sp_test_subset_log, function(x) unlist(x[["radiating_nodes"]])))
tmp_radiating_sp <- as.data.frame(sapply(output_all_sp_test_subset_log, function(x) unlist(x[["radiating_sp"]])))

# Creating a list of the tempoary dataframes
tmp_data <- list(tmp_colonisation_sp,tmp_radiating_nodes,tmp_radiating_sp)

#Setting rownames
tmp_names <- output_all_sp_test_subset_log[[1]][,1]
for (i in 1:length(tmp_data)){
  row.names(tmp_data[[i]]) <- tmp_names
}
# Creating a dataframe where i can congregate the mean, min and max of all the valeiues
output_mean<-output_all_sp_test_subset_log[[1]]

for (i in 1:length(tmp_data)){
  loop_df <- data.frame()
  vars <- c("colonization_sp", "radiating_nodes", "radiating_sp")
  ranges <- c("mean","min","max")
  loop_mean <- tibble::rownames_to_column(as.data.frame(apply(tmp_data[[i]],1,FUN = mean)), "LEVEL3_NAM")
  loop_min <- tibble::rownames_to_column(as.data.frame(apply(tmp_data[[i]],1,FUN = min)), "LEVEL3_NAM")
  loop_max <- tibble::rownames_to_column(as.data.frame(apply(tmp_data[[i]],1,FUN = max)), "LEVEL3_NAM")
  loop_df <- cbind(loop_mean)
  loop_df <- merge(loop_df, loop_min, by="LEVEL3_NAM")
  loop_df <- merge(loop_df, loop_max, by="LEVEL3_NAM")
  names(loop_df) <- c("LEVEL3_NAM",paste(vars[i],ranges[1], sep = "_"),paste(vars[i],ranges[2], sep = "_"),paste(vars[i],ranges[3], sep = "_"))
    
  output_mean <- merge(output_mean, loop_df, by="LEVEL3_NAM")
}


output_mean





## All Seed plants all islands
# Vars
# area     dist    nearest_neighbour_distance_border_scaled    GeologicalOrigin    max30_elev

# Anagenesis
# Ana model 2
glm_ana_pois_all_sp_2 <- glm(
                              colonization_sp_mean+radiating_nodes_mean~area+poly(dist,2,raw=TRUE)+max30_elev,
                              data=output_mean,
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_ana_pois_all_sp_2)

# Ana model 3
glm_ana_pois_all_sp_3 <- glm(
                              colonization_sp_mean+radiating_nodes_mean~area+poly(dist,2,raw=TRUE)+max30_elev+nearest_neighbour_distance_border_scaled+relevel(factor(GeologicalOrigin), ref = "volcanic"),
                              data=output_mean,
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_ana_pois_all_sp_3)

# Ana model 4
glm_ana_pois_all_sp_4 <- glm(
                              colonization_sp_mean+radiating_nodes_mean~area+poly(dist,2,raw=TRUE)+nearest_neighbour_distance_border_scaled,
                              data=output_mean,
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_ana_pois_all_sp_4)

#### Cladogenesis
glm_clad_pois_all_sp_2 <- glm(
                              radiating_sp_mean-radiating_nodes_mean~area+poly(dist,1,raw=TRUE)+max30_elev+nearest_neighbour_distance_border_scaled+relevel(factor(GeologicalOrigin), ref = "volcanic"),
                              data=output_mean[which(output_mean$GeologicalOrigin != "atoll"),],
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_clad_pois_all_sp_2)

glm_clad_pois_all_sp_3 <- glm(
                              radiating_sp_mean-radiating_nodes_mean~area+poly(dist,1,raw=TRUE)+nearest_neighbour_distance_border_scaled+relevel(factor(GeologicalOrigin), ref = "volcanic"),
                              data=output_mean[which(output_mean$GeologicalOrigin != "atoll"),],
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_clad_pois_all_sp_3)

glm_clad_pois_all_sp_4 <- glm(
                              radiating_sp_mean-radiating_nodes_mean~area+max30_elev+nearest_neighbour_distance_border_scaled+relevel(factor(GeologicalOrigin), ref = "volcanic"),
                              data=output_mean[which(output_mean$GeologicalOrigin != "atoll"),],
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_clad_pois_all_sp_4)



# Endemics
# endems model 2
glm_endems_pois_all_sp_2 <- glm(
                              colonization_sp_mean+radiating_nodes_mean+radiating_sp_mean~area+poly(dist,2,raw=TRUE)+max30_elev,
                              data=output_mean,
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_endems_pois_all_sp_2)

# Endems model 3
glm_endems_pois_all_sp_3 <- glm(
                              colonization_sp_mean+radiating_nodes_mean+radiating_sp_mean~area+poly(dist,2,raw=TRUE)+max30_elev+nearest_neighbour_distance_border_scaled+relevel(factor(GeologicalOrigin), ref = "volcanic"),
                              data=output_mean,
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_endems_pois_all_sp_3)

# Endems model 4
glm_endems_pois_all_sp_4 <- glm(
                              colonization_sp_mean+radiating_nodes_mean+radiating_sp_mean~area+poly(dist,2,raw=TRUE)+nearest_neighbour_distance_border_scaled,
                              data=output_mean,
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_endems_pois_all_sp_4)





## All Seed plants only volcanic islands
# Anagenesis
# Ana model 2
glm_ana_pois_all_sp_volc_2 <- glm(
                              colonization_sp_mean+radiating_nodes_mean~area+poly(dist,2,raw=TRUE)+max30_elev+nearest_neighbour_distance_border_scaled,
                              data=output_mean[which(output_mean$GeologicalOrigin == "volcanic"),],
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_ana_pois_all_sp_volc_2)

# Ana model 3
glm_ana_pois_all_sp_volc_3 <- glm(
                              colonization_sp_mean+radiating_nodes_mean~area+poly(dist,2,raw=TRUE),
                              data=output_mean[which(output_mean$GeologicalOrigin == "volcanic"),],
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_ana_pois_all_sp_volc_3)

# Ana model 4
glm_ana_pois_all_sp_volc_4 <- glm(
                              colonization_sp_mean+radiating_nodes_mean~area,
                              data=output_mean[which(output_mean$GeologicalOrigin == "volcanic"),],
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_ana_pois_all_sp_volc_4)

#### Cladogenesis
#Clad model 2
glm_clad_pois_all_sp_volc_2 <- glm(
                              radiating_sp_mean-radiating_nodes_mean~area+poly(dist,1,raw=TRUE)+max30_elev+nearest_neighbour_distance_border_scaled,
                              data=output_mean[which(output_mean$GeologicalOrigin == "volcanic"),],
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_clad_pois_all_sp_volc_2)

#Clad model 3
glm_clad_pois_all_sp_volc_3 <- glm(
                              radiating_sp_mean-radiating_nodes_mean~poly(dist,1,raw=TRUE)+max30_elev+nearest_neighbour_distance_border_scaled,
                              data=output_mean[which(output_mean$GeologicalOrigin == "volcanic"),],
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_clad_pois_all_sp_volc_3)

#Clad model 4
glm_clad_pois_all_sp_volc_4 <- glm(
                              radiating_sp_mean-radiating_nodes_mean~poly(dist,1,raw=TRUE)+max30_elev,
                              data=output_mean[which(output_mean$GeologicalOrigin == "volcanic"),],
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_clad_pois_all_sp_volc_4)


# Endemics
# endems model 2
glm_endems_pois_all_sp_volc_2 <- glm(
                              colonization_sp_mean+radiating_nodes_mean+radiating_sp_mean~area+poly(dist,2,raw=TRUE)+max30_elev+nearest_neighbour_distance_border_scaled,
                              data=output_mean[which(output_mean$GeologicalOrigin == "volcanic"),],
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_endems_pois_all_sp_volc_2)

# Endems model 3
glm_endems_pois_all_sp_volc_3 <- glm(
                              colonization_sp_mean+radiating_nodes_mean+radiating_sp_mean~area+poly(dist,2,raw=TRUE),
                              data=output_mean[which(output_mean$GeologicalOrigin == "volcanic"),],
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_endems_pois_all_sp_volc_3)

# Endems model 4
glm_endems_pois_all_sp_volc_4 <- glm(
                              colonization_sp_mean+radiating_nodes_mean+radiating_sp_mean~area,
                              data=output_mean[which(output_mean$GeologicalOrigin == "volcanic"),],
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_endems_pois_all_sp_volc_4)


## Coryphoideae
# Anagenesis
# Ana model 2
glm_ana_pois_cory_2 <- glm(
                              `No. anagenesis`+`No. edges`~area,
                              data=test_dat_log_cory,
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_ana_pois_cory_2)

# Ana model 3
glm_ana_pois_cory_3 <- glm(
                              `No. anagenesis`+`No. edges`~area+max30_elev+nearest_neighbour_distance_border_scaled,
                              data=test_dat_log_cory,
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_ana_pois_cory_3)

# Ana model 4
glm_ana_pois_cory_4 <- glm(
                              `No. anagenesis`+`No. edges`~area+nearest_neighbour_distance_border_scaled+relevel(factor(GeologicalOrigin), ref = "volcanic"),
                              data=test_dat_log_cory,
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_ana_pois_cory_4)

#### Cladogenesis
#Clad model 2
glm_clad_pois_cory_2 <- glm(
                              `No. cladogenesis`-`No. edges`~area+poly(dist,1,raw=TRUE),
                              data=test_dat_log_cory,
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_clad_pois_cory_2)

#Clad model 3
glm_clad_pois_cory_3 <- glm(
                              `No. cladogenesis`-`No. edges`~area+nearest_neighbour_distance_border_scaled,
                              data=test_dat_log_cory,
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_clad_pois_cory_3)

#Clad model 4
glm_clad_pois_cory_4 <- glm(
                              `No. cladogenesis`-`No. edges`~area+max30_elev,
                              data=test_dat_log_cory,
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_clad_pois_cory_4)


# Endemics
# Endems model 2
glm_endems_pois_cory_2 <- glm(
                              `No. cladogenesis`+`No. anagenesis`~area+nearest_neighbour_distance_border_scaled,
                              data=test_dat_log_cory,
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_endems_pois_cory_2)

# Endems model 3
glm_endems_pois_cory_3 <- glm(
                              `No. cladogenesis`+`No. anagenesis`~area+max30_elev,
                              data=test_dat_log_cory,
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_endems_pois_cory_3)

# Endems model 4
glm_endems_pois_cory_4 <- glm(
                              `No. cladogenesis`+`No. anagenesis`~area+relevel(factor(GeologicalOrigin), ref = "volcanic"),
                              data=test_dat_log_cory,
                              family = quasipoisson(link="log"),na.action = na.fail)

summary(glm_endems_pois_cory_4)
