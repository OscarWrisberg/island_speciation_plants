# Loading the data
output_all_sp_test_subset_log <- readRDS(file.path(datadir,"output_for_test_all_sp_subset_log"))
test_dat_log_cory <- readRDS(file.path(datadir,"output_for_test_coryphoideae_subset_log"))



# Finding the best model for the cladogenesis
best_model_vars_clad <- data.frame()
second_best_model_clad  <- data.frame()
third_best_model_clad  <- data.frame()
fourth_best_model_clad  <- data.frame()

defaultW <- getOption("warn")
options(warn = -1)

for (i in 1:100){
test_dat_log <- output_all_sp_test_subset_log[[i]]

dredge_output <- qdredge(glm(radiating_sp-radiating_nodes~max30_elev+area+poly(dist,1, raw=TRUE)+(nearest_neighbour_distance_border_scaled+1)+GeologicalOrigin,data=test_dat_log,family = quasipoisson,na.action = na.fail))

best_model_vars_clad <- rbind(best_model_vars_clad, dredge_output[which(dredge_output$delta == 0),1:9])

second_best_model_clad  <- rbind(second_best_model_clad,dredge_output[2,1:10])
third_best_model_clad  <- rbind(third_best_model_clad, dredge_output[3,1:10])
fourth_best_model_clad <- rbind(fourth_best_model_clad,dredge_output[4,1:10])

}

options(warn = defaultW)

#best_model_vars_test <- is.na(best_model_vars)
best_model_vars_count_clad <- is.na(best_model_vars_clad) %>% as.data.frame() %>% group_by_all() %>% count()
best_model_vars_count_clad$mean_delta <- 0


second_best_model_clad_count  <- is.na(second_best_model_clad[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()
third_best_model_clad_count  <- is.na(third_best_model_clad[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()
fourth_best_model_clad_count <-is.na(fourth_best_model_clad[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()

second_best_model_clad_count$mean_delta  <- mean(second_best_model_clad[,10])
third_best_model_clad_count$mean_delta <- mean(third_best_model_clad[,10])
fourth_best_model_clad_count$mean_delta <- mean(fourth_best_model_clad[,10])

####

clad_table <- rbind(best_model_vars_count_clad,second_best_model_clad_count,third_best_model_clad_count,fourth_best_model_clad_count)

clad_table <- as.data.frame(clad_table)

clad_table[,1:6][clad_table[,1:6]==TRUE] <- "Drop"
clad_table[,1:6][clad_table[,1:6]==FALSE] <- "Keep"

# -------- Cladogenesis ----------------
glm_clad_pois_all_sp_list_tidy <- list()
glm_clad_pois_all_sp_list_models <- list()
glm_clad_pois_all_sp_list_dispersion <- list()

# Loop which fit the model to each of the 100 datasets
for (i in 1:100){
  glm_clad_pois_all_sp <- glm(radiating_sp-radiating_nodes~area+poly(dist,1,raw=TRUE)+nearest_neighbour_distance_border_scaled,data=output_all_sp_test_subset_log[[i]][which(output_all_sp_test_subset_log[[i]]$GeologicalOrigin != "atoll"),],family = quasipoisson(link="log"),na.action = na.fail)
  glm_clad_pois_all_sp_list_models[[i]] <- glm_clad_pois_all_sp
  glm_clad_pois_all_sp_list_tidy[[i]] <- broom::tidy(glm_clad_pois_all_sp)
 glm_clad_pois_all_sp_list_dispersion[[i]] <- sum(residuals(glm_clad_pois_all_sp,type ="pearson")^2)/glm_clad_pois_all_sp$df.residual
}

# Removing datasets with the one cladogenetic speciation event on the bahamas
glm_clad_pois_all_sp_list_models <- glm_clad_pois_all_sp_list_models[-c(8,26,31,41)]

#Now lets find the confidence intervals for these models.
fit_clad_all_95_list <- list()
for (i in 1:96){
  fit_clad_all_95 <- confint(glm_clad_pois_all_sp_list_models[[i]], level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")
  fit_clad_all_95_list[[i]] <- fit_clad_all_95
  }

# Combining the 95 intervals on all the models
results_clad_all_list <- list()
for (i in 1:96){
  results_clad_all_list[[i]] <- dplyr::bind_cols(glm_clad_pois_all_sp_list_tidy[[i]],fit_clad_all_95_list[[i]]) 
}

# Quickly checking the VIF parameters
car::vif(glm_clad_pois_all_sp_list_models[[1]])

# Finding the best model for the anagenesis
best_model_vars_ana <- data.frame()
second_best_model_ana  <- data.frame()
third_best_model_ana  <- data.frame()
fourth_best_model_ana  <- data.frame()

for (i in 1:100){
test_dat_log <- output_all_sp_test_subset_log[[i]]

dredge_output <- qdredge(glm(colonization_sp+radiating_nodes~max30_elev+area+poly(dist,2, raw = TRUE)+nearest_neighbour_distance_border_scaled+GeologicalOrigin,data=test_dat_log,family = quasipoisson,na.action = na.fail))

row_to_add <- dredge_output[which(dredge_output$delta == 0),1:9]

best_model_vars_ana <- rbind(best_model_vars_ana, row_to_add)

second_best_model_ana  <- rbind(second_best_model_ana,dredge_output[2,1:10])
third_best_model_ana  <- rbind(third_best_model_ana, dredge_output[3,1:10])
fourth_best_model_ana <- rbind(fourth_best_model_ana,dredge_output[4,1:10])

}

#best_model_vars_test <- is.na(best_model_vars)
best_model_vars_count_ana <- is.na(best_model_vars_ana[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()
best_model_vars_count_ana$mean_delta <- 0

second_best_model_count_ana  <- is.na(second_best_model_ana[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()
third_best_model_count_ana  <- is.na(third_best_model_ana[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()
fourth_best_model_count_ana <-is.na(fourth_best_model_ana[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()

second_best_model_count_ana$mean_delta  <- mean(second_best_model_ana[,10])
third_best_model_count_ana$mean_delta <- mean(third_best_model_ana[,10])
fourth_best_model_count_ana$mean_delta <- mean(fourth_best_model_ana[,10])

ana_table <- rbind(best_model_vars_count_ana,second_best_model_count_ana,third_best_model_count_ana,fourth_best_model_count_ana)

ana_table <- as.data.frame(ana_table)

ana_table[,1:6][ana_table[,1:6]==TRUE] <- "Drop"
ana_table[,1:6][ana_table[,1:6]==FALSE] <- "Keep"

ana_table

# Creating a for loop which fits the model using all 100 different datasets 
# -------- Anagenesis ----------------
glm_ana_pois_all_sp_list_tidy <- list()
glm_ana_pois_all_sp_list_models <- list()
glm_ana_pois_all_sp_list_dispersion <- list()

for (i in 1:100){
  glm_ana_all_sp_loop <- glm(colonization_sp+radiating_nodes~area+poly(dist,2, raw=TRUE)+nearest_neighbour_distance_border_scaled+max30_elev,data=output_all_sp_test_subset_log[[i]],family = quasipoisson)
  glm_ana_pois_all_sp_list_models[[i]] <- glm_ana_all_sp_loop
  glm_ana_pois_all_sp_list_tidy[[i]] <- broom::tidy(glm_ana_all_sp_loop)
  glm_ana_pois_all_sp_list_dispersion[[i]] <- sum(residuals(glm_ana_all_sp_loop,type ="pearson")^2)/glm_ana_all_sp_loop$df.residual
}

#Now lets find the confidence intervals for these models.
fit_ana_all_95_list <- list()
for (i in 1:100){
  fit_ana_all_95 <- confint(glm_ana_pois_all_sp_list_models[[i]], level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")
  fit_ana_all_95_list[[i]] <- fit_ana_all_95
  }

# Combining the 95 intervals on all the models
results_ana_all_list <- list()
for (i in 1:100){
  results_ana_all_list[[i]] <- dplyr::bind_cols(glm_ana_pois_all_sp_list_tidy[[i]],fit_ana_all_95_list[[i]]) 
}

car::vif(glm_ana_pois_all_sp_list_models[[1]])


# Finding the best model for endemics

# Finding the best model for the cladogenesis
best_model_vars_endems <- data.frame()
second_best_model_endems  <- data.frame()
third_best_model_endems  <- data.frame()
fourth_best_model_endems <- data.frame()

for (i in 1:100){
test_dat_log <- output_all_sp_test_subset_log[[i]]

dredge_output <- qdredge(glm(colonization_sp+radiating_sp~max30_elev+area+poly(dist,2, raw = TRUE)+nearest_neighbour_distance_border_scaled+GeologicalOrigin,data=test_dat_log,family = quasipoisson,na.action = na.fail))

row_to_add <- dredge_output[which(dredge_output$delta == 0),1:10]

best_model_vars_endems <- rbind(best_model_vars_endems, row_to_add)

second_best_model_endems  <- rbind(second_best_model_endems,dredge_output[2,1:10])
third_best_model_endems  <- rbind(third_best_model_endems, dredge_output[3,1:10])
fourth_best_model_endems <- rbind(fourth_best_model_endems,dredge_output[4,1:10])

}

#best_model_vars_test <- is.na(best_model_vars)
best_model_vars_count_endems <- is.na(best_model_vars_endems[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()
best_model_vars_count_endems$mean_delta <- 0

second_best_model_endems_count  <- is.na(second_best_model_endems[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()
third_best_model_endems_count  <- is.na(third_best_model_endems[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()
fourth_best_model_endems_count <-is.na(fourth_best_model_endems[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()

second_best_model_endems_count$mean_delta  <- mean(second_best_model_endems[,10])
third_best_model_endems_count$mean_delta <- mean(third_best_model_endems[,10])
fourth_best_model_endems_count$mean_delta <- mean(fourth_best_model_endems[,10])

endems_table <- rbind(best_model_vars_count_endems,second_best_model_endems_count,third_best_model_endems_count,fourth_best_model_endems_count)

endems_table <- as.data.frame(endems_table)

endems_table[,1:6][endems_table[,1:6]==TRUE] <- "Drop"
endems_table[,1:6][endems_table[,1:6]==FALSE] <- "Keep"

endems_table



# Endemics

# -------- Endemics ----------------
glm_endems_pois_all_sp_list_tidy <- list()
glm_endems_pois_all_sp_list_models <- list()
best_model_vars_endems <- data.frame()
second_best_model_endems  <- data.frame()
third_best_model_endems  <- data.frame()
fourth_best_model_endems  <- data.frame()

defaultW <- getOption("warn")
options(warn = -1)

for (i in 1:100){
  glm_endems_all_sp_loop <- glm(colonization_sp+radiating_sp~area+poly(dist,2, raw=TRUE)+nearest_neighbour_distance_border_scaled+max30_elev,data=output_all_sp_test_subset_log[[i]],family = quasipoisson)
  glm_endems_pois_all_sp_list_models[[i]] <- glm_endems_all_sp_loop
  glm_endems_pois_all_sp_list_tidy[[i]] <- broom::tidy(glm_endems_all_sp_loop)
}

#Now lets find the confidence intervals for these models.
fit_endems_all_95_list <- list()
for (i in 1:100){
  fit_endems_all_95 <- confint(glm_endems_pois_all_sp_list_models[[i]], level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")
  fit_endems_all_95_list[[i]] <- fit_endems_all_95
  }

# Combining the 95 intervals on all the models
results_endems_all_list <- list()
for (i in 1:100){
  results_endems_all_list[[i]] <- dplyr::bind_cols(glm_endems_pois_all_sp_list_tidy[[i]],fit_endems_all_95_list[[i]]) 
}




#Finding the best model for proportion endemics

# Finding the best model for the cladogenesis
# best_model_vars_prop_endems <- data.frame()
# 
# for (i in 1:100){
# test_dat_log <- output_all_sp_test_subset_log[[i]]
# 
# dredge_output <- qdredge_bin(glm(((colonization_sp+radiating_sp)/Total_sp)~mean_mx30_grd+area+dist+nearest_neighbour_distance_border_scaled+GeologicalOrigin-1,data=test_dat_log,family = quasibinomial, na.action = na.fail))
# 
# #print(dredge_output[which(dredge_output$delta == min(dredge_output[which(dredge_output$delta != 0)]$delta)), 1:10])
# 
# #best_model_vars_prop_endems <- rbind(best_model_vars_prop_endems, dredge_output[which(dredge_output$delta == min(dredge_output[which(dredge_output$delta != 0)]$delta)), 1:10])
# best_model_vars_prop_endems <- rbind(best_model_vars_prop_endems, dredge_output[which(dredge_output$delta == 0),1:10])
# }
# #best_model_vars_test <- is.na(best_model_vars)
# best_model_vars_count_prop_endems <- is.na(best_model_vars_prop_endems) %>% as.data.frame() %>% group_by_all() %>% count()
# best_model_vars_count_prop_endems



# Proportion endemics model

# -------- Endemic proportion ----------------
# glm_prop_endems_all_sp_list_tidy <- list()
# glm_prop_endems_all_sp_list_models <- list()
# 
# for (i in 1:100) {
#   #print((output_all_sp_test[[i]]$colonization_sp+output_all_sp_test[[i]]$radiating_sp)/output_all_sp_test[[i]]$Total_sp )
#   glm_prop_endems_all_sp_loop <- glm(((colonization_sp+radiating_sp)/Total_sp)~nearest_neighbour_distance_border_scaled,data=output_all_sp_test_subset_log[[i]],family = quasibinomial)
#   glm_prop_endems_all_sp_list_models[[i]] <- glm_prop_endems_all_sp_loop
#   glm_prop_endems_all_sp_list_tidy[[i]] <- broom::tidy(glm_prop_endems_all_sp_loop)
# }
# 
# #Now lets find the confidence intervals for these models.
# fit_prop_endems_all_95_list <- list()
# for (i in 1:100){
#   fit_prop_endems_all_95 <- confint(glm_prop_endems_all_sp_list_models[[i]], level = 0.95) %>% 
#   data.frame() %>% 
#   dplyr::rename("conf.low_95" = "X2.5..",
#          "conf.high_95" = "X97.5..")
#   fit_prop_endems_all_95_list[[i]] <- fit_prop_endems_all_95
#   }
# 
# # Combining the 95 intervals on all the models
# results_prop_endems_all_list <- list()
# for (i in 1:100){
#   results_prop_endems_all_list[[i]] <- dplyr::bind_cols(glm_prop_endems_all_sp_list_tidy[[i]],fit_prop_endems_all_95_list[[i]]) 
# }




# Finding the best model for anagenesis proportion


# # Finding the best model for the cladogenesis
# best_model_vars_ana_prop <- data.frame()
# 
# for (i in 1:100){
# test_dat_log <- output_all_sp_test_subset_log[[i]]
# 
# dredge_output <- qdredge_bin(glm(colonization_sp/(colonization_sp+radiating_sp)~mean_mx30_grd+area+dist+nearest_neighbour_distance_border_scaled+GeologicalOrigin-1,data=test_dat_log,family = quasipoisson, na.action = na.fail))
# 
# best_model_vars_ana_prop <- rbind(best_model_vars_ana_prop, dredge_output[which(dredge_output$delta == 0),1:9])
# 
# }
# 
# # dredge_output[which(dredge_output$delta == min(dredge_output[which(dredge_output$delta != 0)]$delta))
# #best_model_vars_test <- is.na(best_model_vars)
# best_model_vars_count_ana_prop <- is.na(best_model_vars_ana_prop) %>% as.data.frame() %>% group_by_all() %>% count()
# best_model_vars_count_ana_prop



# Proportion endemics from anagenesis

# -------- Proportion Regional allopatry ----------------
# glm_regional_allo_prop_all_sp_list_tidy <- list()
# glm_regional_allo_prop_all_sp_list_models <- list()
# 
# for (i in 1:100){
# #print((output_all_sp_test[[i]]$colonization_sp)/(output_all_sp_test[[i]]$colonization_sp+output_all_sp_test[[i]]$radiating_sp))
#   glm_regional_allo_prop_all_sp_loop <- glm(colonization_sp/(colonization_sp+radiating_sp)~dist-1, data=output_all_sp_test[[i]],family = quasibinomial)
#   glm_regional_allo_prop_all_sp_list_models[[i]] <- glm_regional_allo_prop_all_sp_loop
#   glm_regional_allo_prop_all_sp_list_tidy[[i]] <- broom::tidy(glm_regional_allo_prop_all_sp_loop)
# }
# 
# #Now lets find the confidence intervals for these models.
# fit_regional_allo_prop_all_95_list <- list()
# for (i in 1:100){
#   
#   fit_regional_allo_prop_all_95 <- confint(glm_regional_allo_prop_all_sp_list_models[[i]], level = 0.95) %>% 
#   data.frame() 
#     
#   names(fit_regional_allo_prop_all_95) <- "dist"
#   
#   fit_regional_allo_prop_all_95 <- fit_regional_allo_prop_all_95%>%
#     t()  %>%
#     data.frame()
#   
#   
#   fit_regional_allo_prop_all_95 <- dplyr::rename(fit_regional_allo_prop_all_95,"conf.low_95" = "X2.5..",
#          "conf.high_95" = "X97.5..")
#   fit_regional_allo_prop_all_95_list[[i]] <- fit_regional_allo_prop_all_95
# }
# 
# # Combining the 95 intervals on all the models
# results_regional_allo_prop_all_list <- list()
# for (i in 1:100){
#   results_regional_allo_prop_all_list[[i]] <- dplyr::bind_cols(glm_regional_allo_prop_all_sp_list_tidy[[i]],fit_regional_allo_prop_all_95_list[[i]]) 
# }

#pR2(glm_regional_allo_prop_all_sp_list_models[[1]])





# Starting the plotting

# Adding the speciation mode to the results
for (i in 1:100) {
  results_ana_all_list[[i]]$type <- "Anagenesis"
  results_endems_all_list[[i]]$type <- "Number Endemics"
  results_prop_endems_all_list[[i]]$type <- "Proportion Endemics"
  #results_regional_allo_prop_all_list[[i]]$type <- "Proportion Endemics from Regional Allopatry"
}

for (i in 1:96) {
results_clad_all_list[[i]]$type <- "Cladogenesis"
}
# I have some trees which cause problems because they have more parameter,s these are (8,26,31,41)
# these 4 models do not have any parameters for the geological origin because we have a single radiating species pair on the bahamas atolls
# because of this the scripts cannot find the upper confidence interval for any of the geological origins 
# I have therefore decided to exclude these 4 models in my script¨
#results_clad_all_list_final <- results_clad_all_list[-c(8,26,31,41)]

# So now I have 2 x 100 results tables
# I need to combine these into a massive data.frame that I can use to create a combined violin / bar plot for cladogenesis/Anagenesis 
#Creating dataframes with all estimates for each of the coefficients

#Anagenesis
glm_coef_area_ana <- as.data.frame(t(sapply(results_ana_all_list, function(x) unlist(x[2,]))))
glm_coef_dist_ana <- as.data.frame(t(sapply(results_ana_all_list, function(x) unlist(x[3,]))))
glm_coef_dist_squared_ana <- as.data.frame(t(sapply(results_ana_all_list, function(x) unlist(x[4,]))))
glm_coef_fragmentation_ana <- as.data.frame(t(sapply(results_ana_all_list, function(x) unlist(x[5,]))))
glm_coef_max_elev_ana <- as.data.frame(t(sapply(results_ana_all_list, function(x) unlist(x[6,]))))


#Cladogenesis
glm_coef_area_clad <- as.data.frame(t(sapply(results_clad_all_list, function(x) unlist(x[2,]))))
glm_coef_dist_clad <- as.data.frame(t(sapply(results_clad_all_list, function(x) unlist(x[3,]))))
glm_coef_fragmentation_clad <- as.data.frame(t(sapply(results_clad_all_list, function(x) unlist(x[4,]))))


# Endems
glm_coef_area_endems <- as.data.frame(t(sapply(results_endems_all_list, function(x) unlist(x[2,]))))
glm_coef_dist_endems <- as.data.frame(t(sapply(results_endems_all_list, function(x) unlist(x[3,]))))
glm_coef_dist_squared_endems <- as.data.frame(t(sapply(results_endems_all_list, function(x) unlist(x[4,]))))
glm_coef_fragmentation_endems <- as.data.frame(t(sapply(results_endems_all_list, function(x) unlist(x[5,]))))
glm_coef_max_elev_endems <- as.data.frame(t(sapply(results_endems_all_list, function(x) unlist(x[6,]))))


# # Endems proportion
# glm_coef_hight_prop_endems <- as.data.frame(t(sapply(results_prop_endems_all_list, function(x) unlist(x[2,]))))
# glm_coef_area_prop_endems <- as.data.frame(t(sapply(results_prop_endems_all_list, function(x) unlist(x[3,]))))
# glm_coef_dist_prop_endems <- as.data.frame(t(sapply(results_prop_endems_all_list, function(x) unlist(x[4,]))))
# glm_coef_fragmentation_prop_endems <- as.data.frame(t(sapply(results_prop_endems_all_list, function(x) unlist(x[5,]))))
# glm_coef_Geoatoll_prop_endems <- as.data.frame(t(sapply(results_prop_endems_all_list, function(x) unlist(x[6,]))))
# glm_coef_Geocontinental_prop_endems <- as.data.frame(t(sapply(results_prop_endems_all_list, function(x) unlist(x[7,]))))
# glm_coef_Geomixed_prop_endems <- as.data.frame(t(sapply(results_prop_endems_all_list, function(x) unlist(x[8,]))))
# 
# # Proportion of endems from regional allopatric speciation
# glm_coef_hight_regional_allo_prop <- as.data.frame(t(sapply(results_regional_allo_prop_all_list, function(x) unlist(x[2,]))))
# glm_coef_area_regional_allo_prop <- as.data.frame(t(sapply(results_regional_allo_prop_all_list, function(x) unlist(x[3,]))))
# glm_coef_dist_regional_allo_prop <- as.data.frame(t(sapply(results_regional_allo_prop_all_list, function(x) unlist(x[4,]))))
# glm_coef_fragmentation_regional_allo_prop <- as.data.frame(t(sapply(results_regional_allo_prop_all_list, function(x) unlist(x[5,]))))
# glm_coef_Geoatoll_regional_allo_prop <- as.data.frame(t(sapply(results_regional_allo_prop_all_list, function(x) unlist(x[6,]))))
# glm_coef_Geocontinental_regional_allo_prop <- as.data.frame(t(sapply(results_regional_allo_prop_all_list, function(x) unlist(x[7,]))))
# glm_coef_Geomixed_regional_allo_prop <- as.data.frame(t(sapply(results_regional_allo_prop_all_list, function(x) unlist(x[8,]))))

# Stacking the dataframes
glm_coefs <- rbind(glm_coef_area_ana,
                   glm_coef_dist_ana,
                   glm_coef_dist_squared_ana,
                   glm_coef_fragmentation_ana,
                   glm_coef_max_elev_ana,
                   glm_coef_area_clad,
                   glm_coef_dist_clad,
                   glm_coef_fragmentation_clad,
                   glm_coef_area_endems,
                   glm_coef_dist_endems,
                   glm_coef_dist_squared_endems,
                   glm_coef_fragmentation_endems,
                   glm_coef_max_elev_endems)

glm_coefs_list_clad <- list(glm_coef_area_clad,
                   glm_coef_dist_clad,
                   glm_coef_fragmentation_clad)

glm_coefs_list_ana <- list(glm_coef_area_ana,
                   glm_coef_dist_ana,
                   glm_coef_dist_squared_ana,
                   glm_coef_fragmentation_ana,
                   glm_coef_max_elev_ana)

glm_coefs_list_endems <- list(glm_coef_area_endems,
                   glm_coef_dist_endems,
                   glm_coef_dist_squared_endems,
                   glm_coef_fragmentation_endems,
                   glm_coef_max_elev_endems)

# glm_coefs_list_prop_endems <- list(glm_coef_area_prop_endems,
#                    glm_coef_hight_prop_endems,
#                    glm_coef_dist_prop_endems,
#                    glm_coef_fragmentation_prop_endems,
#                    glm_coef_Geocontinental_prop_endems,
#                    glm_coef_Geomixed_prop_endems,
#                    glm_coef_Geoatoll_prop_endems)
# 
# glm_coefs_list_regional_allo_prop <- list(glm_coef_area_regional_allo_prop,
#                    glm_coef_hight_regional_allo_prop,
#                    glm_coef_dist_regional_allo_prop,
#                    glm_coef_fragmentation_regional_allo_prop,
#                    glm_coef_Geocontinental_regional_allo_prop,
#                    glm_coef_Geomixed_regional_allo_prop
#                    ) # glm_coef_Geoatoll_regional_allo_prop


# Renaming some variables
#Writing up the names
#glm_coefs_names <- c("Area","Max elevation","Isolation","Fragmentation","Continental","Mixed Origin", "Atoll")
glm_coefs_names_old <- unique(glm_coefs$term)
glm_coefs_names <- c("Area","Isolation", "Isolation Squared","Fragmentation","Max Elevevation","Isolation")
glm_coefs_names_df <- cbind(glm_coefs_names, glm_coefs_names_old)

# For loop for renaming
for (i in 1:6){
  glm_coefs[which(glm_coefs$term==glm_coefs_names_df[i,2]),1] <- glm_coefs_names_df[i,1]
}

# Now I only need to create the dataframe with the means and then I think I can plot it.
#Anagenesis means
glm_coefs_means_ana <- data.frame()
for (i in 1:length(glm_coefs_list_ana)) {
  temp_means <- vector()
  temp_means <- apply(glm_coefs_list_ana[[i]], 2, function(x) mean(as.numeric(unlist(x))))
  temp_means["term"] <- glm_coefs_list_ana[[i]][1,1]
  temp_means["type"] <- glm_coefs_list_ana[[i]][1,8]
  #print(temp_means)
  glm_coefs_means_ana <- rbind(glm_coefs_means_ana,temp_means)
  print(glm_coefs_means_ana)
}
names(glm_coefs_means_ana) <- names(glm_coef_dist_ana)
# For loop for renaming
for (i in 1:6){
  glm_coefs_means_ana[which(glm_coefs_means_ana$term==glm_coefs_names_df[i,2]),1] <- glm_coefs_names_df[i,1]
}


# making a loop to convert character list rows to numeric 
for(i in 2:7){
  glm_coefs_means_ana[,i] <-as.numeric(unlist(glm_coefs_means_ana[,i]))
}


# Cladogenesis means
glm_coefs_means_clad <- data.frame()
for (i in 1:length(glm_coefs_list_clad)) {
  print(i)
  temp_means <- vector()
  temp_means <- apply(glm_coefs_list_clad[[i]], 2, function(x) mean(as.numeric(unlist(x))))
  temp_means["term"] <- glm_coefs_list_clad[[i]][1,1]
  temp_means["type"] <- glm_coefs_list_clad[[i]][1,8]
  glm_coefs_means_clad <- rbind(glm_coefs_means_clad,temp_means)
}

names(glm_coefs_means_clad) <- names(glm_coef_dist_clad)
# For loop for renaming
for (i in 1:6){
  glm_coefs_means_clad[which(glm_coefs_means_clad$term==glm_coefs_names_df[i,2]),1] <- glm_coefs_names_df[i,1]
}

# Endemics means
glm_coefs_means_endems <- data.frame()
for (i in 1:length(glm_coefs_list_endems)) {
  temp_means <- vector()
  temp_means <- apply(glm_coefs_list_endems[[i]], 2, function(x) mean(as.numeric(unlist(x))))
  temp_means["term"] <- glm_coefs_list_endems[[i]][1,1]
  temp_means["type"] <- glm_coefs_list_endems[[i]][1,8]
  glm_coefs_means_endems <- rbind(glm_coefs_means_endems,temp_means)
}

names(glm_coefs_means_endems) <- names(glm_coef_dist_endems)
# For loop for renaming
for (i in 1:6){
  glm_coefs_means_endems[which(glm_coefs_means_endems$term==glm_coefs_names_df[i,2]),1] <- glm_coefs_names_df[i,1]
}

# proportion endemic means
# glm_coefs_means_prop_endems <- data.frame()
# for (i in 1:7) {
#   temp_means <- vector()
#   temp_means <- apply(glm_coefs_list_prop_endems[[i]], 2, function(x) mean(as.numeric(unlist(x))))
#   temp_means["term"] <- glm_coefs_list_prop_endems[[i]][1,1]
#   temp_means["type"] <- glm_coefs_list_prop_endems[[i]][1,8]
#   glm_coefs_means_prop_endems <- rbind(glm_coefs_means_prop_endems,temp_means)
# }
# names(glm_coefs_means_prop_endems) <- names(glm_coef_dist_prop_endems)

# Regional allopatry proportion means
# glm_coefs_means_regional_allo_prop <- data.frame()
# for (i in 1:6) {
#   temp_means <- vector()
#   temp_means <- apply(glm_coefs_list_regional_allo_prop[[i]], 2, function(x) mean(as.numeric(unlist(x))))
#   temp_means["term"] <- glm_coefs_list_regional_allo_prop[[i]][1,1]
#   temp_means["type"] <- glm_coefs_list_regional_allo_prop[[i]][1,8]
#   glm_coefs_means_regional_allo_prop <- rbind(glm_coefs_means_regional_allo_prop,temp_means)
# }
# names(glm_coefs_means_regional_allo_prop) <- names(glm_coef_dist_regional_allo_prop)

# Stacking the means dataframes for the means
coefs_means <- NULL
coefs_means <- rbind(glm_coefs_means_clad,
                     glm_coefs_means_ana,
                     glm_coefs_means_endems)

coefs_means <- AddVars_mean(coefs_means, terms_list, type_list) # this function is located in Data_setup_and_functions
glm_coefs <- AddVars_violin(glm_coefs, terms_list, type_list)  # this function is located in Data_setup_and_functions

# making a loop to convert character list rows to numeric 
for(i in 2:7){
  coefs_means[,i] <-as.numeric(unlist(coefs_means[,i]))
}

# Changing names in dataframe of estimates
for (i in 1:6){
  coefs_means[which(coefs_means$term==glm_coefs_names_df[i,2]),1] <- glm_coefs_names_df[i,1]
}

# Creating plotting order
plot_order <- c("Area", "Isolation","Isolation Squared","Max Elevation","Fragmentation","Continental","Mixed Origin")

#Now I should be ready to make my violin + box plot for all species on all the islands.
gvil_all_all <- ggplot(glm_coefs, aes(color=type, fill = type)) +
  scale_colour_manual(values = c("#efc86e","#6f9969","#454a74")) +
  scale_fill_manual(values = c("#efc86e","#6f9969","#454a74")) +
  scale_x_discrete(limits = plot_order) +
  geom_violin(aes(x = term, y = as.numeric(estimate)),width = 0.60,position=position_dodge(width=1), scale = "width", alpha = 0.8,lwd=0.2) +
  geom_violin(aes(x = term, y = as.numeric(conf.high_95)),width = 0.60,position=position_dodge(width=1), scale = "width", alpha = 0.8,lwd=0.2) +
  geom_violin(aes(x = term, y = as.numeric(conf.low_95)),width = 0.60,position=position_dodge(width=1), scale = "width", alpha = 0.8,lwd=0.2) +
  geom_crossbar(data = coefs_means, size=0.1, alpha=0.7, aes(x= term, y = estimate,
                    ymin=conf.low_95,
                    ymax=conf.high_95),
                position = position_dodge(width = 1)) +
  geom_text(data = coefs_means,
            aes(x = term, y = conf.high_95 + 0.3, label = ifelse(p.value < 0.05, "*", " ")),
            position = position_dodge(width = 1),
            hjust = -0.7, size = 5, show.legend = FALSE, na.rm = TRUE)+
  geom_hline(yintercept = 0, colour = gray(0), lty = 3) +
  xlab("Coefficients") +
  ylab("Estimate") +
  scale_y_continuous(breaks = seq(-10, 15, by = 5), limits = c(-10,15)) +
  theme_classic() +
  theme(axis.title.y = element_blank(),legend.position = "none") +
  ggtitle("All Plants: All Geological Origins")

gvil_all_all


# Finding the best model for the cladogenesis
dredge_output_clad_cory <- qdredge(glm(`No. cladogenesis`-`No. edges`~max30_elev+area+dist+nearest_neighbour_distance_border_scaled+GeologicalOrigin,data=test_dat_log_cory,family = quasipoisson, na.action = na.fail))

dredge_output_clad_cory


# Looking at the results from the cladogenesis model
glm_cory_clad <- glm(`No. cladogenesis`-`No. edges`~area,data=test_dat_log_cory,family = quasipoisson, na.action = na.fail)

summary(glm_cory_clad)



# Finding the best model for the Anagenesis
dredge_output_ana_cory <- qdredge(glm(`No. anagenesis`+`No. edges`~max30_elev+area+poly(dist,2,raw = TRUE)+nearest_neighbour_distance_border_scaled+GeologicalOrigin,data=test_dat_log_cory,family = quasipoisson, na.action = na.fail))

dredge_output_ana_cory

glm_cory_ana <- glm(`No. anagenesis`+`No. edges`~area+nearest_neighbour_distance_border_scaled, data=test_dat_log_cory,family = quasipoisson, na.action = na.fail)

summary(glm_cory_ana)


# Finding the best model for the cladogenesis
# dredge_output_radiating_nodes_cory <- qdredge(glm(`No. edges`~mean_mx30_grd+area+dist^2+nearest_neighbour_distance_border_scaled+GeologicalOrigin,data=test_dat_log_cory,family = quasipoisson, na.action = na.fail))
# 
# dredge_output_radiating_nodes_cory





# Finding the best model for the Endems
dredge_output_endems_cory <- qdredge(glm(`No. anagenesis`+`No. cladogenesis`~max30_elev+area+poly(dist,2)+nearest_neighbour_distance_border_scaled+GeologicalOrigin,data=test_dat_log_cory,family = quasipoisson, na.action = na.fail))

dredge_output_endems_cory



glm_cory_endems <- glm(`No. anagenesis`+`No. cladogenesis`~area,data=test_dat_log_cory,family = quasipoisson, na.action = na.fail)

summary(glm_cory_endems)


car::vif(glm_cory_ana)
#car::vif(glm_cory_clad) Fewer than two terms so we cant do VIF
#car::vif(glm_cory_endems) Fewer than two terms so we cant do VIF

#Getting the coefficient estimates for cladogenesis and anagenesis

#Anagenesis
results_ana_cory <- broom::tidy(glm_cory_ana)

#Now lets find the confidence intervals for these models.
fit_ana_cory <- confint(glm_cory_ana, level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")

results_ana_cory <- dplyr::bind_cols(results_ana_cory,fit_ana_cory) 
results_ana_cory <- results_ana_cory[-1,]

#cladogenesis
results_clad_cory <- broom::tidy(glm_cory_clad)

#Now lets find the confidence intervals for these models.
# "mean_mx30_grd"
fit_clad_cory <- confint(glm_cory_clad, level = 0.95) %>%  
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")

results_clad_cory <- dplyr::bind_cols(results_clad_cory,fit_clad_cory)
results_clad_cory <- results_clad_cory[-1,]

# No endems
results_endems_cory <- broom::tidy(glm_cory_endems)

#Now lets find the confidence intervals for these models.
fit_endems_cory <- confint(glm_cory_endems, level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")

results_endems_cory <- dplyr::bind_cols(results_endems_cory,fit_endems_cory)
results_endems_cory <- results_endems_cory[-1,]
# #Endemics proportion
# results_prop_endems_cory <- broom::tidy(glm_prop_endems_cory)
# 
# #Now lets find the confidence intervals for these models.
# fit_prop_endems_cory <- confint(glm_prop_endems_cory, level = 0.95) %>% 
#   data.frame() %>% 
#   dplyr::rename("conf.low_95" = "X2.5..",
#          "conf.high_95" = "X97.5..")
# 
# results_prop_endems_cory <- dplyr::bind_cols(results_prop_endems_cory,fit_prop_endems_cory)
# 
# 
# #Proportion endemics from regional allopatry
# results_regional_allo_prop_cory <- broom::tidy(glm_regional_allo_prop_cory)
# 
# #Now lets find the confidence intervals for these models.
# fit_regional_allo_prop_cory <- confint(glm_regional_allo_prop_cory, level = 0.95) %>% 
#   data.frame() %>% 
#   dplyr::rename("conf.low_95" = "X2.5..",
#          "conf.high_95" = "X97.5..")
# 
# results_regional_allo_prop_cory <- dplyr::bind_cols(results_regional_allo_prop_cory,fit_regional_allo_prop_cory)





results_ana_cory$type <-"Anagenesis"
results_clad_cory$type <- "Cladogenesis"
results_endems_cory$type <- "Number Endemics"
#results_prop_endems_cory$type <- "Proportion Endemics"
#results_regional_allo_prop_cory$type <- "Proportion Endemics from Regional Allopatry"

results_cory <- dplyr::bind_rows(results_ana_cory,
                                 results_clad_cory,
                                 results_endems_cory,
                                 #results_prop_endems_cory,
                                 #results_regional_allo_prop_cory
                                 )


#glm_coefs_names <- c("Area","Max elevation","Isolation","Fragmentation","Continental","Mixed Origin", "Atoll")
glm_coefs_names_cory <- c("Area","Fragmentation")
glm_coefs_names_old_cory <- unique(results_cory$term)
glm_coefs_names_df_cory <- cbind(glm_coefs_names_cory, glm_coefs_names_old_cory)

# For loop for renaming
for (i in 1:2){
  results_cory[which(results_cory$term==glm_coefs_names_df_cory[i,2]),1] <- glm_coefs_names_df_cory[i,1]
}


results_cory <- AddVars_mean(results_cory, terms_list, type_list)

for(i in 2:7){
  results_cory[,i] <-as.numeric(unlist(results_cory[,i]))
}

#Now I should be ready to make my violin + box plot for all species on all the islands.
gvil_coryphoideae <- ggplot(results_cory, aes(color=type, fill = type)) +
  scale_colour_manual(values = c("#efc86e","#6f9969","#454a74")) +
  scale_fill_manual(values = c("#efc86e","#6f9969","#454a74")) +
  scale_x_discrete(limits = plot_order) +
  geom_crossbar(data = results_cory, size=0.1, alpha=0.7, aes(x= term, y = estimate,
                    ymin=conf.low_95,
                    ymax=conf.high_95),
                position = position_dodge2(width = 1.5)) +
  geom_text(data = results_cory,
            aes(x = term, y = conf.high_95 + 0.3, label = ifelse(p.value < 0.05, "*", " ")),
            position = position_dodge(width = 1),
            hjust = -0.7, size = 5, show.legend = FALSE, na.rm = TRUE)+
  geom_hline(yintercept = 0, colour = gray(0), lty = 3) +
  xlab("Coefficients") +
  ylab("Estimate") +
  scale_y_continuous(breaks = seq(-10, 15, by = 5), limits = c(-10,15)) +
  theme_classic() +
  labs(fill = element_blank(), colour = element_blank()) +
  theme(axis.title.y = element_blank(), legend.position = "none") +
  ggtitle("Coryphoideae subfamily")

gvil_coryphoideae






prow_2 <- cowplot::plot_grid(gvil_all_all, gvil_coryphoideae+theme(legend.position = "none"), labels = c("a.", "b."), align = "v", rel_widths = c(1,1))

legend_bottom <- get_legend(gvil_coryphoideae + guides(color = guide_legend(nrow=1))+theme(legend.position = "bottom", legend.text=element_text(size=6)))

prow_2 <- cowplot::plot_grid(prow_2,legend_bottom, ncol = 1, rel_heights = c(1,.1))

prow_2


test_dat_log_volc <- list()
for (i in 1:100) {
  test_dat_log_volc[[i]] <- output_all_sp_test_subset_log[[i]][which(output_all_sp_test_subset_log[[i]]$GeologicalOrigin == "volcanic"),]
}




# Finding the best model for the cladogenesis
best_model_vars_clad_volc <- data.frame()
second_best_model_clad_volc  <- data.frame()
third_best_model_clad_volc  <- data.frame()
fourth_best_model_clad_volc <- data.frame()

for (i in 1:100){
test_data_volcanic <- test_dat_log_volc[[i]]

dredge_output <- qdredge(glm(radiating_sp-radiating_nodes~max30_elev+area+poly(dist,1)+nearest_neighbour_distance_border_scaled,data=test_data_volcanic,family = quasipoisson,na.action = na.fail))

row_to_add <- dredge_output[which(dredge_output$delta == 0),1:9]

best_model_vars_clad_volc <- rbind(best_model_vars_clad_volc, row_to_add)

second_best_model_clad_volc  <- rbind(second_best_model_clad_volc,dredge_output[2,1:9])
third_best_model_clad_volc  <- rbind(third_best_model_clad_volc, dredge_output[3,1:9])
fourth_best_model_clad_volc <- rbind(fourth_best_model_clad_volc,dredge_output[4,1:9])
 
}

#best_model_vars_test <- is.na(best_model_vars)
best_model_vars_count_clad_volc <- is.na(best_model_vars_clad_volc[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()
best_model_vars_count_clad_volc$mean_delta <- 0

second_best_model_clad_volc_count  <- is.na(second_best_model_clad_volc[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()
third_best_model_clad_volc_count  <- is.na(third_best_model_clad_volc[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()
fourth_best_model_clad_volc_count <-is.na(fourth_best_model_clad_volc[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()

second_best_model_clad_volc_count$mean_delta  <- mean(second_best_model_clad_volc[,9])
third_best_model_clad_volc_count$mean_delta <- mean(third_best_model_clad_volc[,9])
fourth_best_model_clad_volc_count$mean_delta <- mean(fourth_best_model_clad_volc[,9])

clad_volc_table <- rbind(best_model_vars_count_clad_volc,second_best_model_clad_volc_count,third_best_model_clad_volc_count,fourth_best_model_clad_volc_count)

clad_volc_table <- as.data.frame(clad_volc_table)

clad_volc_table[,1:6][clad_volc_table[,1:6]==TRUE] <- "Drop"
clad_volc_table[,1:6][clad_volc_table[,1:6]==FALSE] <- "Keep"

clad_volc_table









glm_clad_pois_all_sp_list_models_volc <- list()
glm_clad_pois_all_sp_list_tidy_volc <- list()
glm_clad_pois_all_sp_list_dispersion_volc <- list()

for (i in 1:100){
  glm_clad_all_sp_volc_loop <- glm(radiating_sp-radiating_nodes~area+poly(dist,1)+max30_elev,data=test_dat_log_volc[[i]],family = quasipoisson)
  glm_clad_pois_all_sp_list_models_volc[[i]] <- glm_clad_all_sp_volc_loop
  glm_clad_pois_all_sp_list_tidy_volc[[i]] <- broom::tidy(glm_clad_all_sp_volc_loop)
  glm_clad_pois_all_sp_list_dispersion_volc[[i]] <- sum(residuals(glm_clad_all_sp_volc_loop,type ="pearson")^2)/glm_clad_all_sp_volc_loop$df.residual
}

#Now lets find the confidence intervals for these models.
# fit_clad_all_95_list_volc <- list()
# for (i in 1:100){
#   fit_clad_all_95_volc <- confint(glm_clad_pois_all_sp_list_volc_models[[i]], level = 0.95) %>% 
#   data.frame()
#   
#   names(fit_clad_all_95_volc) <- "Area" 
#   
#   fit_clad_all_95_volc <- data.frame(fit_clad_all_95_volc) %>%
#     t()  %>%
#     data.frame()
#   
#     fit_clad_all_95_volc <- dplyr::rename(fit_clad_all_95_volc,"conf.low_95" = "X2.5..",
#          "conf.high_95" = "X97.5..")
#   
#   fit_clad_all_95_list_volc[[i]] <- fit_clad_all_95_volc
# }

#Now lets find the confidence intervals for these models.
fit_clad_all_95_list_volc <- list()
for (i in 1:100){
  fit_clad_all_95_volc <- confint(glm_clad_pois_all_sp_list_models_volc[[i]], level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")
  fit_clad_all_95_list_volc[[i]] <- fit_clad_all_95_volc
  }


# Combining the 95 intervals on all the models
results_clad_all_list_volc <- list()
for (i in 1:100){
  results_clad_all_list_volc[[i]] <- dplyr::bind_cols(glm_clad_pois_all_sp_list_tidy_volc[[i]],fit_clad_all_95_list_volc[[i]]) 
}



# Finding the best model for the Anagenesis on volcanic islands
best_model_vars_ana_volc <- data.frame()
second_best_model_ana_volc  <- data.frame()
third_best_model_ana_volc  <- data.frame()
fourth_best_model_ana_volc <- data.frame()

for (i in 1:100){
test_data_volcanic <- test_dat_log_volc[[i]]

dredge_output <- qdredge(glm(colonization_sp+radiating_nodes~max30_elev+area+poly(dist,2)+nearest_neighbour_distance_border_scaled,data=test_data_volcanic,family = quasipoisson,na.action = na.fail))

row_to_add <- dredge_output[which(dredge_output$delta == 0),1:9]

best_model_vars_ana_volc <- rbind(best_model_vars_ana_volc, row_to_add)

second_best_model_ana_volc  <- rbind(second_best_model_ana_volc,dredge_output[2,1:9])
third_best_model_ana_volc  <- rbind(third_best_model_ana_volc, dredge_output[3,1:9])
fourth_best_model_ana_volc <- rbind(fourth_best_model_ana_volc,dredge_output[4,1:9])
 
}

#best_model_vars_test <- is.na(best_model_vars)
best_model_vars_count_ana_volc <- is.na(best_model_vars_ana_volc[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()
best_model_vars_count_ana_volc$mean_delta <- 0

second_best_model_ana_volc_count  <- is.na(second_best_model_ana_volc[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()
third_best_model_ana_volc_count  <- is.na(third_best_model_ana_volc[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()
fourth_best_model_ana_volc_count <-is.na(fourth_best_model_ana_volc[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()

second_best_model_ana_volc_count$mean_delta  <- mean(second_best_model_ana_volc[,9])
third_best_model_ana_volc_count$mean_delta <- mean(third_best_model_ana_volc[,9])
fourth_best_model_ana_volc_count$mean_delta <- mean(fourth_best_model_ana_volc[,9])

ana_volc_table <- rbind(best_model_vars_count_ana_volc,second_best_model_ana_volc_count,third_best_model_ana_volc_count,fourth_best_model_ana_volc_count)

ana_volc_table <- as.data.frame(ana_volc_table)

ana_volc_table[,1:6][ana_volc_table[,1:6]==TRUE] <- "Drop"
ana_volc_table[,1:6][ana_volc_table[,1:6]==FALSE] <- "Keep"

ana_volc_table




# Creating a for loop which fits the model using all 100 different datasets 
# -------- Anagenesis ----------------
glm_ana_pois_all_sp_list_tidy_volc <- list()
glm_ana_pois_all_sp_list_models_volc <- list()

for (i in 1:100){
  glm_ana_all_sp_loop_volc <- glm(colonization_sp+radiating_nodes~area+poly(dist,2, raw=TRUE)+nearest_neighbour_distance_border_scaled,data=test_dat_log_volc[[i]],family = quasipoisson)
  glm_ana_pois_all_sp_list_models_volc[[i]] <- glm_ana_all_sp_loop_volc
  glm_ana_pois_all_sp_list_tidy_volc[[i]] <- broom::tidy(glm_ana_all_sp_loop_volc)
}

#Now lets find the confidence intervals for these models.
fit_ana_all_95_list_volc <- list()
for (i in 1:100){
  fit_ana_all_95_volc <- confint(glm_ana_pois_all_sp_list_models_volc[[i]], level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")
  fit_ana_all_95_list_volc[[i]] <- fit_ana_all_95_volc
  }

# Combining the 95 intervals on all the models
results_ana_all_list_volc <- list()
for (i in 1:100){
  results_ana_all_list_volc[[i]] <- dplyr::bind_cols(glm_ana_pois_all_sp_list_tidy_volc[[i]],fit_ana_all_95_list_volc[[i]]) 
}






# Finding the best model for the endemics on volcanic islands
best_model_vars_endems_volc <- data.frame()
second_best_model_endems_volc  <- data.frame()
third_best_model_endems_volc  <- data.frame()
fourth_best_model_endems_volc <- data.frame()

for (i in 1:100){
test_data_volcanic <- test_dat_log_volc[[i]]

dredge_output <- qdredge(glm(colonization_sp+radiating_sp~max30_elev+area+poly(dist,2)+nearest_neighbour_distance_border_scaled,data=test_data_volcanic,family = quasipoisson,na.action = na.fail))

row_to_add <- dredge_output[which(dredge_output$delta == 0),1:9]

best_model_vars_endems_volc <- rbind(best_model_vars_endems_volc, row_to_add)

second_best_model_endems_volc  <- rbind(second_best_model_endems_volc,dredge_output[2,1:9])
third_best_model_endems_volc  <- rbind(third_best_model_endems_volc, dredge_output[3,1:9])
fourth_best_model_endems_volc <- rbind(fourth_best_model_endems_volc,dredge_output[4,1:9])
 
}

#best_model_vars_test <- is.na(best_model_vars)
best_model_vars_count_endems_volc <- is.na(best_model_vars_endems_volc[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()
best_model_vars_count_endems_volc$mean_delta <- 0

second_best_model_endems_volc_count  <- is.na(second_best_model_endems_volc[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()
third_best_model_endems_volc_count  <- is.na(third_best_model_endems_volc[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()
fourth_best_model_endems_volc_count <-is.na(fourth_best_model_endems_volc[1:6]) %>% as.data.frame() %>% group_by_all() %>% count()

second_best_model_endems_volc_count$mean_delta  <- mean(second_best_model_endems_volc[,9])
third_best_model_endems_volc_count$mean_delta <- mean(third_best_model_endems_volc[,9])
fourth_best_model_endems_volc_count$mean_delta <- mean(fourth_best_model_endems_volc[,9])

endems_volc_table <- rbind(best_model_vars_count_endems_volc,second_best_model_endems_volc_count,third_best_model_endems_volc_count,fourth_best_model_endems_volc_count)

endems_volc_table <- as.data.frame(endems_volc_table)

endems_volc_table[,1:6][endems_volc_table[,1:6]==TRUE] <- "Drop"
endems_volc_table[,1:6][endems_volc_table[,1:6]==FALSE] <- "Keep"

endems_volc_table



# Creating a for loop which fits the model using all 100 different datasets 
# -------- Endems ----------------
glm_endems_pois_all_sp_list_tidy_volc <- list()
glm_endems_pois_all_sp_list_models_volc <- list()

for (i in 1:100){
  glm_endems_all_sp_loop_volc <- glm(colonization_sp+radiating_sp~area+poly(dist,2, raw=TRUE)+nearest_neighbour_distance_border_scaled,data=test_dat_log_volc[[i]],family = quasipoisson)
  glm_endems_pois_all_sp_list_models_volc[[i]] <- glm_endems_all_sp_loop_volc
  glm_endems_pois_all_sp_list_tidy_volc[[i]] <- broom::tidy(glm_endems_all_sp_loop_volc)
}

#Now lets find the confidence intervals for these models.
fit_endems_all_95_list_volc <- list()
for (i in 1:100){
  fit_endems_all_95_volc <- confint(glm_endems_pois_all_sp_list_models_volc[[i]], level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")
  fit_endems_all_95_list_volc[[i]] <- fit_endems_all_95_volc
  }

# Combining the 95 intervals on all the models
results_endems_all_list_volc <- list()
for (i in 1:100){
  results_endems_all_list_volc[[i]] <- dplyr::bind_cols(glm_endems_pois_all_sp_list_tidy_volc[[i]],fit_endems_all_95_list_volc[[i]]) 
}





# # Checking variance inflation factors
# vif_res <- data.frame(matrix(ncol=7,nrow=0))
# colnames(vif_res) <- c("area","poly(dist, 1, raw = TRUE)","poly(dist, 2, raw = TRUE)","max30_elev","GeologicalOrigin","model","nearest_neighbour_distance_border_scaled")
# 
# for (i in 1:100){
#   vif_ana <- as.data.frame(t(car::vif(glm_ana_pois_all_sp_list_models[[i]])))
#   vif_ana$GeologicalOrigin <- NA
#   vif_ana$"poly(dist, 1, raw = TRUE)" <- NA
#   vif_ana$model <- "Anagenesis all islands"
#   
#   vif_endems <- as.data.frame(t(car::vif(glm_endems_pois_all_sp_list_models[[i]])))
#   vif_endems$mean_mx30_grd <- NA
#   vif_endems$GeologicalOrigin <- NA
#   vif_endems$"poly(dist, 1, raw = TRUE)" <- NA
#   #vif_endems$"I(dist * dist)" <- NA
#   vif_endems$model <- "Endemics all islands"
#   
#   vif_clad_volc <- as.data.frame(t(car::vif(glm_clad_pois_all_sp_list_models_volc[[i]])))
#   vif_clad_volc$GeologicalOrigin <- NA
#   vif_clad_volc$"poly(dist, 1, raw = TRUE)" <- NA
#   vif_clad_volc$"poly(dist, 2, raw = TRUE)" <- NA
#   vif_clad_volc$model <- "Cladogenesis volcanic islands"
# 
#   
#   vif_ana_volc <- as.data.frame(t(car::vif(glm_ana_pois_all_sp_list_models_volc[[i]])))
#   vif_ana_volc$mean_mx30_grd <- NA
#   vif_ana_volc$GeologicalOrigin <- NA
#   vif_ana_volc$"poly(dist, 1, raw = TRUE)" <- NA
#   #vif_ana_volc$dist <- NA
#   vif_ana_volc$model <- "Anagenesis volcanic islands"
#   
#   vif_endems_volc <- as.data.frame(t(car::vif(glm_endems_pois_all_sp_list_models_volc[[i]])))
#   vif_endems_volc$GeologicalOrigin <- NA
#   vif_endems_volc$"poly(dist, 1, raw = TRUE)" <- NA
#   #vif_endems_volc$"I(dist * dist)" <- NA
#   vif_endems_volc$model <- "Endems volcanic islands"
#   
#   vif_res <- rbind(vif_res,
#                    vif_ana[1,],
#                    vif_endems[1,],
#                    vif_ana_volc[1,],
#                    vif_endems_volc[1,],
#                    vif_clad_volc[1,])
# }
# 
# 
# 
# vif_clad_total <- data_frame()
# for (i in 1:96){
#   vif_clad <- as.data.frame(t(car::vif(glm_clad_pois_all_sp_list_models[[i]])[,1]))
#   vif_clad$mean_mx30_grd <- NA
#   vif_clad$"poly(dist, 2, raw = TRUE)" <- NA
#   names(vif_clad)[names(vif_clad) == 'relevel(factor(GeologicalOrigin), ref = "volcanic")'] <- "GeologicalOrigin"
#   vif_clad$model <- "Cladogenesis all islands"
#   
#   vif_clad_total <- rbind(vif_clad_total,vif_clad[1,])
# }
# 
# vif_res <- rbind(vif_res,vif_clad_total)
# 
# vif_res2 <- vif_res[order(vif_res$model),]
# vif_res2



# Adding the speciation mode to the results
for (i in 1:100) {
  results_ana_all_list_volc[[i]]$type <- "Anagenesis"
  results_endems_all_list_volc[[i]]$type <- "Number Endemics"
  results_clad_all_list_volc[[i]]$type <- "Cladogenesis"
}


# I have some trees which cause problems because they have more parameter,s these are (8,26,31,41)
# these 4 models do not have any parameters for the geological origin because we have a single radiating species pair on the bahamas atolls
# because of this the scripts cannot find the upper confidence interval for any of the geological origins 
# I have therefore decided to exclude these 4 models in my script¨
#results_clad_all_list_final <- results_clad_all_list[-c(8,26,31,41)]

# So now I have 2 x 100 results tables
# I need to combine these into a massive data.frame that I can use to create a combined violin / bar plot for cladogenesis/Anagenesis 
#Creating dataframes with all estimates for each of the coefficients

#Anagenesis
glm_coef_area_ana_volc <- as.data.frame(t(sapply(results_ana_all_list_volc, function(x) unlist(x[2,]))))
glm_coef_dist_ana_volc <- as.data.frame(t(sapply(results_ana_all_list_volc, function(x) unlist(x[3,]))))
glm_coef_dist_squared_ana_volc <- as.data.frame(t(sapply(results_ana_all_list_volc, function(x) unlist(x[4,]))))
glm_coef_fragmentation_ana_volc <- as.data.frame(t(sapply(results_ana_all_list_volc, function(x) unlist(x[5,]))))

#Cladogenesis
glm_coef_max_elev_clad_volc <- as.data.frame(t(sapply(results_clad_all_list_volc, function(x) unlist(x[2,]))))
glm_coef_area_clad_volc <- as.data.frame(t(sapply(results_clad_all_list_volc, function(x) unlist(x[3,]))))
glm_coef_dist_clad_volc <- as.data.frame(t(sapply(results_clad_all_list_volc, function(x) unlist(x[4,]))))

# Endems
glm_coef_area_endems_volc <- as.data.frame(t(sapply(results_endems_all_list_volc, function(x) unlist(x[2,]))))
glm_coef_dist_endems_volc <- as.data.frame(t(sapply(results_endems_all_list_volc, function(x) unlist(x[3,]))))
glm_coef_dist_squared_endems_volc <- as.data.frame(t(sapply(results_endems_all_list_volc, function(x) unlist(x[4,]))))
glm_coef_fragmentation_endems_volc <- as.data.frame(t(sapply(results_endems_all_list_volc, function(x) unlist(x[5,]))))

# Stacking the dataframes
glm_coefs_volc <- rbind(glm_coef_area_ana_volc,
                   glm_coef_dist_ana_volc,
                   glm_coef_dist_squared_ana_volc,
                   glm_coef_fragmentation_ana_volc,
                   glm_coef_max_elev_clad_volc,
                   glm_coef_area_clad_volc,
                   glm_coef_dist_clad_volc,
                   glm_coef_area_endems_volc,
                   glm_coef_dist_endems_volc,
                   glm_coef_dist_squared_endems_volc,
                   glm_coef_fragmentation_endems_volc)

glm_coefs_list_clad_volc <- list(glm_coef_area_clad_volc,
                                 glm_coef_max_elev_clad_volc,
                                 glm_coef_dist_clad_volc)

glm_coefs_list_ana_volc <- list(glm_coef_area_ana_volc,
                   glm_coef_dist_ana_volc,
                   glm_coef_dist_squared_ana_volc,
                   glm_coef_fragmentation_ana_volc)

glm_coefs_list_endems_volc <- list(glm_coef_area_endems_volc,
                   glm_coef_dist_endems_volc,
                   glm_coef_dist_squared_endems_volc,
                   glm_coef_fragmentation_endems_volc)

# Renaming some variables
#Writing up the names
#glm_coefs_names <- c("Area","Max elevation","Isolation","Fragmentation","Continental","Mixed Origin", "Atoll")
glm_coefs_names_old_volc <- unique(glm_coefs_volc$term)
glm_coefs_names_volc <- c("Area","Isolation", "Isolation Squared", "Fragmentation","Max Elevation", "Isolation")
glm_coefs_names_df_volc <- cbind(glm_coefs_names_volc, glm_coefs_names_old_volc)

# For loop for renaming
for (i in 1:length(glm_coefs_names_volc)){
  glm_coefs_volc[which(glm_coefs_volc$term==glm_coefs_names_df_volc[i,2]),1] <- glm_coefs_names_df_volc[i,1]
}

# Now I only need to create the dataframe with the means and then I think I can plot it.
#Anagenesis means
glm_coefs_means_ana_volc <- data.frame()
for (i in 1:3) {
  temp_means <- vector()
  temp_means <- apply(glm_coefs_list_ana_volc[[i]], 2, function(x) mean(as.numeric(unlist(x))))
  temp_means["term"] <- glm_coefs_list_ana_volc[[i]][1,1]
  temp_means["type"] <- glm_coefs_list_ana_volc[[i]][1,8]
  #print(temp_means)
  glm_coefs_means_ana_volc <- rbind(glm_coefs_means_ana_volc,temp_means)
}
names(glm_coefs_means_ana_volc) <- names(glm_coef_dist_ana_volc)

#making a loop to convert character list rows to numeric
for(i in 2:7){
  glm_coefs_means_ana_volc[,i] <-as.numeric(unlist(glm_coefs_means_ana_volc[,i]))
}


# Cladogenesis means
glm_coefs_means_clad_volc <- data.frame()
for (i in 1:3) {
  temp_means <- vector()
  temp_means <- apply(glm_coefs_list_clad_volc[[i]], 2, function(x) mean(as.numeric(unlist(x))))
  temp_means["term"] <- glm_coefs_list_clad_volc[[i]][1,1]
  temp_means["type"] <- glm_coefs_list_clad_volc[[i]][1,8]
  glm_coefs_means_clad_volc <- rbind(glm_coefs_means_clad_volc,temp_means)
}

names(glm_coefs_means_clad_volc) <- names(glm_coef_area_clad_volc)

for(i in 2:7){
  glm_coefs_means_clad_volc[i] <- as.numeric(unlist(glm_coefs_means_clad_volc[i]))
}


# Endemics means
glm_coefs_means_endems_volc <- data.frame()
for (i in 1:4) {
  temp_means <- vector()
  temp_means <- apply(glm_coefs_list_endems_volc[[i]], 2, function(x) mean(as.numeric(unlist(x))))
  temp_means["term"] <- glm_coefs_list_endems_volc[[i]][1,1]
  temp_means["type"] <- glm_coefs_list_endems_volc[[i]][1,8]
  glm_coefs_means_endems_volc <- rbind(glm_coefs_means_endems_volc,temp_means)
}
names(glm_coefs_means_endems_volc) <- names(glm_coef_dist_endems_volc)



# Stacking the means dataframes for the means
coefs_means_volc <- NULL
coefs_means_volc <- rbind(glm_coefs_means_clad_volc,
                     glm_coefs_means_ana_volc,
                     glm_coefs_means_endems_volc
                     #glm_coefs_means_prop_endems,
                    #glm_coefs_means_regional_allo_prop
                     )
# I need to add all combinations of models and variables to have an even plot
# So all in all I need 18 "means"

# Changing names in dataframe of estimates
for (i in 1:length(glm_coefs_names_volc)){
  glm_coefs_volc[which(glm_coefs_volc$term==glm_coefs_names_df_volc[i,2]),1] <- glm_coefs_names_df_volc[i,1]
}

# Changing names in dataframe of means
for (i in 1:length(glm_coefs_names_volc)){
  coefs_means_volc[which(coefs_means_volc$term==glm_coefs_names_df_volc[i,2]),1] <- glm_coefs_names_df_volc[i,1]
}


glm_coefs_volc <- AddVars_violin(glm_coefs_volc,terms_list,type_list)
coefs_means_volc <- AddVars_mean(coefs_means_volc,terms_list,type_list)

# making a loop to convert character list rows to numeric 
for(i in 2:7){
  coefs_means_volc[,i] <-as.numeric(unlist(coefs_means_volc[,i]))
}

for(i in 2:7){
  glm_coefs_volc[,i] <-as.numeric(unlist(glm_coefs_volc[,i]))
}


# Creating plotting order
plot_order <- c("Area", "Isolation","Isolation Squared","Max Elevation","Fragmentation","Continental","Mixed Origin")

#Now I should be ready to make my violin + box plot for all species on all the islands.
gvil_all_all_volc <- ggplot(glm_coefs_volc, aes(color=type, fill = type)) +
  scale_colour_manual(values = c("#efc86e","#6f9969","#454a74")) +
  scale_fill_manual(values = c("#efc86e","#6f9969","#454a74")) +
  scale_x_discrete(limits = plot_order) +
  geom_violin(aes(x = term, y = estimate),width = 0.8,
             position=position_dodge(width=1),
              scale = "width", alpha = 0.8,lwd=0.2) +
  geom_violin(aes(x = term, y = conf.high_95),width = 0.8,
             position=position_dodge(width=1),
              scale = "width", alpha = 0.8,lwd=0.2) +
  geom_violin(aes(x = term, y = conf.low_95),width = 0.8,
             position=position_dodge(width=1),
              scale = "width", alpha = 0.8,lwd=0.2) +
  geom_crossbar(data = coefs_means_volc, size=0.1, alpha=0.7, aes(x= term, y = estimate,
                   ymin=conf.low_95,
                   ymax=conf.high_95),
                na.rm = TRUE,
             position = position_dodge(width = 1)) +
  geom_text(data = coefs_means_volc,
            aes(x = term, y = conf.high_95 + 0.3, label = ifelse(p.value < 0.05, "*", " ")),
            position = position_dodge(width = 1),
            hjust = -0.7, size = 5, show.legend = FALSE, na.rm = TRUE)+
  geom_hline(yintercept = 0, colour = gray(0), lty = 3) +
  xlab("Coefficients") +
  ylab("Estimate") +
  scale_y_continuous(breaks = seq(-10, 15, by = 5), limits = c(-10,15)) +
  theme_classic() +
  theme(axis.title.y = element_blank(),legend.position = "none") +
  ggtitle("All Plants: only Volcanic Islands")

gvil_all_all_volc

prow_3 <- cowplot::plot_grid(gvil_coryphoideae,gvil_all_all,gvil_all_all_volc, labels = c("a.", "b.", "c."), align = "hv", scale = c(1,1,1,1), nrow = 3, ncol = 1)
legend_bottom <- get_legend(gvil_coryphoideae + guides(color = guide_legend(nrow=1))+theme(legend.position = "bottom", legend.text=element_text(size=6)))

prow_3 <- cowplot::plot_grid(prow_3,legend_bottom, ncol = 1, rel_heights = c(1,.05))

prow_3





prow4 <- align_plots(gvil_all_all, gvil_coryphoideae,gvil_all_all_volc, align = "hv")

legend_right <- get_legend(gvil_coryphoideae + guides(color = guide_legend(ncol=1, title="Model"), fill = guide_legend(ncol=1, title="Model"))+theme(legend.position = "right",legend.text=element_text(size=9),panel.border = element_rect(colour = "black", fill=NA)))

bottom_row <- plot_grid(
  prow4[[3]],legend_right,
  labels = c("c."),
  rel_widths = c(1,1),
  nrow = 1
)

top_row <- plot_grid(
  prow4[[1]],prow4[[2]],
  labels = c("a.","b."),
  rel_widths = c(1,1),
  nrow = 1
)

prow_5 <- plot_grid(top_row,bottom_row, nrow = 2)
prow_5


# Writing tables for the supplement
write_csv(glm_coefs,file.path("./Figures/supplement","glm_coefs.csv")) # Estimates for all plants all islands
write_csv(coefs_means,file.path("./Figures/supplement","glm_coefs_means.csv")) # Means for all plants all islands 
write_csv(results_cory,file.path("./Figures/supplement","glm_coefs_coryphoideae.csv"))
write_csv(glm_coefs_volc,file.path("./Figures/supplement","glm_coefs_volc.csv")) # EStimates for all plants volcanic origin
write_csv(coefs_means_volc,file.path("./Figures/supplement","glm_coefs_means_volc.csv")) # Means for all plants volcanic origin
#write_csv(vif_res2,file.path("./Figures/supplement","Vif_all_models.csv")) # Variance inflation factors for all models 



(coefs_means[1,2]-coefs_means[5,2])/sqrt((coefs_means[1,3]^2)+(coefs_means[5,3]^2))


plot_col_hump <- ggplot(data = output_all_sp_congregated,aes(dist,(colonization_sp_mean+radiating_nodes_mean))) +
  geom_point(aes(color=GeologicalOrigin, ), size = 3.5) +
  xlab("Isolation") +
  ylab("Anagetic speciation events") +
  ggrepel::geom_text_repel(aes(label=ifelse(colonization_edges_mean>500,as.character(LEVEL3_NAM),'')),nudge_y=0.03,nudge_x=-0.10, size=3) +
  scale_colour_manual(values = met.brewer("Pillement", 4)) +
  scale_y_continuous(breaks = seq(0, 7000, by = 1000), limits = c(0,7000)) +
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)), alpha = 0.3, ) +
  theme_classic() +
  theme(legend.position = "bottom")

plot_col_hump




# plot_clad_volc_height <- ggplot(data = output_all_sp_congregated[which(output_all_sp_congregated$GeologicalOrigin == "volcanic"),],aes((max_elev_30m),(radiating_sp_mean-radiating_nodes_mean))) +
#   geom_point(aes(color=GeologicalOrigin, ), size = 3.5) +
#   xlab("Maximum Elevation") +
#   ylab("Cladogenetic speciation events") +
#   ggrepel::geom_text_repel(aes(label=ifelse(radiating_sp_mean-radiating_nodes_mean>10,as.character(LEVEL3_NAM),'')),nudge_y=0.03,nudge_x=-0.10, size=3) +
#   scale_colour_manual(values = met.brewer("Pillement", 4)) +
#   scale_x_continuous(breaks = c(10,20,50,100,250,500,1000,2000,4500), limits = c(0,5000)) +
#   scale_y_continuous(breaks = c(10,25,50,100,250,400), limits = c(0,350))

# plot_clad_volc_height



# Assuming output_all_sp_test_subset_log is a list of 100 datasets

# Initialize vectors to store the coefficient estimates, p-values, and R-squared values
coefficients <- numeric(100)
p_values <- numeric(100)
r_squared <- numeric(100)

# Loop over all datasets and fit the model
for (i in 1:100) {
  model <- lm(colonization_sp + radiating_nodes ~ radiating_sp - radiating_nodes, data = output_all_sp_test_subset_log[[i]])
  summary_model <- summary(model)
  
  # Extract the coefficient estimate and p-value for radiating_sp
  coefficients[i] <- summary_model$coefficients["radiating_sp", "Estimate"]
  p_values[i] <- summary_model$coefficients["radiating_sp", "Pr(>|t|)"]
  
  # Extract the R-squared value
  r_squared[i] <- summary_model$r.squared
}

# Calculate the average coefficient estimate, average p-value, and average R-squared
average_coefficient <- mean(coefficients)
average_p_value <- mean(p_values)
average_r_squared <- mean(r_squared)

# Print the results
cat("Average coefficient estimate for radiating_sp:", average_coefficient, "\n")
cat("Average p-value for radiating_sp:", average_p_value, "\n")
cat("Average R-squared:", average_r_squared, "\n")







