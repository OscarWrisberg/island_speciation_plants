# Script which prints the VIF values of all the models used
###############
# All Islands #
###############

###############################
# Cladogenesis on all islands #
###############################
#Models
clad_model1 <- (radiating_sp-radiating_nodes~area+poly(dist,1,raw=TRUE)+nearest_neighbour_distance_border_scaled)
clad_model2 <- (radiating_sp-radiating_nodes~area+poly(dist,1,raw=TRUE)+nearest_neighbour_distance_border_scaled+max30_elev+GeologicalOrigin)
clad_model3 <- (radiating_sp-radiating_nodes~area+poly(dist,1,raw=TRUE)+nearest_neighbour_distance_border_scaled+GeologicalOrigin)
clad_model4 <- (radiating_sp-radiating_nodes~area+nearest_neighbour_distance_border_scaled+max30_elev+GeologicalOrigin)
#Model 1

glm_clad_pois_all_sp_list_models_vif1 <- list()

for (i in 1:100){
  glm_clad_pois_all_sp_vif1 <- glm(clad_model1,data=output_all_sp_test_subset_log[[i]][which(output_all_sp_test_subset_log[[i]]$GeologicalOrigin != "atoll"),],family = quasipoisson(link="log"),na.action = na.fail)
  glm_clad_pois_all_sp_list_models_vif1[[i]] <- glm_clad_pois_all_sp_vif1
}

clad_vif_1 <- list()
for (i in 1:100){
  clad_vif_1 <- max(car::vif(glm_clad_pois_all_sp_list_models_vif1[[i]]))
  }
clad_max_vif_1 <- max(clad_vif_1)

#Model 2
glm_clad_pois_all_sp_list_models_vif2 <- list()
for (i in 1:100){
  glm_clad_pois_all_sp <- glm(clad_model2,data=output_all_sp_test_subset_log[[i]][which(output_all_sp_test_subset_log[[i]]$GeologicalOrigin != "atoll"),],family = quasipoisson(link="log"),na.action = na.fail)
  glm_clad_pois_all_sp_list_models_vif2[[i]] <- glm_clad_pois_all_sp
}

clad_vif_2 <- list()
for (i in 1:100){
  print(car::vif(glm_clad_pois_all_sp_list_models_vif2[[i]]))
  clad_vif_2 <- max(car::vif(glm_clad_pois_all_sp_list_models_vif2[[i]])[,1])
}
clad_max_vif_2 <- max(clad_vif_2)

#Model 3
glm_clad_pois_all_sp_list_models_vif3 <- list()
for (i in 1:100){
  glm_clad_pois_all_sp <- glm(clad_model3,data=output_all_sp_test_subset_log[[i]][which(output_all_sp_test_subset_log[[i]]$GeologicalOrigin != "atoll"),],family = quasipoisson(link="log"),na.action = na.fail)
  glm_clad_pois_all_sp_list_models_vif3[[i]] <- glm_clad_pois_all_sp
}

clad_vif_3 <- list()
for (i in 1:100){
  print(car::vif(glm_clad_pois_all_sp_list_models_vif3[[i]]))
  clad_vif_3 <- max(car::vif(glm_clad_pois_all_sp_list_models_vif3[[i]])[,1])
}
clad_max_vif_3 <- max(clad_vif_3)

#Model 4
glm_clad_pois_all_sp_list_models_vif4 <- list()
for (i in 1:100){
  glm_clad_pois_all_sp <- glm(clad_model4,data=output_all_sp_test_subset_log[[i]][which(output_all_sp_test_subset_log[[i]]$GeologicalOrigin != "atoll"),],family = quasipoisson(link="log"),na.action = na.fail)
  glm_clad_pois_all_sp_list_models_vif4[[i]] <- glm_clad_pois_all_sp
}

clad_vif_4 <- list()
for (i in 1:100){
  print(car::vif(glm_clad_pois_all_sp_list_models_vif4[[i]]))
  clad_vif_4 <- max(car::vif(glm_clad_pois_all_sp_list_models_vif4[[i]])[,1])
}
clad_max_vif_4 <- max(clad_vif_4)

clad_max_vif_seed <- c(clad_max_vif_1,clad_max_vif_2,clad_max_vif_3,clad_max_vif_4)


#############################
# Anagenesis on all islands #
#############################
ana_model1 <- (colonization_sp+radiating_nodes~area+poly(dist,2,raw=TRUE)+max30_elev+nearest_neighbour_distance_border_scaled)
ana_model2 <- (colonization_sp+radiating_nodes~area+poly(dist,2,raw=TRUE)+max30_elev)
ana_model3 <- (colonization_sp+radiating_nodes~area+poly(dist,1,raw=TRUE)+GeologicalOrigin+max30_elev+nearest_neighbour_distance_border_scaled+GeologicalOrigin)
ana_model4 <- (colonization_sp+radiating_nodes~area+poly(dist,2,raw=TRUE)+nearest_neighbour_distance_border_scaled)


#Model 1
glm_ana_pois_all_sp_list_models_vif1 <- list()

for (i in 1:100){
  glm_ana_pois_all_sp_vif1 <- glm(ana_model1,data=output_all_sp_test_subset_log[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_ana_pois_all_sp_list_models_vif1[[i]] <- glm_ana_pois_all_sp_vif1
}

ana_vif_1 <- list()
for (i in 1:100){
  ana_vif_1 <- max(car::vif(glm_ana_pois_all_sp_list_models_vif1[[i]])[,1])
}
ana_max_vif_1 <- max(ana_vif_1)

#Model 2
glm_ana_pois_all_sp_list_models_vif2 <- list()
for (i in 1:100){
  glm_ana_pois_all_sp <- glm(ana_model2,data=output_all_sp_test_subset_log[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_ana_pois_all_sp_list_models_vif2[[i]] <- glm_ana_pois_all_sp
}

ana_vif_2 <- list()
for (i in 1:100){
  #print(car::vif(glm_ana_pois_all_sp_list_models_vif2[[i]]))
  ana_vif_2 <- max(car::vif(glm_ana_pois_all_sp_list_models_vif2[[i]])[,1])
}
ana_max_vif_2 <- max(ana_vif_2)

#Model 3
glm_ana_pois_all_sp_list_models_vif3 <- list()
for (i in 1:100){
  glm_ana_pois_all_sp <- glm(ana_model3,data=output_all_sp_test_subset_log[[i]][which(output_all_sp_test_subset_log[[i]]$GeologicalOrigin != "atoll"),],family = quasipoisson(link="log"),na.action = na.fail)
  glm_ana_pois_all_sp_list_models_vif3[[i]] <- glm_ana_pois_all_sp
}

ana_vif_3 <- list()
for (i in 1:100){
  #print(car::vif(glm_ana_pois_all_sp_list_models_vif3[[i]]))
  ana_vif_3 <- max(car::vif(glm_ana_pois_all_sp_list_models_vif3[[i]])[,1])
}
ana_max_vif_3 <- max(ana_vif_3)

#Model 4
glm_ana_pois_all_sp_list_models_vif4 <- list()
for (i in 1:100){
  glm_ana_pois_all_sp <- glm(ana_model4,data=output_all_sp_test_subset_log[[i]][which(output_all_sp_test_subset_log[[i]]$GeologicalOrigin != "atoll"),],family = quasipoisson(link="log"),na.action = na.fail)
  glm_ana_pois_all_sp_list_models_vif4[[i]] <- glm_ana_pois_all_sp
}

ana_vif_4 <- list()
for (i in 1:100){
  print(car::vif(glm_ana_pois_all_sp_list_models_vif4[[i]]))
  ana_vif_4 <- max(car::vif(glm_ana_pois_all_sp_list_models_vif4[[i]])[,1])
}
ana_max_vif_4 <- max(ana_vif_4)

ana_max_vif_seed <- c(ana_max_vif_1,ana_max_vif_2,ana_max_vif_3,ana_max_vif_4)

#############################
#  Endemics on all islands  #
#############################
endem_model1 <- (colonization_sp+radiating_sp~area+poly(dist,2,raw=TRUE)+max30_elev+nearest_neighbour_distance_border_scaled)
endem_model2 <- (colonization_sp+radiating_sp~area+poly(dist,2,raw=TRUE)+max30_elev)
endem_model3 <- (colonization_sp+radiating_sp~area+poly(dist,2,raw=TRUE)+GeologicalOrigin+max30_elev+nearest_neighbour_distance_border_scaled+GeologicalOrigin)
endem_model4 <- (colonization_sp+radiating_sp~area+poly(dist,2,raw=TRUE)+nearest_neighbour_distance_border_scaled)


#Model 1
glm_endems_pois_all_sp_list_models_vif1 <- list()

for (i in 1:100){
  glm_endems_pois_all_sp_vif1 <- glm(endem_model1,data=output_all_sp_test_subset_log[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_endems_pois_all_sp_list_models_vif1[[i]] <- glm_endems_pois_all_sp_vif1
}

endems_vif_1 <- list()
for (i in 1:100){
  endems_vif_1 <- max(car::vif(glm_endems_pois_all_sp_list_models_vif1[[i]])[,1])
}
endems_max_vif_1 <- max(endems_vif_1)

#Model 2
glm_endems_pois_all_sp_list_models_vif2 <- list()
for (i in 1:100){
  glm_endems_pois_all_sp <- glm(endem_model2,data=output_all_sp_test_subset_log[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_endems_pois_all_sp_list_models_vif2[[i]] <- glm_endems_pois_all_sp
}

endems_vif_2 <- list()
for (i in 1:100){
  print(car::vif(glm_endems_pois_all_sp_list_models_vif2[[i]]))
  endems_vif_2 <- max(car::vif(glm_endems_pois_all_sp_list_models_vif2[[i]])[,1])
}
endems_max_vif_2 <- max(endems_vif_2)

#Model 3
glm_endems_pois_all_sp_list_models_vif3 <- list()
for (i in 1:100){
  glm_endems_pois_all_sp <- glm(endem_model3,data=output_all_sp_test_subset_log[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_endems_pois_all_sp_list_models_vif3[[i]] <- glm_endems_pois_all_sp
}

endems_vif_3 <- list()
for (i in 1:100){
  print(car::vif(glm_endems_pois_all_sp_list_models_vif3[[i]]))
  endems_vif_3 <- max(car::vif(glm_endems_pois_all_sp_list_models_vif3[[i]])[,1])
}
endems_max_vif_3 <- max(endems_vif_3)

#Model 4
glm_endems_pois_all_sp_list_models_vif4 <- list()
for (i in 1:100){
  glm_endems_pois_all_sp <- glm(endem_model4,data=output_all_sp_test_subset_log[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_endems_pois_all_sp_list_models_vif4[[i]] <- glm_endems_pois_all_sp
}

endems_vif_4 <- list()
for (i in 1:100){
  print(car::vif(glm_endems_pois_all_sp_list_models_vif4[[i]]))
  endems_vif_4 <- max(car::vif(glm_endems_pois_all_sp_list_models_vif4[[i]])[,1])
}
endems_max_vif_4 <- max(endems_vif_4)

endems_max_vif_seed <- c(endems_max_vif_1,endems_max_vif_2,endems_max_vif_3,endems_max_vif_4)

vif_all_seed <- c(ana_max_vif_seed,clad_max_vif_seed,endems_max_vif_seed)
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################

#####################################
# Only islands with Volcanic origin #
#####################################
####################################
# Cladogenesis on volcanic islands #
####################################
#Models
clad_volc_model1 <- (radiating_sp-radiating_nodes~area+poly(dist,1,raw=TRUE)+max30_elev)
clad_volc_model2 <- (radiating_sp-radiating_nodes~area+poly(dist,1,raw=TRUE)+max30_elev+nearest_neighbour_distance_border_scaled)
clad_volc_model3 <- (radiating_sp-radiating_nodes~poly(dist,1,raw=TRUE)+max30_elev+nearest_neighbour_distance_border_scaled)
clad_volc_model4 <- (radiating_sp-radiating_nodes~poly(dist,1,raw=TRUE)+max30_elev)

#Model 1
glm_clad_volc_pois_all_sp_list_models_vif1 <- list()

for (i in 1:100){
  glm_clad_volc_pois_all_sp_vif1 <- glm(clad_volc_model1,data=test_dat_log_volc[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_clad_volc_pois_all_sp_list_models_vif1[[i]] <- glm_clad_pois_all_sp_vif1
}

clad_volc_vif_1 <- list()
for (i in 1:100){
  print(car::vif(glm_clad_volc_pois_all_sp_list_models_vif1[[i]]))
  clad_volc_vif_1 <- max(car::vif(glm_clad_volc_pois_all_sp_list_models_vif1[[i]]))
}
clad_volc_max_vif_1 <- max(clad_volc_vif_1)

#Model 2
glm_clad_volc_pois_all_sp_list_models_vif2 <- list()
for (i in 1:100){
  glm_clad_volc_pois_all_sp <- glm(clad_volc_model2,data=test_dat_log_volc[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_clad_volc_pois_all_sp_list_models_vif2[[i]] <- glm_clad_volc_pois_all_sp
}

clad_volc_vif_2 <- list()
for (i in 1:100){
  print(car::vif(glm_clad_volc_pois_all_sp_list_models_vif2[[i]]))
  clad_volc_vif_2 <- max(car::vif(glm_clad_volc_pois_all_sp_list_models_vif2[[i]]))
}
clad_volc_max_vif_2 <- max(clad_volc_vif_2)

#Model 3
glm_clad_volc_pois_all_sp_list_models_vif3 <- list()
for (i in 1:100){
  glm_clad_volc_pois_all_sp <- glm(clad_volc_model3,data=test_dat_log_volc[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_clad_volc_pois_all_sp_list_models_vif3[[i]] <- glm_clad_volc_pois_all_sp
}

clad_volc_vif_3 <- list()
for (i in 1:100){
  print(car::vif(glm_clad_volc_pois_all_sp_list_models_vif3[[i]]))
  clad_volc_vif_3 <- max(car::vif(glm_clad_volc_pois_all_sp_list_models_vif3[[i]]))
}
clad_volc_max_vif_3 <- max(clad_volc_vif_3)

#Model 4
glm_clad_volc_pois_all_sp_list_models_vif4 <- list()
for (i in 1:100){
  glm_clad_volc_pois_all_sp <- glm(clad_volc_model4,data=test_dat_log_volc[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_clad_volc_pois_all_sp_list_models_vif4[[i]] <- glm_clad_volc_pois_all_sp
}

clad_volc_vif_4 <- list()
for (i in 1:100){
  print(car::vif(glm_clad_volc_pois_all_sp_list_models_vif4[[i]]))
  clad_volc_vif_4 <- max(car::vif(glm_clad_volc_pois_all_sp_list_models_vif4[[i]]))
}
clad_volc_max_vif_4 <- max(clad_volc_vif_4)

clad_volc_max_vif_seed <- c(clad_volc_max_vif_1,clad_volc_max_vif_2,clad_volc_max_vif_3,clad_volc_max_vif_4)


##################################
# Anagenesis on volcanic islands #
##################################
ana_volc_model1 <- (colonization_sp+radiating_nodes~area+poly(dist,2,raw=TRUE)+nearest_neighbour_distance_border_scaled)
ana_volc_model2 <- (colonization_sp+radiating_nodes~area+poly(dist,2,raw=TRUE)+max30_elev+nearest_neighbour_distance_border_scaled)
ana_volc_model3 <- (colonization_sp+radiating_nodes~area+poly(dist,2,raw=TRUE))
ana_volc_model4 <- (colonization_sp+radiating_nodes~area)


#Model 1
glm_ana_volc_pois_all_sp_list_models_vif1 <- list()

for (i in 1:100){
  glm_ana_volc_pois_all_sp_vif1 <- glm(ana_volc_model1,data=test_dat_log_volc[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_ana_volc_pois_all_sp_list_models_vif1[[i]] <- glm_ana_volc_pois_all_sp_vif1
}

ana_volc_vif_1 <- list()
for (i in 1:100){
  ana_volc_vif_1 <- max(car::vif(glm_ana_volc_pois_all_sp_list_models_vif1[[i]])[,1])
}
ana_volc_max_vif_1 <- max(ana_volc_vif_1)

#Model 2
glm_ana_volc_pois_all_sp_list_models_vif2 <- list()
for (i in 1:100){
  glm_ana_volc_pois_all_sp <- glm(ana_volc_model2,data=test_dat_log_volc[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_ana_volc_pois_all_sp_list_models_vif2[[i]] <- glm_ana_volc_pois_all_sp
}

ana_volc_vif_2 <- list()
for (i in 1:100){
  #print(car::vif(glm_ana_pois_all_sp_list_models_vif2[[i]]))
  ana_volc_vif_2 <- max(car::vif(glm_ana_pois_all_sp_list_models_vif2[[i]])[,1])
}
ana_volc_max_vif_2 <- max(ana_volc_vif_2)

#Model 3
glm_ana_volc_pois_all_sp_list_models_vif3 <- list()
for (i in 1:100){
  glm_ana_volc_pois_all_sp <- glm(ana_volc_model3,data=test_dat_log_volc[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_ana_volc_pois_all_sp_list_models_vif3[[i]] <- glm_ana_volc_pois_all_sp
}

ana_vif_3 <- list()
for (i in 1:100){
  #print(car::vif(glm_ana_pois_all_sp_list_models_vif3[[i]]))
  ana_volc_vif_3 <- max(car::vif(glm_ana_pois_all_sp_list_models_vif3[[i]])[,1])
}
ana_volc_max_vif_3 <- max(ana_volc_vif_3)

#Model 4
# Cant do a Variance inflation investigation on a model with only 1 variable
# glm_ana_volc_pois_all_sp_list_models_vif4 <- list()
# for (i in 1:100){
#   glm_ana_volc_pois_all_sp <- glm(ana_volc_model4,data=test_dat_log_volc[[i]],family = quasipoisson(link="log"),na.action = na.fail)
#   glm_ana_volc_pois_all_sp_list_models_vif4[[i]] <- glm_ana_volc_pois_all_sp
# }
# 
# ana_volc_vif_4 <- list()
# for (i in 1:100){
#   print(car::vif(glm_ana_volc_pois_all_sp_list_models_vif4[[i]]))
#   ana_vif_4 <- max(car::vif(glm_ana_volc_pois_all_sp_list_models_vif4[[i]])[,1])
# }
# ana_volc_max_vif_4 <- max(ana_volc_vif_4)
# 
ana_volc_max_vif_seed <- c(ana_volc_max_vif_1,ana_volc_max_vif_2,ana_volc_max_vif_3) # ana_volc_max_vif_4

##################################
#  Endemics on volcanic islands  #
##################################
endem_volc_model1 <- (colonization_sp+radiating_sp~area+poly(dist,2,raw=TRUE)+nearest_neighbour_distance_border_scaled)
endem_volc_model2 <- (colonization_sp+radiating_sp~area+poly(dist,2,raw=TRUE)+max30_elev+nearest_neighbour_distance_border_scaled)
endem_volc_model3 <- (colonization_sp+radiating_sp~area+poly(dist,2,raw=TRUE))
endem_volc_model4 <- (colonization_sp+radiating_sp~area)


#Model 1
glm_endems_volc_pois_all_sp_list_models_vif1 <- list()

for (i in 1:100){
  glm_endems_volc_pois_all_sp_vif1 <- glm(endem_volc_model1,data=test_dat_log_volc[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_endems_volc_pois_all_sp_list_models_vif1[[i]] <- glm_endems_volc_pois_all_sp_vif1
}

endems_volc_vif_1 <- list()
for (i in 1:100){
  print(car::vif(glm_endems_volc_pois_all_sp_list_models_vif1[[i]]))
  endems_volc_vif_1 <- max(car::vif(glm_endems_volc_pois_all_sp_list_models_vif1[[i]])[,1])
}
endems_volc_max_vif_1 <- max(endems_volc_vif_1)

#Model 2
glm_endems_volc_pois_all_sp_list_models_vif2 <- list()
for (i in 1:100){
  glm_endems_volc_pois_all_sp <- glm(endem_volc_model2,data=test_dat_log_volc[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_endems_volc_pois_all_sp_list_models_vif2[[i]] <- glm_endems_pois_all_sp
}

endems_volc_vif_2 <- list()
for (i in 1:100){
  print(car::vif(glm_endems_volc_pois_all_sp_list_models_vif2[[i]]))
  endems_volc_vif_2 <- max(car::vif(glm_endems_volc_pois_all_sp_list_models_vif2[[i]])[,1])
}
endems_volc_max_vif_2 <- max(endems_volc_vif_2)

#Model 3
glm_endems_volc_pois_all_sp_list_models_vif3 <- list()
for (i in 1:100){
  glm_endems_volc_pois_all_sp <- glm(endem_volc_model3,data=test_dat_log_volc[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_endems_volc_pois_all_sp_list_models_vif3[[i]] <- glm_endems_volc_pois_all_sp
}

endems_volc_vif_3 <- list()
for (i in 1:100){
  print(car::vif(glm_endems_volc_pois_all_sp_list_models_vif3[[i]]))
  endems_volc_vif_3 <- max(car::vif(glm_endems_volc_pois_all_sp_list_models_vif3[[i]])[,1])
}
endems_volc_max_vif_3 <- max(endems_volc_vif_3)

#Model 4
glm_endems_pois_all_sp_list_models_vif4 <- list()
for (i in 1:100){
  glm_endems_pois_all_sp <- glm(endem_volc_model4,data=test_dat_log_volc[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_endems_pois_all_sp_list_models_vif4[[i]] <- glm_endems_pois_all_sp
}

endems_vif_4 <- list()
for (i in 1:100){
  print(car::vif(glm_endems_pois_all_sp_list_models_vif4[[i]]))
  endems_vif_4 <- max(car::vif(glm_endems_pois_all_sp_list_models_vif4[[i]])[,1])
}
endems_volc_max_vif_4 <- max(endems_vif_4)

endems_volc_max_vif_seed <- c(endems_volc_max_vif_1,endems_volc_max_vif_2,endems_volc_max_vif_3,endems_volc_max_vif_4)

vif_all_seed_volc <- c()
vif_all_seed_volc <- c(ana_volc_max_vif_seed,clad_volc_max_vif_seed,endems_volc_max_vif_seed)



#########################################################################################################################################################
#########################################################################################################################################################
########################################################################################################################################################

#########################
#  Coryphoideae Models  #
#########################


#############################
# Anagenesis on Coryphoids #
#############################
ana_cory_model1 <- (`No. anagenesis`+`No. edges`~area+nearest_neighbour_distance_border_scaled)
ana_cory_model2 <- (`No. anagenesis`+`No. edges`~area)
ana_cory_model3 <- (`No. anagenesis`+`No. edges`~area+max30_elev+nearest_neighbour_distance_border_scaled)
ana_cory_model4 <- (`No. anagenesis`+`No. edges`~area+GeologicalOrigin+nearest_neighbour_distance_border_scaled)


#Model 1
glm_ana_cory_pois_all_sp_vif1 <- glm(ana_cory_model1,data=test_dat_log_cory,family = quasipoisson(link="log"),na.action = na.fail)
ana_cory_vif_1 <- max(car::vif(glm_ana_cory_pois_all_sp_vif1))
ana_cory_max_vif_1 <- max(ana_cory_vif_1)

#Model 2
# glm_ana_cory_pois_all_sp <- glm(ana_cory_model2,data=test_dat_log_cory,family = quasipoisson(link="log"),na.action = na.fail)
# ana_cory_vif_2 <- max(car::vif(glm_ana_cory_pois_all_sp))
# ana_cory_max_vif_2 <- max(ana_cory_vif_2)


#Model 3
glm_ana_cory_pois_all_sp_vif3 <- glm(ana_cory_model3,data=test_dat_log_cory,family = quasipoisson(link="log"),na.action = na.fail)
ana_cory_vif_3 <- max(car::vif(glm_ana_cory_pois_all_sp_vif3))
ana_cory_max_vif_3 <- max(ana_cory_vif_3)

#Model 4
glm_ana_cory_pois_all_sp_vif4 <- glm(ana_cory_model4,data=test_dat_log_cory,family = quasipoisson(link="log"),na.action = na.fail)
ana_cory_vif_4 <- max(car::vif(glm_ana_cory_pois_all_sp_vif4)[,1])
ana_cory_max_vif_4 <- max(ana_cory_vif_4)

ana_cory_max_vif_seed <- c(ana_cory_max_vif_1,NA,ana_cory_max_vif_3,ana_cory_max_vif_4) #ana_cory_max_vif_2


#############################
# Cladogenesis on Coryphoids #
#############################
clad_cory_model1 <- (`No. cladogenesis`-`No. edges`~area)
clad_cory_model2 <- (`No. cladogenesis`-`No. edges`~area+poly(dist,1,raw=TRUE))
clad_cory_model3 <- (`No. cladogenesis`-`No. edges`~area+nearest_neighbour_distance_border_scaled)
clad_cory_model4 <- (`No. cladogenesis`-`No. edges`~area+max30_elev)


#Model 1
# glm_clad_cory_pois_all_sp_vif1 <- glm(clad_cory_model1,data=test_dat_log_cory,family = quasipoisson(link="log"),na.action = na.fail)
# clad_cory_vif_1 <- max(car::vif(glm_clad_cory_pois_all_sp_vif1))
# clad_cory_max_vif_1 <- max(clad_cory_vif_1)

#Model 2
glm_clad_cory_pois_all_sp <- glm(clad_cory_model2,data=test_dat_log_cory,family = quasipoisson(link="log"),na.action = na.fail)
clad_cory_vif_2 <- max(car::vif(glm_clad_cory_pois_all_sp))
clad_cory_max_vif_2 <- max(clad_cory_vif_2)


#Model 3
glm_clad_cory_pois_all_sp_vif3 <- glm(clad_cory_model3,data=test_dat_log_cory,family = quasipoisson(link="log"),na.action = na.fail)
clad_cory_vif_3 <- max(car::vif(glm_clad_cory_pois_all_sp_vif3))
clad_cory_max_vif_3 <- max(clad_cory_vif_3)

#Model 4
glm_clad_cory_pois_all_sp_vif4 <- glm(clad_cory_model4,data=test_dat_log_cory,family = quasipoisson(link="log"),na.action = na.fail)
clad_cory_vif_4 <- max(car::vif(glm_clad_cory_pois_all_sp_vif4))
clad_cory_max_vif_4 <- max(clad_cory_vif_4)

clad_cory_max_vif_seed <- c(NA,clad_cory_max_vif_2,clad_cory_max_vif_3,clad_cory_max_vif_4) #clad_cory_max_vif_1


#############################
# Endemics on Coryphoids #
#############################
endems_cory_model1 <- (`No. cladogenesis`-`No. anagenesis`~area)
endems_cory_model2 <- (`No. cladogenesis`+`No. anagenesis`~area+nearest_neighbour_distance_border_scaled)
endems_cory_model3 <- (`No. cladogenesis`+`No. anagenesis`~area+max30_elev)
endems_cory_model4 <- (`No. cladogenesis`+`No. anagenesis`~area+GeologicalOrigin)

#Model 1
# glm_endems_cory_pois_all_sp_vif1 <- glm(endems_cory_model1,data=test_dat_log_cory,family = quasipoisson(link="log"),na.action = na.fail)
# endems_cory_vif_1 <- max(car::vif(glm_endems_cory_pois_all_sp_vif1))
# endems_cory_max_vif_1 <- max(endems_cory_vif_1)

#Model 2
glm_endems_cory_pois_all_sp <- glm(endems_cory_model2,data=test_dat_log_cory,family = quasipoisson(link="log"),na.action = na.fail)
endems_cory_vif_2 <- max(car::vif(glm_endems_cory_pois_all_sp))
endems_cory_max_vif_2 <- max(endems_cory_vif_2)


#Model 3
glm_endems_cory_pois_all_sp_vif3 <- glm(endems_cory_model3,data=test_dat_log_cory,family = quasipoisson(link="log"),na.action = na.fail)
endems_cory_vif_3 <- max(car::vif(glm_endems_cory_pois_all_sp_vif3))
endems_cory_max_vif_3 <- max(endems_cory_vif_3)

#Model 4
glm_endems_cory_pois_all_sp_vif4 <- glm(endems_cory_model4,data=test_dat_log_cory,family = quasipoisson(link="log"),na.action = na.fail)
endems_cory_vif_4 <- max(car::vif(glm_endems_cory_pois_all_sp_vif4))
endems_cory_max_vif_4 <- max(endems_cory_vif_4)

endems_cory_max_vif_seed <- c(NA,endems_cory_max_vif_2,endems_cory_max_vif_3,endems_cory_max_vif_4) #endems_cory_max_vif_1
