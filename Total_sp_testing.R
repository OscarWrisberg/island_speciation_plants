# Finding the best model for the cladogenesis
best_model_vars_total<- data.frame()

for (i in 1:100){
  test_dat_log <- output_all_sp_test_subset_log[[i]]
  
  dredge_output <- qdredge(glm(radiating_nodes~mean_mx30_grd+area+dist+I(dist*dist)+nearest_neighbour_distance_border_scaled+GeologicalOrigin,data=test_dat_log,family = quasipoisson,na.action = na.fail))
  
  best_model_vars_total <- rbind(best_model_vars_total, dredge_output[which(dredge_output$delta == 0),1:9])
}

#best_model_vars_test <- is.na(best_model_vars)
best_model_vars_count_total <- is.na(best_model_vars_total) %>% as.data.frame() %>% group_by_all() %>% count()
best_model_vars_count_total


glm_total_pois_all_sp_list_tidy <- list()
glm_total_pois_all_sp_list_models <- list()
glm_total_pois_all_sp_list_dispersion <- list()

# Loop which fit the model to each of the 100 datasets
for (i in 1:100){
  glm_total_pois_all_sp <- glm(Total_sp~area+dist+I(dist*dist)+GeologicalOrigin,data=output_all_sp_test_subset_log[[i]],family = quasipoisson(link="log"),na.action = na.fail)
  glm_total_pois_all_sp_list_models[[i]] <- glm_total_pois_all_sp
  glm_total_pois_all_sp_list_tidy[[i]] <- broom::tidy(glm_total_pois_all_sp)
  glm_total_pois_all_sp_list_dispersion[[i]] <- sum(residuals(glm_total_pois_all_sp,type ="pearson")^2)/glm_total_pois_all_sp$df.residual
}
glm_total_pois_all_sp_list_tidy

plot(output_all_sp_test_subset_log[[1]]$mean_mx30_grd, output_all_sp_test_subset_log[[1]]$Total_sp)
plot(output_all_sp_test_subset_log[[1]]$area, output_all_sp_test_subset_log[[1]]$radiating_nodes)
plot(output_all_sp_test_subset_log[[1]]$dist, output_all_sp_test_subset_log[[1]]$radiating_nodes)
plot(output_all_sp_test_subset_log[[1]]$dist*output_all_sp_test_subset_log[[1]]$area, output_all_sp_test_subset_log[[1]]$Total_sp)


poly_3 <-ggplot(output_all_sp_test_subset_log[[1]], aes(dist,Total_sp)) +
  geom_point() +
  stat_smooth(method = "glm",formula = y ~ poly(x, 4),se = TRUE, data = output_all_sp_test_subset_log[[1]])
poly_3


poly_4 <- ggplot(output_all_sp_test_subset_log[[1]], aes()) +
  geom_point(aes(x=dist,y=colonization_sp, colour = "Red")) +
  geom_point(aes(dist,radiating_sp)) +
  stat_smooth(method = "glm",formula = y ~ poly(x, 4),se = TRUE, data = output_all_sp_test_subset_log[[1]], aes(x=dist,y=colonization_sp, colour = "Red")) +
  stat_smooth(method = "glm",formula = y ~ x,se = TRUE, data = output_all_sp_test_subset_log[[1]][which(output_all_sp_test_subset_log[[1]]$radiating_sp >= 2),], aes(x=dist,y=radiating_sp))
poly_4




ggplot(test_dat_log_volc[[1]]) + 
  geom_point(aes(x=area,y=radiating_sp)) +
  stat_smooth(method = "glm",formula = y ~ poly(x,1),se = TRUE, data = output_all_sp_test_subset_log[[1]], aes(x=area,y=radiating_sp))
