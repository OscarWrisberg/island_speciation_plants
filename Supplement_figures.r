

# loading environment with data from Melanie
load(file.path(datadir,"/workspace_TACTed_species.RData"))
output_with_island_no_log <- readRDS(file.path(datadir,"output_for_test_all_sp"))
names(dat)[names(dat) =="bot_coutry"] <- "LEVEL3_COD"
#output_with_island <- merge(output_with_island, dat, by ="LEVEL3_COD")

# Finding the island botanical countries with the biggest and smalles size.
output_with_island_no_log[[1]][which.min(output_with_island_no_log[[1]]$area),]
output_with_island_no_log[[1]][which.max(output_with_island_no_log[[1]]$area),]

output_with_island_no_log[[1]][which.min(output_with_island_no_log[[1]]$dist),]
output_with_island_no_log[[1]][which.max(output_with_island_no_log[[1]]$dist),]


#Loading the one with only endemics
load(file.path(datadir,"/workspace_TACTed_species2.RData"))
names(dat2)[names(dat2) =="bot_coutry"] <- "LEVEL3_COD"
names(dat2)[names(dat2) =="non_TACTed"] <- "non_TACTed_endems"
names(dat2)[names(dat2) =="TACTed"] <- "TACTed_endems"
output_all_sp_congregated <- merge(output_all_sp_congregated, dat2, by ="LEVEL3_COD")


# Coefficient estimates for the models fitted to the 100 datasets based on the tacted phylogenies and plants on all islands
write_csv(glm_coefs,file.path("./Figures/supplement","glm_coefs.csv"))

# Means for the models fitted to the 100 datasets based on the tacted phylogenies and plants on all islands
write_csv(coefs_means,file.path("./Figures/supplement","glm_coefs_means.csv"))

# Results from the models fitted to the dataset created by the coryphoid phylogeny.
write_csv(results_cory,file.path("./Figures/supplement","glm_coefs_coryphoideae.csv"))

# Coefficient estimates for the models fitted to the 100 datasets based on the tacted phylogenies and plants on only volcanic islands
write_csv(glm_coefs_volc,file.path("./Figures/supplement","glm_coefs_volc.csv"))

# Means for the models fitted to the 100 datasets based on the tacted phylogenies and plants on only volcanic islands
write_csv(coefs_means_volc,file.path("./Figures/supplement","glm_coefs_means_volc.csv"))

# Variance inflation factors for all models 
#write_csv(vif_res2,file.path("./Figures/supplement","Vif_all_models.csv"))


# Plot showing correlation between Island Isolation and Anagenesis
plot_col_hump <- ggplot(data = output_all_sp_congregated,aes(dist,(colonization_sp_mean+radiating_nodes_mean))) +
  geom_point(aes(color=GeologicalOrigin, ), size = 3.5) +
  xlab("Isolation") +
  ylab("Anagetic speciation events") +
  ggrepel::geom_text_repel(aes(label=ifelse(colonization_edges_mean>500,as.character(LEVEL3_NAM),'')),nudge_y=0.03,nudge_x=-0.10, size=3) +
  scale_colour_manual(values = met.brewer("Pillement", 4)) +
  #scale_y_continuous(breaks = seq(0, 7000, by = 1000), limits = c(0,7000)) +
  scale_y_log10(breaks = seq(0, 7000), limits = c(0,7000)) +
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)), alpha = 0.3) +
  theme_classic() +
  theme(legend.position = "bottom")

plot_col_hump

plot_col_hump <- ggplot(data = output_all_sp_congregated, aes(x = dist, y = log(colonization_sp_mean + radiating_nodes_mean))) +
  geom_point(aes(color = GeologicalOrigin), size = 3.5) +
  xlab("Isolation") +
  ylab("Anagetic speciation events") +
  ggrepel::geom_text_repel(aes(label = ifelse(colonization_edges_mean > 500, as.character(LEVEL3_NAM), '')), 
                           nudge_y = 0.03, nudge_x = -0.10, size = 3) +
  scale_colour_manual(values = met.brewer("Pillement", 4)) +
  #scale_y_continuous(breaks = seq(0, 7000, by = 1000), limits = c(0, 7000)) +
  #scale_y_log10(breaks = c(1, 10, 100, 1000,2500,5000, 10000), limits = c(1, 7000)) +
  scale_x_continuous(breaks = seq(0, 6000, by = 500), limits = c(0,6000)) +
  geom_smooth(method = "loess", aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)), alpha = 0.3) +
  theme_classic() +
  theme(legend.position = "bottom")

plot_col_hump



plot_col_anaprop <- ggplot(data = output_all_sp_congregated,aes((TACTed_endems/(TACTed_endems+non_TACTed_endems)),(colonization_sp_mean/(colonization_sp_mean+radiating_sp_mean)))) +
  geom_point(aes(color=GeologicalOrigin ), size = 3.5) +
  ylab("Proportion anagenesis") +
  xlab("Proportion TACTed all sp") +
  ggrepel::geom_text_repel(aes(label=ifelse(colonization_edges_mean>500,as.character(LEVEL3_NAM),'')),nudge_y=0.03,nudge_x=-0.03, size=3) +
  scale_colour_manual(values = met.brewer("Pillement", 4)) +
  scale_y_continuous(limits = c(0,1,15)) +
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)), alpha = 0.3) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("Proportion of Anagenesis plotted against the number of tacted species in the island botanical country")

plot_col_anaprop


# Can we also test if there is a linear correlation ?
Tact_x_endems_model_linear <- lm(I(TACTed_endems/(TACTed_endems+non_TACTed_endems))~poly(I(colonization_sp_mean/(colonization_sp_mean+radiating_sp_mean)),1,raw = TRUE), data = output_all_sp_congregated)
summary(Tact_x_endems_model_linear)


# Can we also test if there is a hump shaped correlation between proportion tacted endems and isolation
Tact_endems_x_dist_model_hump <- lm(I(TACTed_endems/(TACTed_endems+non_TACTed_endems))~poly(dist,2,raw = TRUE), data = output_all_sp_congregated)
summary(Tact_endems_x_dist_model_hump)

# Propotion tacted X dist linear
Tact_endems_x_dist_model_linear <- lm(I(TACTed_endems/(TACTed_endems+non_TACTed_endems))~poly(dist,1,raw = TRUE), data = output_all_sp_congregated)
summary(Tact_endems_x_dist_model_linear)

# Propotion tacted X area
Tact_endems_x_area_model_linear <- lm(I(TACTed_endems/(TACTed_endems+non_TACTed_endems))~area, data = output_all_sp_congregated)
summary(Tact_endems_x_area_model_linear) # This is significant

# Proportion TACTED x fragmentation
Tact_endems_x_fragmentation_model_linear <- lm(I(TACTed_endems/(TACTed_endems+non_TACTed_endems))~nearest_neighbour_distance_border_scaled, data = output_all_sp_congregated)
summary(Tact_endems_x_fragmentation_model_linear) 

# Proportion TACTED x Geological origin
Tact_endems_x_geo_model_linear <- lm(I(TACTed_endems/(TACTed_endems+non_TACTed_endems))~GeologicalOrigin, data = output_all_sp_congregated)
summary(Tact_endems_x_geo_model_linear)

# Proportion TACTED x Elevation
Tact_endems_x_elev_model_linear <- lm(I(TACTed_endems/(TACTed_endems+non_TACTed_endems))~max_elev_30m, data = output_all_sp_congregated) # This is currently broken
summary(Tact_endems_x_elev_model_linear)  


# Is there a hump shaped relationship between area and isolation
Tact_endems_x_dist_hump <- lm(I(TACTed_endems/(TACTed_endems+non_TACTed_endems))~poly(dist,2,raw = TRUE), data = output_all_sp_congregated)
summary(Tact_endems_x_dist_hump) # this is NOT significant. 


plot_col_anaprop_endems <- ggplot(data = output_all_sp_congregated, aes((TACTed_endems/(TACTed_endems+non_TACTed_endems)), (colonization_sp_mean/(colonization_sp_mean+radiating_sp_mean)))) +
  geom_point(aes(color = GeologicalOrigin), size = 3.5) +
  ylab("Proportion anagenesis") +
  xlab("Proportion TACTed endems") +
  coord_cartesian(ylim = c(0.5, 1)) +
  ggrepel::geom_text_repel(aes(label = ifelse(colonization_edges_mean > 500, as.character(LEVEL3_NAM), '')), nudge_y = 0.03, nudge_x = -0.03, size = 3) +
  scale_colour_manual(values = met.brewer("Pillement", 4)) +
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)), alpha = 0.3) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("Proportion of Anagenesis x number of tacted endemic species")

plot_col_anaprop_endems



plot_col_area_tact <- ggplot(data = output_all_sp_congregated, aes((TACTed_endems/(TACTed_endems+non_TACTed_endems)), log10(area))) +
  geom_point(aes(color = GeologicalOrigin), size = 3.5) +
  ylab("Area") +
  xlab("Proportion TACTed endems") +
  ggrepel::geom_text_repel(aes(label = ifelse(colonization_edges_mean > 500, as.character(LEVEL3_NAM), '')), nudge_y = 0.03, nudge_x = -0.03, size = 3) +
  scale_colour_manual(values = met.brewer("Pillement", 4)) +
  geom_smooth(aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)), alpha = 0.3) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("Proportion of tacted endemics X Area of islands")

plot_col_area_tact




plot_clad_volc_height <- ggplot(data = output_all_sp_congregated[which(output_all_sp_congregated$GeologicalOrigin == "volcanic"),],aes(log10(mean_mx30_grd),(radiating_sp_mean-radiating_nodes_mean))) +
  geom_point(aes(color=GeologicalOrigin, ), size = 3.5) +
  xlab("Isolation") +
  ylab("Anagetic speciation events") +
  ggrepel::geom_text_repel(aes(label=ifelse(colonization_edges_mean>500,as.character(LEVEL3_NAM),'')),nudge_y=0.03,nudge_x=-0.10, size=3) +
  scale_colour_manual(values = met.brewer("Pillement", 4)) +
  geom_abline(intercept = 0 , slope = 8.058323902, size = 3)

plot_clad_volc_height





#Loading shapefile of botanical countries
botanical_countries <- st_read(dsn = "/home/au543206/Documents/Coryphoideae/tree_analysis/NRI/wgsrpd-master/level3", layer = "level3")

islands <- c("Borneo", "Sumatera", "Hainan", "Taiwan", "Christmas I.", "Maluku", "New Guinea", "Philippines", "Jawa", "Sulawesi", "Nansei-shoto", "Madagascar", "Gulf of Guinea Is.", "Lesser Sunda Is.", "Sri Lanka", "Mexican Pacific Is.", "Andaman Is.", "Nicobar Is.", "Vanuatu", "Bismarck Archipelago", "Solomon Is.",  "Sicilia", "Cuba", "Leeward Is.", "Bahamas", "Turks-Caicos Is.",  "Dominican Republic", "Trinidad-Tobago", "Venezuelan Antilles", "Windward Is.", "Jamaica", "Cayman Is.","Mozambique Channel I", "Hawaii" )

island_2 <- output_all_sp_test_subset_log[[1]][,1]

island_3 <- output_all_sp_test_subset_log[[1]]$LEVEL3_NAM # Finding all the islands where there is anagenetic or cladogenetic speciation. 

island_4 <- botanical_countries_islands

#Islands that are not classified according to Melanies script..
#bouvet I.
# Chagos Archipelago
# Cocos (Keeling) I.
# Crozet Is.
# Føroyar
# Galápagos
# Greenland
# Great Britain
# Howland-Baker Is.
# Heard-McDonald Is.
# Ireland
# Kazan-retto
# Laccadive Is.
# Mozambique Channel Is.
# Marcus I.
# Maldives
# Marshall Is.
# Nauru
# Prince Edward I.
# Phoenix Is.
# Réunion
# South China Sea
# South Sandwich Is.
# Tuvalu

#Selecting all botanical countries which are in 
botanical_countries$island <- 0
botanical_countries$island[which(botanical_countries$LEVEL3_NAM %in% island_3)] <- 1

# Loading Dataset with Geological Origin
island_data_origin <- read_csv(file.path(datadir,"/Islands_TDWG_AllData.txt"))

#Renaming column to be the same as in the gift database
names(island_data_origin)[1] <- "LEVEL3_COD"

#selecting only the geological origin and the name
island_data_origin_merge <- island_data_origin[,c(1,6)]

# geological origin with island counts data.
botanical_countries$GeologicalOrigin <- "Mainland"
for (i in 1:nrow(botanical_countries)){
  cod <- botanical_countries[[i,2]]
  #print(cod)
  if (cod %in% island_data_origin_merge$LEVEL3_COD){
    print(which(island_data_origin_merge$LEVEL3_COD == cod))
  botanical_countries[which(botanical_countries$LEVEL3_COD == cod),]$GeologicalOrigin <-island_data_origin_merge[which(island_data_origin_merge$LEVEL3_COD == cod),]$GeologicalOrigin
  }
}

# Setting Nunavut as Mainland:
botanical_countries[which(botanical_countries$LEVEL3_NAM == "Nunavut"),]$GeologicalOrigin <- "Mainland"

#Selecting all botanical countries which are in 
botanical_countries$island <- 0
botanical_countries$island[which(botanical_countries$GeologicalOrigin != "Mainland")] <- 1

# Adding islands with no speciation
other_isles <- c("Greenland","Great Britain","Ireland","Marshall Is.","South China Sea","Laccadive Is.","Maldives","Chagos Archipelago","Phoenix Is.","Tuvalu","Kazan-retto","Heard-McDonald Is.","Crozet Is.","South Sandwich Is.","Howland-Baker Is.","Marcus I.")
botanical_countries$oisland <- 0
botanical_countries$oisland[which(botanical_countries$LEVEL3_NAM %in% other_isles)] <- 1

# Adding species with no species recorded there 
miss_isles<- c("bouvet I.","Cocos (Keeling) I.","Mozambique Channel Is.","Nauru, Prince Edward I.")
botanical_countries$misland <- 0
botanical_countries$misland[which(botanical_countries$LEVEL3_NAM %in% miss_isles)] <- 1

# Create a column in botanical_countries with a 1 if any of island, oisland or misland contains a 1
botanical_countries$any_island <- ifelse(botanical_countries$island == "1" | botanical_countries$oisland == "1" | botanical_countries$misland == "1", 1, 0)





# Create an sf object from the polygon data
polygons <- st_as_sf(botanical_countries)
  
  # Calculate the areas of the polygons
sf::sf_use_s2(FALSE)

# for (x in 1:length(polygons$LEVEL3_NAM)) {
#   print(polygons$LEVEL3_NAM[x])
#   print(st_area(polygons[x,]))
# }
  polygons$area <- st_area(polygons)

  polygons$area <- units::set_units(polygons$area, "km^2")
  
  # Filter polygons smaller than the size threshold
  small_polygons <- polygons[polygons$area < units::set_units(20000, "km^2"), ]
  
  # Create circles around the small polygons
  circles <- st_buffer(small_polygons, dist = sqrt(1/pi))
  convex_hull <- st_convex_hull(circles) # Dont use it on 3 and 46
  convex_hull[which(convex_hull$LEVEL3_NAM == "Aleutian Is."),] <- circles[which(circles$LEVEL3_NAM == "Aleutian Is."),] # Aleutian Is.
  convex_hull[which(convex_hull$LEVEL3_NAM == "Fiji"),] <- circles[which(circles$LEVEL3_NAM == "Fiji"),] # Fiji
  
  # Plot the polygons and circles
  ggplot(data = polygons) +
    geom_sf(aes(fill = ifelse(polygons$any_island == 1, GeologicalOrigin, NA), color = "black")) +
    geom_sf(data = convex_hull[which(convex_hull$any_island ==1),], aes(fill = NA, color=polygons[which(convex_hull$any_island ==1),]$GeologicalOrigin, linetype = "solid")) +
    scale_colour_manual(values = met.brewer("Pillement", 6)) +
    scale_fill_manual(values = met.brewer("Pillement", 6)) +
    coord_sf() +
    theme_void()


  

  botanical_map <- ggplot(data = polygons) +
  scale_colour_manual(values = c("#a9845b","#697852","#738e8e", "#2b4655","#E0218A")) +
  scale_fill_manual(values = c("#a9845b","#697852","#738e8e", "#2b4655","#E0218A")) +
  geom_sf(aes(fill = ifelse(polygons$any_island == 1, GeologicalOrigin, "white"))) +
  #geom_sf(data = convex_hull[which(convex_hull$any_island == 1),], fill = NA, aes(color = GeologicalOrigin, linetype = "solid"), linewidth = 0.6) +
  geom_sf(data = convex_hull[which(convex_hull$oisland == 1),], fill = NA, aes(color = geo_origin), linetype = "solid") +
  coord_sf() +
  theme_void() +
  theme(legend.position = "none") +
  ggtitle("Island Botanical Countries")

botanical_map

# pdf(file = file.path("./Figures","Botanical_countries_map.pdf"),
#     width = 10,
#     height = 7)
# 
# plot(botanical_map)
# 
# dev.off()

# Blue are islands from which there apparently is no endemic species according to World checklist of Vascular plants.
# Mustard coloured islands have no species at all according to World checklist of vascular plants.
# Red islands does have endemic species according to the world checklist of vascular plants.

for (i in 1:length(circles[which(circles$island==1),]$island)) {
 print(i)
 print(ggplot() +
    geom_sf(data = polygons) +    
    geom_sf(data = circles[i,], fill = NA, color = "red", linetype = "solid") +
    ggtitle(circles[which(circles$island==1),]$island[i]))
  }



length(circles[which(circles$island==1),])



# Calculating how much of the total land area in the world is made up of island botanical countries.
#Now I should be able to attach all the island characteristics from the GIFT database to the database 
gift_data <- read.csv(file.path(datadir,"env_tdwg_GIFT.csv"), sep = ",")

# Calculating the proportion of botanical countries area which is made up of island botanical countries. 
all_mainland_names <- botanical_countries$LEVEL3_NAM[which(botanical_countries$any_island == 0)]
all_island_names <- botanical_countries$LEVEL3_NAM[which(botanical_countries$any_island == 1)]

Gift_area_all_islands <- gift_data$area[which(gift_data$LEVEL3_NAM %in% all_island_names)]
Gift_area_mainland <- gift_data$area[which(gift_data$LEVEL3_NAM %in% all_mainland_names)]

Island_botanical_country_area_proportion <- sum(Gift_area_all_islands)/(sum(Gift_area_all_islands)+sum(Gift_area_mainland))*100

cat("The proportion of the earths area which is made of of Island Botanical Countries is: ", Island_botanical_country_area_proportion, "%")

# Now I want to calculate how big a proportion of all species are endemic on island botanical countries.

#Creating a list of all accepted names within Coryphoideae
data <- as.data.frame(fread(file.path(datadir,"wcvp_names.csv"), sep = "|"))
dist_data <- read.csv(file.path(datadir,"wcvp_distribution.csv"), sep = "|") 

# First of lets find out how many accepted species there are according to the wcvp names
data_sp_accepted <- data[which(data$taxon_status == "Accepted" & data$taxon_rank == "Species"),]

# Maybe I should remove ferns
fern_families <- c("Equisetaceae","Psilotaceae","Ophioglossidae","Marattiaceae","Osmundaceae","Hymenophyllaceae","Matoniaceae","Dipteridaceae","Gleicheniaceae","Lygodiaceae","Schizaeaceae","Anemiaceae","Salviniaceae","Marsileaceae","Thyrsopteridaceae","Loxsomataceae","Culcitaceae","Plagiogyriaceae","Cibotiaceae","Metaxyaceae","Dicksoniaceae","Cyatheaceae","Saccolomataceae","Cystodiaceae","Lonchitidaceae","Lindsaeaceae","Pteridaceae","Dennstaedtiaceae","Cystopteridaceae","Rhachidosoraceae","Diplaziopsidaceae","Desmophlebiaceae","Hemidictyaceae","Aspleniaceae","Woodsiaceae","Onocleaceae","Blechnaceae","Athyriaceae","Thelypteridaceae","Didymochlaenaceae","Hypodematiaceae","Dryopteridaceae","Nephrolepidaceae","Lomariopsidaceae","Tectariaceae","Oleandraceae","Davalliaaceae","Polypodiaceae")

other_fams <- data_sp_accepted$family

other_fams <- unique(other_fams[which(!other_fams %in% fern_families)])

data_sp_accepted <- data_sp_accepted[which(data_sp_accepted$family %in% other_fams),]
total_sp <- length(data_sp_accepted$plant_name_id) # 345365



all_fles <- list.files(file.path("./TACT_trees"), pattern = ".rds")
dat_list <- lapply(all_fles, function(x) readRDS(file.path("./TACT_trees",x)))
output_with_island <- readRDS(file.path("./","output_with_island_data.rds"))

# Using the Tacted trees to find the number of species endemic to island botanical countries.
# Calculating the lengths of each list
output_all_sp_uniq <- lapply(dat_list, function(x) lapply(x, function(y) lapply(y, unique)))
output_all_sp_lengths <- lapply(output_all_sp_uniq, function(x) lapply(x, function(y) lapply(y, length)))

# Adding 0's to the places where there is no radiation nodes or edges
output_all_sp_lengths <- lapply(output_all_sp_lengths, function(x) lapply(x,function(y) if(length(y)==3) {append(y,c(0,0))}else{y}))

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

sum(unlist(output_as_df_all_sp[[1]]$colonization_sp))
sum(unlist(output_as_df_all_sp[[1]]$radiating_sp))

total_endemics <- sum(unlist(output_as_df_all_sp[[1]]$radiating_sp)) + sum(unlist(output_as_df_all_sp[[1]]$colonization_sp))

total_endemics/total_sp



# I want to create a figure which compares the number of species from anagenesis and cladogenesis of the coryphoid phylogeny and the all seed plants phylogeny.
head(output_with_island_congregated)

head(output_with_island_cory)

head(output_all_sp_congregated)

output_with_island_cory[,c(51,52)]

# I need to combine datasets to get both the number for the coryphoid phylogeny and the all seed plants phylogeny.
# First i create a subset of the cory phylogeny.
cory_subset_for_merge <- output_with_island_cory[,c(1,3,4,5,51,52)]

# Now I need to merge these two datasets together based on LEVEL3_NAM. The LEVEL3_NAM which are in output_all_sp_congregated but not in cory_subset_for_merge should be added as 0's.
output_all_sp_congregated <- merge(output_all_sp_congregated, cory_subset_for_merge, by = "LEVEL3_NAM", all.x = TRUE)
output_all_sp_congregated[is.na(output_all_sp_congregated)] <- 0


# Fit a linear model
fit_clad <- lm(I(`No. cladogenesis`-`No. edges`) ~ I(radiating_sp_mean-radiating_nodes_mean), data = output_all_sp_congregated)
summary(fit_clad)

# Calculate the R^2 value
r_squared_clad <- summary(fit_clad)$r.squared


p_clad_comparison <- ggplot(output_all_sp_congregated, aes((`No. cladogenesis`-`No. edges`),(radiating_sp_mean-radiating_nodes_mean))) +
  geom_point() +
  labs(title = "Number of Coryphoideae within-region speciation events compared to all plants") +
  xlab("No. Coryphoideae") +
  ylab("No. All plants") +
  geom_text(aes(label=ifelse(`No. anagenesis`+`No. cladogenesis`>1,as.character(LEVEL3_NAM),'')),hjust=0.8 ,vjust=1.8, size=3) +
  geom_smooth(method = lm) +
  annotate("text", x = Inf, y = Inf, label = paste("R^2 = ", round(r_squared_clad, 2)), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  theme_classic()

p_clad_comparison


# Fit a linear model
fit_ana <- lm(I(`No. anagenesis`+`No. edges`) ~ I(colonization_sp_mean+radiating_nodes_mean), data = output_all_sp_congregated)
summary(fit_ana)

# Calculate the R^2 value
r_squared_ana <- summary(fit_ana)$r.squared


p_ana_comparison <- ggplot(output_all_sp_congregated, aes((`No. anagenesis`+`No. edges`),(colonization_sp_mean+radiating_nodes_mean))) +
  geom_point() +
  #labs(title = "Number of Coryphoideae between-region speciation events compared to all plants") +
  xlab("No. Coryphoideae") +
  ylab("No. All plants") +
  geom_text(aes(label=ifelse(`No. anagenesis`+`No. cladogenesis`>2,as.character(LEVEL3_NAM),'')),hjust=0.8 ,vjust=1.8, size=3) +
  geom_smooth(method = lm) +
    annotate("text", x = Inf, y = Inf, label = paste("R^2 = ", round(r_squared_ana, 2)), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  theme_classic()

p_ana_comparison
