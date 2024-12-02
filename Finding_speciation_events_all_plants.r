#Loading data for plotting
all_fles <- list.files(file.path(datadir,"TACT_trees"), pattern = ".rds")
dat_list <- lapply(all_fles, function(x) readRDS(file.path(datadir,"TACT_trees",x))) # Results from Melanies Tacted trees
output_with_island <- readRDS(file.path(datadir,"output_with_island_data.rds"))
#data.table::setnames(output_with_island,'GeologicalOrigin.y','GeologicalOrigin')



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


# Loading Dataset with Geological Origin
island_data_origin <- read_csv(file.path(datadir,"/Islands_TDWG_AllData.txt"))

#Renaming column to be the same as in the gift database
names(island_data_origin)[1] <- "LEVEL3_COD"

#selecting only the geological origin and the name
island_data_origin_merge <- island_data_origin[,c(1,6)]

#Merging geological origin with island counts data.
output_all_sp <- lapply(output_as_df_all_sp, function(x) merge(x, island_data_origin_merge, by="LEVEL3_COD"))

#Loading the botanical country dataframe
botanical_countries <- st_read(dsn = file.path(datadir,"wgsrpd-master/level3"), layer = "level3")


# There are some problems with Melanies data where some of the island botanical countries have bad names

# Faulty island names
# "F\xf8royar"            "Gal\xe1pagos"          "Juan Fern\xe1ndez Is." "R\xe9union"
#Should be Faroe, Galapagos, Juan Fernandes Islands, Reunion
# LEVEL3_COD names = FOR GAL JNF REU
bad_names <- c("FOR","GAL","JNF","REU")
nnames <- botanical_countries[which(botanical_countries$LEVEL3_COD %in% bad_names),]

# Renaming fields with bad names 
# This function might crash from time to time because the order of the columns shift ??
# I have no idea why this happens but that is life... If it happens, change the number after the output_all_sp[[i]][,x] to the column named LEVEL3_NAM
for (i in 1:100) {
  problems <- output_all_sp[[i]][which(!(output_all_sp[[i]][,8] %in% botanical_countries$LEVEL3_NAM)),1]
  print(problems)
  output_all_sp[[i]][which(output_all_sp[[i]][,1] %in% problems),8] <- nnames[[1]][which(nnames$LEVEL3_COD %in% problems)]
  #output_all_sp[[i]][which(!(output_all_sp[[i]][,1] %in% botanical_countries$LEVEL3_COD)),1] <- nnames[[1]][which(nnames$LEVEL3_COD %in% problems)]
}

#Now I need to create a subset of the botanical countries which are the islands.
botanical_countries_islands <- botanical_countries[which(botanical_countries$LEVEL3_NAM %in% output_all_sp[[1]]$LEVEL3_NAM),]

botanical_countries_islands$LEVEL3_COD[which(botanical_countries_islands$LEVEL3_COD == "HAI")] <- "DOM"

sf::sf_use_s2(FALSE)


#Combining Haiti and the Dominican Republic
botanical_countries_islands %>% 
    group_by(LEVEL3_COD) %>%
    summarise(geometry = sf::st_union(geometry)) %>%
    ungroup()


# test_dom <- st_union(botanical_countries_islands[which(botanical_countries_islands$LEVEL3_COD == "DOM"),],botanical_countries_islands[which(botanical_countries_islands$LEVEL3_COD == "HAI"),])
# 
# botanical_countries_islands[which(botanical_countries_islands$LEVEL3_COD == "HAI"),]
# 
# plot(botanical_countries_islands[which(botanical_countries_islands$LEVEL3_COD == "DOM"),])


#Setting sf to NOT use S2 distances
sf::sf_use_s2(FALSE)

near_neighbour_dist <- data.frame()

print("Starting on the near neighbour distance calculation")
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
#which(!(output_all_sp[[1]][,8] %in% near_neighbour_dist$LEVEL3_NAM))
#length(output_all_sp[[1]][,1])

# Converting some columns in my dataframe from lists to be a numeric. 
for (i in 1:100){
  for (k in 3:7)
  output_all_sp[[i]][,k] <- as.numeric(unlist(output_all_sp[[i]][,k]))
}


# Making tempoary dataframes for the different variables.
tmp_colonisation_nodes <- as.data.frame(sapply(output_all_sp, function(x) unlist(x[["colonisation_nodes"]])))
tmp_colonisation_sp <- as.data.frame(sapply(output_all_sp, function(x) unlist(x[["colonization_sp"]])))
tmp_colonisation_edges <- as.data.frame(sapply(output_all_sp, function(x) unlist(x[["colonization_edges"]])))
tmp_radiating_nodes <- as.data.frame(sapply(output_all_sp, function(x) unlist(x[["radiating_nodes"]])))
tmp_radiating_sp <- as.data.frame(sapply(output_all_sp, function(x) unlist(x[["radiating_sp"]])))

# Creating a list of the tempoary dataframes
tmp_data <- list(tmp_colonisation_edges,tmp_colonisation_nodes,tmp_colonisation_sp,tmp_radiating_nodes,tmp_radiating_sp)

#Setting rownames
tmp_names <- output_all_sp[[1]][,2]
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



p1 <- ggplot(output_all_sp_congregated, aes(area, dist+1)) +
  geom_point(aes(size=colonization_edges_mean, color=GeologicalOrigin)) +
  scale_colour_manual(values = met.brewer("Pillement", 4)) +
  labs(size ="Species", title = "Between region") +
  xlab("Island Area (Km²)") +
  ylab("Distance to Mainland (Km) ") +
  scale_size(range=c(0.1,8), limits = c(1,7000)) +
  #ggrepel::geom_text_repel(aes(label=ifelse(colonization_edges_mean>1500,as.character(LEVEL3_NAM),'')),nudge_y=0.03,nudge_x=-0.10, size=4) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(100,1000,10000,100000,1000000), limits = c(100,856017.9)) +
  scale_y_continuous(breaks = c(10,50,100,250,500,1000,2000,4000,8000), limits = c(1,8000)) +
   guides(size = guide_legend("Endemic species."), colour="none") +
  coord_trans(x="log10", y="log10")

p1



p2 <- ggplot(output_all_sp_congregated, aes(area,dist+1)) +
  geom_point(aes(size=radiating_sp_mean, color=GeologicalOrigin)) +
  scale_colour_manual(values = met.brewer("Pillement", 4)) +
  labs(size ="Endemic Sp", title = "Within region") +
  xlab("Island Area (Km²)") +
  ylab("Distance to Mainland (Km)") +
  scale_size(range=c(0.1,8), limits = c(1,7000)) +
  #ggrepel::geom_text_repel(aes(label=ifelse(radiating_sp_mean>150,as.character(LEVEL3_NAM),'')),nudge_y=0.03,nudge_x=-0.10, size=4) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(100,1000,10000,100000,1000000), limits = c(100,856017.9)) +
  scale_y_continuous(breaks = c(10,50,100,250,500,1000,2000,4000,8000), limits = c(1,8000)) +
  guides(size = "none", colour=guide_legend("Geological Origin")) +
  coord_trans(x="log10", y="log10")

p2


# Creaint combined plot
all_vasc_scatter <- cowplot::plot_grid(p1,p2)
legend_bottom <- get_legend(p1+theme(legend.position = "bottom", legend.text=element_text(size=12)))
all_vasc_scatter <- cowplot::plot_grid(all_vasc_scatter,legend_bottom, ncol=1,nrow = 2, rel_heights = c(1,.1))

# now add the title
title_all_vasc <- ggdraw() + 
  draw_label(
    "All vascular plants",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 40)
  )
all_vasc_scatter <- plot_grid(
  title_all_vasc, all_vasc_scatter,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

all_vasc_scatter


# Creating a boxplot to complement the scatterplot
# I want the boxplot to show the number of species resulting from anagenesis and cladogenesis for each origin of island botanical country.
# Then I would also like to have a boxplot showing the isolation, area fragmentation and max elevation for each island type aswell.

# Boxplot for Anagenetic Speciation Events
p1_boxplot <- ggplot(output_all_sp_congregated, aes(x = GeologicalOrigin, y = log(colonization_edges_mean+radiating_nodes_mean))) +
  geom_boxplot(aes(fill = GeologicalOrigin)) +
  scale_fill_manual(values = met.brewer("Pillement", 4)) +
  labs(title = "Between-region speciation") +
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.title.x = element_blank()   # Remove x-axis title
  ) +
  theme(legend.position = "none")

p1_boxplot

# Boxplot for Cladogenetic Speciation Events
p2_boxplot <- ggplot(output_all_sp_congregated, aes(x = GeologicalOrigin, y = log(radiating_sp_mean-radiating_nodes_mean))) +
  geom_boxplot(aes(fill = GeologicalOrigin)) +
  scale_fill_manual(values = met.brewer("Pillement", 4)) +
  labs(title = "Within-region speciation") +
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.title.x = element_blank()   # Remove x-axis title
  ) +
  theme(legend.position = "none")

p2_boxplot

# Use cowplot to combine the two boxplots
all_vasc_boxplot <- cowplot::plot_grid(p1_boxplot,p2_boxplot)
all_vasc_boxplot

# now add the title
title_all <- cowplot::ggdraw() + 
  draw_label(
    "All seed plants",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 40)
  )
all_vasc_boxplot <- plot_grid(
  title_all, all_vasc_boxplot,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)
all_vasc_boxplot



# Boxplot for area
p3_boxplot <- ggplot(output_all_sp_congregated, aes(x = GeologicalOrigin, y = log(area))) +
  geom_boxplot(aes(fill = GeologicalOrigin)) +
  scale_fill_manual(values = met.brewer("Pillement", 4)) +
  labs(title = "Area") +
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(
  axis.text.x = element_blank(),  # Remove x-axis text
  axis.ticks.x = element_blank(),  # Remove x-axis ticks
  axis.title.x = element_blank()   # Remove x-axis title
  ) +
  theme(legend.position = "none")

p3_boxplot

# Boxplot for Isolation
p4_boxplot <- ggplot(output_all_sp_congregated, aes(x = GeologicalOrigin, y = log(dist))) +
  geom_boxplot(aes(fill = GeologicalOrigin)) +
  scale_fill_manual(values = met.brewer("Pillement", 4)) +
  labs(title = "Isolation") +
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(
  axis.text.x = element_blank(),  # Remove x-axis text
  axis.ticks.x = element_blank(),  # Remove x-axis ticks
  axis.title.x = element_blank()   # Remove x-axis title
  ) +
  theme(legend.position = "none")

p4_boxplot


# Boxplot for Fragmentation
p5_boxplot <- ggplot(output_all_sp_congregated, aes(x = GeologicalOrigin, y = log(nearest_neighbour_distance_border_scaled))) +
  geom_boxplot(aes(fill = GeologicalOrigin)) +
  scale_fill_manual(values = met.brewer("Pillement", 4)) +
  labs(title = "Fragmentation") +
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(
  axis.text.x = element_blank(),  # Remove x-axis text
  axis.ticks.x = element_blank(),  # Remove x-axis ticks
  axis.title.x = element_blank()   # Remove x-axis title
  ) +
  theme(legend.position = "none")

p5_boxplot


# Boxplot for Maximum elevation
p6_boxplot <- ggplot(output_all_sp_test_subset_log[[1]], aes(x = GeologicalOrigin, y = log(max30_elev))) +
  geom_boxplot(aes(fill = GeologicalOrigin)) +
  scale_fill_manual(values = met.brewer("Pillement", 4)) +
  labs(title = "Maximum elevation") +
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.title.x = element_blank()   # Remove x-axis title
  ) +
  theme(legend.position = "none")

p6_boxplot


# Create a cowplott of p3,p4,p5 and p6 in one row
boxplot_phys_vars <- cowplot::plot_grid(p3_boxplot,p4_boxplot,p5_boxplot,p6_boxplot, nrow = 1)

# Adding the legend to the boxplot
legend_bottom <- get_legend(p4_boxplot+theme(legend.position = "bottom", legend.text=element_text(size=12)))
boxplot_phys_vars_with_legend <- cowplot::plot_grid(boxplot_phys_vars,legend_bottom, ncol=1,nrow = 2, rel_heights = c(1,.1))

# now add the title
title_phys <- ggdraw() + 
  draw_label(
    "Island physical characteristics",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 40)
  )
boxplot_phys_vars_with_legend <- plot_grid(
  title_phys, boxplot_phys_vars_with_legend,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

boxplot_phys_vars_with_legend




output_with_island <- merge(output_with_island, output_all_sp_congregated[,c(1,57:71)], by="LEVEL3_COD",all = FALSE, sort = FALSE, no.dups = TRUE)

p3 <- ggplot(output_with_island, aes(log10(`No. cladogenesis`+`No. anagenesis`), log10(radiating_sp_mean+colonization_sp_mean)))
p3 + geom_point(aes(color=GeologicalOrigin)) +
scale_colour_manual(values = met.brewer("Pillement", 4)) +
labs(title = "Number of endemic Coryphoideae species compared to endemics of all plants") +
  xlab("Endemic Coryphoids") +
  ylab("Endemic Plants") +
  geom_text(aes(label=ifelse(`No. anagenesis`>2,as.character(LEVEL3_NAM),'')),hjust=0.8 ,vjust=1.8, size=3) +
  geom_smooth(method = lm) +
  theme_classic()



p4 <- ggplot(output_with_island, aes(`No. cladogenesis`, radiating_sp_mean))
p4 + geom_point(aes(color=GeologicalOrigin)) +
scale_colour_manual(values = met.brewer("Pillement", 4)) +
labs(title = "Number of species resulting from cladogenesis: \n Coryphoideae compared to all plants") +
  xlab("Endemics from Radiation: Coryphoids") +
  ylab("Endemics from Radiation: all plants") +
  geom_text(aes(label=ifelse(`No. anagenesis`+`No. cladogenesis`>2,as.character(LEVEL3_NAM),'')),hjust=0.8 ,vjust=1.8, size=3) +
  geom_smooth(method = lm) +
  theme_classic()



p5 <- ggplot(output_with_island, aes(`No. anagenesis`, colonization_sp_mean))
p5 + geom_point(aes(color=GeologicalOrigin)) +
scale_colour_manual(values = met.brewer("Pillement", 4)) +
labs(title = "Number of endemic species with no sister species found on the island: \n Coryphoideae compared to all plants") +
  xlab("Endemic Coryphoids") +
  ylab("Endemic Plants") +
  ggrepel::geom_text_repel(aes(label=ifelse(`No. anagenesis`>1,as.character(LEVEL3_NAM),'')),size=3) +
  geom_smooth(method = lm) +
  theme_classic()



p6 <- ggplot(output_with_island, aes((`No. cladogenesis`+`No. anagenesis`)/Total_sp_coryphoids, (radiating_sp_mean+colonization_sp_mean)/Total_sp))
p6 + geom_point(aes(color=GeologicalOrigin)) +
scale_colour_manual(values = met.brewer("Pillement", 4)) +
labs(title = "Proportion of endemic Coryphoideae species compared to \n Proportion of endemics of all plants") +
  xlab("Proportion Endemic Coryphoids") +
  ylab("Proportion Endemic Plants") +
  geom_text_repel(aes(label=ifelse(colonization_sp_mean>10 & `No. cladogenesis`+`No. anagenesis`>=1 ,as.character(LEVEL3_NAM),'')),size=3) +
  geom_smooth(method = lm) +
  theme_classic()






p7 <- ggplot(output_with_island, aes((`No. cladogenesis`)/Total_sp_coryphoids, (radiating_sp_mean)/Total_sp))
p7 + geom_point(aes(color=GeologicalOrigin)) +
labs(title = "Proportion of endemic Coryphoideae species compared to \n Proportion of endemics of all plants") +
  scale_colour_manual(values = met.brewer("Pillement", 4)) +
  xlab("Proportion Endemic Coryphoids") +
  ylab("Proportion Endemic Plants") +
  geom_text(aes(label=ifelse(radiating_sp_mean>2 & `No. cladogenesis`>=1 ,as.character(LEVEL3_NAM),'')),hjust=0.8 ,vjust=1.8, size=3) +
  geom_smooth(method = lm) +
  theme_classic()


# Cladogenesis model on mean number of radiating sp
glm_clad_pois_all_sp <- glm(radiating_sp_mean-radiating_nodes_mean~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled)+GeologicalOrigin,data=output_all_sp_congregated[which(output_all_sp_congregated$radiating_sp_mean >= 1),],family = quasipoisson(link='log'))

# Investigating the dispersion of my dataset
dp_glm_clad_pois_all_sp = sum(residuals(glm_clad_pois_all_sp,type ="pearson")^2)/glm_clad_pois_all_sp$df.residual # Calculating dispersion

summary(glm_clad_pois_all_sp, dispersion = dp_glm_clad_pois_all_sp) # Only area is of significant importance and they both have a positive effect on the proportion of endemics



#Changing NA to 0 in the dataframes as island botanical countries consisting of a single island should have a fragmentation (near neighbour distance) of 0 and not NA
output_all_sp_test <- lapply(output_all_sp, function(x) x <- x[!is.na(x$nearest_neighbour_distance_border_scaled),])

#  This does tho create some problems because I use LOG10 in the fragmentation scaling and the Log10 of 0 is infinite
output_all_sp_test <-lapply(output_all_sp_test, function(x) x[which(!x$dist == 0),])

# First lets find the number of species present on each island botanical country
# Trying to load WCVP distribution data to see if it is better
dist_data <- read.csv(file.path(datadir,"checklist_distribution.txt"), sep = "|") 
dist_data <- dist_data[which(dist_data$introduced != 1 & dist_data$location_doubtful != 1 ),]

dist_data[which(dist_data$area == "Amsterdam-St.Paul Is"),7] <- "Amsterdam-St.Paul Is."
dist_data[which(dist_data$area == "Central American Pac"),7] <- "C. American Pacific Is."
dist_data[which(dist_data$area == "Leeward Is."),7] <- "Leeward Is. AB Ant"
dist_data[which(dist_data$area == "Mozambique Channel I"),7] <- "Mozambique Channel Is."
dist_data[which(dist_data$area == "Marion-Prince Edward"),7] <- "Marion-Prince Edward Is."


#So now I need to look through the islands that I have and count the total number of species on the islands
total_island_sp <- data.frame()



for (i in 1:length(botanical_countries$LEVEL3_NAM)){
  total_island_sp[i,1] <- botanical_countries$LEVEL3_NAM[i]
  total_island_sp[i,2] <- length(which(dist_data$area == botanical_countries$LEVEL3_NAM[i]))
}

names(total_island_sp) <- c("LEVEL3_NAM","Total_sp")

output_all_sp_test <- lapply(output_all_sp_test, function(x) merge(x, total_island_sp, by="LEVEL3_NAM"))

saveRDS(output_all_sp_test, file.path(datadir,"Output_all_sp_test.rds"))

# Creating a for loop which fits the model using all 100 different datasets 
# -------- Cladogenesis ----------------
glm_clad_pois_all_sp_list_tidy <- list()
glm_clad_pois_all_sp_list_models <- list()
for (i in 1:100){
  glm_clad_pois_all_sp <- glm(radiating_sp-radiating_nodes~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled)+relevel(factor(GeologicalOrigin), ref = "volcanic"),data=output_all_sp_test[[i]][which(output_all_sp_test[[i]]$radiating_sp >= 1),],family = quasipoisson(link='log'))
  glm_clad_pois_all_sp_list_models[[i]] <- glm_clad_pois_all_sp
  glm_clad_pois_all_sp_list_tidy[[i]] <- broom::tidy(glm_clad_pois_all_sp)
}


#Now lets find the confidence intervals for these models.
fit_clad_all_95_list <- list()
for (i in 1:100){
  fit_clad_all_95 <- confint(glm_clad_pois_all_sp_list_models[[i]], level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")
  fit_clad_all_95_list[[i]] <- fit_clad_all_95
  }

# Combining the 95 intervals on all the models
results_clad_all_list <- list()
for (i in 1:100){
  results_clad_all_list[[i]] <- dplyr::bind_cols(glm_clad_pois_all_sp_list_tidy[[i]],fit_clad_all_95_list[[i]]) 
}


#Trying to plot the cladogenetic coefficients table
# for (i in 1:100) {
#   clad_plot <- ggplot(results_clad_all_list[[i]], 
#        aes(x = term, y = estimate)) +
#         geom_hline(yintercept = 0, 
#                    colour = gray(1/2), lty = 2) +
#         geom_point(aes(x = term, 
#                     y = estimate)) + 
#         geom_linerange(aes(x = term, 
#                      ymin = conf.low_95,
#                      ymax = conf.high_95),
#                    lwd = 1/2) + 
#         ggtitle("Cladogenesis coefficient table all sp") +
#         coord_flip()
#   print(clad_plot)
# }



######### This below is different models I have not converted to the list method yet. 
# Removing all islands which has no fragmentation index Ie. island botanical countries that only consist of a single island
output_all_sp_congregated_no_single_islands <- output_all_sp_congregated[!is.na(output_all_sp_congregated$nearest_neighbour_distance_border_scaled),]

# Removing Island botanical countries with an isolation index of 0 (Aegean islands in the mediterranean)
output_all_sp_congregated_no_single_islands <- output_all_sp_congregated_no_single_islands[which(!output_all_sp_congregated_no_single_islands$dist == 0),]

#Fitting the model on the mean values across all trees
glm_ana_pois_all_sp <- glm(colonization_sp_mean+radiating_nodes_mean ~
                      log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled)+GeologicalOrigin,
                    data=output_all_sp_congregated_no_single_islands,
                    family = quasipoisson(link = "log"))

# Calculating dispersion
dp_glm_ana_pois_all_sp = sum(residuals(glm_ana_pois_all_sp,type ="pearson")^2)/glm_ana_pois_all_sp$df.residual

#Results
summary(glm_ana_pois_all_sp,  dispersion = dp_glm_ana_pois_all_sp) # Area and isolation are significant


# Creating a for loop which fits the model using all 100 different datasets 
# -------- Anagenesis ----------------
glm_ana_pois_all_sp_list_tidy <- list()
glm_ana_pois_all_sp_list_models <- list()
for (i in 1:100){
  glm_ana_all_sp_loop <- glm(colonization_sp+radiating_nodes~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled)+relevel(factor(GeologicalOrigin), ref = "volcanic"),data=output_all_sp_test[[i]],family = quasipoisson(link='log'))
  glm_ana_pois_all_sp_list_models[[i]] <- glm_ana_all_sp_loop
  glm_ana_pois_all_sp_list_tidy[[i]] <- broom::tidy(glm_ana_all_sp_loop)
}




glm_ana_pois_all_sp_list_tidy_test <- list()
glm_ana_pois_all_sp_list_models_test <- list()
for (i in 1:100){
  glm_ana_all_sp_loop_test <- glm(log10(colonization_sp+radiating_nodes)~log10(mean_mx30_grd)+log10(area)+log10(dist^2+1)+log10(nearest_neighbour_distance_border_scaled+1)+relevel(factor(GeologicalOrigin), ref = "volcanic"),data=output_all_sp[[i]],family = quasipoisson(link='log'))
  glm_ana_pois_all_sp_list_models_test[[i]] <- glm_ana_all_sp_loop_test
  glm_ana_pois_all_sp_list_tidy_test[[i]] <- broom::tidy(glm_ana_all_sp_loop_test)
}

for (i in 1:100){
print(sum(residuals(glm_ana_pois_all_sp_list_models_test[[i]],type ="pearson")^2)/glm_ana_pois_all_sp_list_models_test[[i]]$df.residual)
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


#Trying to plot the anagenetic coefficients table
# for (i in 1:100) {
#   ana_plot <- ggplot(results_ana_all_list[[i]], 
#        aes(x = term, y = estimate)) +
#         geom_hline(yintercept = 0, 
#                    colour = gray(1/2), lty = 2) +
#         geom_point(aes(x = term, 
#                     y = estimate)) + 
#         geom_linerange(aes(x = term, 
#                      ymin = conf.low_95,
#                      ymax = conf.high_95),
#                    lwd = 1/2) + 
#         ggtitle(paste0("Anagenesis coefficient table all sp",i, sep = "  ")) +
#         coord_flip()
#   print(ana_plot)
# }

# Selecting only volcanic islands
output_all_sp_test_volcanic <- lapply(output_all_sp_test, function(x) x[which(x$GeologicalOrigin == "volcanic"),])

# Fitting cladogenetic model for volcanic islands.
glm_clad_pois_all_sp_volcanic_list_tidy <- list()
glm_clad_pois_all_sp_volcanic_list_models <- list()
for (i in 1:100){
  glm_clad_pois_all_sp_volcanic <- glm(radiating_sp-radiating_nodes~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled),data=output_all_sp_test_volcanic[[i]][which(output_all_sp_test_volcanic[[i]]$radiating_sp >= 1),],family = quasipoisson(link='log'))
  glm_clad_pois_all_sp_volcanic_list_models[[i]] <- glm_clad_pois_all_sp_volcanic
  glm_clad_pois_all_sp_volcanic_list_tidy[[i]] <- broom::tidy(glm_clad_pois_all_sp_volcanic)
}

#-------------- preparing the coefficient plot for these models --------------
#Now lets find the confidence intervals for these models.
fit_clad_all_volcanic_95_list <- list()
for (i in 1:100){
  fit_clad_all_volcanic_95 <- confint(glm_clad_pois_all_sp_volcanic_list_models[[i]], level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")
  fit_clad_all_volcanic_95_list[[i]] <- fit_clad_all_volcanic_95
  }

# Combining the 95 intervals on all the models
results_clad_all_volcanic_list <- list()
for (i in 1:100){
  results_clad_all_volcanic_list[[i]] <- dplyr::bind_cols(glm_clad_pois_all_sp_volcanic_list_tidy[[i]],fit_clad_all_volcanic_95_list[[i]]) 
}


#Trying to plot the cladogenetic coefficients table
# for (i in 1:100) {
#   clad_volcanic_plot <- ggplot(results_clad_all_volcanic_list[[i]], 
#        aes(x = term, y = estimate)) +
#         geom_hline(yintercept = 0, 
#                    colour = gray(1/2), lty = 2) +
#         geom_point(aes(x = term, 
#                     y = estimate)) + 
#         geom_linerange(aes(x = term, 
#                      ymin = conf.low_95,
#                      ymax = conf.high_95),
#                    lwd = 1/2) + 
#         ggtitle(paste0("Cladogenetic coefficient table all sp /n only volcanic",i, sep = "  ")) +
#         coord_flip()
#   print(clad_volcanic_plot)
# }


# Before I continue down this strange road of finding significant values, let us just explore some plotting options.
# I will try to create a frankenstein of a violin plot where I will add 3 violin plots to the same coefficient one for the mean 
glm_coef_hight_clad_volc <- as.data.frame(t(sapply(results_clad_all_volcanic_list, function(x) unlist(x[2,]))))
glm_coef_area_clad_volc <- as.data.frame(t(sapply(results_clad_all_volcanic_list, function(x) unlist(x[3,]))))
glm_coef_dist_clad_volc <- as.data.frame(t(sapply(results_clad_all_volcanic_list, function(x) unlist(x[4,]))))
glm_coef_fragmentation_clad_volc <- as.data.frame(t(sapply(results_clad_all_volcanic_list, function(x) unlist(x[5,]))))

# Stacking the dataframes
glm_coefs_clad_volc <- rbind(glm_coef_area_clad_volc, glm_coef_hight_clad_volc, glm_coef_dist_clad_volc, glm_coef_fragmentation_clad_volc)
glm_coefs_list_clad_volc <- list(glm_coef_hight_clad_volc,glm_coef_area_clad_volc,glm_coef_dist_clad_volc,glm_coef_fragmentation_clad_volc)

# #
# glm_coef_area[,2:7] <- as.numeric(unlist(glm_coef_area[,2:7]))
# as.numeric(glm_coef_area$estimate)


#I have to create a dataframe with the means of all the 
glm_coefs_mean_clad_volc <- data.frame()
for (i in 1:4) {
  temp_means <- vector()
  temp_means <- apply(glm_coefs_list_clad_volc[[i]], 2, function(x) mean(as.numeric(unlist(x))))
  temp_means["term"] <- glm_coefs_list_clad_volc[[i]][1,1]
  print(temp_means)
  glm_coefs_mean_clad_volc <- rbind(glm_coefs_mean_clad_volc,temp_means)
}
names(glm_coefs_mean_clad_volc) <- names(glm_coef_dist_clad_volc)

# making a loop to convert character list rows to numeric 
for(i in 2:7){
  glm_coefs_mean_clad_volc[,i] <-as.numeric(unlist(glm_coefs_mean_clad_volc[,i]))
}

# Changeing some names
glm_coefs_mean_names_clad_volc <- c("Max elevation","Area","Isolation","Fragmentation")
glm_coefs_mean_names_clad_volc_old <- c("log10(mean_mx30_grd)","log10(area)","log10(dist)","log10(nearest_neighbour_distance_border_scaled)")
glm_coefs_names_df_clad_volc <- cbind(glm_coefs_mean_names_clad_volc, glm_coefs_mean_names_clad_volc_old)
glm_coefs_mean_clad_volc[,1] <- glm_coefs_mean_names_clad_volc

for (i in 1:4){
  glm_coefs_clad_volc[which(glm_coefs_clad_volc$term==glm_coefs_names_df_clad_volc[i,2]),1] <- glm_coefs_names_df_clad_volc[i,1]
}


# Creating a plot using violins to show the distribution of the mean and the 95 confidence interval
gvil_volcanic_clad <- ggplot(glm_coefs_clad_volc, aes(color=term, fill = term)) +
  scale_colour_manual(values = met.brewer("Hiroshige", 4)) +
  scale_fill_manual(values = met.brewer("Hiroshige", 4)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 10) +
  geom_violin(aes(x = term, y = as.numeric(estimate),)) +
  geom_violin(aes(x = term, y = as.numeric(conf.high_95),)) +
  geom_violin(aes(x = term, y = as.numeric(conf.low_95),)) +
  geom_crossbar(data = glm_coefs_mean_clad_volc, size=0.1, alpha=0.5, aes(x= term, y = estimate,
                    ymin=conf.low_95,
                    ymax=conf.high_95)) +
  xlab("Coefficients") +
  ylab("Estimate") +
  theme_classic() +
  theme(legend.position = "none") +
  
  coord_flip()

gvil_volcanic_clad

# Fitting anagenetic model for volcanic islands.
glm_ana_pois_all_sp_volcanic_list_tidy <- list()
glm_ana_pois_all_sp_volcanic_list_models <- list()
for (i in 1:100){
  glm_ana_pois_all_sp_volcanic <- glm(colonization_sp+radiating_nodes~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled),data=output_all_sp_test_volcanic[[i]],family = quasipoisson(link='log'))
  glm_ana_pois_all_sp_volcanic_list_models[[i]] <- glm_ana_pois_all_sp_volcanic
  glm_ana_pois_all_sp_volcanic_list_tidy[[i]] <- broom::tidy(glm_ana_pois_all_sp_volcanic)
}

#-------------- preparing the coefficient plot for these models --------------
#Now lets find the confidence intervals for these models.
fit_ana_all_volcanic_95_list <- list()
for (i in 1:100){
  fit_ana_all_volcanic_95 <- confint(glm_ana_pois_all_sp_volcanic_list_models[[i]], level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")
  fit_ana_all_volcanic_95_list[[i]] <- fit_ana_all_volcanic_95
  }

# Combining the 95 intervals on all the models
results_ana_all_volcanic_list <- list()
for (i in 1:100){
  results_ana_all_volcanic_list[[i]] <- dplyr::bind_cols(glm_ana_pois_all_sp_volcanic_list_tidy[[i]],fit_ana_all_volcanic_95_list[[i]]) 
}

#Trying to plot the Anagenetic coefficients table
# for (i in 1:100) {
#   ana_volcanic_plot <- ggplot(results_ana_all_volcanic_list[[i]],
#        aes(x = term, y = estimate)) +
#         geom_hline(yintercept = 0,
#                    colour = gray(1/2), lty = 2) +
#         geom_point(aes(x = term,
#                     y = estimate)) +
#         geom_linerange(aes(x = term,
#                      ymin = conf.low_95,
#                      ymax = conf.high_95),
#                    lwd = 1/2) +
#         ggtitle(paste0("Anagenetic coefficient table all sp /n only volcanic",i, sep = "  ")) +
#         coord_flip()
#   print(ana_volcanic_plot)
# }






###########################################################################################
#---------------------------------- Violin and bars plot ---------------------------------#
###########################################################################################

#Creating dataframes with all estimates for each of the coefficients
glm_coef_hight_ana_volcanic <- as.data.frame(t(sapply(results_ana_all_volcanic_list, function(x) unlist(x[2,]))))
glm_coef_area_ana_volcanic <- as.data.frame(t(sapply(results_ana_all_volcanic_list, function(x) unlist(x[3,]))))
glm_coef_dist_ana_volcanic <- as.data.frame(t(sapply(results_ana_all_volcanic_list, function(x) unlist(x[4,]))))
glm_coef_fragmentation_ana_volcanic <- as.data.frame(t(sapply(results_ana_all_volcanic_list, function(x) unlist(x[5,]))))

# Stacking the dataframes
glm_coefs_ana_volc <- rbind(glm_coef_area_ana_volcanic, glm_coef_hight_ana_volcanic, glm_coef_dist_ana_volcanic, glm_coef_fragmentation_ana_volcanic)
glm_coefs_list_ana_volc <- list(glm_coef_hight_ana_volcanic,glm_coef_area_ana_volcanic,glm_coef_dist_ana_volcanic,glm_coef_fragmentation_ana_volcanic)


#I have to create a dataframe with the means of all the 
glm_coefs_mean_ana_volc <- data.frame()
for (i in 1:4) {
  temp_means <- vector()
  temp_means <- apply(glm_coefs_list_ana_volc[[i]], 2, function(x) mean(as.numeric(unlist(x))))
  temp_means["term"] <- glm_coefs_list_ana_volc[[i]][1,1]
  print(temp_means)
  glm_coefs_mean_ana_volc <- rbind(glm_coefs_mean_ana_volc,temp_means)
}
names(glm_coefs_mean_ana_volc) <- names(glm_coef_dist_ana_volcanic)

# making a loop to convert character list rows to numeric 
for(i in 2:7){
  glm_coefs_mean_ana_volc[,i] <-as.numeric(unlist(glm_coefs_mean_ana_volc[,i]))
}

# Changeing some names
glm_coefs_mean_names_ana_volc <- c("Max elevation","Area","Isolation","Fragmentation")
glm_coefs_mean_names_ana_volc_old <- c("log10(mean_mx30_grd)","log10(area)","log10(dist)","log10(nearest_neighbour_distance_border_scaled)")
glm_coefs_names_df_ana_volc <- cbind(glm_coefs_mean_names_ana_volc, glm_coefs_mean_names_ana_volc_old)

# Changing names in dataframe of means
glm_coefs_mean_ana_volc[,1] <- glm_coefs_mean_names_ana_volc

# Changing names in dataframe of estimates
for (i in 1:4){
  glm_coefs_ana_volc[which(glm_coefs_ana_volc$term==glm_coefs_names_df_ana_volc[i,2]),1] <- glm_coefs_names_df_ana_volc[i,1]
}


# Creating a plot using violins to show the distribution of the mean and the 95 confidence interval
gvil_volcanic_ana <- ggplot(glm_coefs_ana_volc, aes(color=term, fill = term)) +
  scale_colour_manual(values = met.brewer("Hiroshige", 4)) +
  scale_fill_manual(values = met.brewer("Hiroshige", 4)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 10) +
  geom_violin(aes(x = term, y = as.numeric(estimate),)) +
  geom_violin(aes(x = term, y = as.numeric(conf.high_95),)) +
  geom_violin(aes(x = term, y = as.numeric(conf.low_95),)) +
  geom_crossbar(data = glm_coefs_mean_ana_volc, size=0.1, alpha=0.5, aes(x= term, y = estimate,
                    ymin=conf.low_95,
                    ymax=conf.high_95)) +
  xlab("Coefficients") +
  ylab("Estimate") +
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()

gvil_volcanic_ana






# Fitting endemics model for volcanic islands.
glm_endems_all_sp_volcanic_list_tidy <- list()
glm_endems_all_sp_volcanic_list_models <- list()
for (i in 1:100){
  glm_endems_all_sp_volcanic <- glm(colonization_sp+radiating_sp~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled),data=output_all_sp_test_volcanic[[i]],family = quasipoisson(link='log'))
  glm_endems_all_sp_volcanic_list_models[[i]] <- glm_endems_all_sp_volcanic
  glm_endems_all_sp_volcanic_list_tidy[[i]] <- broom::tidy(glm_endems_all_sp_volcanic)
}

#-------------- preparing the coefficient plot for these models --------------
#Now lets find the confidence intervals for these models.
fit_endems_all_volcanic_95_list <- list()
for (i in 1:100){
  fit_endems_all_volcanic_95 <- confint(glm_endems_all_sp_volcanic_list_models[[i]], level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")
  fit_endems_all_volcanic_95_list[[i]] <- fit_endems_all_volcanic_95
  }

# Combining the 95 intervals on all the models
results_endems_all_volcanic_list <- list()
for (i in 1:100){
  results_endems_all_volcanic_list[[i]] <- dplyr::bind_cols(glm_endems_all_sp_volcanic_list_tidy[[i]],fit_endems_all_volcanic_95_list[[i]]) 
}

# Endems
glm_coef_hight_endems_volcanic <- as.data.frame(t(sapply(results_endems_all_volcanic_list, function(x) unlist(x[2,]))))
glm_coef_area_endems_volcanic <- as.data.frame(t(sapply(results_endems_all_volcanic_list, function(x) unlist(x[3,]))))
glm_coef_dist_endems_volcanic <- as.data.frame(t(sapply(results_endems_all_volcanic_list, function(x) unlist(x[4,]))))
glm_coef_fragmentation_endems_volcanic <- as.data.frame(t(sapply(results_endems_all_volcanic_list, function(x) unlist(x[5,]))))

# Stacking the dataframes
glm_coefs_endems_volc <- rbind(glm_coef_hight_endems_volcanic,
                            glm_coef_area_endems_volcanic,
                            glm_coef_dist_endems_volcanic,
                            glm_coef_fragmentation_endems_volcanic)

# List of dataframes
glm_coefs_list_endems_volc <- list(glm_coef_hight_endems_volcanic,
                            glm_coef_area_endems_volcanic,
                            glm_coef_dist_endems_volcanic,
                            glm_coef_fragmentation_endems_volcanic)

# Endemic means
glm_coefs_means_endems_volc <- data.frame()
for (i in 1:4) {
  temp_means <- vector()
  temp_means <- apply(glm_coefs_list_endems_volc[[i]], 2, function(x) mean(as.numeric(unlist(x))))
  temp_means["term"] <- glm_coefs_list_endems_volc[[i]][1,1]
  glm_coefs_means_endems_volc <- rbind(glm_coefs_means_endems_volc,temp_means)
}
names(glm_coefs_means_endems_volc) <- names(glm_coef_dist_endems_volcanic)





# Fitting endemic proportopm model for volcanic islands.
glm_prop_endems_all_sp_volcanic_list_tidy <- list()
glm_prop_endems_all_sp_volcanic_list_models <- list()
for (i in 1:100){
  glm_prop_endems_all_sp_volcanic <- glm((colonization_sp+radiating_sp)/Total_sp~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled),data=output_all_sp_test_volcanic[[i]],family = quasibinomial(link='log'))
  glm_prop_endems_all_sp_volcanic_list_models[[i]] <- glm_prop_endems_all_sp_volcanic
  glm_prop_endems_all_sp_volcanic_list_tidy[[i]] <- broom::tidy(glm_prop_endems_all_sp_volcanic)
}

#-------------- preparing the coefficient plot for these models --------------
#Now lets find the confidence intervals for these models.
fit_prop_endems_all_volcanic_95_list <- list()
for (i in 1:100){
  fit_prop_endems_all_volcanic_95 <- confint(glm_prop_endems_all_sp_volcanic_list_models[[i]], level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")
  fit_prop_endems_all_volcanic_95_list[[i]] <- fit_prop_endems_all_volcanic_95
  }

# Combining the 95 intervals on all the models
results_prop_endems_all_volcanic_list <- list()
for (i in 1:100){
  results_prop_endems_all_volcanic_list[[i]] <- dplyr::bind_cols(glm_prop_endems_all_sp_volcanic_list_tidy[[i]],fit_prop_endems_all_volcanic_95_list[[i]]) 
}

# Endemc proportion
glm_coef_hight_prop_endems_volcanic <- as.data.frame(t(sapply(results_prop_endems_all_volcanic_list, function(x) unlist(x[2,]))))
glm_coef_area_prop_endems_volcanic <- as.data.frame(t(sapply(results_prop_endems_all_volcanic_list, function(x) unlist(x[3,]))))
glm_coef_dist_prop_endems_volcanic <- as.data.frame(t(sapply(results_prop_endems_all_volcanic_list, function(x) unlist(x[4,]))))
glm_coef_fragmentation_prop_endems_volcanic <- as.data.frame(t(sapply(results_prop_endems_all_volcanic_list, function(x) unlist(x[5,]))))


# Stacking the dataframes
glm_coefs_prop_endems_volc <- rbind(glm_coef_hight_prop_endems_volcanic,
                            glm_coef_area_prop_endems_volcanic,
                            glm_coef_dist_prop_endems_volcanic,
                            glm_coef_fragmentation_prop_endems_volcanic)

# List of dataframes
glm_coefs_list_prop_endems_volc <- list(glm_coef_hight_prop_endems_volcanic,
                            glm_coef_area_prop_endems_volcanic,
                            glm_coef_dist_prop_endems_volcanic,
                            glm_coef_fragmentation_prop_endems_volcanic)

# Endemic means
glm_coefs_means_prop_endems_volc <- data.frame()
for (i in 1:4) {
  temp_means <- vector()
  temp_means <- apply(glm_coefs_list_prop_endems_volc[[i]], 2, function(x) mean(as.numeric(unlist(x))))
  temp_means["term"] <- glm_coefs_list_prop_endems_volc[[i]][1,1]
  glm_coefs_means_prop_endems_volc <- rbind(glm_coefs_means_prop_endems_volc,temp_means)
}
names(glm_coefs_means_prop_endems_volc) <- names(glm_coef_dist_prop_endems_volcanic)




# Fitting proportion endemics from regional allopatric model for volcanic islands.
glm_regional_allo_prop_all_sp_volcanic_list_tidy <- list()
glm_regional_allo_prop_all_sp_volcanic_list_models <- list()
for (i in 1:100){
  glm_regional_allo_prop_all_sp_volcanic <- glm((colonization_sp+radiating_nodes)/(colonization_sp+radiating_sp)~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled),data=output_all_sp_test_volcanic[[i]],family = quasibinomial)
  glm_regional_allo_prop_all_sp_volcanic_list_models[[i]] <- glm_regional_allo_prop_all_sp_volcanic
  glm_regional_allo_prop_all_sp_volcanic_list_tidy[[i]] <- broom::tidy(glm_regional_allo_prop_all_sp_volcanic)
}

#-------------- preparing the coefficient plot for these models --------------
#Now lets find the confidence intervals for these models.
fit_regional_allo_prop_all_volcanic_95_list <- list()
for (i in 1:100){
  fit_regional_allo_prop_all_volcanic_95 <- confint(glm_regional_allo_prop_all_sp_volcanic_list_models[[i]], level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")
  fit_regional_allo_prop_all_volcanic_95_list[[i]] <- fit_regional_allo_prop_all_volcanic_95
  }

# Combining the 95 intervals on all the models
results_regional_allo_prop_all_volcanic_list <- list()
for (i in 1:100){
  results_regional_allo_prop_all_volcanic_list[[i]] <- dplyr::bind_cols(glm_regional_allo_prop_all_sp_volcanic_list_tidy[[i]],fit_regional_allo_prop_all_volcanic_95_list[[i]]) 
}

# Proportion endemics from regional allopatric speciation
glm_coef_hight_regional_allo_prop_volcanic <- as.data.frame(t(sapply(results_regional_allo_prop_all_volcanic_list, function(x) unlist(x[2,]))))
glm_coef_area_regional_allo_prop_volcanic <- as.data.frame(t(sapply(results_regional_allo_prop_all_volcanic_list, function(x) unlist(x[3,]))))
glm_coef_dist_regional_allo_prop_volcanic <- as.data.frame(t(sapply(results_regional_allo_prop_all_volcanic_list, function(x) unlist(x[4,]))))
glm_coef_fragmentation_regional_allo_prop_volcanic <- as.data.frame(t(sapply(results_regional_allo_prop_all_volcanic_list, function(x) unlist(x[5,]))))

# Stacking the dataframes
glm_coefs_regional_allo_prop_volc <- rbind(glm_coef_hight_regional_allo_prop_volcanic,
                            glm_coef_area_regional_allo_prop_volcanic,
                            glm_coef_dist_regional_allo_prop_volcanic,
                            glm_coef_fragmentation_regional_allo_prop_volcanic)

# List of dataframes
glm_coefs_list_regional_allo_prop_volc <- list(glm_coef_hight_regional_allo_prop_volcanic,
                            glm_coef_area_regional_allo_prop_volcanic,
                            glm_coef_dist_regional_allo_prop_volcanic,
                            glm_coef_fragmentation_regional_allo_prop_volcanic)

# regional allopatric proportion means
glm_coefs_means_regional_allo_prop_volc <- data.frame()
for (i in 1:4) {
  temp_means <- vector()
  temp_means <- apply(glm_coefs_list_regional_allo_prop_volc[[i]], 2, function(x) mean(as.numeric(unlist(x))))
  temp_means["term"] <- glm_coefs_list_regional_allo_prop_volc[[i]][1,1]
  glm_coefs_means_regional_allo_prop_volc <- rbind(glm_coefs_means_regional_allo_prop_volc,temp_means)
}
names(glm_coefs_means_regional_allo_prop_volc) <- names(glm_coef_dist_regional_allo_prop_volcanic)







# Make a combined model for only volcanic islands
#Stacking all Coef estimates for violins
glm_coefs_ana_volc$type <- "Regional Allopatry"
glm_coefs_clad_volc$type <- "Regional Sympatry"
glm_coefs_endems_volc$type <- "Number Endemics"
glm_coefs_prop_endems_volc$type <- "Proportion Endemics"
glm_coefs_regional_allo_prop_volc$type <- "Proportion Endemics from Regional Allopatry"

glm_coefs_volc <- data.frame()
glm_coefs_volc <- rbind(glm_coefs_ana_volc,
                        glm_coefs_clad_volc,
                        glm_coefs_endems_volc,
                        glm_coefs_prop_endems_volc,
                        glm_coefs_regional_allo_prop_volc)

#stacking all means for barplots
glm_coefs_mean_ana_volc$type <- "Regional Allopatry"
glm_coefs_mean_clad_volc$type <- "Regional Sympatry"
glm_coefs_means_endems_volc$type <- "Number Endemics"
glm_coefs_means_prop_endems_volc$type <- "Proportion Endemics"
glm_coefs_means_regional_allo_prop_volc$type <- "Proportion Endemics from Regional Allopatry"

glm_coefs_mean_volc <- data.frame()
glm_coefs_mean_volc <- rbind(glm_coefs_mean_ana_volc,
                             glm_coefs_mean_clad_volc,
                             glm_coefs_means_endems_volc,
                             glm_coefs_means_prop_endems_volc,
                             glm_coefs_means_regional_allo_prop_volc)




# making a loop to convert character list rows to numeric 
for(i in 2:7){
  glm_coefs_mean_volc[,i] <-as.numeric(unlist(glm_coefs_mean_volc[,i]))
}

for(i in 2:7){
  glm_coefs_volc[,i] <-as.numeric(unlist(glm_coefs_volc[,i]))
}



# Changeing some names of the means
glm_coefs_mean_names_volc <- c("Max elevation","Area","Isolation","Fragmentation")
glm_coefs_mean_names_volc_old <- c("log10(mean_mx30_grd)","log10(area)","log10(dist)","log10(nearest_neighbour_distance_border_scaled)")
glm_coefs_names_df_volc <- cbind(glm_coefs_mean_names_volc, glm_coefs_mean_names_volc_old)
c

for (i in 1:4){
  glm_coefs_mean_volc[which(glm_coefs_mean_volc$term==glm_coefs_names_df_volc[i,2]),1] <- glm_coefs_names_df_volc[i,1]
}

for (i in 1:4){
  glm_coefs_volc[which(glm_coefs_volc$term==glm_coefs_names_df_volc[i,2]),1] <- glm_coefs_names_df_volc[i,1]
}

# Creating plotting order
plot_order_volcanic <- c("Area", "Isolation", "Max elevation", "Fragmentation")


# Creating a plot using violins to show the distribution of the mean and the 95 confidence interval
gvil_volcanic <- ggplot(glm_coefs_volc, aes(color=type, fill = type)) +
  scale_colour_manual(values = met.brewer("Hiroshige", 5)) +
  scale_fill_manual(values = met.brewer("Hiroshige", 5)) +
  scale_x_discrete(limits = plot_order_volcanic) +
  geom_violin(aes(x = term, y = as.numeric(estimate)),width = 0.68,position=position_dodge(width=0.9), scale = "width", alpha = 0.7,lwd=0.2) +
  geom_violin(aes(x = term, y = as.numeric(conf.high_95)),width = 0.68,position=position_dodge(width=0.9), scale = "width", alpha = 0.7,lwd=0.2) +
  geom_violin(aes(x = term, y = as.numeric(conf.low_95)),width = 0.68,position=position_dodge(width=0.9), scale = "width", alpha = 0.7,lwd=0.2) +
  geom_crossbar(data = glm_coefs_mean_volc, size=0.1, alpha=0.5, aes(x= term, y = estimate,
                    ymin=conf.low_95,
                    ymax=conf.high_95),
                position = position_dodge2(width = 1)) +
    geom_hline(yintercept = 0, colour = gray(0), lty = 3) +
  xlab("Coefficients") +
  ylab("Estimate") +
  theme_classic() +
  labs(fill = element_blank(), colour = element_blank()) +
  theme(legend.position = "None",legend.text=element_text(size=6)) +
  labs(title = "Volcanic Island Botanical Countries") +
  coord_flip()

gvil_volcanic



legend_bottom_volc <- get_legend(gvil_volcanic + guides(color = guide_legend(nrow=1))+theme(legend.position = "bottom", legend.text=element_text(size=6)))

prow_volc <- cowplot::plot_grid(gvil_volcanic,legend_bottom_volc, ncol = 1, rel_heights = c(1,.1))
prow_volc


# Lets quickly do the fit of the some additional models
# Now I can fit the additional models

# -------- Endemics ----------------
glm_endems_pois_all_sp_list_tidy <- list()
glm_endems_pois_all_sp_list_models <- list()
for (i in 1:100){
  glm_endems_all_sp_loop <- glm(colonization_sp+radiating_sp~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled)+relevel(factor(GeologicalOrigin), ref = "volcanic"),data=output_all_sp_test[[i]],family = quasipoisson(link='log'))
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

# -------- Endemic proportion ----------------
glm_prop_endems_all_sp_list_tidy <- list()
glm_prop_endems_all_sp_list_models <- list()
for (i in 1:100) {
  #print((output_all_sp_test[[i]]$colonization_sp+output_all_sp_test[[i]]$radiating_sp)/output_all_sp_test[[i]]$Total_sp )
  glm_prop_endems_all_sp_loop <- glm(((colonization_sp+radiating_sp)/Total_sp)~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled)+relevel(factor(GeologicalOrigin), ref = "volcanic"),data=output_all_sp_test[[i]],family = quasibinomial)
  glm_prop_endems_all_sp_list_models[[i]] <- glm_prop_endems_all_sp_loop
  glm_prop_endems_all_sp_list_tidy[[i]] <- broom::tidy(glm_prop_endems_all_sp_loop)
}

#Now lets find the confidence intervals for these models.
fit_prop_endems_all_95_list <- list()
for (i in 1:100){
  fit_prop_endems_all_95 <- confint(glm_prop_endems_all_sp_list_models[[i]], level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")
  fit_prop_endems_all_95_list[[i]] <- fit_prop_endems_all_95
  }

# Combining the 95 intervals on all the models
results_prop_endems_all_list <- list()
for (i in 1:100){
  results_prop_endems_all_list[[i]] <- dplyr::bind_cols(glm_prop_endems_all_sp_list_tidy[[i]],fit_prop_endems_all_95_list[[i]]) 
}

# -------- Proportion Regional allopatry ----------------
glm_regional_allo_prop_all_sp_list_tidy <- list()
glm_regional_allo_prop_all_sp_list_models <- list()
for (i in 1:100){
#print((output_all_sp_test[[i]]$colonization_sp)/(output_all_sp_test[[i]]$colonization_sp+output_all_sp_test[[i]]$radiating_sp))
  glm_regional_allo_prop_all_sp_loop <- glm(colonization_sp/(colonization_sp+radiating_sp)~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled)+relevel(factor(GeologicalOrigin), ref = "volcanic"), data=output_all_sp_test[[i]],family = quasibinomial)
  glm_regional_allo_prop_all_sp_list_models[[i]] <- glm_regional_allo_prop_all_sp_loop
  glm_regional_allo_prop_all_sp_list_tidy[[i]] <- broom::tidy(glm_regional_allo_prop_all_sp_loop)
}

#Now lets find the confidence intervals for these models.
fit_regional_allo_prop_all_95_list <- list()
for (i in 1:100){
  fit_regional_allo_prop_all_95 <- confint(glm_regional_allo_prop_all_sp_list_models[[i]], level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")
  fit_regional_allo_prop_all_95_list[[i]] <- fit_regional_allo_prop_all_95
  }

# Combining the 95 intervals on all the models
results_regional_allo_prop_all_list <- list()
for (i in 1:100){
  results_regional_allo_prop_all_list[[i]] <- dplyr::bind_cols(glm_regional_allo_prop_all_sp_list_tidy[[i]],fit_regional_allo_prop_all_95_list[[i]]) 
}


# Adding the speciation mode to the results
for (i in 1:100) {
  results_clad_all_list[[i]]$type <- "Regional Sympatry"
  results_ana_all_list[[i]]$type <- "Regional Allopatry"
  results_endems_all_list[[i]]$type <- "Number Endemics"
  results_prop_endems_all_list[[i]]$type <- "Proportion Endemics"
  results_regional_allo_prop_all_list[[i]]$type <- "Proportion Endemics from Regional Allopatry"
}

# I have some trees which cause problems because they have more parameter,s these are (8,26,31,41)
# these 4 models do not have any parameters for the geological origin because we have a single radiating species pair on the bahamas atolls
# because of this the scripts cannot find the upper confidence interval for any of the geological origins 
# I have therefore decided to exclude these 4 models in my script¨
results_clad_all_list_final <- results_clad_all_list[-c(8,26,31,41)]



# So now I have 2 x 100 results tables
# I need to combine these into a massive data.frame that I can use to create a combined violin / bar plot for cladogenesis/Anagenesis 
#Creating dataframes with all estimates for each of the coefficients
#Anagenesis
glm_coef_hight_ana <- as.data.frame(t(sapply(results_ana_all_list, function(x) unlist(x[2,]))))
glm_coef_area_ana <- as.data.frame(t(sapply(results_ana_all_list, function(x) unlist(x[3,]))))
glm_coef_dist_ana <- as.data.frame(t(sapply(results_ana_all_list, function(x) unlist(x[4,]))))
glm_coef_fragmentation_ana <- as.data.frame(t(sapply(results_ana_all_list, function(x) unlist(x[5,]))))
glm_coef_Geoatoll_ana <- as.data.frame(t(sapply(results_ana_all_list, function(x) unlist(x[6,]))))
glm_coef_Geocontinental_ana <- as.data.frame(t(sapply(results_ana_all_list, function(x) unlist(x[7,]))))
glm_coef_Geomixed_ana <- as.data.frame(t(sapply(results_ana_all_list, function(x) unlist(x[8,]))))

#Cladogenesis
glm_coef_hight_clad <- as.data.frame(t(sapply(results_clad_all_list_final, function(x) unlist(x[2,]))))
glm_coef_area_clad <- as.data.frame(t(sapply(results_clad_all_list_final, function(x) unlist(x[3,]))))
glm_coef_dist_clad <- as.data.frame(t(sapply(results_clad_all_list_final, function(x) unlist(x[4,]))))
glm_coef_fragmentation_clad <- as.data.frame(t(sapply(results_clad_all_list_final, function(x) unlist(x[5,]))))
glm_coef_Geocontinental_clad <- as.data.frame(t(sapply(results_clad_all_list_final, function(x) unlist(x[6,]))))
glm_coef_Geomixed_clad <- as.data.frame(t(sapply(results_clad_all_list_final, function(x) unlist(x[7,]))))

# Endems
glm_coef_hight_endems <- as.data.frame(t(sapply(results_endems_all_list, function(x) unlist(x[2,]))))
glm_coef_area_endems <- as.data.frame(t(sapply(results_endems_all_list, function(x) unlist(x[3,]))))
glm_coef_dist_endems <- as.data.frame(t(sapply(results_endems_all_list, function(x) unlist(x[4,]))))
glm_coef_fragmentation_endems <- as.data.frame(t(sapply(results_endems_all_list, function(x) unlist(x[5,]))))
glm_coef_Geoatoll_endems <- as.data.frame(t(sapply(results_endems_all_list, function(x) unlist(x[6,]))))
glm_coef_Geocontinental_endems <- as.data.frame(t(sapply(results_endems_all_list, function(x) unlist(x[7,]))))
glm_coef_Geomixed_endems <- as.data.frame(t(sapply(results_endems_all_list, function(x) unlist(x[8,]))))

# Endems proportion
glm_coef_hight_prop_endems <- as.data.frame(t(sapply(results_prop_endems_all_list, function(x) unlist(x[2,]))))
glm_coef_area_prop_endems <- as.data.frame(t(sapply(results_prop_endems_all_list, function(x) unlist(x[3,]))))
glm_coef_dist_prop_endems <- as.data.frame(t(sapply(results_prop_endems_all_list, function(x) unlist(x[4,]))))
glm_coef_fragmentation_prop_endems <- as.data.frame(t(sapply(results_prop_endems_all_list, function(x) unlist(x[5,]))))
glm_coef_Geoatoll_prop_endems <- as.data.frame(t(sapply(results_prop_endems_all_list, function(x) unlist(x[6,]))))
glm_coef_Geocontinental_prop_endems <- as.data.frame(t(sapply(results_prop_endems_all_list, function(x) unlist(x[7,]))))
glm_coef_Geomixed_prop_endems <- as.data.frame(t(sapply(results_prop_endems_all_list, function(x) unlist(x[8,]))))

# Proportion of endems from regional allopatric speciation
glm_coef_hight_regional_allo_prop <- as.data.frame(t(sapply(results_regional_allo_prop_all_list, function(x) unlist(x[2,]))))
glm_coef_area_regional_allo_prop <- as.data.frame(t(sapply(results_regional_allo_prop_all_list, function(x) unlist(x[3,]))))
glm_coef_dist_regional_allo_prop <- as.data.frame(t(sapply(results_regional_allo_prop_all_list, function(x) unlist(x[4,]))))
glm_coef_fragmentation_regional_allo_prop <- as.data.frame(t(sapply(results_regional_allo_prop_all_list, function(x) unlist(x[5,]))))
glm_coef_Geoatoll_regional_allo_prop <- as.data.frame(t(sapply(results_regional_allo_prop_all_list, function(x) unlist(x[6,]))))
glm_coef_Geocontinental_regional_allo_prop <- as.data.frame(t(sapply(results_regional_allo_prop_all_list, function(x) unlist(x[7,]))))
glm_coef_Geomixed_regional_allo_prop <- as.data.frame(t(sapply(results_regional_allo_prop_all_list, function(x) unlist(x[8,]))))

# Stacking the dataframes
glm_coefs <- rbind(glm_coef_area_ana,
                   glm_coef_hight_ana,
                   glm_coef_dist_ana,
                   glm_coef_fragmentation_ana,
                   glm_coef_Geocontinental_ana,
                   glm_coef_Geomixed_ana,
                   glm_coef_Geoatoll_ana,
                   glm_coef_hight_clad,
                   glm_coef_area_clad,
                   glm_coef_dist_clad,
                   glm_coef_fragmentation_clad,
                   glm_coef_Geomixed_clad,
                   glm_coef_Geocontinental_clad,
                   glm_coef_hight_endems,
                   glm_coef_area_endems,
                   glm_coef_dist_endems,
                   glm_coef_fragmentation_endems,
                   glm_coef_Geoatoll_endems,
                   glm_coef_Geocontinental_endems,
                   glm_coef_Geomixed_endems,
                   glm_coef_hight_prop_endems,
                   glm_coef_area_prop_endems,
                   glm_coef_dist_prop_endems,
                   glm_coef_fragmentation_prop_endems,
                   glm_coef_Geoatoll_prop_endems,
                   glm_coef_Geocontinental_prop_endems,
                   glm_coef_Geomixed_prop_endems,
                   glm_coef_hight_regional_allo_prop,
                   glm_coef_area_regional_allo_prop,
                   glm_coef_dist_regional_allo_prop,
                   glm_coef_fragmentation_regional_allo_prop,
                   #glm_coef_Geoatoll_regional_allo_prop,
                   glm_coef_Geocontinental_regional_allo_prop,
                   glm_coef_Geomixed_regional_allo_prop)

glm_coefs_list_clad <- list(glm_coef_hight_clad,
                   glm_coef_area_clad,
                   glm_coef_dist_clad,
                   glm_coef_fragmentation_clad,
                   glm_coef_Geomixed_clad,
                   glm_coef_Geocontinental_clad)

glm_coefs_list_ana <- list(glm_coef_area_ana,
                   glm_coef_hight_ana,
                   glm_coef_dist_ana,
                   glm_coef_fragmentation_ana,
                   glm_coef_Geocontinental_ana,
                   glm_coef_Geomixed_ana,
                   glm_coef_Geoatoll_ana)

glm_coefs_list_endems <- list(glm_coef_area_endems,
                   glm_coef_hight_endems,
                   glm_coef_dist_endems,
                   glm_coef_fragmentation_endems,
                   glm_coef_Geocontinental_endems,
                   glm_coef_Geomixed_endems,
                   glm_coef_Geoatoll_endems)

glm_coefs_list_prop_endems <- list(glm_coef_area_prop_endems,
                   glm_coef_hight_prop_endems,
                   glm_coef_dist_prop_endems,
                   glm_coef_fragmentation_prop_endems,
                   glm_coef_Geocontinental_prop_endems,
                   glm_coef_Geomixed_prop_endems,
                   glm_coef_Geoatoll_prop_endems)

glm_coefs_list_regional_allo_prop <- list(glm_coef_area_regional_allo_prop,
                   glm_coef_hight_regional_allo_prop,
                   glm_coef_dist_regional_allo_prop,
                   glm_coef_fragmentation_regional_allo_prop,
                   glm_coef_Geocontinental_regional_allo_prop,
                   glm_coef_Geomixed_regional_allo_prop
                   ) # glm_coef_Geoatoll_regional_allo_prop


# Renaming some variables
#Writing up the names
glm_coefs_names <- c("Area","Max elevation","Isolation","Fragmentation","Continental","Mixed Origin", "Atoll")
glm_coefs_names_old <- unique(glm_coefs$term)
glm_coefs_names_df <- cbind(glm_coefs_names, glm_coefs_names_old)

# For loop for renaming
for (i in 1:7){
  glm_coefs[which(glm_coefs$term==glm_coefs_names_df[i,2]),1] <- glm_coefs_names_df[i,1]
}

# Now I only need to create the dataframe with the means and then I think I can plot it.
#Anagenesis means
glm_coefs_means_ana <- data.frame()
for (i in 1:7) {
  temp_means <- vector()
  temp_means <- apply(glm_coefs_list_ana[[i]], 2, function(x) mean(as.numeric(unlist(x))))
  temp_means["term"] <- glm_coefs_list_ana[[i]][1,1]
  temp_means["type"] <- glm_coefs_list_ana[[i]][1,8]
  #print(temp_means)
  glm_coefs_means_ana <- rbind(glm_coefs_means_ana,temp_means)
  print(glm_coefs_means_ana)
}
names(glm_coefs_means_ana) <- names(glm_coef_dist_ana)

# making a loop to convert character list rows to numeric 
for(i in 2:7){
  glm_coefs_means_ana[,i] <-as.numeric(unlist(glm_coefs_means_ana[,i]))
}


# Cladogenesis means
glm_coefs_means_clad <- data.frame()
for (i in 1:6) {
  temp_means <- vector()
  temp_means <- apply(glm_coefs_list_clad[[i]], 2, function(x) mean(as.numeric(unlist(x))))
  temp_means["term"] <- glm_coefs_list_clad[[i]][1,1]
  temp_means["type"] <- glm_coefs_list_clad[[i]][1,8]
  glm_coefs_means_clad <- rbind(glm_coefs_means_clad,temp_means)
}
names(glm_coefs_means_clad) <- names(glm_coef_dist_clad)

# Endemics means
glm_coefs_means_endems <- data.frame()
for (i in 1:7) {
  temp_means <- vector()
  temp_means <- apply(glm_coefs_list_endems[[i]], 2, function(x) mean(as.numeric(unlist(x))))
  temp_means["term"] <- glm_coefs_list_endems[[i]][1,1]
  temp_means["type"] <- glm_coefs_list_endems[[i]][1,8]
  glm_coefs_means_endems <- rbind(glm_coefs_means_endems,temp_means)
}
names(glm_coefs_means_endems) <- names(glm_coef_dist_endems)

# proportion endemic means
glm_coefs_means_prop_endems <- data.frame()
for (i in 1:7) {
  temp_means <- vector()
  temp_means <- apply(glm_coefs_list_prop_endems[[i]], 2, function(x) mean(as.numeric(unlist(x))))
  temp_means["term"] <- glm_coefs_list_prop_endems[[i]][1,1]
  temp_means["type"] <- glm_coefs_list_prop_endems[[i]][1,8]
  glm_coefs_means_prop_endems <- rbind(glm_coefs_means_prop_endems,temp_means)
}
names(glm_coefs_means_prop_endems) <- names(glm_coef_dist_prop_endems)

# Regional allopatry proportion means
glm_coefs_means_regional_allo_prop <- data.frame()
for (i in 1:6) {
  temp_means <- vector()
  temp_means <- apply(glm_coefs_list_regional_allo_prop[[i]], 2, function(x) mean(as.numeric(unlist(x))))
  temp_means["term"] <- glm_coefs_list_regional_allo_prop[[i]][1,1]
  temp_means["type"] <- glm_coefs_list_regional_allo_prop[[i]][1,8]
  glm_coefs_means_regional_allo_prop <- rbind(glm_coefs_means_regional_allo_prop,temp_means)
}
names(glm_coefs_means_regional_allo_prop) <- names(glm_coef_dist_regional_allo_prop)

# Stacking the means dataframes for the means
coefs_means <- NULL
coefs_means <- rbind(glm_coefs_means_clad, glm_coefs_means_ana, glm_coefs_means_endems, glm_coefs_means_prop_endems, glm_coefs_means_regional_allo_prop)

# making a loop to convert character list rows to numeric 
for(i in 2:7){
  coefs_means[,i] <-as.numeric(unlist(coefs_means[,i]))
}

# Changing names in dataframe of estimates
for (i in 1:7){
  coefs_means[which(coefs_means$term==glm_coefs_names_df[i,2]),1] <- glm_coefs_names_df[i,1]
}

# Creating plotting order
plot_order <- c("Area", "Isolation", "Max elevation", "Fragmentation", "Continental", "Mixed Origin", "Atoll")

#Now I should be ready to make my violin + box plot for all species on all the islands.
gvil_all_all <- ggplot(glm_coefs, aes(color=type, fill = type)) +
  scale_colour_manual(values = met.brewer("Hiroshige", 5)) +
  scale_fill_manual(values = met.brewer("Hiroshige", 5)) +
  scale_x_discrete(limits = plot_order) +
  geom_violin(aes(x = term, y = as.numeric(estimate)),width = 0.60,position=position_dodge(width=0.9), scale = "width", alpha = 0.5,lwd=0.2) +
  geom_violin(aes(x = term, y = as.numeric(conf.high_95)),width = 0.60,position=position_dodge(width=0.9), scale = "width", alpha = 0.5,lwd=0.2) +
  geom_violin(aes(x = term, y = as.numeric(conf.low_95)),width = 0.60,position=position_dodge(width=0.9), scale = "width", alpha = 0.5,lwd=0.2) +
  geom_crossbar(data = coefs_means, size=0.1, alpha=0.7, aes(x= term, y = estimate,
                    ymin=conf.low_95,
                    ymax=conf.high_95),
                position = position_dodge2(width = 1)) +
  geom_hline(yintercept = 0, colour = gray(0), lty = 3) +
  xlab("Coefficients") +
  ylab("Estimate") +
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()

gvil_all_all


# 
# 
# logo_phyl <- readPNG(file.path("./","Phylo_for_plot.png")) %>% 
#   image_read()
# bitmap <- logo_phyl[[1]]
# bitmap[4,,] <- as.raw(as.integer(bitmap[4,,])*0.3)
# logo_phyl <- image_read(bitmap)
# 
# px <- ggdraw() +
#   draw_plot(px) +
#   draw_image(logo_phyl, scale = 0.25, hjust = 0.15, vjust = 0.23)
# 
# px



# Now I want to run the best model for only the coryphoidea and then add them through the cowplot
glm_ana_cory <- glm(`No. anagenesis`+`No. edges` ~
                      log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled)+
                    relevel(factor(GeologicalOrigin), ref = "volcanic"),
                    data=output_with_island[which(output_with_island$nearest_neighbour_distance_border_scaled != 0),],
                    na.action = na.fail,
                    family = quasipoisson(link='log'))

dp_glm_ana_cory = sum(residuals(glm_ana_cory,type ="pearson")^2)/glm_ana_cory$df.residual

summary(glm_ana_cory, dispersion = dp_glm_ana_cory)



#Cladogenesis
glm_clad_cory <- glm(`No. cladogenesis`-`No. edges`~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled), data=output_with_island[which(output_with_island$nearest_neighbour_distance_border_scaled != 0),], na.action = na.fail,family = quasipoisson(link='log'))

summary(glm_clad_cory)




#Endems
glm_endems_cory <- glm(`No. cladogenesis`+`No. anagenesis`~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled)+
                    relevel(factor(GeologicalOrigin), ref = "volcanic"), data=output_with_island[which(output_with_island$nearest_neighbour_distance_border_scaled != 0),], na.action = na.fail,family = quasipoisson(link='log'))

summary(glm_endems_cory)


#proportion Endems
glm_prop_endems_cory <- glm((`No. cladogenesis`+`No. anagenesis`)/Total_sp_coryphoids~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled)+
                    relevel(factor(GeologicalOrigin), ref = "volcanic"), data=output_with_island[which(output_with_island$nearest_neighbour_distance_border_scaled != 0),], na.action = na.fail,family = quasibinomial)

summary(glm_prop_endems_cory)






# Proportion of endemics from regional allopatric speciation
glm_regional_allo_prop_cory <- glm(`No. anagenesis`/(`No. anagenesis`+`No. cladogenesis`)~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled)+
                    relevel(factor(GeologicalOrigin), ref = "volcanic"), data=output_with_island[which(output_with_island$nearest_neighbour_distance_border_scaled != 0 & output_with_island$`No. anagenesis` > 0),], na.action = na.fail,family = quasipoisson(link='log'))

summary(glm_regional_allo_prop_cory)




#Getting the coefficient estimates for cladogenesis and anagenesis

#Anagenesis
results_ana_cory <- broom::tidy(glm_ana_cory)

#Now lets find the confidence intervals for these models.
fit_ana_cory <- confint(glm_ana_cory, level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")

results_ana_cory <- dplyr::bind_cols(results_ana_cory,fit_ana_cory) 

#cladogenesis
results_clad_cory <- broom::tidy(glm_clad_cory)

#Now lets find the confidence intervals for these models.
fit_clad_cory <- confint(glm_clad_cory, level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")

results_clad_cory <- dplyr::bind_cols(results_clad_cory,fit_clad_cory)

# No endems
results_endems_cory <- broom::tidy(glm_endems_cory)

#Now lets find the confidence intervals for these models.
fit_endems_cory <- confint(glm_endems_cory, level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")

results_endems_cory <- dplyr::bind_cols(results_endems_cory,fit_endems_cory)
results_endems_cory <- results_endems_cory[c(2,3,4,5,7,8),]

#Endemics proportion
results_prop_endems_cory <- broom::tidy(glm_prop_endems_cory)

#Now lets find the confidence intervals for these models.
fit_prop_endems_cory <- confint(glm_prop_endems_cory, level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")

results_prop_endems_cory <- dplyr::bind_cols(results_prop_endems_cory,fit_prop_endems_cory)


#Proportion endemics from regional allopatry
results_regional_allo_prop_cory <- broom::tidy(glm_regional_allo_prop_cory)

#Now lets find the confidence intervals for these models.
fit_regional_allo_prop_cory <- confint(glm_regional_allo_prop_cory, level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")

results_regional_allo_prop_cory <- dplyr::bind_cols(results_regional_allo_prop_cory,fit_regional_allo_prop_cory)


results_ana_cory$type <-"Regional Allopatry"
results_clad_cory$type <- "Regionnal Sympatry"
results_endems_cory$type <- "Number Endemics"
results_prop_endems_cory$type <- "Proportion Endemics"
results_regional_allo_prop_cory$type <- "Proportion Endemics from Regional Allopatry"

results_cory <- dplyr::bind_rows(results_ana_cory,
                                 results_clad_cory,
                                 results_endems_cory,
                                 results_prop_endems_cory,
                                 results_regional_allo_prop_cory)


# For loop for renaming
for (i in 1:7){
  results_cory[which(results_cory$term==glm_coefs_names_df[i,2]),1] <- glm_coefs_names_df[i,1]
}


#Now I should be ready to make my violin + box plot for all species on all the islands.
gvil_coryphoideae <- ggplot(results_cory, aes(color=type, fill = type)) +
  scale_colour_manual(values = met.brewer("Hiroshige", 5)) +
  scale_fill_manual(values = met.brewer("Hiroshige", 5)) +
  scale_x_discrete(limits = plot_order) +
  geom_crossbar(data = results_cory, size=0.2, alpha=0.7, aes(x= term, y = estimate,
                    ymin=conf.low_95,
                    ymax=conf.high_95),
                position = position_dodge2(width = 1.5)) +
  geom_hline(yintercept = 0, colour = gray(0), lty = 3) +
  xlab("Coefficients") +
  ylab("Estimate") +
  scale_y_continuous(breaks = seq(-8, 8, by = 4), limits = c(-8,8)) +
  theme_classic() +
  labs(fill = element_blank(), colour = element_blank()) +
  theme(axis.title.y = element_blank(), legend.position = "none") +
  coord_flip()

gvil_coryphoideae



prow <- cowplot::plot_grid(gvil_all_all, gvil_coryphoideae+theme(legend.position = "none"), labels = c("a.", "b."), align = "v", rel_widths = c(1,1))
prow

legend_bottom <- get_legend(gvil_coryphoideae + guides(color = guide_legend(nrow=1))+theme(legend.position = "bottom", legend.text=element_text(size=6)))

prow <- cowplot::plot_grid(prow,legend_bottom, ncol = 1, rel_heights = c(1,.1))
prow

