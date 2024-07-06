# Read tree
cory_phylo <- read.tree(file.path(datadir,"astral_tree_orthologs.tre"))

# Read renaming file
figurename <- read.table(file.path(datadir,"names_for_tips.csv"), sep=",", colClasses = "character")

#figurename <- figurename[]
figurename_idx <- figurename$V2
names(figurename_idx) <- figurename$V1

#Renaming tips
cory_phylo$tip.label = figurename_idx[cory_phylo$tip.label]

# Species which I have multiple of
# Brahea dulcis # drop 1 and 3 keep 2
# Brahea edulis # drop 2 keep 1
# Coccothrinax acuminata # drop 2 keep 1
# Coccothrinax argentata # drop 1, 2 and 4 keep 3
# Coccothrinax argentea # drop 2 keep 1
# Coccothrinax barbadensis # drop 1 keep 2
# Coccothrinax garciana # drop 2 keep 1
# Coccothrinax yuraguana # there is no second species or maybe there is but it is not in the phylogeny
# Copernicia cowellii # drop 2 keep 1
# Corypha lecomtei # drop 2 keep 1 
# Licuala lauterbachii # drop 1 keep 2
# Phoenix paludosa # drop 1 keep 2
# Pritchardia flynnii # drop 2 and 3 keep 1
# Pritchardia glabrata # drop 2 keep 1
# Pritchardia hardyi # drop 2 and 3 keep 1
# Pritchardia kahukuensis # drop 2 keep 1
# Pritchardia munroi # drop 1 keep 2
# Pritchardia perlmanii # drop 1 keep 2
# Pritchadia viscosa # drop 2 keep 1
# Rhapis robusta # drop 2 keep 1
# Rhapis puhuongensis # drop 2 keep 1
# Trithrinax brasiliensis # drop 2 and 3 keep 1
# Trithrinax schizophylla # drop 2 and 3 keep 1
# Zombia antillarum # drop 2

# Tips which I need to drop because I have duplicates
tips_to_drop <- c("Brahea dulcis 1", "Brahea dulcis 3", "Brahea edulis 2", "Coccothrinax acuminata 2", "Coccothrinax argentata 1", "Coccothrinax argentata 2", "Coccothrinax argentata 4", "Coccothrinax argentea 2", "Coccothrinax barbadensis 1", "Coccothrinax garciana 2", "Copernicia cowellii 2", "Corypha lecomtei 2", "Licuala lautherbachii 1", "Phoenix paludosa 1", "Pritchardua flynnii 2", "Pritchardua flynnii 3", "Pritchardia glabrata 2", "Pritchardia hardyi 2","Pritchardia hardyi 3", "Pritchardia kahukuensis 2", "Pritchardia munroi 1", "Pritchardia perlmanii 2", "Pritchardia viscosa 2", "Rhapis robusta 2", "Rhapus puhuongensis 2", "Trithrinax brasiliensis 2", "Trithrinax brasiliensis 3", "Trithrinax schizophylla 2", "Trithrinax schizophylla 3", "Zombia antillarum 2") 
length(tips_to_drop) # 30

#Dropping outgroup, novo species and a missing subsp.
root <- c("Aphandra natalia","Eugeissona tristis","Nypa fruticans","Chrysalidocarpus ambositrae") # 4

# Species which are not published or identified
unpublished <- c("Chelyocarpus sp. José De Gracia 826","Chelyocarpus sp. José De Gracia 892","Coccothrinax aff. Jimenezii T. A. Zanoni 33375","Coccothrinax bermudezi","Coccothrinax bonettiana","Coccothrinax montgomeryana","Cryosophila sp.  José De Gracia 824","Lanonia sp. Henderson A. 3590","Lanonia sp. Henderson A. 3727","Licuala heatubunii","Licuala vanroyensis","Livistona sp. Henderson A. 3730","Sabal miamiensis","Licuala polybracteae")
length(unpublished) #14

#Sub sp and Variants
subspANDvars <- c("Cryosophila kalbreyeri subsp. cogolloi","Cryosophila kalbreyeri subsp. kalbreyeri","Licuala glabra var. selangorensis","Licuala lauterbachii subsp. Peekelii","Licuala lauterbachii var. bougainvillensis","Licuala lauterbachii var. lauterbachii","Licuala malajana var. humilis","Licuala mattanensis var. paucisecta","Licuala peltata var. sumawongii","Phoenix loureiroi var. pedunculata","Rhapis subtilis subsp. Siamensis","Rhapis subtilis subsp. Subtilis","Washingtonia filifera var. robusta","Rhapis laosensis subsp. macrantha","Licuala ramsayi var. tuckeri")
length(subspANDvars) #15

# Dropping tips
cory_phylo <- drop.tip(cory_phylo,tips_to_drop) # tips to drop
cory_phylo <- drop.tip(cory_phylo,root) # root species
cory_phylo <- drop.tip(cory_phylo,unpublished) # unpublished or unidentified species
cory_phylo <- drop.tip(cory_phylo, subspANDvars) # sub species and vars

# Remove the trailing numbers
cory_phylo$tip.label <- gsub(" \\d+$", "", cory_phylo$tip.label)

# Tips where I need to convert a sub sp to the sp
cory_phylo$tip.label["3130"][1] <- "Licuala glabra"
cory_phylo$tip.label["3176"][1] <- "Licuala malajana"
cory_phylo$tip.label["3412"][1] <- "Licuala mattanensis"
cory_phylo$tip.label["3334"][1] <- "Licuala peltata"
cory_phylo$tip.label["3390"][1] <- "Licuala ramsayi"
cory_phylo$tip.label["1022"][1] <- "Phoenix loureiroi"



#Making tree ultrametric
cory_phylo$edge.length[is.na(cory_phylo$edge.length)] <- 0.001



# Loading all genes tree, renaming tips to species names and dropping duplicates and bad species.


#Loading tree file
cory_phylo_all_genes<- read.tree(file.path(datadir,"astral_tree.tre"))

# Dropping tips
cory_phylo_all_genes <- drop.tip(cory_phylo_all_genes,tips_to_drop)
cory_phylo_all_genes <- drop.tip(cory_phylo_all_genes,unpublished_and_root)


# Renaming tips with problematic names
#cory_phylo$tip.label["1042"][1] <- "Brahea edulis"
cory_phylo_all_genes$tip.label["1193"][1] <- "Coccothrinax barbadensis"
cory_phylo_all_genes$tip.label["2004"][1] <- "Arenga undulatifolia"
cory_phylo_all_genes$tip.label["1044"][1] <- "Chamaerops humilis"
cory_phylo_all_genes$tip.label["1194"][1] <- "Chelyocarpus chuco"
cory_phylo_all_genes$tip.label["1249"][1] <- "Coccothrinax argentata"
cory_phylo_all_genes$tip.label["1193"][1] <- "Coccothrinax barbadensis"
cory_phylo_all_genes$tip.label["1049"][1] <- "Copernicia prunifera"
cory_phylo_all_genes$tip.label["1251"][1] <- "Hemithrinax ekmaniana"
cory_phylo_all_genes$tip.label["2007"][1] <- "Johannesteijsmannia altifrons"
cory_phylo_all_genes$tip.label["1224"][1] <- "Trithrinax brasiliensis"
cory_phylo_all_genes$tip.label["1437"][1] <- "Trithrinax schizophylla"
cory_phylo_all_genes$tip.label["1078"][1] <- "Zombia antillarum"
#cory_phylo$tip.label["1061"][1] <- "Licuala telifera"
cory_phylo_all_genes$tip.label["1074"][1] <- "Serenoa repens"
cory_phylo_all_genes$tip.label["1434"][1] <- "Schippia concolor"
cory_phylo_all_genes$tip.label["1022"][1] <- "Phoenix loureiroi var. loureiroi"
cory_phylo_all_genes$tip.label["1023"][1] <- "Phoenix loureiroi var. pedunculata"
cory_phylo_all_genes$tip.label["3158"][1] <- "Licuala lauterbachii var. bougainvillensis"
cory_phylo_all_genes$tip.label["4012"][1] <- "Coccothrinax acuminata"
cory_phylo_all_genes$tip.label["4015"][1] <- "Coccothrinax garciana"
cory_phylo_all_genes$tip.label["4020"][1] <- "Coccothrinax yuraguana"
cory_phylo_all_genes$tip.label["4022"][1] <- "Copernicia brittonorum"
cory_phylo_all_genes$tip.label["4023"][1] <- "Copernicia cowellii"
cory_phylo_all_genes$tip.label["4011"][1] <- "Cryosophila kalbreyeri subsp. kalbreyeri"
cory_phylo_all_genes$tip.label["4061"][1] <- "Pritchardia glabrata"
cory_phylo_all_genes$tip.label["4057"][1] <- "Pritchardia flynnii"
cory_phylo_all_genes$tip.label["4063"][1] <- "Pritchardia hardyi"
cory_phylo_all_genes$tip.label["4068"][1] <- "Pritchardia kahukuensis"
cory_phylo_all_genes$tip.label["4074"][1] <- "Pritchardia munroi"
cory_phylo_all_genes$tip.label["4077"][1] <- "Pritchardia pacifica"
cory_phylo_all_genes$tip.label["4078"][1] <- "Pritchardia perlmanii"
cory_phylo_all_genes$tip.label["4081"][1] <- "Pritchardia viscosa"
cory_phylo_all_genes$tip.label["2058"][1] <- "Rhapis laosensis subsp. macrantha"

#Making tree ultrametric
cory_phylo_all_genes$edge.length[is.na(cory_phylo_all_genes$edge.length)] <- 0.001


#Now were loading the WCSP name dataset in order to create a list of all accepted names.´
#Creating a list of all accepted names within Coryphoideae
data <- as.data.frame(fread(file.path(datadir,"wcvp_names.csv"), sep = "|"))

species <- data[which(data$taxon_rank == "Species" | data$taxon_rank == "Variety" | data$taxon_rank == "Subspecies"),] # Selecting only species and varieties/subspecies
palms <- species[which(species$family=="Arecaceae"),] # Selecting only Palms
apalms <- palms[which(palms$taxon_status=="Accepted"),] # Selecting only Accepted species

# List of Genera in Coryphoideae
genera <- c("Sabal","Schippia","Trithrinax","Zombia","Coccothrinax","Hemithrinax","Leucothrinax","Thrinax","Chelyocarpus",
            "Cryosophila","Itaya","Phoenix","Chamaerops","Guihaia","Trachycarpus","Rhapidophyllum","Maxburretia","Rhapis",
            "Livistona","Licuala","Johannesteijsmannia","Pholidocarpus","Pritchardiopsis","Acoelorraphe","Serenoa","Brahea",
            "Colpothrinax","Copernicia","Pritchardia","Washingtonia","Chuniophoenix","Kerriodoxa","Nannorrhops","Tahina",
            "Caryota","Arenga","Wallichia","Corypha","Bismarckia","Satranala","Hyphaene","Medemia","Latania","Lodoicea",
            "Borassodendron","Borassus","Lanonia","Saribus","Sabinaria","Acoelorraphe")

# Selecting only palms in Coryphoideae
cory_data <- apalms[which(apalms$genus %in% genera),]

#Dropping synonym species names rows in data.frame
cory_data_accepted <- cory_data[which(cory_data$taxon_status =="Accepted"),] # (cory_data$taxon_rank=="Variety"))
unique(cory_data_accepted$genus)

# Find out how many accepted subsp and vars there are in Coryphoidear
subsp_var <- cory_data_accepted[which( cory_data_accepted$taxon_rank == "Variety" | cory_data_accepted$taxon_rank == "Subspecies"),]
dim(subsp_var) #  there are 56 sub sp and vars

dim(cory_data_accepted)

#Dropping synonym species names rows in data.frame
cory_data_accepted_only_sp <- cory_data_accepted[which(cory_data_accepted$infraspecific_rank != "var."),] 
cory_data_accepted_only_sp <- cory_data_accepted_only_sp[which(cory_data_accepted_only_sp$infraspecific_rank != "subsp."),]
cory_data_accepted_only_sp <- cory_data_accepted_only_sp[which(cory_data_accepted_only_sp$species_hybrid == ""),]
unique(cory_data_accepted_only_sp$genus)

# Do i just need to split it
# splitting the plant_name_id column
cory_data_accepted_only_sp$plant_name_id <- gsub("-wcs", "", cory_data_accepted_only_sp$plant_name_id)
cory_data_accepted_only_sp
unique(cory_data_accepted_only_sp$genus)


# Trying to load WCVP distribution data to see if it is better
dist_data <- read.csv(file.path(datadir,"wcvp_distribution.csv"), sep = "|") 
cory_dist_data <- dist_data[which(dist_data$plant_name_id %in% cory_data_accepted_only_sp$plant_name_id),]
cory_dist_data

#Renaming WCS ID's to species names
uniqe_names_cory <- unique(cory_dist_data$plant_name_id)
for (i in 1:length(unique(cory_dist_data$plant_name_id))){
  id <- uniqe_names_cory[i]
  
  if (id %in% cory_dist_data$plant_name_id){
    cory_dist_data[which(cory_dist_data$plant_name_id == id),2] <-cory_data_accepted_only_sp[which(cory_data_accepted_only_sp$plant_name_id ==id),22]
  } else {
    print("not in")
  }
}

#We probably need to remove the introduced locations.
cory_dist_data_no_introduced <- cory_dist_data[which(cory_dist_data$introduced != 1 | cory_dist_data$extinct != 1),]


#Soo how many different locations are there in the new database??
areas <- unique(cory_dist_data_no_introduced$area) # 203. Seems doable.
cory_dist_binary <- as.data.frame(cory_dist_data_no_introduced[,1])

#Converting data to wide format and 
community_matrixes <- cory_dist_data_no_introduced %>%
  dplyr::count(plant_name_id, area) %>%
  tidyr::pivot_wider(names_from = plant_name_id, values_from = n, values_fill = 0) %>%
  as.matrix()

# Setting Column names
names(community_matrixes) <- community_matrixes[1,]
#community_matrixes <- community_matrixes[-1,]

#Setting row names 
row.names(community_matrixes) <- community_matrixes[,1]
community_matrixes <- community_matrixes[,-1]

#Checking data overlap between tree and distrubution data. 
community_matrixes <- t(community_matrixes)
chk <- name.check(cory_phylo, community_matrixes)
chk

#Removing tips which are not in data
cory_phylo <- drop.tip(cory_phylo, chk$tree_not_data)

#Removing species from data which are not in tree
community_matrixes_orthologs <- community_matrixes[!(rownames(community_matrixes) %in% chk$data_not_tree),]

community_matrixes_orthologs <- t(community_matrixes_orthologs)

#Renaming countries which do not match with the botanical countries names.
rownames(community_matrixes_orthologs)[which(rownames(community_matrixes_orthologs) == "Gambia")] <- "Gambia, The"
rownames(community_matrixes_orthologs)[which(rownames(community_matrixes_orthologs) == "Panamá")] <- "Panama"
rownames(community_matrixes_orthologs)[which(rownames(community_matrixes_orthologs) == "Central African Repu")] <- "Central African Republic"
rownames(community_matrixes_orthologs)[which(rownames(community_matrixes_orthologs) == "Leeward Is. AB Ant")] <- "Leeward Is."
#rownames(community_matrixes_orthologs)[which(rownames(community_matrixes_orthologs) == "Mozambique Channel I")] <- "Mozambique Channel Is."
community_matrixes_orthologs <- t(community_matrixes_orthologs)


islands_cory <- c("Borneo", "Sumatera", "Hainan", "Taiwan", "Christmas I.", "Maluku", "New Guinea", "Philippines", "Jawa", "Sulawesi", "Nansei-shoto", "Madagascar", "Gulf of Guinea Is.", "Lesser Sunda Is.", "Sri Lanka", "Mexican Pacific Is.", "Andaman Is.", "Nicobar Is.", "Vanuatu", "Bismarck Archipelago", "Solomon Is.",  "Sicilia", "Cuba", "Leeward Is.", "Bahamas", "Turks-Caicos Is.",  "Dominican Republic","Haiti", "Trinidad-Tobago", "Venezuelan Antilles", "Windward Is.", "Jamaica", "Cayman Is.","Mozambique Channel I", "Hawaii" )

# Find number of single island endemics in Coryphoideae
community_matrixes_numbs <- community_matrixes

#community_matrixes_numbs <- t(community_matrixes)
#community_matrixes_numbs <- community_matrixes_numbs[,islands]

# Defining the class of the data of the matrixes
class(community_matrixes_numbs) <- "numeric"

#Adding a Total column to find single island endemics
community_matrixes_numbs_total <- cbind(community_matrixes_numbs, rowSums(community_matrixes_numbs[]))

#names(community_matrixes_numbs_total)[names]

colnames(community_matrixes_numbs_total)[184] <- "Total"
community_matrixes_numbs_endems <- community_matrixes_numbs[which(community_matrixes_numbs_total[,184] == 1),]



count_species_areas <- function(community_matrix, areas_list) {
  # Find species exclusive to one area
  exclusive_species <- rowSums(community_matrix[, areas_list] == 1) == 1
  
  # Find species found in multiple areas
  multiple_areas_species <- rowSums(community_matrix[, areas_list] == 1) > 1
  
  # Count the number of species exclusive to one area
  num_exclusive_species <- sum(exclusive_species)
  
  # Count the number of species found in multiple areas
  num_multiple_areas_species <- sum(multiple_areas_species)
  
  # Return the results as a named list
  result <- list(
    exclusive_species = num_exclusive_species,
    multiple_areas_species = num_multiple_areas_species
  )
  
  return(result)
}

count_species_areas(community_matrixes_numbs,islands_cory)

edges_per_island <- list()
species_per_island <- list()
comb <- list()

  for (i in 1:length(islands_cory)) {
  edges_with_only_island_decen <- c()
  sp_on_island <- c()
  decen_edges_island <- c()
  decendants_of_edges <- c()
  for(node in (1:length(cory_phylo$node.label) + Ntip(cory_phylo)+2)){ # Looping through internal branches minus the root
    decendants <- getDescendants(cory_phylo,node) # Finding descendants from each node
    decen_tips <- cory_phylo$tip.label[which(cory_phylo$tip.label %in% cory_phylo$tip.label[decendants])] #Getting the name of the descendant tips
    decen_edges <- decendants[which(!(decendants %in% decen_tips))] # Finding all descendant edges
    found_on_island <- (community_matrixes_orthologs[decen_tips,islands_cory[i]] == "1") 
    island_chk <- all(found_on_island == TRUE) # Checking if all descendants are found on the island
    if (island_chk == TRUE) {
      if (!(is.null(getParent(cory_phylo,node))) && !(node %in% decen_edges_island)) { #This line checks if the parent node is the root, then it checks if the node is a descendant of a node already added to the list.
      edges_with_only_island_decen <- append(edges_with_only_island_decen,node) #Adding edge to a list of edges endemic to the island 
      decen_edges_island <- append(decen_edges_island,decen_edges) # Adding all edges descending from this node to a list
      decendants_of_edges <- append(decendants_of_edges,decen_tips)
      }
    }
  }
  for(ntip in (1:Ntip(cory_phylo))) {
    if (!(cory_phylo$tip.label[ntip] %in% decendants_of_edges)){ # Checking if species is a descendant of already saved edges
      if (community_matrixes_orthologs[cory_phylo$tip.label[ntip],islands_cory[i]] == "1"){ # Checking if species is found on island
        sp_on_island <- append(sp_on_island,cory_phylo$tip.label[ntip] ) # Saving tips
      }
    }
  }
  edges_per_island[[i]] <- list(edges_with_only_island_decen) # Saving all edges where all descendants are endemic to an island
  species_per_island[[i]] <- list(sp_on_island) # Saving all other species found on island
}

names(edges_per_island) <- islands_cory
names(species_per_island) <- islands_cory

comb <- list()
for(i in 1:length(islands_cory)){
  print(islands_cory[i])
  comb[[i]] <- list(edges_per_island[islands_cory[i]], species_per_island[islands_cory[i]])
  names(comb[[i]][[1]]) <- "Edges"
  names(comb[[i]][[2]]) <- "Species"
}
names(comb) <- islands_cory


FindEndemicLineages <- function(tree,matrix,areas){

  edges_per_island <- list()
  species_per_island <- list()
  endemics_from_radiation <- list()
  comb <- list()
  
    for (i in 1:length(areas)) {
    edges_with_only_island_decen <- c()
    sp_on_island <- c()
    decen_edges_island <- c()
    decendants_of_edges <- c()
    
      for(node in (1:length(tree$node.label) + Ntip(tree))){ # Looping through internal branches
        decendants <- getDescendants(tree,node) # Finding descendants from each node
        decen_tips <- tree$tip.label[which(tree$tip.label %in% tree$tip.label[decendants])] #Getting the name of the descendant tips
        decen_edges <- decendants[which(!(decendants %in% decen_tips))] # Finding all descendant edges

        found_on_island <- (matrix[decen_tips,areas[i]] == "1")
        
        only_on_island <- (matrix[decen_tips, 
                                  which(colnames(matrix[decen_tips,]) != areas[i])]
                           == "0") # This line should check if all the descendants of edge is only found on this island
        
        
        island_chk <- all(found_on_island == TRUE) # Checking if all descendants are found on the island
        endemic_chk <- all(only_on_island == TRUE) # Checking if all descendants are endemic to the island
        
        if (island_chk == TRUE && endemic_chk ==TRUE) {
          if (!(is.null(getParent(cory_phylo,node))) && !(node %in% decen_edges_island)) { #This line checks if the parent node is the root, then it checks if the node is a descendant of a node already added to the list.
          edges_with_only_island_decen <- append(edges_with_only_island_decen,node) #Adding edge to a list of edges endemic to the island 
          decen_edges_island <- append(decen_edges_island,decen_edges) # Adding all edges descending from this node to a list
          decendants_of_edges <- append(decendants_of_edges,decen_tips)
          }
        }
      }
    for(ntip in (1:Ntip(tree))) {
      if (!(tree$tip.label[ntip] %in% decendants_of_edges)){ # Checking if species is a descendant of already saved edges
        endemic_sp <- all(matrix[tree$tip.label[ntip], # Selecting Row
                                 which(names(matrix[tree$tip.label[[ntip]],]) != areas[i])] # Selecting Columns
                                 == "0")
        if (matrix[tree$tip.label[ntip],areas[i]] == "1" && endemic_sp ==TRUE ){ # Checking if species is found on island
          sp_on_island <- append(sp_on_island,tree$tip.label[ntip]) # Saving tips
        }
      }
    }
    edges_per_island[[i]] <- edges_with_only_island_decen # Saving all edges where all descendants are endemic to an island
    species_per_island[[i]] <- sp_on_island # Saving all colonizing Endemic species found on island
    endemics_from_radiation[[i]] <- decendants_of_edges # Saving all species from edges that radiated
    
  }
  
  names(edges_per_island) <- areas
  names(species_per_island) <- areas
  names(endemics_from_radiation) <- areas
  
  comb <- list()
  for(i in 1:length(areas)){
    comb[[i]] <- list(edges_per_island[areas[i]], endemics_from_radiation[areas[i]], species_per_island[areas[i]])
    names(comb[[i]][[1]]) <- "Radiating Edges"
    names(comb[[i]][[2]]) <- "Endemics from Radiating edges"
    names(comb[[i]][[3]]) <- "Endemics from Colonization"
  }
  names(comb) <- areas
  return(comb)
}


# Endemic lineages and species on islands based on orthologs tree and community matrix with all species on islands
output_community_matrix_islands <- FindEndemicLineages(cory_phylo,community_matrixes,islands_cory)

#Can I create a database consisting of the islands and then the characteristics of each island.
#Potential Columns I want for my database.
  # Number of colonization edges (Single species and Edges)
    #Edges should ideally only be for endemic species.(This needs to be written into the function, currently it uses all species present on the island)
    # Maximum number and Min number? 
  # Number of edges leading to a radiation
  # Number of endemic species resulting from radiation
  # Number of species colonizing but not radiating

output_as_df <- as.data.frame(t(as.data.frame(do.call(cbind, output_community_matrix_islands))))
output_as_df_lengths <- output_as_df

for (i in 1:length(output_as_df[[1]])){
  for (k in 1:ncol(output_as_df)){
    output_as_df_lengths[[k]][[i]] <- (length(output_as_df[[k]][[i]][[1]]))
  }
}

colnames(output_as_df_lengths) <- c("No. edges", "No. cladogenesis", "No. anagenesis")
output_as_df_lengths <- rownames_to_column(output_as_df_lengths, var="LEVEL3_NAM")


#Editing Gift dataset to make Haiti and Dominican republic one botanical country.
# Adding the Areas together
gift_data[which(gift_data$LEVEL3_COD == "DOM"),8] <- gift_data[which(gift_data$LEVEL3_COD == "DOM"),8] + gift_data[which(gift_data$LEVEL3_COD == "HAI"),8]
# Selecting the lowest distance to mainland
gift_data[which(gift_data$LEVEL3_COD == "DOM"),9] <- min(gift_data[which(gift_data$LEVEL3_COD == "HAI" | gift_data$LEVEL3_COD == "DOM"),9]) 
# Calculating the mean maximum elevation for the entire island.
gift_data[which(gift_data$LEVEL3_COD == "DOM"),42] <- gift_data[which(gift_data$LEVEL3_COD == "HAI"),42]*(gift_data[which(gift_data$LEVEL3_COD == "HAI"),8]/gift_data[which(gift_data$LEVEL3_COD == "DOM"),8]) + gift_data[which(gift_data$LEVEL3_COD == "DOM"),42]*(1-gift_data[which(gift_data$LEVEL3_COD == "HAI"),8]/gift_data[which(gift_data$LEVEL3_COD == "DOM"),8])


output_with_island_cory <- merge(output_as_df_lengths, gift_data, by="LEVEL3_NAM")
output_with_island_cory[,2:4] <- lapply(output_with_island_cory[,2:4], unlist)

# Adding the geological origin
# Loading Dataset with Geological Origin
island_data_origin <- read_csv(file.path(datadir,"/Islands_TDWG_AllData.txt"))

#Renaming column to be the same as in the gift database
names(island_data_origin)[1] <- "LEVEL3_COD"

#selecting only the geological origin and the name
island_data_origin_merge <- island_data_origin[,c(1,6)]

#Merging geological origin with island counts data. 
output_with_island_cory <- merge(output_with_island_cory,island_data_origin_merge, by="LEVEL3_COD")


# The venezuelan antilles are not in the Geological origin dataset
#which(!(output_as_df_all_sp[[1]][,1] %in% island_data_origin_merge$LEVEL3_COD))
#print(output_as_df_all_sp[[1]][101,1])


p <- ggplot(output_with_island_cory, aes(area, dist))
p + geom_point(aes(size=`No. cladogenesis`+`No. anagenesis`), ) +
  labs(size ="Endemic Sp", title = "Cladogenesis + Anagenesis") +
  xlab("Island Area") +
  ylab("Distance to Mainland") +
  scale_size_area(limits = c(1,50),
                  breaks = c(1,2,5,10,25,50),labels =c("1","2","<5","<10","<25","<50")) +
  ggrepel::geom_text_repel(aes(label=ifelse((`No. cladogenesis`+`No. anagenesis`)>2,as.character(LEVEL3_NAM),'')),hjust=1.2, size=3) +
  theme_classic()


p1_cory <- ggplot(output_with_island_cory, aes(area, dist)) +
  geom_point(aes(size = ifelse(`No. cladogenesis` >= 2, `No. cladogenesis`, 1.2), color=GeologicalOrigin, shape = ifelse(`No. cladogenesis` >= 2, "Shape1", "Shape2"))) +
  scale_colour_manual(values = met.brewer("Pillement", 4)) +
  scale_shape_manual(values = c("Shape1" = 16, "Shape2" = 4)) +
  labs(size ="Endemic Sp", title = "Cladogenesis") +
  xlab("Island Area (Km²)") +
  ylab("Distance to Mainland (Km)") +
  scale_size(range=c(1,8), limits = c(1,33), breaks = c(1, 10, 20, 30)) +
  #ggrepel::geom_text_repel(aes(label=ifelse(`No. cladogenesis`>2,as.character(LEVEL3_NAM),'')), size=4) +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c(100,1000,10000,100000,1000000), limits = c(100,856017.9)) +
  scale_y_continuous(breaks = c(10,50,100,250,500,1000,2000,4000,8000), limits = c(1,8000)) +
  guides(size = guide_legend("Endemic species."), colour="none", shape="none") +
  coord_trans(x="log10", y="log10")

p1_cory


p2_cory <- ggplot(output_with_island_cory, aes(area, dist)) +
  geom_point(aes(size = ifelse(`No. anagenesis`+`No. edges` >= 1, `No. anagenesis`+`No. edges`, 1.2), color=GeologicalOrigin, shape = ifelse(`No. anagenesis`+`No. edges` >= 1, "Shape1", "Shape2"))) +
  scale_colour_manual(values = met.brewer("Pillement", 4)) +
  labs(size ="Endemic Sp", title = "Anagenesis") +
  scale_shape_manual(values = c("Shape1" = 16, "Shape2" = 4)) +
  xlab("Island Area (Km²)") +
  ylab("Distance to Mainland (Km)") +
  scale_size(range=c(1,8), limits = c(1,33)) +
  #ggrepel::geom_text_repel(aes(label=ifelse(`No. anagenesis`>5,as.character(LEVEL3_NAM),'')),hjust=1.2, size=4) +
  theme_classic() +
  scale_x_continuous(breaks = c(100,1000,10000,100000,1000000), limits = c(100,856017.9)) +
  scale_y_continuous(breaks = c(10,50,100,250,500,1000,2000,4000,8000), limits = c(1,8000)) +
  guides(size = "none", colour="none", shape = "none") +
  coord_trans(x="log10", y="log10")

#theme(legend.position = "none") +
#guides(size = "none", colour=guide_legend("Geological Origin")) +

p2_cory

cory_scatter <- cowplot::plot_grid(p2_cory,p1_cory+theme(legend.position = "none"), rel_widths = c(1,1))
legend_bottom_cory <- get_legend(p1_cory)
cory_scatter <- cowplot::plot_grid(cory_scatter,legend_bottom_cory, ncol=1,nrow = 2, rel_heights = c(1,.1))

# now add the title
title_cory <- cowplot::ggdraw() + 
  draw_label(
    "Coryphoideae",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 40)
  )
cory_scatter <- plot_grid(
  title_cory, cory_scatter,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

cory_scatter










p4 <- ggplot(output_with_island_cory, aes(log10(area), log10(dist)))
p4 + geom_point(aes(size=`No. anagenesis`+`No. edges`, color=GeologicalOrigin)) +
  labs(size ="No. Colonizations", title = "Possible number of Colonizations from Phylogeny") +
  xlab("Island Area") +
  ylab("Distance to Mainland") +
  scale_size_area(limits = c(1,35),
                  breaks = c(1,2,6,10,30)) +
  geom_text(aes(label=ifelse(`No. anagenesis`+`No. edges`>3,as.character(LEVEL3_NAM),'')),hjust=0.8 ,vjust=1.8, size=3) +
  theme_classic2()



p3 <- ggplot(output_with_island_cory, aes(area, `No. cladogenesis` + `No. anagenesis`))
p3 + geom_point() +
  labs(size ="Endemic Sp", title = "Species area Relationship") +
  xlab("Island Area") +
  ylab("Species Number")


variables <- colnames(output_with_island_cory)
variables <- variables[11:49]

for (i in variables){
  px <- ggplot(output_with_island_cory, aes(output_with_island_cory[,i],`No. cladogenesis` + `No. anagenesis`)) +
    geom_point() +
    ylab("Species Richness") +
    xlab(i) +
    geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
    theme_classic2()
  print(px)
}

model <- glm(`No. cladogenesis`/(`No. cladogenesis`+`No. anagenesis`) ~ dist + log10(area)+mean_mx30_grd,family = gaussian, output_with_island_cory)

summary(model)


px1 <- ggplot(output_with_island_cory, aes(area,`No. cladogenesis`/(`No. anagenesis`+`No. cladogenesis`))) +
  geom_point() +
  ylab("Proportion of Endemics from Cladogenesis") +
  xlab("Area") +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
    ggrepel::geom_text_repel(aes(label = ifelse(`No. cladogenesis`/(`No. anagenesis`+`No. cladogenesis`)>0,as.character(LEVEL3_NAM),'')),hjust=1.2, size=3) +
  ylim(0,1) +
  theme_classic2()
px1


Now we can do the same but with the Cladogenesis and Anagenesis


# for (i in variables){
#   plot_dat <- output_with_island_cory %>% select('No. cladogenesis','No. anagenesis',!!rlang::sym(i)) %>% pivot_longer(-!!rlang::sym(i))
#   COLS = c("red","turquoise")
#   names(COLS) = c("No. cladogenesis","No. anagenesis")
# 
#   py <- ggplot(plot_dat)+
#     geom_point(aes(!!rlang::sym(i),value,colour=name,fill=name))+
#     geom_smooth(aes(!!rlang::sym(i),value,colour=name,fill=name), method = lm)+
#     scale_fill_manual(name="Speciation Process", values = COLS) +
#     scale_color_manual(name="Speciation Process", values = COLS) +
#     ylab("Number of Species") +
#     theme_classic2()
#   print(py)
# }




Okay So there is not really any of these that show any nice correlation between number of species and different climatic variables.
Can we instead investigate the phylogeny a little deeper.
I will now use GGtree to produce some better pictures of the phylogenies to show how the dispersal has been onto the islands


# cuba_edges <- output_as_df["Cuba",][[1]]
# cuba_sp <- output_as_df["Cuba",][[3]]
# 
# # 
# # test_tree <- groupClade(cory_phylo,cuba_edges[[1]][[1]],group_name = "cladogenesis edges") 
# # 
# # 
# # p_tree <- ggtree(cory_phylo)
# # p_tree <- viewClade(p_tree, MRCA(p_tree,"Coccothrinax clarensis", "Copernicia hospita")) +
# #   collapse(p_tree,node = 610)
# # p_tree
# # 
# # 
# # MRCA(p_tree,"Arenga ryukyuensis", "Tahina spectabilis") #420
# # MRCA(p_tree, "Copernicia tectorum", "Copernicia fallaensis") #571
# # MRCA(p_tree, "Coccothrinax macroglossa", "Coccothrinax gracilis") #532
# # MRCA(p_tree, "Hemithrinax ekmaniana", "Hemithrinax compacta") #530
# # 
# cory_phylo_group <- groupClade(cory_phylo,cuba_edges[[1]][[1]], group_name = "cladogenesis_edges") 
# # cory_phylo_group <- groupClade(cory_phylo_group, cuba_sp, group_name = "anagenesis")
# # ?groupClade()
# # ptree2
# # 
# # cuba_sp
# 
# #################################
# # I should find all edges and just color them using the "simpler" method
# all_edges_cuba <- c()
# for (i in 1:length(cuba_edges$Cuba$`Radiating Edges`)){
#   all_edges_cuba<- append(all_edges_cuba,getDescendants(cory_phylo,cuba_edges$Cuba$`Radiating Edges`[i])) #111 112
# }
# 
# for (i in 1:length(cuba_sp$Cuba$`Endemics from Colonization`)){
#   print(cuba_sp$Cuba$`Endemics from Colonization`[i])
#   all_edges_cuba <- append(all_edges_cuba,which(cory_phylo$tip.label == cuba_sp$Cuba$`Endemics from Colonization`[i]))
# }
# 
# which(cory_phylo$tip.label == "Coccothrinax miraguama") # I need to do this as this will find the edge which needs to be coloured
# 
# edgecol <- data.frame(node=1:(cory_phylo$Nnode+length(cory_phylo$tip.label)), color="Black")
# 
# edgecol[which(edgecol$node %in% all_edges_cuba),2] <- "Red"
# 
# 
# ###############################
# 
# ptree3 <- ggtree(cory_phylo_group) %<+% edgecol + aes(color=I(color))
# ptree3 <- collapse(ptree3,node = 610) + geom_point2(aes(subset=(node==610)), shape=18,size=2)
# #ptree3 <- collapse(ptree3,node = 420) + geom_point2(aes(subset=(node==420)), shape=18,size=2)
# 
# ptree3 <- viewClade(ptree3, MRCA(ptree3,"Copernicia glabrescens", "Coccothrinax macroglossa")) +
#   geom_cladelab(node = 571, label = "Copernicia") +
#   geom_cladelab(node=532, label = "Coccothrinax") +
#   geom_cladelab(node=530, label = "Hemithrinax") +
#   labs(FALSE)
#   #scale_color_manual(group=c("Black,","Orange"))
# 
# #ptree3  




Now I need to find the total number of palm species on each of the islands that I am investigating.
I should be able to do this using the world checklist and the checklist distribution.
I should start by making a list of all the accepted palm species names 


apalms # a list of all the accepted palm species (subspecies and varieties)

dist_data # the data for the distribution of all species

palm_dist_data <- dist_data[which(dist_data$plant_name_id %in% apalms$plant_name_id),] # dist data for all palm species.

#So now I need to look through the islands that I have and count the total number of species on the islands
output_as_df_lengths

total_palm_sp <- data.frame()

for (i in 1:length(output_as_df_lengths$LEVEL3_NAM)){
  print(output_as_df_lengths$LEVEL3_NAM[i])
  print(length(which(palm_dist_data$area == output_as_df_lengths$LEVEL3_NAM[i])))
  total_palm_sp[i,1] <- output_as_df_lengths$LEVEL3_NAM[i]
  total_palm_sp[i,2] <- length(which(palm_dist_data$area == output_as_df_lengths$LEVEL3_NAM[i]))
}

names(total_palm_sp) <- c("LEVEL3_NAM", "Total_sp")

output_with_island_cory <- merge(output_with_island_cory, total_palm_sp, by="LEVEL3_NAM")


# I also want to add the total number of Coryphoid species just to see if an island has Coryphoid species present
dist_data_coryphoids <- dist_data[which(dist_data$plant_name_id %in% cory_data_accepted$plant_name_id),]
total_coryphoid_sp <- data.frame()

for (i in 1:length(output_as_df_lengths$LEVEL3_NAM)){
  print(output_as_df_lengths$LEVEL3_NAM[i])
  print(length(which(dist_data_coryphoids$area == output_as_df_lengths$LEVEL3_NAM[i])))
  total_coryphoid_sp[i,1] <- output_as_df_lengths$LEVEL3_NAM[i]
  total_coryphoid_sp[i,2] <- length(which(dist_data_coryphoids$area == output_as_df_lengths$LEVEL3_NAM[i]))
}

names(total_coryphoid_sp) <- c("LEVEL3_NAM", "Total_sp_coryphoids")

output_with_island_cory <- merge(output_with_island_cory, total_coryphoid_sp, by="LEVEL3_NAM")


p7 <- ggplot(output_with_island_cory, aes(`No. cladogenesis`+`No. anagenesis`, Total_sp))
p7 + geom_point() +
  labs(title = "Number of Coryphoideae species compared to total Arecaceae numbers") +
  xlab("No. Coryphoideae") +
  ylab("No. Arecaceae") +
  geom_text(aes(label=ifelse(`No. anagenesis`+`No. cladogenesis`>2,as.character(LEVEL3_NAM),'')),hjust=0.8 ,vjust=1.8, size=3) +
  geom_smooth(method = lm) +
  theme_classic()



#sp_area_dist_ap <- output_with_island_cory %>%  pivot_longer(cols = ,)

p8 <- ggplot(output_with_island_cory, aes(area, dist))
p8 + geom_point(aes(size=Total_sp), ) + 
  geom_point(aes(size=`No. cladogenesis` + `No. anagenesis`), color = "red") +
  labs(size ="Species", title = "Species Number compared to Area and Isolation") +
  xlab("Island Area") +
  ylab("Distance to Mainland") +
  scale_size_area(limits = c(1,317),
                  breaks = c(10,50,100,150,200,300)) +
  ggrepel::geom_text_repel(aes(label = ifelse(Total_sp>50,as.character(LEVEL3_NAM),'')),hjust=1.2, size=3) +
  theme_classic()



p9 <- ggplot(output_with_island_cory, aes(area, dist))
p9 + geom_point(aes(size=Total_sp_coryphoids), ) + 
  geom_point(aes(size=`No. cladogenesis` + `No. anagenesis`), color = "red") +
  labs(size ="Species", title = "Species Number compared to Area and Isolation") +
  xlab("Island Area") +
  ylab("Distance to Mainland") +
  scale_size_area(limits = c(1,93),
                  breaks = c(10,20,30,40,50,90)) +
  ggrepel::geom_text_repel(aes(label = ifelse(Total_sp>50,as.character(LEVEL3_NAM),'')),hjust=1.2, size=3) +
  theme_classic()




p10 <- ggplot(output_with_island_cory, aes(`No. cladogenesis`+`No. anagenesis`, Total_sp_coryphoids))
p10 + geom_point() +
  labs(title = "Number of Endemic Coryphoideae species compared to total Coryphoid species numbers") +
  xlab("No. Endemic Coryphoideae") +
  ylab("No. Coryphoideae") +
  geom_text(aes(label=ifelse(`No. anagenesis`+`No. cladogenesis`>2,as.character(LEVEL3_NAM),'')),hjust=0.8 ,vjust=1.8, size=3) +
  geom_smooth(method = lm) +
  theme_classic()




p11 <- ggplot(output_with_island_cory, aes(Total_sp_coryphoids,Total_sp))
p11 + geom_point() +
  labs(title = "Number of Coryphoideae species compared to total Arecaceae species numbers") +
  xlab("No. Coryphoideae") +
  ylab("No. Arecaceae") +
  geom_text(aes(label=ifelse(`No. anagenesis`+`No. cladogenesis`>2,as.character(LEVEL3_NAM),'')),hjust=0.8 ,vjust=1.8, size=3) +
  geom_smooth(method = lm) +
  theme_classic()

I think I have all the species data I need for now. But I would really like to get some data on how disjunct the island botanical countries are in order to do some tests which should show the impact of connectivity within the island botanical countries

# st_convex_hull() This is the function I want to use
# For each of the island botanical countries I should loop through them and calculate the proportion of the convex hull which is land area

#Loading the botanical country dataframe
botanical_countries <- st_read(dsn = file.path("../data","wgsrpd-master/level3"), layer = "level3")

#Now I need to create a subset of the botanical countries which are the islands.
botanical_countries_islands_cory <- botanical_countries[which(botanical_countries$LEVEL3_NAM %in% islands_cory),]

# #Now I need to create a for loop and a merge.
# island_area_chull <- data_frame()
# 
# for (i in 1:length(botanical_countries_islands[[1]])){
#   chul <- st_convex_hull(botanical_countries_islands[i,])
#   area <- st_area(chul)
#   area <- units::set_units(area,km^2)
#   print(area)
#   island_area_chull[i,1] <- botanical_countries_islands[[1]][i]
#   island_area_chull[i,2] <- area
# }
# #Setting the names of the tibble
# names(island_area_chull) <- c("LEVEL3_NAM", "convex_hull_area")
# 
# #Now I should be able to merge the two datasets
# output_with_island_cory <- merge(output_with_island_cory, island_area_chull, by="LEVEL3_NAM")
# 
# output_with_island_cory$area_prop <- output_with_island_cory$area/output_with_island_cory$convex_hull_area


Question is still how does this area proportion scale?
Some islands might have a really big proportion because the island only consists of a single small island while others might have a small proportion because it consists of a single big island and then 3 or 4 really distant small islands.

Therefore it would probably be a good idea to scale this proportion to something with the area
I think we could just scale it to the overall area, but this does not discern between 4 semi big islands which are medium isolated and a really big island with a few small satellite islands.

I think a better approach would be to use the area of the biggest island in the Island Botanical Country

sf::sf_use_s2(FALSE)

biggest_island <- data.frame()

for (k in 1:length(botanical_countries_islands_cory[[1]])) {
  k_area <- c()
  test_cast <- st_cast(botanical_countries_islands_cory[k,1], "POLYGON")
  for (i in 1:length(test_cast[[1]])) {
    k_area <- append(k_area,st_area(test_cast[i,]))
    k_area <- units::set_units(k_area,km^2)
  }
  biggest_island[k,1] <- botanical_countries_islands_cory[[1]][k]
  biggest_island[k,2] <- max(k_area)
  print(botanical_countries_islands_cory[[1]][k])
  print(max(k_area))
}

# #Setting the names of the tibble
# names(biggest_island) <- c("LEVEL3_NAM", "biggest_island_area")
# 
# output_with_island_cory <- merge(output_with_island_cory, biggest_island, by="LEVEL3_NAM")
# output_with_island_cory <- units::drop_units(output_with_island_cory)
# 
# convex_hull_area_standardized_to_big_island <- (output_with_island_cory$area_prop*biggest_island$biggest_island_area)

# p12 <- ggplot(output_with_island_cory, aes(log10(convex_hull_area_standardized_to_big_island),Total_sp))
# p12 + geom_point(aes(log10(convex_hull_area_standardized_to_big_island),Total_sp)) +
#   #labs(title = "Number of Total sp compared to (proportion land of chull)/Area of Biggest Island") +
#   #xlab("-Log10(proportion land of chull/Area of Biggest Island)") +
#   #ylab("No. Arecaceae") +
#   #geom_text_repel(aes(label=ifelse(Total_sp>1,as.character(LEVEL3_NAM),'')), size=3) +
#   theme_classic()


Through a conversation with Erick Kush, I have come to the conclusion that I should just use the mean centroid distance between all the islands, and IF i want to do something with the size of the island at the same time I should probably multiply each distance with the size of the source island as this would mean that you capture large distances with large islands better (Dividing by size would mean you could get same result from small not so isolated islands as you would with large very isolated islands)

So goal is now to calculate the mean distance for a distance matrix of all island botanical countries in order to capture the disjunct nature of these Island botanical countries. 


#Calculate distances between all island from island border and centroid
island_distances <- data.frame()

for (k in 1:length(botanical_countries_islands_cory[[1]])) {
  test_cast <- st_cast(botanical_countries_islands_cory[k,1], "POLYGON")
  dist_bord <- st_distance(test_cast)
  test_cast_cent <- st_centroid(test_cast)
  dist_cent <- st_distance(test_cast_cent)
  island_distances[k,1] <- botanical_countries_islands_cory[[1]][k]
  island_distances[k,2] <- mean(dist_bord[upper.tri(dist_bord)])
  island_distances[k,3] <- mean(dist_cent[upper.tri(dist_cent)])
}

#Setting the names of the tibble
names(island_distances) <- c("LEVEL3_NAM", "mean_border_dist", "mean_centroid_dist") # the unit is in meters so maybe I should convert it to km

island_distances[,2:3] <- island_distances[,2:3]/1000

output_with_island_cory <- merge(output_with_island_cory, island_distances, by="LEVEL3_NAM")


These numbers have a problem because the mean distance is really high for some island botanical countries due to their total size and extent, this is the example for Borneo and Hawaii. A quick look at the different island botanical countries shows how these numbers might be affected by size of islands.

# for (k in 1:length(botanical_countries_islands[[1]])) {
#   plot(botanical_countries_islands[k,1], main = botanical_countries_islands[[1]][k])
# }


So the idea is, maybe I should standardize the distance between all islands based on the size of the island being measured from?
This would mean that small far away islands will have less of an impact.

Problem is that distance is pairwise and in the previous method i just took the distances in the upper triangle as there were no differences if the distance was from Island A to B or B to A. But if you somehow standardize for the size of the island were measuring from, suddenly the distance from A to B can be very different to the distance from B to A.

I therefore need to figure out if I should always measure the distance between two islands using the smallest island or the biggest island.
I could probably do both and see how it affects the results

Unfortunately using the previous approach we either lose the area of the smallest or biggest island.
Can we somehow figure out a method that also incorporates the size of the largest island. maybe if we standardize the distance based on the size of the island compared to the total size of the island botanical country.


#First we need to calculate the distances between the islands of the island botanical country and then divide each row in the distance matrix with the area of the island of that row.

# island_scaled_distances <- data.frame()
# 
# for (k in 1:length(botanical_countries_islands[[1]])) { # Looping through island botanical countries
#   print(c("starting on ",botanical_countries_islands[[1]][k]))
#   k_area <- c()
#   test_cast <- st_cast(botanical_countries_islands[k,1], "POLYGON") # Changing the shape file of a island botanical country to several polygons
#   for (i in 1:length(test_cast[[1]])) { # Looping through polygons
#     k_area <- append(k_area,as.numeric(st_area(test_cast[i,]))/1000000) # Calculating the area of each single polygon
#   }
#   
#   dist_bord <- st_distance(test_cast) #Calculating the distances between the borders of each polygon for each island botanical country
#   test_cast_cent <- st_centroid(test_cast) # Finding the centroid for each polygon
#   dist_cent <- st_distance(test_cast_cent) # Finding the distances between all centroids
#   
#   #Now I want to multiply the values in each row with the areas of the island in the row
#   for (l in 1:length(k_area)){
#     dist_bord[l,] <- dist_bord[l,]*k_area[l]
#     dist_cent[l,] <- dist_cent[l,]*k_area[l]
#   }
# 
#   print(c("Area vector for islands of ", botanical_countries_islands[[1]][k]))
#   print(k_area)
#   if (length(k_area)==1){
#   island_scaled_distances[k,1] <- botanical_countries_islands[[1]][k]
#   island_scaled_distances[k,2] <- NaN
#   island_scaled_distances[k,3] <- NaN
#   island_scaled_distances[k,4] <- NaN
#   island_scaled_distances[k,5] <- NaN
#   } else {
#     #Now that all the rows in the distance matrix is standardized to island size I will need to find out which of the islands in each A-B comparison is the large and which one is the smallest. I can do this if i transpose the matrix as this mean I can compare cell [1,2] in one data set with [2,1] in the other data set by only looking at the [1,2] cells
#   dist_bord_t <- t(dist_bord)
#   dist_cent_t <- t(dist_cent)
#   
#   ##########################################################
#   ## Adding the distances based on area of biggest island ##
#   ##########################################################
#   
#   island_dist_scaled_bord_big <- as.data.frame(matrix(0,length(dist_bord[1,]),length(dist_bord[,1])))
#   island_dist_scaled_cent_big <- data.frame(matrix(0,length(dist_bord[1,]),length(dist_bord[,1])))
#   
#   
#   for (g in 1:(length(dist_bord[,1])-1)){ # looping through rows
#     for (h in (g+1):length(dist_bord[1,])){ # looping through columns in the upper triangle
#       if (k_area[g] > k_area[h]) { # Checking if the area of island g is bigger than island h
#         island_dist_scaled_bord_big[g,h] <- dist_bord[g,h] # if island measured from is bigger add the scaled distance to the distance matrix
#       } else {
#         island_dist_scaled_bord_big[g,h] <- dist_bord_t[g,h] # else if the other island is bigger add the distance scaled to that islands size.
#       }
#     }
#   }
# 
#   print(c("Done dist_bord_big for ",botanical_countries_islands[[1]][k]))
#   
#   for (g in 1:(length(dist_cent[,1])-1)){ # Looping through rows
#     for (h in (g+1):length(dist_cent[1,])){ # Looping through columns in the upper triangle
#       if (k_area[g] > k_area[h]) { # Checking which of the islands is the biggest 
#         island_dist_scaled_cent_big[g,h] <- dist_cent[g,h] # If island measured from is biggest add the scaled distance to the distance matrix
#       } else {
#         island_dist_scaled_cent_big[g,h] <- dist_cent_t[g,h] # Else add the scaled distance from the other island to the measured island.
#       }
#     }
#   }
#   
#   print(c("Done dist_cent_big for ",botanical_countries_islands[[1]][k]))
#   ###########################################################
#   ## Adding the distances based on area of smallest island ##
#   ###########################################################
#   
#   island_dist_scaled_bord_small <- as.data.frame(matrix(0,length(dist_bord[1,]),length(dist_bord[,1])))
#   island_dist_scaled_cent_small <- as.data.frame(matrix(0,length(dist_bord[1,]),length(dist_bord[,1])))
#   
#   for (g in 1:(length(dist_bord[,1])-1)){ # looping through rows
#     for (h in (g+1):length(dist_bord[1,])){ # looping through columns in the upper triangle 
#       if (k_area[g] < k_area[h]) { # Checking if the area of island g is bigger than island h
#         island_dist_scaled_bord_small[g,h] <- dist_bord[g,h] # if island measured from is bigger add the scaled distance to the distance matrix
#       } else {
#         island_dist_scaled_bord_small[g,h] <- dist_bord_t[g,h] # else if the other island is bigger add the distance scaled to that islands size.
#       }
#     }
#   }
#   
#   for (g in 1:(length(dist_cent[,1])-1)){ # Looping through rows
#     for (h in (g+1):length(dist_cent[1,])){ # Looping through columns in the upper triangle
#       if (k_area[g] < k_area[h]) { # Checking which of the islands is the biggest 
#         island_dist_scaled_cent_small[g,h] <- dist_cent[g,h] # If island measured from is biggest add the scaled distance to the distance matrix
#       } else {
#         island_dist_scaled_cent_small[g,h] <- dist_cent_t[g,h] # Else add the scaled distance from the other island to the measured island.
#       }
#     }
#   }
#   
#   # So now I should have 4 dataframes only filled in the upper triangle with the scaled distances between island borders and centroids
#   # I can then calculate a mean from these and add them to my dataframe. 
#   
#   island_scaled_distances[k,1] <- botanical_countries_islands[[1]][k]
#   island_scaled_distances[k,2] <- mean(island_dist_scaled_bord_big[upper.tri(island_dist_scaled_bord_big)])
#   island_scaled_distances[k,3] <- mean(island_dist_scaled_bord_small[upper.tri(island_dist_scaled_bord_small)])
#   island_scaled_distances[k,4] <- mean(island_dist_scaled_cent_big[upper.tri(island_dist_scaled_cent_big)])
#   island_scaled_distances[k,5] <- mean(island_dist_scaled_cent_small[upper.tri(island_dist_scaled_cent_small)])
#   }
# }
# 
# 
# #Setting the names of the tibble
# names(island_scaled_distances) <- c("LEVEL3_NAM", "mean_border_dist_scaled_big", "mean_border_scaled_small","mean_centroid_dist_scaled_big","mean_centroid_scaled_small") # the unit is in meters so maybe I should convert it to km
# 
# island_scaled_distances[,2:5] <- island_scaled_distances[,2:5]/1000
# 
# output_with_island_cory <- merge(output_with_island_cory, island_scaled_distances, by="LEVEL3_NAM")





After talking with Wolf I think that I should also calculate the mean nearest neighbour distance for all the islands and then scale all those distances based on the size of both islands.

# This could be done using the distance matrix for all the islands
# Then I should just select the lowest value for each row which is larger than 0. as this would be the nearest neighbour
# In order to scale this nearest neighbor distance with the area of the islands you should just multiply the distance with the area of island[row nr] and island[column nr]

near_neighbour_dist <- data.frame()

for (k in 1:length(botanical_countries_islands_cory[[1]])) { # Looping through island botanical countries
  print(c("starting on ",botanical_countries_islands_cory[[1]][k]))
  k_area <- c()
  test_cast <- st_cast(botanical_countries_islands_cory[k,1], "POLYGON") # Changing the shape file of a island botanical country to several polygons
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
  
  print(mean(unlist(near_neighbour_dist_bord_scaled)))
  print(mean(unlist(near_neighbour_dist_cent_scaled)))
  
  near_neighbour_dist[k,1] <- botanical_countries_islands_cory[[1]][k]
  near_neighbour_dist[k,2] <- mean(unlist(near_neighbour_dist_bord_scaled))
  near_neighbour_dist[k,3] <- mean(unlist(near_neighbour_dist_cent_scaled))
  near_neighbour_dist[k,4] <- mean(near_neighbour_dist_bord)
  near_neighbour_dist[k,5] <- mean(near_neighbour_dist_cent)
}


names(near_neighbour_dist) <- c("LEVEL3_NAM", "nearest_neighbour_distance_border_scaled", "nearest_neighbour_distance_centroid_scaled", "neareast_neighbour_bord", "neareast_neighbour_cent") # the unit is in meters so maybe I should convert it to km

output_with_island_cory <- merge(output_with_island_cory, near_neighbour_dist, by="LEVEL3_NAM")
near_neighbour_dist
output_with_island_cory
order(output_with_island_cory$area)

We also need to load the dataset from Melanie about Island Origin

# # # Loading Dataset with Geological Origin
# island_data_origin <- read_csv(file.path("../data","/Islands_TDWG_AllData.txt"))
# 
# #Renaming column to be the same as in the gift database
# names(island_data_origin)[1] <- "LEVEL3_COD"
# 
# #Merging some names from gift to geological origin
# island_data_origin <- merge(island_data_origin, gift_data[,2:3], by="LEVEL3_COD")
# 
# #selecting only the geological origin and the name
# island_data_origin_merge <- island_data_origin[,c(6,25)]
# 
# #Merging geological origin with island counts data.
# output_with_island_cory <- merge(output_with_island_cory, island_data_origin_merge, by="LEVEL3_NAM")
# output_with_island_cory


output_with_island_cory_quasi_subset <- output_with_island_cory[which(output_with_island_cory$`No. cladogenesis`+output_with_island_cory$`No. anagenesis` != 0),]
output_with_island_cory_quasi_subset$GeologicalOrigin <- as.factor(output_with_island_cory_quasi_subset$GeologicalOrigin)



# So I am investigating the effect of overdispersion on the statistical models.
# Furthermore I need to change all my models so we try to investigate how the different island characteristics affect cladogenesis and anagenesis.
# The models I want to investigate are cladogenesis-edges ~ island characteristics & anagenesis+edges ~island characteristics

# Lets try to run a poisson model for both of these things.
# Because Dredge fails if we have NA's in our dataset, and we do have NA's in the nearest neighbour stat, I will create a different dataframe where the 
# NA's are changed to 0's because these botanical countries consists of a single island and therefore setting the fragmentation to 0 is justified. 



output_with_island_cory$nearest_neighbour_distance_border_scaled[is.na(output_with_island_cory$nearest_neighbour_distance_border_scaled)] <- 0
output_with_island_cory_no_christmas_or_jamaica <- output_with_island_cory[-(which(output_with_island_cory$LEVEL3_NAM == "Jamaica" | output_with_island_cory$LEVEL3_NAM =="Christmas I.")),]
  
# Cladogenesis model
glm_clad_pois_all <- glm(`No. cladogenesis`-`No. edges`~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled)+GeologicalOrigin,
                    data=output_with_island_cory_no_christmas_or_jamaica, na.action = na.fail,family = quasipoisson(link='log'))

dp_glm_clad_pois_all = sum(residuals(glm_clad_pois_all,type ="pearson")^2)/glm_clad_pois_all$df.residual # Calculating dispersion

summary(glm_clad_pois_all, dispersion = dp_glm_clad_pois_all) # Only area and isolation is of significant importance and they both have a positive effect on the proportion of endemics

# Anagenesis model

glm_ana_pois_all <- glm(`No. anagenesis`+`No. edges` ~
                      mean_mx30_grd+area+dist+nearest_neighbour_distance_border_scaled+GeologicalOrigin,
                    data=output_with_island_cory,
                    na.action = na.fail,
                    family = poisson(link='log'))

dp_glm_ana_pois_all = sum(residuals(glm_ana_pois_all,type ="pearson")^2)/glm_ana_pois_all$df.residual

summary(glm_ana_pois_all, dispersion = dp_glm_ana_pois_all) 

#Log10
glm_ana_pois_all_log10 <- glm(`No. anagenesis`+`No. edges` ~
                      log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled)+GeologicalOrigin,
                    data=output_with_island_cory_no_christmas_or_jamaica,
                    na.action = na.fail,
                    family = quasipoisson(link='log'))

dp_glm_ana_pois_all_log10 = sum(residuals(glm_ana_pois_all,type ="pearson")^2)/glm_ana_pois_all$df.residual

summary(glm_ana_pois_all_log10, dispersion = dp_glm_ana_pois_all_log10)







Because our variance is overdispersed for the cladogenesis data we have to resort to a quasipoisson model.
Although quasipoisson models have no AIC because they have no likelihood there has been some workarounds done.
My solution is to use the AIC of the poisson model and then use the rest of the parameters from the quasipoisson model in the dredge function from MuMIn.


This is all done Based on the vignette PDF by Ben Bolker (https://cran.r-project.org/web/packages/bbmle/vignettes/quasi.pdf)
and the further work by Stacy DeRuiter in her  STAT 245, Advanced Data Analysis course at Calvin University. (https://stacyderuiter.github.io/s245-notes-bookdown/regression-for-count-data.html)  section 6.7.3

This method of using the AIC value of the non quasi model requires some setting up and that is what we will be doing here.


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
qdredge <- function(model, family='x.quasipoisson', na.action=na.fail, chat = dfun(model), rank='QAIC', ...){
  model2 <- update(model, family=family, na.action=na.action)
  (dt <- dredge(model2, rank=rank, chat=chat, ...))
}


Now lets run this qdredge function on our cladogenesis models


#Cladogenesis
#dredge(glm_clad_cory)
#qdredge(glm_clad_cory)


Here we can see that the model using all variables and the model only 



And compare them with the results from the dredge function on the anagenesis data.

#Anagenesis
qdredge(glm_ana_pois_all_log10)



dredge(glm_ana_pois_all_log10)




Now we have the best models, and I would therefore like to make a coefficients plot for the coefficients of the best models.
To start with I will use broom::tidy in order to get the estimates.

Something is wrong with my estimates of the importance of geological origin. The standard errors are massive for my cladogenesis model.
This means I cannot get reliable confidence intervals because they approach 





#Getting the coefficient estimates for cladogenesis
results_ana <- broom::tidy(glm_ana_pois_all)

#Now lets find the confidence intervals for these models.
fit_ana_90 <- confint(glm_ana_pois_all, level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")

results_ana <- dplyr::bind_cols(results_ana,fit_ana_90) 

#Trying to plot the anagenetic coefficients table
ggplot(results_ana, 
       aes(x = term, y = estimate)) +
        geom_hline(yintercept = 0, 
                   colour = gray(1/2), lty = 2) +
        geom_point(aes(x = term, 
                    y = estimate)) + 
        geom_linerange(aes(x = term, 
                     ymin = conf.low_95,
                     ymax = conf.high_95),
                   lwd = 1/2) + 
        ggtitle("Anagenesis coefficient table") +
        coord_flip()




#Getting the coefficient estimates for cladogenesis
results_ana_log <- broom::tidy(glm_ana_pois_all_log10)

#Now lets find the confidence intervals for these models.
fit_ana_90_log <- confint(glm_ana_pois_all_log10, level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")

results_ana_log <- dplyr::bind_cols(results_ana_log,fit_ana_90_log) 

#Trying to plot the anagenetic coefficients table
ggplot(results_ana_log, 
       aes(x = term, y = estimate)) +
        geom_hline(yintercept = 0, 
                   colour = gray(1/2), lty = 2) +
        geom_point(aes(x = term, 
                    y = estimate)) + 
        geom_linerange(aes(x = term, 
                     ymin = conf.low_95,
                     ymax = conf.high_95),
                   lwd = 1/2) + 
        ggtitle("Anagenesis coefficient table Log10 and \n Excluded Jamaica and Christmas Island") +
        coord_flip()







# a bit different model
glm_clad_pois_no_geo <- glm(`No. cladogenesis`-`No. edges`~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled),
                    data=output_with_island_cory_no_christmas_or_jamaica, na.action = na.fail,family = quasipoisson(link='log'))

summary(glm_clad_pois_no_geo)


Let us also try to fit a negative binomial model for the cladogenesis

glm_clad_nb_no_geo <- glm.nb(`No. cladogenesis`-`No. edges`~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled)+GeologicalOrigin,
                    data=output_with_island_cory_no_christmas_or_jamaica)

summary(glm_clad_nb_no_geo)







#Now lets find the confidence intervals for these models.
fit_clad_95 <- confint(glm_clad_pois_no_geo, level = 0.95) %>% 
  data.frame() %>% 
  dplyr::rename("conf.low_95" = "X2.5..",
         "conf.high_95" = "X97.5..")


results_clad <- broom::tidy(glm_clad_pois_no_geo) # quasi model
results_clad <- dplyr::bind_cols(results_clad,fit_clad_95) 

#Trying to plot the cladogenetic coefficients table
# Just not possible right now because the standard errors for the geological origin of the cladogenesis dataset is so huge!
ggplot(results_clad, 
       aes(x = term, y = estimate)) +
        geom_hline(yintercept = 0, 
                   colour = gray(1/2), lty = 2) +
        geom_point(aes(x = term, 
                    y = estimate)) + 
        geom_linerange(aes(x = term, 
                     ymin = conf.low_95,
                     ymax = conf.high_95),
                   lwd = 1/2) + 
        ggtitle("Cladogenetic coefficient table") +
        coord_flip()

I have an idea that this huge standard error is caused by the amount of islands having no cladogenesis at all, and that I therefore have to incorporate some sort of zeroinflation model into my dataset.




model1 <- zeroinfl(
  `No. cladogenesis`-`No. edges`~log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled)+GeologicalOrigin |
    log10(mean_mx30_grd)+log10(area)+log10(dist)+log10(nearest_neighbour_distance_border_scaled)+GeologicalOrigin,
                   data=output_with_island_cory_no_christmas_or_jamaica,
                   dist = "negbin",)

summary(model1)





ggplot(output_with_island_cory_no_christmas_or_jamaica, aes(x=`No. cladogenesis`-`No. edges`)) +
  geom_bar()


ggplot(output_with_island_cory_no_christmas_or_jamaica, aes(x=`No. anagenesis`+`No. edges`)) +
  geom_bar()


Creating some boxplots for the predictor variables  for the islands which have cladogenesis and the islands which does not have cladogenesis. 

First I need to find the islands where there is Cladogenesis


clad_islands <- output_with_island_cory$LEVEL3_NAM[which(output_with_island_cory$`No. cladogenesis` >= 1)]
no_clad <- output_with_island_cory$LEVEL3_NAM[-(which(output_with_island_cory$`No. cladogenesis` >= 1))]

boxplot_data <- output_with_island_cory[c("LEVEL3_NAM","mean_mx30_grd","area","dist","nearest_neighbour_distance_border_scaled")]
boxplot_data$clad <- ifelse(boxplot_data$LEVEL3_NAM %in% clad_islands, "Yes", "No")

variabl <- c("Max Height", "Area", "Isolation", "Fragmentation")
for (l in 2:5) {
  print(l)
  xyz <- ggplot(boxplot_data, aes(x=clad, y=log10(boxplot_data[,l]), fill=clad)) + 
    geom_boxplot() +
    geom_jitter(color="black", size=0.4, alpha=0.8) +
    ggtitle(variabl[l-1])
  print(xyz)
}

variabl <- c("Max Height", "Area", "Isolation", "Fragmentation")



anova(glm_ana_pois_all)
glm_ana_pois_all



# Old code which I am a bit afraid to just delete.

Old_tips_to_drop <- c("Brahea edulis - Wolf Eiserhardt", "Copernicia prunifera angela", "Coccothrinax argentata angela 3", "Coccothrinax argentata angela 2","Coccothrinax argentata MSL26 S8", "Coccothrinax barbadensis angela 2", "Zombia antillarum angela", "Hemithrinax ekmaniana SBL209","Thrinax morissii", "Schippia concolor MSL47", "Chelyocarpus chuco SBL223", "Trithrinax brasiliensis angela 2","Trithrinax brasiliensis angela 3", "Trithrinax brasiliensis angela 4", "Trithrinax brasiliensis angela 5","Trithrinax brasiliensis MSL21", "Trithrinax schizophylla angela 2", "Trithrinax schizophylla angela 3", "Serenoa repens angela","Sabal etonia miamiensis according to Patrick", "Undetermined S0 L001 R1 001", "Arenga undulatifolia SBL593","Arenga undulatifolia MSL72 S42", "Hyphaene coriacea - Wolf Eiserhardt", "Johannesteijsmannia altifrons SBL576","Undetermined_S0_L001_R1_001", "Brahea bella","Licuala telifera - Wolf Eiserhardt","Pritchardia flynnii_2", "Pritchardia flynnii_3", "Pritchardia viscosa_2", "Pritchardia viscosa_3", "Pritchardia hardyi_2", "Pritchardia hardyi_3","Pritchardia perlmanii_2", "Pritchardia perlmanii_3", "Pritchardia glabrata_2","pritchardia munroii_2", "Pritchardia kahukuensis_2","Arenga distincta 2","Rhapis puhuongensis 2","Rhapis robusta_2", "Brahea dulcis 2", "Phoenix paludosa 2", "Corypha lecomtei 2", "coccothrinax acuminata 2","Licuala lauterbachii 2", "Rhapis robusta 2","Rhapis subtilis subsp subtilis")

Old_unpublished_and_root <- c("Licuala bakeri", "Licuala bankae","Chuniophoenix suoitiensis", "Licuala heatubunii", "Licuala polybracteae", "Licuala suprafolia","Licuala vanroyensis", "Aphandra natalia", "Dypsis ambositrae", "Eugeissona tristis", "Nypa fructicans","Livistona sp nov hybrid or chinensis","Lanonia sp. Nov “ngoclinhensis”","Lanonia sp. Nov “kontumensis”","Licuala lauterbachii subsp. bismarck", "Coccothrinax sp. nov_T._A._Zanoni_15741","Coccothrinax sp. nov_T._A._Zanoni_33375", "Coccothrinax sp. nov_T._A._Zanoni_34387","Coccothrinax sp.nov _ M._M._Mejía_Pimentel_23684","Licuala sp nov hybrid or chinensis", "Arecaceae chylyocarpus 1", "Arecaceae chylyocarpus 2", "Arecaceae cryophilla")

old_subspANDvars <- c("Phoenix loureiroi var pedunculata","Cryosophila kalbreyeri subsp kalbreyeri","Rhapis laoensis subsp. macrantha","Rhapis subtilis subsp. siamensis","Rhapis subtilis subsp. subtilis","Licuala mattanensis var. paucisecta","Licuala malajana var. humilis","Licuala peltata var. sumawongii","Licuala glabra var. selangorensis","Licuala ramsayi var. tuckeri","Licuala lauterbachii var. lauterbachii","Licuala lauterbachii subsp. bougainvillensis","1023","3158")

# # Renaming tips with problematic names
# #cory_phylo$tip.label["1042"][1] <- "Brahea edulis"
# cory_phylo$tip.label["1193"][1] <- "Coccothrinax barbadensis"
# cory_phylo$tip.label["2004"][1] <- "Arenga undulatifolia"
# cory_phylo$tip.label["1044"][1] <- "Chamaerops humilis"
# cory_phylo$tip.label["1194"][1] <- "Chelyocarpus chuco"
# cory_phylo$tip.label["1249"][1] <- "Coccothrinax argentata"
# cory_phylo$tip.label["1193"][1] <- "Coccothrinax barbadensis"
# cory_phylo$tip.label["1049"][1] <- "Copernicia prunifera"
# cory_phylo$tip.label["1251"][1] <- "Hemithrinax ekmaniana"
# cory_phylo$tip.label["2007"][1] <- "Johannesteijsmannia altifrons"
# cory_phylo$tip.label["1224"][1] <- "Trithrinax brasiliensis"
# cory_phylo$tip.label["1437"][1] <- "Trithrinax schizophylla"
# cory_phylo$tip.label["1078"][1] <- "Zombia antillarum"
# #cory_phylo$tip.label["1061"][1] <- "Licuala telifera"
# cory_phylo$tip.label["1074"][1] <- "Serenoa repens"
# cory_phylo$tip.label["1434"][1] <- "Schippia concolor"
# cory_phylo$tip.label["1022"][1] <- "Phoenix loureiroi var. loureiroi" # Just make this the species for this study
# #cory_phylo$tip.label["1023"][1] <- "Phoenix loureiroi var. pedunculata"
# #cory_phylo$tip.label["3158"][1] <- "Licuala lauterbachii var. bougainvillensis"
# cory_phylo$tip.label["4012"][1] <- "Coccothrinax acuminata"
# cory_phylo$tip.label["4015"][1] <- "Coccothrinax garciana"
# cory_phylo$tip.label["4020"][1] <- "Coccothrinax yuraguana"
# cory_phylo$tip.label["4022"][1] <- "Copernicia brittonorum"
# cory_phylo$tip.label["4023"][1] <- "Copernicia cowellii"
# #cory_phylo$tip.label["4011"][1] <- "Cryosophila kalbreyeri subsp. kalbreyeri"
# cory_phylo$tip.label["4061"][1] <- "Pritchardia glabrata"
# cory_phylo$tip.label["4057"][1] <- "Pritchardia flynnii"
# cory_phylo$tip.label["4063"][1] <- "Pritchardia hardyi"
# cory_phylo$tip.label["4068"][1] <- "Pritchardia kahukuensis"
# cory_phylo$tip.label["4074"][1] <- "Pritchardia munroi"
# cory_phylo$tip.label["4077"][1] <- "Pritchardia pacifica"
# cory_phylo$tip.label["4078"][1] <- "Pritchardia perlmanii"
# cory_phylo$tip.label["4081"][1] <- "Pritchardia viscosa"
# #cory_phylo$tip.label["2058"][1] <- "Rhapis laosensis subsp. macrantha"



df <- data.frame(pet = c("cat","dog","snake"),
  count = c(20,5,94))

plot <- ggplot2::ggplot(df, ggplot2::aes(pet, count, fill = pet)) +
  ggplot2::geom_col() +
  ggplot2::theme(legend.position = "top") # <- specify position of legend

legend = cowplot::get_plot_component(plot, 'guide-box-top', return_all = TRUE)
cowplot::ggdraw(legend)


saveRDS(output_with_island_cory, file.path(datadir,"output_with_island_data.rds"))

