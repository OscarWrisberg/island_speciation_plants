
library(ggtree)
library(tidytree)


cory_phylo <- read.tree(file.path(datadir,"astral_tree_orthologs.tre"))

figurename <- read.table(file.path(datadir,"names_for_tips.csv"), sep=",", colClasses = "character")
#figurename <- figurename[]
figurename_idx <- figurename$V2
names(figurename_idx) <- figurename$V1

#Renaming tips based on the figurename_idx
cory_phylo$tip.label = figurename_idx[cory_phylo$tip.label]

# Tips which I need to drop because I have duplicates
tips_to_drop <- c("Brahea dulcis 1", "Brahea dulcis 3", "Brahea edulis 2", "Coccothrinax acuminata 2", "Coccothrinax argentata 1", "Coccothrinax argentata 2", "Coccothrinax argentata 4", "Coccothrinax argentea 2", "Coccothrinax barbadensis 1", "Coccothrinax garciana 2", "Copernicia cowellii 2", "Corypha lecomtei 2", "Licuala lauterbachii 1", "Phoenix paludosa 1", "Pritchardia flynnii 2", "Pritchardia flynnii 3", "Pritchardia glabrata 2", "Pritchardia hardyi 2","Pritchardia hardyi 3", "Pritchardia kahukuensis 2", "Pritchardia munroi 1", "Pritchardia perlmanii 2","Pritchardia perlmanii 3", "Pritchardia viscosa 2", "Rhapis robusta 2", "Rhapis puhuongensis 2", "Trithrinax brasiliensis 2", "Trithrinax brasiliensis 3", "Trithrinax schizophylla 2", "Trithrinax schizophylla 3", "Zombia antillarum 2")

#Dropping outgroup, novo species and a missing subsp.
root <- c("Aphandra natalia","Eugeissona tristis","Nypa fruticans","Chrysalidocarpus ambositrae")

# Species which are not published or identified
unpublished <- c("Chelyocarpus sp. José De Gracia 826","Chelyocarpus sp. José De Gracia 892","Coccothrinax aff. Jimenezii T. A. Zanoni 33375","Coccothrinax bermudezi","Coccothrinax bonettiana","Coccothrinax montgomeryana","Cryosophila sp.  José De Gracia 824","Lanonia sp. Henderson A. 3590","Lanonia sp. Henderson A. 3727","Licuala heatubunii","Licuala vanroyensis","Livistona sp. Henderson A. 3730","Sabal miamiensis","Licuala polybracteae")

#Sub sp and Variants
subspANDvars <- c("Cryosophila kalbreyeri subsp. cogolloi","Cryosophila kalbreyeri subsp. kalbreyeri","Licuala glabra var. selangorensis","Licuala lauterbachii subsp. Peekelii","Licuala lauterbachii var. bougainvillensis","Licuala lauterbachii var. lauterbachii","Licuala malajana var. humilis","Licuala mattanensis var. paucisecta","Licuala peltata var. sumawongii","Phoenix loureiroi var. pedunculata","Rhapis subtilis subsp. Siamensis","Rhapis subtilis subsp. Subtilis","Washingtonia filifera var. robusta","Rhapis laosensis subsp. macrantha","Licuala ramsayi var. tuckeri")

# Dropping tips
cory_phylo <- drop.tip(cory_phylo,tips_to_drop) # tips to drop
cory_phylo <- drop.tip(cory_phylo,root) # root species
cory_phylo <- drop.tip(cory_phylo,unpublished) # unpublished or unidentified species
cory_phylo <- drop.tip(cory_phylo, subspANDvars) # sub species and vars

# Remove the trailing numbers
cory_phylo$tip.label <- gsub(" \\d+$", "", cory_phylo$tip.label)

# Check for duplicates
duplicates <- duplicated(cory_phylo$tip.label)


# Tips where I need to convert a sub sp to the sp
cory_phylo$tip.label["3130"][1] <- "Licuala glabra"
cory_phylo$tip.label["3176"][1] <- "Licuala malajana"
cory_phylo$tip.label["3412"][1] <- "Licuala mattanensis"
cory_phylo$tip.label["3334"][1] <- "Licuala peltata"
cory_phylo$tip.label["3390"][1] <- "Licuala ramsayi"
cory_phylo$tip.label["1022"][1] <- "Phoenix loureiroi"


#Making tree ultrametric
cory_phylo$edge.length[is.na(cory_phylo$edge.length)] <- 0.001


# Lets try to plot this phylogeny and see how it looks so far.

ggtree(cory_phylo, layout = "dendrogram", ladderize = FALSE, size = 0.6, )

ggtree(cory_phylo, ladderize = TRUE, size = 0.4, branch.length = "none" ) +
  coord_flip() +
  coord_flip()


#Creating a list of all accepted names within Coryphoideae
data <- read.csv(file.path(datadir,"checklist_names.txt"), sep = "|")

species <- data[which(data$taxon_rank == "Species" | data$taxon_rank == "Variety" | data$taxon_rank == "Subspecies"),] # Selecting only species and varieties/subspecies
palms <- species[which(species$family=="Arecaceae"),] # Selecting only Palms
apalms <- palms[which(palms$taxon_status=="Accepted"),] # Selecting only Accepted species

# List of Genera in Coryphoideae
genera <- c("Sabal","Schippia","Trithrinax","Zombia","Coccothrinax","Hemithrinax","Leucothrinax","Thrinax","Chelyocarpus",
            "Cryosophila","Itaya","Phoenix","Chamaerops","Guihaia","Trachycarpus","Rhapidophyllum","Maxburretia","Rhapis",
            "Livistona","Licuala","Johannesteijsmannia","Pholidocarpus","Pritchardiopsis","Acoelorraphe","Serenoa","Brahea",
            "Colpothrinax","Copernicia","Pritchardia","Washingtonia","Chuniophoenix","Kerriodoxa","Nannorrhops","Tahina",
            "Caryota","Arenga","Wallichia","Corypha","Bismarckia","Satranala","Hyphaene","Medemia","Latania","Lodoicea",
            "Borassodendron","Borassus","Lanonia","Saribus","Sabinaria")


subfamilies <- c("Phoeniceae","Trachycarpeae","Sabaleae","Cryosophileae","Chuniophoeniceae","Caryoteae","Corypheae","Borasseae")


# Selecting only palms in Coryphoideae
cory_data <- apalms[which(apalms$genus %in% genera),]

#Dropping synonym species names rows in data.frame
cory_data_accepted <- cory_data[which(cory_data$taxon_status =="Accepted"),] # | (cory_data$taxon_rank=="Variety"))


# Trying to load WCVP distribution data to see if it is better
dist_data <- read.csv(file.path(datadir,"checklist_distribution.txt"), sep = "|") 
cory_dist_data <- dist_data[which(dist_data$plant_name_id %in% cory_data_accepted$plant_name_id),]

#Renaming WCS ID's to species names
uniqe_names_cory <- unique(cory_dist_data$plant_name_id)
for (i in 1:length(unique(cory_dist_data$plant_name_id))){
  id <- uniqe_names_cory[i]
  
  if (id %in% cory_dist_data$plant_name_id){
    cory_dist_data[which(cory_dist_data$plant_name_id == id),1] <- cory_data_accepted[which(cory_data_accepted$plant_name_id ==id),22]
  } else {
    print("not in")
  }
}

#We probably need to remove the introduced locations.
cory_dist_data_no_introduced <- cory_dist_data[which(cory_dist_data$introduced != 1 & cory_dist_data$extinct != 1 & cory_dist_data$location_doubtful != 1),]


#Soo how many different locations are there in the new database??
areas <- unique(cory_dist_data_no_introduced$area) # 183. Seems doable.
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
# There are currently 4 species in my tree which arent in the World checklist downloaded on the 9th of May
# These are "Chuniophoenix suoitienensis","Licuala bakeri","Licuala bankae","Licuala suprafolia" 
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


islands <- c("Borneo", "Sumatera", "Hainan", "Taiwan", "Christmas I.", "Maluku", "New Guinea", "Philippines", "Jawa", "Sulawesi", "Nansei-shoto", "Madagascar", "Gulf of Guinea Is.", "Lesser Sunda Is.", "Sri Lanka", "Mexican Pacific Is.", "Andaman Is.", "Nicobar Is.", "Vanuatu", "Bismarck Archipelago", "Solomon Is.",  "Sicilia", "Cuba", "Leeward Is.", "Bahamas", "Turks-Caicos Is.",  "Dominican Republic", "Trinidad-Tobago", "Venezuelan Antilles", "Windward Is.", "Jamaica", "Cayman Is.","Mozambique Channel I", "Hawaii" )

# We decide to join Haiti and the Dominican republic because it is mainly a single island hispaniola
community_matrixes_orthologs[which(community_matrixes_orthologs[,"Haiti"] == 1),"Dominican Republic"] <- "1" # Making species found in Haiti also be found on Dominical Republic.
community_matrixes_orthologs <- community_matrixes_orthologs[,-which(colnames(community_matrixes_orthologs)=="Haiti")] # Removing Haiti from the community matrix



FindEndemicLineages<- function(tree,matrix,areas){

  edges_per_island <- list()
  species_per_island <- list()
  endemics_from_radiation <- list()
  comb <- list()
  area_with_sp <- c()
  
  
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
        #Widespread_chk <- all(only_on_islands == TRUE) # Checking if all descendants are endemic to the island
        
        if (island_chk == TRUE && endemic_chk ==TRUE) {
          if (!(is.null(getParent(cory_phylo,node))) && !(node %in% decen_edges_island)) { #This line checks if the parent node is the root, then it checks if the node is a descendant of a node already added to the list.
          edges_with_only_island_decen <- append(edges_with_only_island_decen,node) #Adding edge to a list of edges endemic to the island 
          decen_edges_island <- append(decen_edges_island,decen_edges) # Adding all edges descending from this node to a list
          decendants_of_edges <- append(decendants_of_edges,decen_tips) # Adding species descending from these edges to a list
          if( !(areas[i] %in% area_with_sp)) {
          area_with_sp <- append(area_with_sp,areas[i])
            }
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
          if( !(areas[i] %in% area_with_sp)) {
          area_with_sp <- append(area_with_sp,areas[i])
            }
        }
      }
    }
      if (length(edges_with_only_island_decen) > 0) {
        edges_per_island[[i]] <- edges_with_only_island_decen
      } else {
        edges_per_island[[i]] <- list(NULL)
      }
      
      if (length(sp_on_island) > 0) {
        species_per_island[[i]] <- sp_on_island
      } else {
        species_per_island[[i]] <- list(NULL)
      }
      
      if (length(decendants_of_edges) > 0) {
        endemics_from_radiation[[i]] <- decendants_of_edges
      } else {
        endemics_from_radiation[[i]] <- list(NULL)
      }

    #edges_per_island[[i]] <- edges_with_only_island_decen # Saving all edges where all descendants are endemic to an island
    #species_per_island[[i]] <- sp_on_island # Saving all colonizing Endemic species found on island
    #endemics_from_radiation[[i]] <- decendants_of_edges # Saving all species from edges that radiated
    
    #print(i)
    #print()
    #print(areas[i])
    #print(edges_per_island[[i]])
    #print(species_per_island[[i]])
    #print(endemics_from_radiation[[i]])
    
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






FindEndemicLineagesMultipleAreas <- function(tree, presence_matrix, areas) {
  endemic_species <- c()

  for (ntip in 1:Ntip(tree)) {
    # Counting the number of target areas where the species is present
    num_present_areas <- sum(presence_matrix[tree$tip.label[ntip], areas] == "1")

    # Checking if the species is present in at least 2 target areas
    multiple_endemic_sp <- num_present_areas >= 2
    print(multiple_endemic_sp)

    # If the species is present in at least 2 target areas, add it to the list of endemic species
    if (multiple_endemic_sp) {
      endemic_species <- c(endemic_species, tree$tip.label[ntip])
    }
  }

  return(endemic_species)
}






AddSpeciation <- function(tree, endemic, widespread) {
  
  # Convert the input tree to a tibble (data frame) for easier manipulation
  tbbl_tree <- as_tibble(tree)

  # Initialize a new column 'speciation' with default value
  tbbl_tree$speciation <- "Not endemic on Island"
  
  # Loop through all labels in the tree
  for (i in 1:length(tbbl_tree$label)) {
    # If the label is in the widespread list, mark it as "Widespread on Islands"
    if (tbbl_tree$label[i] %in% widespread) {
      tbbl_tree$speciation[i] <- "Widespread on Islands"
    }
  }

  # Loop through each endemic group
  for (i in 1:length(endemic)) {
    for (j in 1:length(endemic[[i]])) {
      # Check if the element is for Cladogenesis (edges or species)
      if (j == 1 | j == 2) {
        # If it is the first element, mark edges
        if (j == 1) {
          # Check if the current element is null
          if (is.null(endemic[[i]][[j]][[1]][[1]]) == TRUE) {
            next # Skip to the next iteration if it is null
          } else {
            # Mark edges in the tree as cladogenetic edges
            tbbl_tree$speciation[which(tbbl_tree$node %in% endemic[[i]][[j]][[1]])] <- "Within-region speciation"
            
            if (!is.null(endemic[[i]][[j]])) {
              # Find all descendants from each node
              decendants <- Descendants(tree, endemic[[i]][[j]][[1]], type = "all")
              # Find only tip descendants
              decen_tips <- Descendants(tree, endemic[[i]][[j]][[1]], type = "tips")
              # Find all descendant edges that are not tips
              decen_edges <- decendants[which(!(decendants %in% decen_tips))]
            }
            # Mark all descendant edges as cladogenetic
            if (!(purrr::is_empty(decen_edges))) {
              for (k in 1:length(decen_edges)) {
                tbbl_tree$speciation[which(tbbl_tree$node %in% decen_edges[[k]])] <- "Within-region speciation"
              }
            }
            # Marking the basal most edge of each Cladogenetic clade as Anagenesis
            basal_edges <- which(tbbl_tree$node %in% endemic[[i]][[j]][[1]])
            for (edge in basal_edges) {
              parent_node <- tbbl_tree$parent[edge]
              # If the parent node is not already marked as Cladogenesis, mark the edge as Anagenesis
              if (tbbl_tree$speciation[parent_node] != "Within-region speciation") {
                tbbl_tree$speciation[edge] <- "Between-region speciation"
              }
            }
          }
        }
        # If it is the second element, mark species
        if (j == 2) {
          # Check if the current element is null
          if (is.null(endemic[[i]][[j]][[1]][[1]])) {
            next # Skip to the next iteration if it is null
          } else {
            # Mark species in the tree as cladogenetic species
            tbbl_tree$speciation[which(tbbl_tree$label %in% endemic[[i]][[j]][[1]])] <- "Within-region speciation"
          }
        }
      } else {
        # For anagenesis, mark species
        if (is.null(endemic[[i]][[j]][[1]][[1]])) {
          next # Skip to the next iteration if it is null
        } else {
          # Mark species in the tree as anagenetic species
          tbbl_tree$speciation[which(tbbl_tree$label %in% endemic[[i]][[j]][[1]])] <- "Between-region speciation"
        }
      }
    }
  }
  
  # Convert the tibble back to a treedata object and return it
  return(as.treedata(tbbl_tree))
}




#Sabaleae
  #Sabal 
  # Sabal Yapa Sabal Etonio
MRCA_sabal <- MRCA(cory_phylo, c("Sabal yapa", "Sabal etonia")) #502

# Cryosophileae
  # Schippia, Itaya, Cryosophila, Trithrinax, Chelyocarpus, Zombia, Coccothrinax, Thrinax, Hemithrinax, Leucothrinax
  # Cryosophila nana, Coccothrinax macroglossa
MRCA_cryosophileae <- MRCA(cory_phylo, c("Sabinaria magnifica", "Coccothrinax macroglossa")) #495

#Phoeniceae
  # Phoenix
  # Phoenix rupicola, Phoenix dactylifera
MRCA_phoeniceae <- MRCA(cory_phylo, c("Phoenix rupicola", "Phoenix dactylifera")) #447

# Brahea
  # Brahea
  # Brahea aculeata , Brahea calcarea
MRCA_brahea <- MRCA(cory_phylo,c("Brahea aculeata","Brahea calcarea"))

#Trachycarpeae
# Subtribe Rhapidinae
  # Chamaerops, Guihaia, Trachycarpus, Rhapidophyllum, Maxburretia, Rhapis, (Brahea, Colpothrinax)
  # Colpothrinax aphanopetala, Rhapis puhuongensis
MRCA_rhapidinae <- MRCA(cory_phylo, c("Rhapidophyllum hystrix", "Rhapis humilis")) #636

#Subtribe livistoninae 
  #Livistona, Licuala, johannesteijsmannia, Pholidocarpus, Saribus, Acoelorrhaphe, Serenoa
  # Serenoa repens, Licuala simplex
MRCA_livistoninae <-MRCA(cory_phylo, c("Livistona carinensis", "Licuala simplex")) #636

# Chuniophoeniceae
  # Chuniophoenix, Kerriodoxa, Nannorrhops, Tahina
  # Nannorrhops ritchieana, chuniophoenix hainanensis
MRCA_chuniophoeniceae <- MRCA(cory_phylo, c("Nannorrhops ritchieana", "Chuniophoenix hainanensis")) #409

# Caryoteae
  # Caryota, Wallichia, Arenga
  # Caryota obtusa, Wallichia gracilis
MRCA_caryoteae <- MRCA(cory_phylo, c("Caryota obtusa", "Wallichia gracilis")) #436

#Corypheae
  # Corypha
  # Corypha lecomtei, Corypha taliera
MRCA_corypheae <- MRCA(cory_phylo, c("Corypha lecomtei", "Corypha taliera")) #415

# Borasseae
  # Subtribe Hyphaeninae
  # Bismarckia, Satranala, Hyphaene, Medemia
  # Bismarkia nobilis, Hyphaene thebaica
MRCA_borasseae <- MRCA(cory_phylo, c("Bismarckia nobilis", "Hyphaene thebaica")) #419

  # Subtribe Lataniinae
  # Borassodendron, Latania, Borassus, Lodoicea
  # Borassus aethiopum, Latania lontaroides
MRCA_lataniinae <- MRCA(cory_phylo, c("Borassus aethiopum", "Latania lontaroides")) #427

# Serenoa and Acoelorraphe
MRCA_Serenoa_Acoelorraphe <- MRCA(cory_phylo, c("Serenoa repens", "Acoelorraphe wrightii"))

#Colpothrinax
MRCA_colpothrinax <- MRCA(cory_phylo, c("Colpothrinax cookii", "Colpothrinax aphanopetala"))

# Pritchardia
MRCA_pritchardia <- MRCA(cory_phylo, c("Pritchardia mitiaroana", "Pritchardia viscosa"))

# Copernicia
MRCA_copernicia <- MRCA(cory_phylo, c("Copernicia alba", "Copernicia hospita"))

# Washingtonia
MRCA_washingtonia <- MRCA(cory_phylo,"Washingtonia filifera")


dat_subfam <- data.frame(
           node = c(MRCA_sabal,MRCA_cryosophileae,MRCA_phoeniceae,MRCA_rhapidinae,MRCA_livistoninae,MRCA_chuniophoeniceae,MRCA_caryoteae,MRCA_corypheae,MRCA_borasseae,MRCA_lataniinae,MRCA_brahea,MRCA_Serenoa_Acoelorraphe,MRCA_colpothrinax,MRCA_pritchardia,MRCA_copernicia,MRCA_washingtonia),
           name = c("Sabaleae","Cryosophileae","Phoeniceae","Rhapidineae","Livistoninae","Chuniophoeniceae","Caryoteae","Corypheae", "Hyphaeninae","Lataniinae","Brahea", "Serenoa & Acoelorraphe","Colpothrinax","Pritchardia","Copernicia","Washingtonia")
       )




output_community_matrix_islands <- FindEndemicLineages(cory_phylo,community_matrixes_orthologs,islands)

widespread_island_lineages <- FindEndemicLineagesMultipleAreas(cory_phylo, community_matrixes_orthologs, islands)
  
cory_phylo_dat <- AddSpeciation(cory_phylo, output_community_matrix_islands, widespread_island_lineages)


cory_phylo_dat_plot <- ggtree(cory_phylo_dat, aes(color=speciation), layout = "circular", ladderize = FALSE, size = 0.8, branch.length = "none") +
  geom_cladelab(data = dat_subfam,mapping = aes(node = node,label = name),fontsize = 4, align = TRUE, angle = "auto", offset = 2, horizontal = TRUE) +
  scale_colour_manual(values = c("#efc86e","#28292b","#EE4B2B","#6f9969")) +
  theme(legend.position = "bottom") +
  theme(legend.title=element_blank())

print(cory_phylo_dat_plot)


# Testing whether there is a phylogenetic signal in the cory_phylo_dat
# I want to know if there is a signal in being on an island or not and also whether there is a signal between within-region speciation and betweem-region speciation.
# I will use the geiger package with the fitdiscrete function to test this.

# First I create a copy of the data to work with
cory_phylo_dat_tibble <- as.tibble(cory_phylo_dat)
tip_labels <- cory_phylo_dat_tibble$label[which(cory_phylo_dat_tibble$label %in% cory_phylo$tip.label)]

# Now I add the data from the cory_phylo_dat_tibble to the tip labels to create a dataframe for testing phylogenetic signal
cory_phylo_dat_phylosig <- cbind(tip_labels, cory_phylo_dat_tibble[which(cory_phylo_dat_tibble$label %in% cory_phylo$tip.label),])

# can I add the numbers from the tip labels of the phylogeny to the data frame?
tip_nr <- names(cory_phylo$tip.label)
rownames(cory_phylo_dat_phylosig) <- tip_nr
cory_phylo_dat_phylosig[,c(1,6)]

speciation_method <- cory_phylo_dat_phylosig$speciation
names(speciation_method) <- cory_phylo$tip.label

# Then I need to fit the model
phylosig_model_endemic_1 <- fitDiscrete(phy = cory_phylo, dat = speciation_method, model = "ER", nsim = 1000, transform = "lambda")
phylosig_model_endemic_2 <- fitDiscrete(phy = cory_phylo, dat = speciation_method, model = "ER", nsim = 1000, transform = "none")

phylosig_model_endemic_1$opt$lnL
phylosig_model_endemic_2$opt$lnL

logLik_1 <- phylosig_model_endemic_1$opt$lnL
logLik_2 <- phylosig_model_endemic_2$opt$lnL

LR <- 2 * (logLik_1 - logLik_2)
p_value <- pchisq(LR, df = 1, lower.tail = FALSE)

print(paste("Testing wether there is a phylogenetic signal in being Between-region speciation, Within-region speciation, Widespread on Islands and Not endemic on Island"))
print(paste("Likelihood ratio test p-value:", p_value))

# Now I want to do the same but I want to test wether there is a phylogenetic signal in being on an island or not.
# Create the data from the speciation_method.
# if speciation method is "Not endemic on Island" then the value is 0, if it is "Widespread on Islands", "Between-region speciation" or "Within-region speciation" then the value is 1.
speciation_num <- c()
for (i in 1:length(speciation_method)) {
  if (speciation_method[i] == "Not endemic on Island") {
    speciation_num <- append(speciation_num, "Not island")
  } else {
    speciation_num <- append(speciation_num, "Island")
  }
}

names(speciation_num) <- cory_phylo$tip.label

phylosig_model_speciation_1 <- fitDiscrete(phy = cory_phylo, dat = speciation_num, model = "ER", nsim = 1000, transform = "lambda")
phylosig_model_speciation_2 <- fitDiscrete(phy = cory_phylo, dat = speciation_num, model = "ER", nsim = 1000, transform = "none")

logLik_spec_1 <- phylosig_model_speciation_1$opt$lnL
logLik_spec_2 <- phylosig_model_speciation_2$opt$lnL

LR_spec <- 2 * (logLik_spec_1 - logLik_spec_2)
p_value_spec <- pchisq(LR_spec, df = 1, lower.tail = FALSE)

print(paste("Testing wether there is a phylogenetic signal in being found on an island or not"))
print(paste("Likelihood ratio test p-value:", p_value_spec))

# phoenix_subtree <- extract.clade(cory_phylo, MRCA(cory_phylo, c("Phoenix rupicola", "Phoenix dactylifera")))
# 
# ggtree(phoenix_subtree) +
#     geom_tiplab(size = 2.5)
  

# Now I need to create some sub phylogenies of the big island radiations.

# Ill have to create some subplots of
# These should be from the most recent common ancestor from the following species.
#  - Licuala ( Licuala simplex, Licuala longipes)
#  - Pritchardia / Copernicia / Washingtonia ( Washingtonia filifera, Pritchardia pacifica )
#  - CoccoThrinax ( Coccothrinax macroglossa, Thrinax radiata )
#  - Corypha Arenga ( Kerriodixa elegans, Arenga longicarpa )
 
 

# MRCA(cory_phylo_dat, c("Licuala simplex", "Licuala longipes")) # 691
# MRCA(cory_phylo_dat, c("Washingtonia filifera", "Pritchardia pacifica")) #554
# MRCA(cory_phylo_dat, c("Coccothrinax macroglossa", "Thrinax radiata")) #511
# MRCA(cory_phylo_dat, c("Kerriodoxa elegans", "Arenga longicarpa")) # 408

 


licuala_subtree <- extract.clade(cory_phylo, MRCA(cory_phylo_dat, c("Licuala simplex", "Licuala longipes")))

licuala_output_community_matrix_islands <- FindEndemicLineages(licuala_subtree,community_matrixes_orthologs,islands)

licuala_widespread_island_lineages <- FindEndemicLineagesMultipleAreas(licuala_subtree, community_matrixes_orthologs, islands)

licuala_phylo_dat <- AddSpeciation(licuala_subtree, licuala_output_community_matrix_islands, licuala_widespread_island_lineages)

licuala_plot <- ggtree(licuala_phylo_dat, aes(color=speciation), ladderize = TRUE, size = 0.8) +
  scale_colour_manual(values = c("#efc86e","#6f9969","#28292b", "#EE4B2B")) +
  geom_tiplab(size = 3) +
  theme(legend.position = "none") +
  ggplot2::xlim(0, 3.8)



pritchardia_subtree <- extract.clade(cory_phylo, MRCA(cory_phylo_dat, c("Washingtonia filifera", "Pritchardia pacifica")))

pritchardia_output_community_matrix_islands <- FindEndemicLineages(pritchardia_subtree,community_matrixes_orthologs,islands)

pritchardia_widespread_island_lineages <- FindEndemicLineagesMultipleAreas(pritchardia_subtree, community_matrixes_orthologs, islands)

pritchardia_phylo_dat <- AddSpeciation(pritchardia_subtree, pritchardia_output_community_matrix_islands, pritchardia_widespread_island_lineages)

pritchardia_plot <- ggtree(pritchardia_phylo_dat, aes(color=speciation), ladderize = TRUE, size = 0.8) +
  scale_colour_manual(values = c("#efc86e","#6f9969","#28292b", "#EE4B2B")) +
  geom_tiplab(size = 3) +
  theme(legend.position = "none") +
  ggplot2::xlim(0, 7.5)






coccothrinax_subtree <- extract.clade(cory_phylo, MRCA(cory_phylo_dat, c("Coccothrinax macroglossa", "Thrinax radiata")))

coccothrinax_output_community_matrix_islands <- FindEndemicLineages(coccothrinax_subtree,community_matrixes_orthologs,islands)

coccothrinax_widespread_island_lineages <- FindEndemicLineagesMultipleAreas(coccothrinax_subtree, community_matrixes_orthologs, islands)

coccothrinax_phylo_dat <- AddSpeciation(coccothrinax_subtree, coccothrinax_output_community_matrix_islands, coccothrinax_widespread_island_lineages)

coccothrinax_plot <- ggtree(coccothrinax_phylo_dat, aes(color=speciation), ladderize = TRUE, size = 0.8) +
  scale_colour_manual(values = c("#efc86e","#6f9969","#28292b", "#EE4B2B")) +
  geom_tiplab(size = 3) +
  theme(legend.position = "none") +
  ggplot2::xlim(0, 5.5)





arenga_subtree <- extract.clade(cory_phylo, MRCA(cory_phylo_dat, c("Kerriodoxa elegans", "Arenga longicarpa")))

arenga_output_community_matrix_islands <- FindEndemicLineages(arenga_subtree,community_matrixes_orthologs,islands)

arenga_widespread_island_lineages <- FindEndemicLineagesMultipleAreas(arenga_subtree, community_matrixes_orthologs, islands)

arenga_phylo_dat <- AddSpeciation(arenga_subtree, arenga_output_community_matrix_islands, arenga_widespread_island_lineages)

arenga_plot <- ggtree(arenga_phylo_dat, aes(color=speciation), ladderize = TRUE, size = 0.8) +
  scale_colour_manual(values = c("#efc86e","#6f9969","#28292b", "#EE4B2B")) +
  geom_tiplab(size = 3) +
  theme(legend.position = "none") +
  ggplot2::xlim(0, 15)






small_phylos <- cowplot::plot_grid(licuala_plot, pritchardia_plot, coccothrinax_plot, arenga_plot)

cowplot::plot_grid(cory_phylo_dat_plot, small_phylos)



small_phylo_2 <-cowplot::plot_grid(pritchardia_plot, coccothrinax_plot, ncol = 1)

small_phylo_3 <-cowplot::plot_grid(licuala_plot, small_phylo_2, ncol = 2, nrow = 1)

phylo4 <- cowplot::plot_grid(cory_phylo_dat_plot, small_phylo_3, ncol = 2, nrow = 1)

phylo4


