# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## create the following files:
### Traits_compilation_zooxanthellate_cut.csv: 798 species, 22 traits, same as Traits_compilation_zooxanthellate.csv but only with species having available hylogenetic information
### supertree_names_validated_100_1510sp.tre : 100 trees, 1500 + 10 species not present in the original tree but present in Traits_compilation_zooxanthellate_cut.csv
### supertree_names_validated_100_798sp.tre  : 100 trees, 798 species present in Traits_compilation_zooxanthellate_cut.csv
### supertree_names_validated_1000_798sp.tre : same as above but for the 1000 trees

## From:
### supertree_names_validated_100.tre    : updated (names) 100 supertrees from (Huang and Roy 2015), created in phylogeny_coral.R
### Traits_compilation_zooxanthellate.csv: 828 species, 22 traits, created in trait_data_compilation.R

require(here)
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")

require(phytools)

source(paste(wd_Rscripts,"functions.R",sep="/"))

### Import updated and cut supertree for coral from (Huang and Roy 2015) --> names have been updated
supertree <- read.tree(paste(wd_Datasets,"supertree_names_validated_100.tre",sep="/"))
numberTrees <- length(supertree) # 100
numberSpeciesTree <- length(supertree[[1]]$tip.label) # 1500
head(supertree[[1]]$tip.label)

### change tip.labels: "Leptopenus_solidus" --> "Leptopenus solidus"
for(i in 1:numberTrees){
  supertree[[i]]$tip.label <- gsub("_"," ", supertree[[i]]$tip.label)  
}

### Import the coral trait dataset (only zooxanthellae)
Functional_traits_zooxan <- read.csv(paste(wd_Datasets,"Traits_compilation_zooxanthellate.csv",sep="/"),header=T,stringsAsFactors = F)
length(Functional_traits_zooxan$species) # 828
head(Functional_traits_zooxan)

### check if species in Functional_traits_zooxan not in the tree: --> remove the ones with no numerical trait information
speciesNotInTree <- speciesNotPresentFun(Functional_traits_zooxan$species,supertree[[1]]$tip.label)  # species present in the 1st list but not the second:
length(speciesNotInTree) # 40
Functional_traits_zooxan[Functional_traits_zooxan$species %in% speciesNotInTree,]

### remove the species that have no numerical trait values from Functional_traits_zooxan:
speciesToRemove <- character()
numberSpeciesNotInTree <- length(speciesNotInTree) # 40
numericalTraits <- c("age_at_maturity","aggressiveness","bleaching_response_index","chlorophyll_a_concentration","colony_max_diameter","corallite_area","dark_respiration","egg_diameter","fecundity_polyp","growth_rate",
                     "lipid_content","reduced_scattering_coefficient","size_at_maturity","skeletal_density","symbiodinium_density","tissue_thickness")
x <- 1
for(i in 1:numberSpeciesNotInTree){
  FTDF_cut <- Functional_traits_zooxan[Functional_traits_zooxan$species == speciesNotInTree[i],]
  FTDF_cut <- FTDF_cut[,numericalTraits]
  if(sum(is.na(FTDF_cut[1,])) == 16){
    speciesToRemove[x] <- speciesNotInTree[i]
    x <- x + 1
  }
}
length(speciesToRemove) # 29
Functional_traits_zooxan[Functional_traits_zooxan$species %in% speciesToRemove,]

# remove from Functional_traits_zooxan and speciesNotInTree 
Functional_traits_zooxan_cut <- Functional_traits_zooxan
speciesNotInTree_cut <- speciesNotInTree
numberSpeciesToRemove <- length(speciesToRemove) # 29
for(i in 1:numberSpeciesToRemove){
  Functional_traits_zooxan_cut <- Functional_traits_zooxan_cut[Functional_traits_zooxan_cut$species != speciesToRemove[i],]
  speciesNotInTree_cut <- speciesNotInTree_cut[speciesNotInTree_cut != speciesToRemove[i]]
}
length(Functional_traits_zooxan_cut$species) # 799 = 828 - 29
length(speciesNotInTree_cut) # 11 = 40 - 29

### check if the genus is already present in the tree:
genus_species_tree <- genus.species.DF.fun(supertree[[1]]$tip.label)
head(genus_species_tree)
genusTree <- levels(as.factor(genus_species_tree$Genus))
length(genusTree) # 240

genus_species_SpToAdd <- genus.species.DF.fun(speciesNotInTree_cut)
genusSpeciesToAdd <- levels(as.factor(genus_species_SpToAdd$Genus))
length(genusSpeciesToAdd) # 10

# genus present in speciesNotInTree_cut but not in supertree
speciesNotPresentFun(genus_species_SpToAdd$Genus,genus_species_tree$Genus)
# "Sclerophyllia" --> "Sclerophyllia maxima" --> the only available numerical value is for corallite area --> to remove from functional traits DF and speciesNotInTree_cut
Functional_traits_zooxan_cut <- Functional_traits_zooxan_cut[Functional_traits_zooxan_cut$species != "Sclerophyllia maxima",]
length(Functional_traits_zooxan_cut$species) # 798
speciesNotInTree_cut <- speciesNotInTree_cut[speciesNotInTree_cut != "Sclerophyllia maxima"]
length(speciesNotInTree_cut) # 10

# write.csv(Functional_traits_zooxan_cut,paste(wd_Datasets,"Traits_compilation_zooxanthellate_cut.csv",sep="/"),row.names = F)

### check number of species with available data for each traits 
### Number of species available for each trait in traits_compilation_zooxanthellate.csv:
length(Functional_traits_zooxan_cut$species) # 798
numberTraits <- length(colnames(Functional_traits_zooxan_cut))
for(i in 1:numberTraits){
  print(paste(colnames(Functional_traits_zooxan_cut)[i],"",length(na.omit(Functional_traits_zooxan_cut[,i]))))
}

### add the 10 speciesNotInTree_cut species in supertree at the genus level:
# https://rdrr.io/cran/phytools/man/add.species.to.genus.html
# http://blog.phytools.org/2013/11/new-function-to-add-species-to-genus-in.html
# The species are added randomly along the edges to maintain binary nodes in the trees (as opposed to adding them at the root of the edges, which would lead to multifurcating nodes, also called “ploytomies”) 
# because several measures only handle binary topologies

Functional_traits_zooxan_cut <- read.csv(paste(wd_Datasets,"Traits_compilation_zooxanthellate_cut.csv",sep="/"),stringsAsFactors = F,header = T)

supertree_speciesAdded <- supertree
for(i in 1:numberTrees){
  supertree_speciesAdded[[i]]$tip.label <- gsub(" ","_", supertree[[i]]$tip.label)  
}
head(supertree_speciesAdded[[20]]$tip.label)
# speciesNotInTree_underscore <- gsub(" ","_",speciesNotInTree_cut)

numberTrees <- length(supertree_speciesAdded) # 100
numberSpeciesToAdd <- length(speciesNotInTree_cut)  # 10
for(i in 1: numberTrees){
  for(j in 1:numberSpeciesToAdd){
    treeTempo <- add.species.to.genus(supertree_speciesAdded[[i]],speciesNotInTree_cut[j],genus = NULL,where="random") 
    supertree_speciesAdded[[i]] <- treeTempo
  }
  print(i)
}
length(supertree_speciesAdded[[1]]$tip.label) # 1510
head(supertree_speciesAdded[[1]]$tip.labe)
sum(duplicated(supertree_speciesAdded[[1]]$tip.labe)) # 0

speciesPresentFun(speciesNotInTree[1],supertree_speciesAdded[[1]]$tip.labe)

# write tree:
# write.tree(supertree_speciesAdded,paste(wd_Datasets,"supertree_names_validated_100_1510sp.tre"),sep="/")

### Remove in supertree_speciesAdded the species that are not in Functional_traits_zooxan_cut
# Import Traits_compilation_zooxanthellate_cut.csv
Functional_traits_zooxan_cut <- read.csv(paste(wd_Datasets,"Traits_compilation_zooxanthellate_cut.csv",sep="/"),header=T,stringsAsFactors = F)
length(Functional_traits_zooxan_cut$species) # 798
sum(duplicated(Functional_traits_zooxan_cut$species))

# Import supertree_names_validated_100_1510sp.tre
supertree_names_validated_100_1510sp <- read.tree(paste(wd_Datasets,"supertree_names_validated_100_1510sp.tre",sep="/"))

length(supertree_names_validated_100_1510sp[[1]]$tip.label) # 1510
numberTrees <- length(supertree_names_validated_100_1510sp)
for(i in 1:numberTrees){
  supertree_names_validated_100_1510sp[[i]]$tip.label <- gsub("_"," ", supertree_names_validated_100_1510sp[[i]]$tip.label)  
}
head(supertree_names_validated_100_1510sp[[20]]$tip.label)
sum(duplicated(supertree_names_validated_100_1510sp[[20]]$tip.label))

speciesToRemove <- speciesNotPresentFun(supertree_names_validated_100_1510sp[[1]]$tip.label,Functional_traits_zooxan_cut$species) # species present in tree but not in FT dataset
length(speciesToRemove) # 712 = 1510 - 798

speciesPresentinBoth <- speciesPresentBothList(supertree_names_validated_100_1510sp[[1]]$tip.label,Functional_traits_zooxan_cut$species)
length(speciesPresentinBoth) # 798  good

# Remove the species from the tree:
supertree_names_validated_100_798sp <- supertree_names_validated_100_1510sp
numberTrees <- length(supertree_names_validated_100_798sp) # 100
numberSpToRemove <- length(speciesToRemove) # 712
for(i in 1: numberTrees){
  supertree_names_validated_100_798sp[[i]] <- drop.tip(supertree_names_validated_100_1510sp[[i]],speciesToRemove)
  print(i)
}
length(supertree_names_validated_100_798sp[[1]]$tip.label) # 798
head(supertree_names_validated_100_798sp[[1]]$tip.labe)
sum(duplicated(supertree_names_validated_100_798sp[[1]]$tip.labe)) # 0

# write tree:
# write.tree(supertree_names_validated_100_798sp,paste(wd_Datasets,"supertree_names_validated_100_798sp.tre",sep="/"))

#***************************************************************************************
# ******************    DO THE SAME FOR THE 1000 trees **************************
#***************************************************************************************

### Import supertree_names_validated.tre
supertree_1000 <- read.tree(paste(wd_Datasets,"supertree_names_validated.tre",sep="/"))
length(supertree_1000[[1]]$tip.label) # 1500
head(supertree_1000[[1]]$tip.label)

numberTrees <- length(supertree_1000) # 1000
for(i in 1:numberTrees){
  supertree_1000[[i]]$tip.label <- gsub("_"," ", supertree_1000[[i]]$tip.label) 
}
head(supertree_1000[[1]]$tip.label)

### Import Traits_compilation_zooxanthellate_cut.csv
Functional_traits_zooxan_cut <- read.csv(paste(wd_Datasets,"Traits_compilation_zooxanthellate_cut.csv",sep="/"),header=T,stringsAsFactors = F)
length(Functional_traits_zooxan_cut$species) # 798

speciesToAddTrees <- speciesNotPresentFun(Functional_traits_zooxan_cut$species,supertree_1000[[1]]$tip.label) # species present in trait df table but not in tree
length(speciesToAddTrees) # 10

speciesToRemoveTrees <- speciesNotPresentFun(supertree_1000[[1]]$tip.label,Functional_traits_zooxan_cut$species) # species present in tree but not in trait df
length(speciesToRemoveTrees) # 712 = 1500 - 798 + 10

# Add the 10 species:
supertree_1000_add <- supertree_1000
for(i in 1:numberTrees){
  supertree_1000_add[[i]]$tip.label <- gsub(" ","_", supertree_1000[[i]]$tip.label)  
}
head(supertree_1000_add[[20]]$tip.label)
library(phytools)
numberSpeciesToAdd <- length(speciesToAddTrees)  # 10
for(i in 1: numberTrees){
  for(j in 1:numberSpeciesToAdd){
    treeTempo <- add.species.to.genus(supertree_1000_add[[i]],speciesToAddTrees[j],genus = NULL,where="random") 
    supertree_1000_add[[i]] <- treeTempo
  }
  print(i)
}
length(supertree_1000_add[[1]]$tip.label) # 1510
head(supertree_1000_add[[1]]$tip.labe)
sum(duplicated(supertree_1000_add[[1]]$tip.labe)) # 0
for(i in 1:numberTrees){
  supertree_1000_add[[i]]$tip.label <- gsub("_"," ", supertree_1000_add[[i]]$tip.label)  
}
head(supertree_1000_add_cut[[20]]$tip.label)

supertree_1000_add_cut <- supertree_1000_add
for(i in 1: numberTrees){
  supertree_1000_add_cut[[i]] <- drop.tip(supertree_1000_add[[i]],speciesToRemoveTrees)
  print(i)
}
length(supertree_1000_add_cut[[1]]$tip.label) # 798
head(supertree_1000_add_cut[[1]]$tip.labe)
sum(duplicated(supertree_1000_add_cut[[1]]$tip.labe)) # 0

# write tree:
# write.tree(supertree_1000_add_cut,paste(wd_Datasets,"supertree_names_validated_1000_798sp.tre",sep="/"))
