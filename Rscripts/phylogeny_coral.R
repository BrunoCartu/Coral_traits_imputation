# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## check the species names in (Huang and Roy 2015)'s supertree 
## create:
### supertree_names_validated.tre : 1500 species, 1000 trees
### supertree_names_validated_100.tre : 1500 species, 100 1st trees

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_phylogeny_coral <- paste(wd,"/Traits_extra/phylogeny_coral",sep="")

require(phytools)

source(paste(wd_Rscripts,"functions.R",sep="/"))

### import of coral superTree publiched by (Huang and Roy 2015):
# http://dx.doi.org/10.5061/dryad.178n3
# http://rstb.royalsocietypublishing.org/content/370/1662/20140010
supertree_original <- read.nexus(file=paste(wd_Datasets_original,"Huang&Roy_Supertree.tre",sep="/"))
print(supertree_original)
# 1000 phylogenetic trees
length(supertree_original[[1]]$tip.label) # 1547
head(supertree_original[[1]]$tip.label)   # "ACR_Acropora_abrolhosensis" "ACR_Acropora_abrotanoides"  "ACR_Acropora_aculeus"       "ACR_Acropora_acuminata"     "ACR_Acropora_akajimensis"   "ACR_Acropora_anthocercis" 
sum(duplicated(supertree_original[[1]]$tip.label)) # 0
supertree_modif <- supertree_original

### change tip.lables from "ACR_Acropora_abrotanoides" to "Acropora_abrotanoides", and then "Acropora_abrotanoides" to "Acropora abrotanoides"
supertree_modif <- .uncompressTipLabel(supertree_modif) # uncompress tip labels
numberTree <- length(supertree_modif) # 1000
for(i in 1:numberTree){
  supertree_modif[[i]]$tip.label <- gsub("^.*?_","", supertree_modif[[i]]$tip.label) # remove symboles on the left of 1st "_"
  supertree_modif[[i]]$tip.label <- gsub("_"," ", supertree_modif[[i]]$tip.label)    # replace "_" by " "
}
head(supertree_modif[[2]]$tip.label) # "Acropora abrolhosensis" "Acropora abrotanoides"  "Acropora aculeus"       "Acropora acuminata"     "Acropora akajimensis"   "Acropora anthocercis"   

### correction for names:
speciesToCheck <- checkSpeciesNames(supertree_modif[[1]]$tip.label, wd_Datasets = wd_Datasets) #  "All species are in the list"
speciesToCheck

### import list of names to change: speciesNameChecked contains all the coral species names encountered (1st cloumn) and the corresponding correct name (second column)
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp) # 1786
length(levels(as.factor(speciesNameChecked$nameSp_checked))) # 1581

### list of species whose name has to be updated AND whose updated name is not already present in the list:
speciesToUpdate <- character()

### list of species to remove from the tree because they have to be corrected but updated name is already in the tree:
speciesToRemoveTree <- character()

### go through all the species in one tree and fill up speciesToUpdate and speciesToRemoveTree
species <- supertree_modif[[1]]$tip.label
numberspeciesTree <- length(species) # 1547
x <- 1 # position in speciesToRemoveTree
y <- 1 # position in speciesToUpdate
for(i in 1:numberspeciesTree){
  speciesNames <- speciesNameChecked[speciesNameChecked$nameSp == species[i],]
  if(speciesNames$nameSp != speciesNames$nameSp_checked){ # the name has to be updated
    if(speciesPresentFun(speciesNames$nameSp_checked,species)){ # if the updated name is already present in the tree --> to remove
      speciesToRemoveTree[x] <- species[i]
      x <- x + 1 
    }else{             # if the updated name is not already present in the tree --> to change
      speciesToUpdate[y] <- species[i]
      y <- y + 1
    }
  }
}
length(speciesToUpdate) # 115
length(speciesToRemoveTree) # 45

### remove the species whose updated name is already present:
supertree_modif_cut <- supertree_modif
for(i in 1:numberTree){
  supertree_modif_cut[[i]] <- ape::drop.tip(supertree_modif[[i]],speciesToRemoveTree)
}
numberSpeciesTree <- length(supertree_modif_cut[[1]]$tip.label) # 1502 = 1547 - 45
sum(duplicated(supertree_modif_cut[[1]]$tip.label)) # 0

### update name of species:
supertree_modif_cut_modif <- supertree_modif_cut
numberSpeciesToUpdate <- length(speciesToUpdate) # 115

# check if duplicated names in the updated names of the species that have to be updated:
speciesNameChecked_cut <- speciesNameChecked[speciesNameChecked$nameSp %in% speciesToUpdate,]
sum(duplicated(speciesNameChecked_cut$nameSp_checked)) # 2 --> 2x2 species in the list have the same updated name:
speciesDuplicated_corrected <- speciesNameChecked_cut$nameSp_checked[duplicated(speciesNameChecked_cut$nameSp_checked)]
speciesDuplicated <- speciesNameChecked_cut[speciesNameChecked_cut$nameSp_checked %in% speciesDuplicated_corrected,]
#                 nameSp         nameSp_checked
#   Acanthastrea bowerbanki Homophyllia bowerbanki
#       Acanthastrea hillae Homophyllia bowerbanki
#     Dichocoenia stellaris   Dichocoenia stokesii
#       Dichocoenia stokesi   Dichocoenia stokesii
# --> Remove Acanthastrea hillae and Dichocoenia stellaris
speciesToRemove <- c("Acanthastrea hillae","Dichocoenia stellaris")

# remove from numberSpeciesToUpdate:
speciesToUpdate <- speciesToUpdate[! (speciesToUpdate %in% speciesToRemove)]
numberSpeciesToUpdate <- length(speciesToUpdate) # 113

# remove from supertree_modif_cut_modif:
for(i in 1:numberTree){
  supertree_modif_cut_modif[[i]] <- ape::drop.tip(supertree_modif_cut_modif[[i]],speciesToRemove)
}
numberSpeciesTree <- length(supertree_modif_cut_modif[[1]]$tip.label) # 1500 
sum(duplicated(supertree_modif_cut_modif[[1]]$tip.label)) # 0

# update the rest of species names in supertree_modif_cut_modif
for(i in 1:numberTree){
  for(j in 1:numberSpeciesToUpdate){
    speciesNames <- speciesNameChecked[speciesNameChecked$nameSp == speciesToUpdate[j],]
    supertree_modif_cut_modif[[i]]$tip.label[supertree_modif_cut_modif[[i]]$tip.label == speciesNames$nameSp] <- speciesNames$nameSp_checked  # this changes all the names that ar...
  }
}
length(supertree_modif_cut_modif[[1]]$tip.label) # 1500
sum(duplicated(supertree_modif_cut_modif[[1]]$tip.label)) # 0

class(supertree_modif_cut_modif) # multiphylo
class(supertree_modif_cut_modif[[1]]) # phylo

### create a smaller supertree with only the 1st 100 trees:
supertree_modif_cut_modif_100 <- supertree_modif_cut_modif[1:100]
length(supertree_modif_cut_modif_100) # 100
length(supertree_modif_cut_modif_100[[1]]$tip.label) # 1500
head(supertree_modif_cut_modif_100[[1]]$tip.label)

### write the updated supertree:
# write.tree(supertree_modif_cut_modif_100,paste(wd_Datasets,"supertree_names_validated_100.tre",sep="/"))
# write.tree(supertree_modif_cut_modif,paste("supertree_names_validated.tre",sep="/"))

plot(supertree_modif_cut_modif_100[[1]],type = "fan", cex = .1)

