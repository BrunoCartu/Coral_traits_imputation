# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org
## produce SymbiontDensity_original.csv and Symbiodinium_density.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_zooxanthellae <- paste(wd,"Traits_extra/zooxanthellae",sep="/")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### Import dataset from coraltraits.org: (April 5 2017)
# zooxanthellate_original <- read.csv("https://coraltraits.org/traits/41.csv", as.is = TRUE)
# there is one species with a space at the end: "Astrangia lajollaensis "
# zooxanthellate_original$specie_name[zooxanthellate_original$specie_name == "Astrangia lajollaensis "] <- "Astrangia lajollaensis"
# write.csv(zooxanthellate_original,paste(wd_Datasets_original,"zooxanthellate_original.csv",sep="/"),row.names = F)

### Import dataset from folder:
zooxanthellate_original <- read.csv(paste(wd_Datasets_original,"zooxanthellate_original.csv",sep="/"),header = T, stringsAsFactors = F)

nameColumn <- c("specie_name","trait_name","value")

zooxanthellate <- zooxanthellate_original[nameColumn]

# check if different traits:
levels(as.factor(zooxanthellate$trait_name))
# "Zooxanthellate"

# check different values of trait: 
levels(as.factor(zooxanthellate$value))
# "azooxanthellate" "both"  "zooxanthellate"

# check for duplicate species:        
sum(duplicated(zooxanthellate$specie_name)) # 0

zooxanthellate <- zooxanthellate[,c("specie_name","value")]
colnames(zooxanthellate) <- c("species","zooanthellae")

### check for names:
tableNamesChecked <- zooxanthellate
speciesToCheck <- checkSpeciesNames(tableNamesChecked$species,wd_Datasets = wd_Datasets) 
speciesToCheck 

### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep = "/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp)

### correct names in species_colonialityFinal_cut:
species <- tableNamesChecked$species
numberSpecies <- length(species) # 1547
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] # only keep the species in vector species
length(speciesChecked$nameSp) # 1547
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}

# Check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 70
speciesDuplicated <- tableNamesChecked$species[duplicated(tableNamesChecked$species)]
speciesDuplicated <- levels(as.factor((speciesDuplicated)))
length(speciesDuplicated) # 61

subset(tableNamesChecked,tableNamesChecked$species %in% speciesDuplicated)
# remove duplicate rows from main table, and place then species in a new table. The zooxanthellae value is taken from zooxanthellate_original
tableNamesChecked_backup <- tableNamesChecked
nbSpDuplicated <- length(speciesDuplicated)     # 61
speciesCorrected <- data.frame(species = rep(NA,nbSpDuplicated),zooanthellae=rep(NA,nbSpDuplicated),stringsAsFactors = F)
for(i in 1:nbSpDuplicated){
  speciesCorrected$species[i] <- speciesDuplicated[i]
  original_species_name <- subset(speciesChecked$nameSp,speciesChecked$nameSp_checked == speciesDuplicated[i])
  original_species_name <- subset(zooxanthellate_original$specie_name,zooxanthellate_original$specie_name %in% original_species_name)[1] # main species names have had added the subgenus name so they are no longer present in original_species_name, so need to go in speciesChecked
  speciesCorrected$zooanthellae[i] <- subset(zooxanthellate_original$value,zooxanthellate_original$specie_name == original_species_name)
  tableNamesChecked_backup <- subset(tableNamesChecked_backup,tableNamesChecked_backup$species != speciesDuplicated[i])
}
tableNamesChecked_backup <- rbind(tableNamesChecked_backup,speciesCorrected)
length(speciesCorrected$species) # 61
length(tableNamesChecked_backup$species) # 1477
sum(duplicated(tableNamesChecked_backup$species)) # 0

tableNamesChecked_backup$species <- gsub(" ","_",tableNamesChecked_backup$species)
tableNamesChecked_backup <- tableNamesChecked_backup[order(tableNamesChecked_backup$species),]
tableNamesChecked_backup$species <- gsub("_"," ",tableNamesChecked_backup$species)

# write.csv(tableNamesChecked_backup,paste(wd_Datasets,"zooxanthellate.csv",sep="/"),row.names = F)

zooxan <- read.csv(paste(wd_Datasets,"zooxanthellate.csv",sep="/"),header = T, stringsAsFactors = F)
length(zooxan$species) # 1477
levels(as.factor(zooxan$zooanthellae)) # "azooxanthellate" "both"            "zooxanthellate" 
length(zooxan$species[zooxan$zooanthellae == "zooxanthellate"]) # 816
length(zooxan$species[zooxan$zooanthellae == "both"]) # 12


