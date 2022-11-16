# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org
## produce modeLarvalDev_original.csv and mode_larval_development.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_Mode_larval_development <- paste(wd,"/Traits_extra/mode_larval_development",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### Import original dataset from coraltraits.org: (April 8 2017)
# modeLarvalDev_original <- read.csv("https://coraltraits.org/traits/5.csv", as.is = TRUE)
# write.csv(modeLarvalDev_original,paste(wd_Datasets_original,"modeLarvalDev_original.csv",sep="/"),row.names = F)

### Import original dataset from folder:
modeLarvalDev_original <- read.csv(paste(wd_Datasets_original,"modeLarvalDev_original.csv",sep="/"),header = T, stringsAsFactors = F)

nameColumnSelectedlarvalDevModeAll <- c("specie_name","trait_name","value")

modeLarvalDev_original_cut <- modeLarvalDev_original[nameColumnSelectedlarvalDevModeAll]

# check if other traits:
levels(as.factor(modeLarvalDev_original_cut$trait_name))
# "Mode of larval development"

# check if a species appears several times:
sum(duplicated(modeLarvalDev_original_cut$specie_name)) # 0

# ckeck different values:
levels(as.factor(modeLarvalDev_original_cut$value))    # "brooder" "spawner"

modeLarvalDev_original_final <- data.frame(
  species = modeLarvalDev_original_cut$specie_name,
  modeLarvalDevelpment = modeLarvalDev_original_cut$value,
  stringsAsFactors = F
)
length(modeLarvalDev_original_final$species) # 347

### there is one species with a " " in the name:
modeLarvalDev_original_final$species[modeLarvalDev_original_final$species == "Astrangia lajollaensis "] <- "Astrangia lajollaensis"

### check for names:
tableNamesChecked <- modeLarvalDev_original_final

speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets = wd_Datasets) # "All species are present in the list"
speciesToCheck

### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp) # 1767

### correct names in species_colonialityFinal_cut:
species <- tableNamesChecked$species
numberSpecies <- length(species) # 347
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] # only keep the species in vector species
length(speciesChecked$nameSp) # 347
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}

### check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 5
duplicatedSpecies <- tableNamesChecked$species[duplicated(tableNamesChecked$species)]
subset(tableNamesChecked,tableNamesChecked$species %in%  duplicatedSpecies) # same values for all the duplicates --> "spawner"
tableNamesChecked <- subset(tableNamesChecked,!(tableNamesChecked$species %in% duplicatedSpecies))
length(tableNamesChecked$species) # 337
tableNamesChecked <- rbind(tableNamesChecked,data.frame(species=duplicatedSpecies,modeLarvalDevelpment=rep("spawner",length(duplicatedSpecies))))

### Add Pocillopora damicornis as both brooder and spawner:
# In (Castrillón-Cifuentes et al., 2015): "it is believed to be both a brooder and a spawner in the Indo-Pacific, while in the TEP it is only a spawner with synchronous development of gametes (Glynn et al. 1991;Chávez-Romo & Reyes-Bonilla 2007; Carpizo-Ituarte et al. 2011; Rodríguez-Troncoso et al. 2011).
# In (Gleason et al., 2011): "The simultaneous hermaphrodite Pocillopora damicornis is considered to be an internal brooder in most areas but, surprisingly, both broadcasts and broods on Western Australian Reefs (Ward, 1992)"
tableNamesChecked <- rbind(tableNamesChecked,c("Pocillopora damicornis","both"))

### alphabetic order:
tableNamesChecked$species <- gsub(" ","zzz",tableNamesChecked$species)
tableNamesChecked <- tableNamesChecked[order(tableNamesChecked$species),]
tableNamesChecked$species <- gsub("zzz"," ",tableNamesChecked$species)
length(tableNamesChecked$species) # 343
sum(duplicated(tableNamesChecked$species)) # 0

# write.csv(tableNamesChecked,paste(wd_Datasets,"mode_larval_development.csv",sep="/"),row.names = F)

