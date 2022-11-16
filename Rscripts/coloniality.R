# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org
## produce coloniality_original.csv and coloniality.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_Coloniality <- paste(wd,"/Traits_extra/coloniality",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### import dataset from coraltraits.org: (April 11 2017)
# coloniality_original <- read.csv("https://coraltraits.org/traits/104.csv", as.is = TRUE)
# write.csv(coloniality_original,paste(wd_Datasets_original,"coloniality_original.csv",sep="/"),row.names = F)

### Import coloniality_original:
coloniality_original <- read.csv(paste(wd_Datasets_original,"coloniality_original.csv",sep="/"),header=T,stringsAsFactors = F)

nameColuns <- c("specie_name","trait_name","standard_unit","value","value_type","precision","precision_type")

colonialityCut <- coloniality_original[nameColuns]

## check if different trait_names:
levels(as.factor(colonialityCut$trait_name))             # coloniality

## check if different values for same species:
length(levels(as.factor(colonialityCut$specie_name)))    # 779
length((colonialityCut$specie_name))                     # 891

## check if different value type:
levels(as.factor(colonialityCut$value_type))             # "expert_opinion"

## check different values:
levels(as.factor(colonialityCut$value))                  # "both"     "colonial" "solitary"

## remove uncessary columns:
nameColuns <- c("specie_name","value")
colonialityCut <- colonialityCut[nameColuns]
head(colonialityCut)

## show species that appear twice in the table and don't have the same value:
duplicatedSpecies <- colonialityCut$specie_name[duplicated(colonialityCut$specie_name)]
numberDuplicates <- length(duplicatedSpecies) # 112
length(levels(as.factor(duplicatedSpecies))) # 112 --> the is no species appearing more than 2 times
for(i in 1:numberDuplicates){
  print(subset(colonialityCut,colonialityCut$specie_name==duplicatedSpecies[i]))
}

# there are many duplicated names but the values are the same, only select the ones with different values:
numberRow <- length(colonialityCut$specie_name) # 891
for(i in 1:numberRow){
  if(i < numberRow){
    if(colonialityCut$specie_name[i] == colonialityCut$specie_name[i+1]){
      if(colonialityCut$value[i] != colonialityCut$value[i+1]){
        print(colonialityCut[i,])
      }
    }
  }
}
# --> none

## remove duplicates:
colonialityFinal <- data.frame(species = as.character(),coloniality = as.character())
colonialityFinal <- rbind(colonialityFinal,colonialityCut[1,])
numberRow <- length(colonialityCut$specie_name) # 891
colonialityFinal <- colonialityCut[!duplicated(colonialityCut$specie_name),]
sum(duplicated(colonialityFinal$specie_name)) # 0
length(colonialityFinal$specie_name) # 779

colnames(colonialityFinal) <- c("species","coloniality")
head(colonialityFinal)

tableNamesChecked <- colonialityFinal

### check for names:
speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets = wd_Datasets)
speciesToCheck
### import list of names to change:
setwd(wd_Datasets)
speciesNameChecked <- read.csv("speciesNameChecked.csv",header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp)

### correct names in species_colonialityFinal_cut:
sp <- tableNamesChecked$species
numberSpecies <- length(sp) # 779
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% sp,] # only keep the species in vector species
length(speciesChecked$nameSp) # 779
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}

### check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 0

## check duplicates:
speciesDuplicated <- tableNamesChecked$species[duplicated(tableNamesChecked$species)]
length(speciesDuplicated) # 33
subset(tableNamesChecked,tableNamesChecked$species %in% speciesDuplicated) # they are all colonial --> remove them

tableNamesChecked_cut <- tableNamesChecked
tableNamesChecked_cut <- subset(tableNamesChecked_cut,!duplicated(tableNamesChecked_cut$species))
sum(duplicated(tableNamesChecked_cut$species)) # 0
length(tableNamesChecked_cut$species) # 746

# creation csv file in the wd dataCompilation
# write.csv(tableNamesChecked_cut, file = paste(wd_Datasets,"coloniality.csv",sep="/"),row.names=F)




