# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org
## produce age_at_maturity_original.csv and age_at_maturity.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_Age_at_maturity <- paste(wd,"Traits_extra/age_at_maturity",sep="/")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### Import original dataset from coraltraits.org: (08 April 2017)
# AgeAtMaturity_original <- read.csv("https://coraltraits.org/traits/47.csv", as.is=TRUE)
# write.csv(AgeAtMaturity_original,paste(wd_Datasets_original,"age_at_maturity_original.csv",sep="/"),row.names = F)

### Import original dataset from folder:
AgeAtMaturity_original <- read.csv(paste(wd_Datasets_original,"age_at_maturity_original.csv",sep="/"),header = T, stringsAsFactors = F)

nameColumnSelectedSAM <- c("specie_name","trait_name","standard_unit","value","value_type","precision","precision_type","replicates","notes")

AgeAtMaturity_cut <- AgeAtMaturity_original[nameColumnSelectedSAM]

## check if different type of traits :
levels(as.factor(AgeAtMaturity_cut$trait_name))

# "Age at maturity" "Generation time"
AgeAtMaturity_cut <- subset(AgeAtMaturity_cut,AgeAtMaturity_cut$trait_name =="Age at maturity")

## check if different units used:
levels(as.factor(AgeAtMaturity_cut$standard_unit))
#years

## check if several values for a same species:
length(levels(as.factor(AgeAtMaturity_cut$specie_name)))
# 3
length(AgeAtMaturity_cut$specie_name)
# 3

ageAtMaturityFinal <- data.frame(
  species = AgeAtMaturity_cut$specie_name,
  ageAtMaturity_yr = as.numeric(AgeAtMaturity_cut$value),
  stringsAsFactors = F
)

### check for names:
tableNamesChecked <- ageAtMaturityFinal

speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets = wd_Datasets) # "All species are present in the list"
speciesToCheck

### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp)

### correct names in species_colonialityFinal_cut:
species <- tableNamesChecked$species
numberSpecies <- length(species) # 3
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] # only keep the species in vector species
length(speciesChecked$nameSp) # 7
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}

### check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 0

### alphabetic order:
tableNamesChecked$species <- gsub(" ","zzz",tableNamesChecked$species)
tableNamesChecked <- tableNamesChecked[order(tableNamesChecked$species),]
tableNamesChecked$species <- gsub("zzz"," ",tableNamesChecked$species)
length(tableNamesChecked$species) # 7
sum(duplicated(tableNamesChecked$species)) # 0

# write.csv(tableNamesChecked,paste(wd_Datasets,"age_at_maturity.csv",sep="/"),row.names = F)
