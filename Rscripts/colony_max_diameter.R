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
wd_ColMaxDiam <- paste(wd,"/Traits_extra/colony_max_diameter",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

####### Import dataset from coraltraits.org: (April 3 2017)
# colonyMaxDiameter_orginal <- read.csv("https://coraltraits.org/traits/90.csv", as.is = TRUE)
# write.csv(colonyMaxDiameter_orginal,paste(wd_Datasets_original,"colonyMaxDiameter_orginal.csv",sep="/"),row.names = F)

### import orginal dataset from foler:
colonyMaxDiameter_orginal <- read.csv(paste(wd_Datasets_original,"colonyMaxDiameter_orginal.csv",sep="/"),header=T,stringsAsFactors = F)

nameColumnSelectedColonyMaxDiamter <- c("specie_name","trait_name","standard_unit","value","value_type","precision","precision_type")

dataColonyMaxDiameterCut <- colonyMaxDiameter_orginal[nameColumnSelectedColonyMaxDiamter]

## check if there are other traits measured:
levels(as.factor(dataColonyMaxDiameterCut$trait_name))    # "Colony maximum diameter" "Water depth"  
# --> remove "Water depth"  :
dataColonyMaxDiameterCut <- subset(dataColonyMaxDiameterCut,dataColonyMaxDiameterCut$trait_name == "Colony maximum diameter")

## check if several values for a same species:
length(levels(as.factor(dataColonyMaxDiameterCut$specie_name)))    # 320
length(as.factor(dataColonyMaxDiameterCut$specie_name))            # 532

## check if different value type:
levels(as.factor(dataColonyMaxDiameterCut$value_type))             # "maximum"     "raw_value"
# --> now need to consider sample size 

## check if different unit of measurement:
levels(as.factor(dataColonyMaxDiameterCut$standard_unit))      # "cm" "mm"
# --> conversion of mm in cm:
dataColonyMaxDiameterCut$value[dataColonyMaxDiameterCut$standard_unit == "mm"] <- dataColonyMaxDiameterCut$value[dataColonyMaxDiameterCut$standard_unit == "mm"] /10
dataColonyMaxDiameterCut$standard_unit[dataColonyMaxDiameterCut$standard_unit == "mm"] <- "cm"

### for duplicated names, I am only keeping the maximum value:
meanColonyMaxDiameter <- tapply(dataColonyMaxDiameterCut$value,dataColonyMaxDiameterCut$specie_name,max)
length(meanColonyMaxDiameter) # 320
colonyMaxDiameterFinal <- data.frame(
  species = as.character(rownames(meanColonyMaxDiameter)),
  colonyMaxDiameter_cm =as.numeric(meanColonyMaxDiameter),
  stringsAsFactors = F
  )
length(colonyMaxDiameterFinal$species) # 320

tableNamesChecked <- colonyMaxDiameterFinal

### check for names:
speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets = wd_Datasets) # "All species are present in the list"
speciesToCheck

### import list of names to change:
setwd(wd_Datasets)
speciesNameChecked <- read.csv("speciesNameChecked.csv",header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp)

### correct names in species_colonialityFinal_cut:
species <- tableNamesChecked$species
numberSpecies <- length(species) # 320
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] # only keep the species in vector species
length(speciesChecked$nameSp) # 320
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}

### alphabetic order:
tableNamesChecked$species <- gsub(" ","zzz",tableNamesChecked$species)
tableNamesChecked <- tableNamesChecked[order(tableNamesChecked$species),]
tableNamesChecked$species <- gsub("zzz"," ",tableNamesChecked$species)

### check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 11
duplicatedSpecies <- tableNamesChecked$species[duplicated(tableNamesChecked$species)]
subset(tableNamesChecked,tableNamesChecked$species %in%  duplicatedSpecies)

### for duplicated names, I am only keeping the maximum value:
meanColonyMaxDiameter <- tapply(tableNamesChecked$colonyMaxDiameter_cm,tableNamesChecked$species,max)
length(meanColonyMaxDiameter) # 309
colonyMaxDiameterFinal <- data.frame(
  species = as.character(rownames(meanColonyMaxDiameter)),
  colonyMaxDiameter_cm =as.numeric(meanColonyMaxDiameter),
  stringsAsFactors = F
)
length(colonyMaxDiameterFinal$species) # 309
sum(duplicated(colonyMaxDiameterFinal$species)) # 0

# creation csv file
# write.csv(colonyMaxDiameterFinal, file = paste(wd_Datasets,"colony_max_diameter.csv",sep="/"),row.names=F)
