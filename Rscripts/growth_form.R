# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org
## produce GrowthFormTypical_original.csv and growth_form.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_growthForm <- paste(wd,"/Traits_extra/growth_form",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### Import dataset from coraltraits.org: trait: GROWTH FROM TYPICAL (April 4 2017)
# GrowthFormTypical_original <- read.csv("https://coraltraits.org/traits/183.csv", as.is = TRUE)
# write.csv(GrowthFormTypical_original,paste(wd_Datasets_original,"GrowthFormTypical_original.csv",sep="/"),row.names = F)

### Import dataset from folder:
GrowthFormTypical_original <- read.csv(paste(wd_Datasets_original,"GrowthFormTypical_original.csv",sep="/"),header=T,stringsAsFactors = F)
sum(duplicated(GrowthFormTypical_original$specie_name)) # 0
nameColumn <- c("specie_name","trait_name","value")

growthFormTypical <- GrowthFormTypical_original[nameColumn]

## check if different traits:
levels(as.factor(growthFormTypical$trait_name))   # "Growth form typical"

## check the different trait values:
levels(as.factor(growthFormTypical$value))
# "branching_closed"         "branching_open"           "columnar"                 "corymbose"                "digitate"                 "encrusting"               "encrusting_long_uprights" "hispidose"                "laminar"                 
# "massive"                  "submassive"               "tables_or_plates"  

## cleck is duplicates in species names:
sum(duplicated(growthFormTypical$specie_name)) # 0
length(growthFormTypical$specie_name)          # 857

growthFormTypical <- growthFormTypical[,c(1,3)]

colnames(growthFormTypical) <- c("species","growthForm")

tableNamesChecked <- growthFormTypical

### check for names:
speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets = wd_Datasets) #  "All species are present in the list"
speciesToCheck

### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp)

### correct names in species_colonialityFinal_cut:
species <- tableNamesChecked$species
numberSpecies <- length(species) # 857
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] # only keep the species in vector species
length(speciesChecked$nameSp) # 857
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}

# Check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 40
speciesDuplicated <- tableNamesChecked$species[duplicated(tableNamesChecked$species)]
speciesDuplicated <- levels(as.factor((speciesDuplicated)))
length(speciesDuplicated) # 35
subset(tableNamesChecked,tableNamesChecked$species %in% speciesDuplicated)
# remove duplicate rows from main table, and place then species in a new table. The growth form value is taken from GrowthFormTypical_original
tableNamesChecked_backup <- tableNamesChecked
nbSpDuplicated <- length(speciesDuplicated)     # 35
speciesCorrected <- data.frame(species = rep(NA,nbSpDuplicated),growthForm=rep(NA,nbSpDuplicated),stringsAsFactors = F)
for(i in 1:nbSpDuplicated){
  speciesCorrected$species[i] <- speciesDuplicated[i]
  speciesCorrected$growthForm[i] <- subset(GrowthFormTypical_original$value,GrowthFormTypical_original$specie_name == speciesDuplicated[i])
  tableNamesChecked_backup <- subset(tableNamesChecked_backup,tableNamesChecked_backup$species != speciesDuplicated[i])
}
tableNamesChecked_backup <- rbind(tableNamesChecked_backup,speciesCorrected)
length(speciesCorrected$species) # 35
length(tableNamesChecked_backup$species) # 817
sum(duplicated(tableNamesChecked_backup$species)) # 0

tableNamesChecked_backup$species <- gsub(" ","_",tableNamesChecked_backup$species)
tableNamesChecked_backup <- tableNamesChecked_backup[order(tableNamesChecked_backup$species),]
tableNamesChecked_backup$species <- gsub("_"," ",tableNamesChecked_backup$species)

### add a column for growth form:
# branching: branching_closed, branching_open, hispidose
# massive: massive, submassive
# the reste does not change
FT_table_modif <- tableNamesChecked_backup
FT_table_modif$growthFormClasses <- NA
numberSpecies <- length(FT_table_modif$species) # 817
for(i in 1:numberSpecies){
  if(FT_table_modif$growthForm[i] == "branching_closed" | FT_table_modif$growthForm[i] == "branching_open" | FT_table_modif$growthForm[i] == "hispidose"){
    FT_table_modif$growthFormClasses[i] <- "branching"
  }else if(FT_table_modif$growthForm[i] == "massive" | FT_table_modif$growthForm[i] == "submassive"){
    FT_table_modif$growthFormClasses[i] <- "massive"
  }else{
    FT_table_modif$growthFormClasses[i] <- FT_table_modif$growthForm[i]
  }
}

### 
# write.csv(FT_table_modif,paste(wd_Datasets,"growth_form.csv",sep="/"),row.names = F)

