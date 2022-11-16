# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org for (1) corallite width min, (2) corallite max
## produce CoralliteWidthMin_original.csv, CoralliteWidthMax_original.csv and corallite_area.csv

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_CoralliteArea <- paste(wd,"/Traits_extra/corallite_area",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### Use Max and Min coralitte diameter. When both values are present, the mean between the 2 values is given. Otherwise the 
### avalable measure is chosen.
### the trait "corallite width is not used as only 4 species are present and they are present in Max and Min corallite width table.

#### Import original datasets from coraltraits.org : (April 4 2017)
### corallite width min (mm) --> from Veron, only 1 raw value by species  ####
# CoralliteWidthMin_original <- read.csv("https://coraltraits.org/traits/181.csv", as.is = TRUE)
# write.csv(CoralliteWidthMin_original,paste(wd_Datasets_original,"CoralliteWidthMin_original.csv",sep="/"),row.names=F)

#### corallite max (mm) --> from Veron, only 1 raw value by species  ####
# CoralliteWidthMax_original <- read.csv("https://coraltraits.org/traits/182.csv", as.is = TRUE)
# write.csv(CoralliteWidthMax_original,paste(wd_Datasets_original,"CoralliteWidthMax_original.csv",sep="/"),row.names=F)

### Import original datasets form folder:
dataCoralliteWidthMin_original <- read.csv(paste(wd_Datasets_original,"CoralliteWidthMin_original.csv",sep="/"),header=T,stringsAsFactors = F)
dataCoralliteWidthMax_original <- read.csv(paste(wd_Datasets_original,"CoralliteWidthMax_original.csv",sep="/"),header=T,stringsAsFactors = F)

nameColumnCoralliteWidth <- c("specie_name","value","value_type","precision","standard_unit")
dataCoralliteWidthMinCut <- dataCoralliteWidthMin_original[nameColumnCoralliteWidth]
dataCoralliteWidthMaxCut <- dataCoralliteWidthMax_original[nameColumnCoralliteWidth]

levels(as.factor(dataCoralliteWidthMinCut$value_type)) # "raw_value"
levels(as.factor(dataCoralliteWidthMaxCut$value_type)) # "raw_value"

levels(as.factor(dataCoralliteWidthMinCut$precision)) # 
levels(as.factor(dataCoralliteWidthMaxCut$precision)) # 

levels(as.factor(dataCoralliteWidthMinCut$standard_unit)) # "mm"
levels(as.factor(dataCoralliteWidthMaxCut$standard_unit)) # "mm"

### check if duplicate in species names 
sum(duplicated(dataCoralliteWidthMinCut$specie_name)) # 0
sum(duplicated(dataCoralliteWidthMaxCut$specie_name)) # 0

length(dataCoralliteWidthMinCut$specie_name) # 682
length(dataCoralliteWidthMaxCut$specie_name) # 726

### merge the two tables:
dataFrameCWMin <- data.frame(
  species = dataCoralliteWidthMinCut$specie_name,
  coralliteWidthMin = dataCoralliteWidthMinCut$value,
  stringsAsFactors = F
)
dataFrameCWMax <- data.frame(
  species = dataCoralliteWidthMaxCut$specie_name,
  coralliteWidthMax = dataCoralliteWidthMaxCut$value,
  stringsAsFactors = F
)
tableMerged <- merge(dataFrameCWMin,dataFrameCWMax, all.x = TRUE,all.y = TRUE)
head(tableMerged)
length(tableMerged$species) # 727

### check for duplicates:
sum(duplicated(tableMerged$species)) # 0

# put the species in an alphabetoc order:
tableMerged$species <- gsub(" ","_",tableMerged$species)
tableMerged <- tableMerged[order(tableMerged$species),]
tableMerged$species <- gsub("_"," ",tableMerged$species)

# creation final table: create an additional column that take the mean between the min and max values. If values not present for both max and min, only consider the value given:
tableMerged$meanCoralliteDiameter_mm <- NA
numberSpecies <- length(tableMerged$species) # 727
for(i in 1:numberSpecies){
  minNA <- is.na(tableMerged$coralliteWidthMin[i])
  maxNA <- is.na(tableMerged$coralliteWidthMax[i])
  if(minNA & !maxNA){
    tableMerged$meanCoralliteDiameter_mm[i] <- tableMerged$coralliteWidthMax[i]
  }else if(!minNA & maxNA){
    tableMerged$meanCoralliteDiameter_mm[i] <- tableMerged$coralliteWidthMin[i]
  }else if(!minNA & !maxNA){
    tableMerged$meanCoralliteDiameter_mm[i] <- (tableMerged$coralliteWidthMax[i] + tableMerged$coralliteWidthMin[i])/2
  }else {
    print(paste("What the heck",tableMerged$species[i],"!"))
  }
}
sum(is.na(tableMerged$meanCoralliteDiameter_mm)) # 0

tableNamesChecked <- tableMerged

### check for names:
speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets = wd_Datasets) #  "All species are present in the list"
speciesToCheck

### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp)

### correct names in species_colonialityFinal_cut:
species <- tableNamesChecked$species
numberSpecies <- length(species) # 727
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] # only keep the species in vector species
length(speciesChecked$nameSp) # 727
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}

# Check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 7
speciesDuplicated <- tableNamesChecked$species[duplicated(tableNamesChecked$species)]
subset(tableNamesChecked,tableNamesChecked$species %in% speciesDuplicated)

# keep the value corresponding to coralTraits.org:
rowsToRemove <- c(95,106,346,594,633,671,690)
tableNamesChecked_cut <- tableNamesChecked[-rowsToRemove,]
length(tableNamesChecked_cut$species) # 720
sum(duplicated(tableNamesChecked_cut$species)) # 0

### conversion from mm (diameter) to cm2 (area):
coralliteArea <- tableNamesChecked_cut
coralliteArea$coralliteWidthMin <- pi/400*tableNamesChecked_cut$coralliteWidthMin^2
coralliteArea$coralliteWidthMax <- pi/400*tableNamesChecked_cut$coralliteWidthMax^2
coralliteArea$meanCoralliteDiameter_mm <- pi/400*tableNamesChecked_cut$meanCoralliteDiameter_mm^2
colnames(coralliteArea) <- c("species","coralliteAreaMin","coralliteAreaMax","coralliteAreaMean_cm2")
head(coralliteArea)

# put the species in an alphabetoc order:
coralliteArea$species <- gsub(" ","_",coralliteArea$species)
coralliteArea <- coralliteArea[order(coralliteArea$species),]
coralliteArea$species <- gsub("_"," ",coralliteArea$species)

head(coralliteArea)

coralliteArea <- coralliteArea[c("species","coralliteAreaMean_cm2","coralliteAreaMin","coralliteAreaMax")]

# creation csv file in the wd dataCompilation
# (coralliteArea, file = paste(wd_Datasets,"corallite_area.csv",sep="/"),row.names=F)


