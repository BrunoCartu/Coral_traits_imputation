# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org
## produce geographical_Regions_original.csv and geographical_region.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")

source(paste(wd_Rscripts,"functions.R",sep="/"))

#### Import original datasets from coraltraits.org : (July 19 2017)
# geoRegion_original <- read.csv("https://coraltraits.org/traits/35.csv", as.is = TRUE)
# write.csv(geoRegion_original,paste(wd_Datasets_original,"geographical_Regions_original.csv",sep="/"),row.names=F)

### Import original datasets form folder:
geoRegion_original <- read.csv(paste(wd_Datasets_original,"geographical_Regions_original.csv",sep="/"),header=T,stringsAsFactors = F)
# View(geoRegion_original)
colnames(geoRegion_original)

colTeeKeep <- c("specie_name","trait_name","value","value_type")
geoRegion_cut <- geoRegion_original[,colTeeKeep]

levels(as.factor(geoRegion_cut$trait_name)) # "Geographical region"
levels(as.factor(geoRegion_cut$value)) 
# "Eastern Atlantic"            "Eastern Pacific"             "Indian Ocean"                "Subantarctic and Antarctic"  "Western and Central Pacific" "Western Atlantic"
length(levels(as.factor(geoRegion_cut$specie_name))) # 1508
length(geoRegion_cut$specie_name) # 2316

geoRegion_cut <- geoRegion_cut[,c("specie_name","value")]
colnames(geoRegion_cut) <- c("species","geographical_region")

# one of the species name is "Astrangia lajollaensis " --> remove space at the end:
geoRegion_cut$species[geoRegion_cut$species == "Astrangia lajollaensis "] <- "Astrangia lajollaensis"

# some species belong to different regions: make a table of presence / absence by regions
geoRegion_presenceAbs <- data.frame(
  species = levels(as.factor(geoRegion_cut$species)),
  Eastern_Atlantic = rep(0,length(levels(as.factor(geoRegion_cut$species)))),
  Eastern_Pacific = rep(0,length(levels(as.factor(geoRegion_cut$species)))),
  Indian_Ocean = rep(0,length(levels(as.factor(geoRegion_cut$species)))),
  Subantarctic_and_Antarctic = rep(0,length(levels(as.factor(geoRegion_cut$species)))),
  Western_and_Central_Pacific = rep(0,length(levels(as.factor(geoRegion_cut$species)))),
  Western_Atlantic = rep(0,length(levels(as.factor(geoRegion_cut$species)))),
  stringsAsFactors = F
)

### fill up the table:
numberRow <- length(geoRegion_cut$species) # 2316
for(i in 1:numberRow){
  Eastern_Atlantic = geoRegion_cut$geographical_region[i] == "Eastern Atlantic"
  Eastern_Pacific = geoRegion_cut$geographical_region[i] == "Eastern Pacific"
  Indian_Ocean = geoRegion_cut$geographical_region[i] == "Indian Ocean"
  Subantarctic_and_Antarctic = geoRegion_cut$geographical_region[i] == "Subantarctic and Antarctic"
  Western_and_Central_Pacific = geoRegion_cut$geographical_region[i] == "Western and Central Pacific"
  Western_Atlantic =  geoRegion_cut$geographical_region[i] == "Western Atlantic"
  
  speciesHere <- geoRegion_cut$species[i]
  
  if(Eastern_Atlantic){
    geoRegion_presenceAbs$Eastern_Atlantic[geoRegion_presenceAbs$species == speciesHere] <- 1
  }else if(Eastern_Pacific){
    geoRegion_presenceAbs$Eastern_Pacific[geoRegion_presenceAbs$species == speciesHere] <- 1
  }else if(Indian_Ocean){
    geoRegion_presenceAbs$Indian_Ocean[geoRegion_presenceAbs$species == speciesHere] <- 1
  }else if(Subantarctic_and_Antarctic){
    geoRegion_presenceAbs$Subantarctic_and_Antarctic[geoRegion_presenceAbs$species == speciesHere] <- 1
  }else if(Western_and_Central_Pacific){
    geoRegion_presenceAbs$Western_and_Central_Pacific[geoRegion_presenceAbs$species == speciesHere] <- 1
  }else if(Western_Atlantic){
    geoRegion_presenceAbs$Western_Atlantic[geoRegion_presenceAbs$species == speciesHere] <- 1
  }
}
head(geoRegion_presenceAbs)
length(geoRegion_presenceAbs$species) # 1508

### check for names:
tableNamesChecked <- geoRegion_presenceAbs
speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets = wd_Datasets) #  "All species are present in the list"
speciesToCheck

### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp) 

### correct names in species_colonialityFinal_cut:
species <- tableNamesChecked$species
numberSpecies <- length(species) # 1508
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] # only keep the species in vector species
length(speciesChecked$nameSp) # 1508
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}
# Check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 70
speciesDuplicated <- tableNamesChecked$species[duplicated(tableNamesChecked$species)]
subset(tableNamesChecked,tableNamesChecked$species %in% speciesDuplicated)
sum(duplicated(speciesDuplicated)) # some species are treplicated or more
speciesDuplicated <- levels(as.factor(speciesDuplicated))

### replace 0 by 1 when conflect
tableNamesChecked_withoutDuplicatedSp <- subset(tableNamesChecked,!tableNamesChecked$species %in% speciesDuplicated)
length(tableNamesChecked_withoutDuplicatedSp$species) # 1377

numberSpeciesDuplicated <- length(speciesDuplicated) # 61

geoRegion_presenceAbs_duplicates <- data.frame(
  species = speciesDuplicated,
  Eastern_Atlantic = rep(0,numberSpeciesDuplicated),
  Eastern_Pacific = rep(0,numberSpeciesDuplicated),
  Indian_Ocean = rep(0,numberSpeciesDuplicated),
  Subantarctic_and_Antarctic = rep(0,numberSpeciesDuplicated),
  Western_and_Central_Pacific = rep(0,numberSpeciesDuplicated),
  Western_Atlantic = rep(0,numberSpeciesDuplicated),
  stringsAsFactors = F
)

for(i in 1:numberSpeciesDuplicated){
  subsettable <- subset(tableNamesChecked,tableNamesChecked$species == speciesDuplicated[i])
  if(sum(subsettable$Eastern_Atlantic) > 0){
    geoRegion_presenceAbs_duplicates$Eastern_Atlantic[geoRegion_presenceAbs_duplicates$species == speciesDuplicated[i] ] <- 1
  }
  if(sum(subsettable$Eastern_Pacific) > 0){
    geoRegion_presenceAbs_duplicates$Eastern_Pacific[geoRegion_presenceAbs_duplicates$species == speciesDuplicated[i] ] <- 1
  }
  if(sum(subsettable$Indian_Ocean) > 0){
    geoRegion_presenceAbs_duplicates$Indian_Ocean[geoRegion_presenceAbs_duplicates$species == speciesDuplicated[i] ] <- 1
  }
  if(sum(subsettable$Subantarctic_and_Antarctic) > 0){
    geoRegion_presenceAbs_duplicates$Subantarctic_and_Antarctic[geoRegion_presenceAbs_duplicates$species == speciesDuplicated[i] ] <- 1
  }
  if(sum(subsettable$Western_and_Central_Pacific) > 0){
    geoRegion_presenceAbs_duplicates$Western_and_Central_Pacific[geoRegion_presenceAbs_duplicates$species == speciesDuplicated[i] ] <- 1
  }
  if(sum(subsettable$Western_Atlantic) > 0){
    geoRegion_presenceAbs_duplicates$Western_Atlantic[geoRegion_presenceAbs_duplicates$species == speciesDuplicated[i] ] <- 1
  }
}

### combine the 2 tables:
tableNamesChecked <- rbind(tableNamesChecked_withoutDuplicatedSp,geoRegion_presenceAbs_duplicates)
# put the species in an alphabetoc order:
tableNamesChecked$species <- gsub(" ","_",tableNamesChecked$species)
tableNamesChecked <- tableNamesChecked[order(tableNamesChecked$species),]
tableNamesChecked$species <- gsub("_"," ",tableNamesChecked$species)

length(tableNamesChecked$species) # 1438
sum(duplicated(tableNamesChecked$species)) # 0

name <- "Acropora cerealis"
name <- sample(tableNamesChecked$species,1)
subset(tableNamesChecked,tableNamesChecked$species == name )
subset(geoRegion_presenceAbs_duplicates,geoRegion_presenceAbs_duplicates$species == name )

# creation csv file in the wd dataCompilation
# write.csv(tableNamesChecked,paste(wd_Datasets,"geographical_region.csv",sep="/"),row.names=F)

