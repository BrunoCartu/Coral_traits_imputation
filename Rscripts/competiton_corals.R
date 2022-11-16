## Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate Precoda_et_al_2017's species-pair interction probability
## produce coral_pair_competition_probability.csv and coral_pair_competition_probability_names.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### Import output species-pair interction probability of Precoda_et_al_2017
Probability_original <- read.csv(paste(wd_Datasets_original,"Precoda_et_al_2017_model_output.csv",sep="/"),header=T,stringsAsFactors = F)
length(Probability_original$species1) # 299925 = (774 * 773) / 2 + 774 
length(levels(as.factor(Probability_original$species1))) # 774
length(levels(as.factor(Probability_original$species2))) # 774
length(speciesPresentBothList(levels(as.factor(Probability_original$species1)),levels(as.factor(Probability_original$species2)))) # 774 --> same species names in both column

### import my personal list of coral species names updated:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp)
## check if there are species names not in my list yet:
species <- as.character(levels(as.factor(Probability_original$species1)))
speciesToCheck <- checkSpeciesNames(species, wd_Datasets = wd_Datasets)  #  "All species are present in the list"
## only select the species present in the DF:
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] 
numberSpecies <- length(speciesChecked$nameSp) # 774
## create a DF for checking species that need to be corrected and if they are already present in the list:
speciesCheckDF <- data.frame(
  speciesName = speciesChecked$nameSp,
  speciesNameCorrected = speciesChecked$nameSp_checked,
  toChange = rep("NO",length(species)),
  newNameAlreadyPresent = rep(NA,length(species)),
  stringsAsFactors = F
)
head(speciesCheckDF)
for(i in 1:numberSpecies){
  if(speciesCheckDF$speciesName[i] != speciesCheckDF$speciesNameCorrected[i]){
    speciesCheckDF$toChange[i] <- "YES"
    if(speciesCheckDF$speciesNameCorrected[i] %in% speciesCheckDF$speciesName){
      speciesCheckDF$newNameAlreadyPresent[i] <- "YES"
    }else{
      speciesCheckDF$newNameAlreadyPresent[i] <- "NO"
    }
  }
}
head(speciesCheckDF)
## number names to correct:
sum(speciesCheckDF$toChange == "YES") # 63
## number of species name correcred that are already present in the list:
sum(na.omit(speciesCheckDF$newNameAlreadyPresent == "YES")) # 33

### correct for species names whose updated name is not alredy present in the list:
Probability_nameCorrected <- Probability_original
speciesCheckDF_cut <- speciesCheckDF[speciesCheckDF$toChange == "YES",]
speciesCheckDF_cut <- speciesCheckDF_cut[speciesCheckDF_cut$newNameAlreadyPresent == "NO",]
speciesNumber <- length(speciesCheckDF_cut$speciesName) # 30
for(i in 1:speciesNumber){
  Probability_nameCorrected$species1[Probability_nameCorrected$species1 == speciesCheckDF_cut$speciesName[i]] <- speciesCheckDF_cut$speciesNameCorrected[i]
  Probability_nameCorrected$species2[Probability_nameCorrected$species2 == speciesCheckDF_cut$speciesName[i]] <- speciesCheckDF_cut$speciesNameCorrected[i]
}

### remove the species for which the updated name is already present in the table:
Probability_nameCorrected_cut <- Probability_nameCorrected
speciesCheckDF_cut_2 <- speciesCheckDF[speciesCheckDF$toChange == "YES",]
speciesCheckDF_cut_2 <- speciesCheckDF_cut_2[speciesCheckDF_cut_2$newNameAlreadyPresent == "YES",]
speciesNumber <- length(speciesCheckDF_cut_2$speciesName) # 33
## remove rows in Probability_nameCorrected_cut$species1 or species2 with species already present
length(levels(as.factor(Probability_nameCorrected_cut$species1))) # 774
for(i in 1:speciesNumber){
  Probability_nameCorrected_cut <- Probability_nameCorrected_cut[Probability_nameCorrected_cut$species1 != speciesCheckDF_cut_2$speciesName[i], ]
  Probability_nameCorrected_cut <- Probability_nameCorrected_cut[Probability_nameCorrected_cut$species2 != speciesCheckDF_cut_2$speciesName[i], ]
}
length(levels(as.factor(Probability_nameCorrected_cut$species1))) # 741 = 774 - 33 , all good
length(levels(as.factor(Probability_nameCorrected_cut$species2))) # 741 = 774 - 33 , all good

### write csv
# write.csv(Probability_nameCorrected_cut,paste(wd_Datasets,"coral_pair_competition_probability.csv",sep="/"),row.names = F)

### write csv with only the species names present;
species <- data.frame(species = levels(as.factor(Probability_nameCorrected_cut$species1)))
length(species$species)
sum(duplicated(species$species))
# write.csv(species,paste(wd_Datasets,"coral_pair_competition_probability_names.csv",sep="/"),row.names = F)






