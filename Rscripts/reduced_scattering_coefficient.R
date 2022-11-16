# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## currate and combine RSC values from (Marcelino et al., 2013, Table_S1) and (Swain et al., 2016)
## produce reduced_scattering_coefficient.csv 
# using:
## Marcelino_et_al_2013_Table_S1.csv
## 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_RSC <- paste(wd,"/Traits_extra/reduced_scattering_coefficient",sep="")

source(paste(wd_Rscripts,"/functions.R",sep=""))

# Import dataset from Marcelino et al., 2013_Table_S1 on Reduced Scattering Coefficient (?'S,m)
## (averaged over 20 measurements per skeleton)
RSC <- read.csv(paste(wd_Datasets_original,"Marcelino_et_al_2013_Table_S1.csv",sep="/"),header=T,stringsAsFactors = F)
head(RSC)

colnames(RSC)[1]  <- "species"

### correct for certain names under column species: 
# Acropora intermedia/nobilis <- Acropora intermedia  (nobilis is not accepted and should be robusta)
# Acropora sp. -> to remove
# Euphyllia sp. -> to remove
# Montipora sp. -> to remove
# Mycedium elephantotus/mancaoi <- Mycedium elephantotus (both are correct names but I have more information for the 1st name)
# Porites sp.  -> to remove
# Siderastrea sp. -> to remove

RSC$species[RSC$species == "Acropora intermedia/nobilis"] <- "Acropora intermedia"
RSC$species[RSC$species == "Mycedium elephantotus/mancaoi"] <- "Mycedium elephantotus"

length(RSC$species)  # 150

RSC <- subset(RSC,species != "Acropora sp.")
RSC <- subset(RSC,species != "Euphyllia sp.")
RSC <- subset(RSC,species != "Montipora sp.")
RSC <- subset(RSC,species != "Porites sp.")
RSC <- subset(RSC,species != "Siderastrea sp.")

length(RSC$species)  # 143

### There are different values for certain species that correspond to different coral fragment --> I averaged them
### mean and sd:
speciesNames <- levels(as.factor(RSC$species))
numberSpecies <- length(speciesNames)  # 89
RSC_mean <- as.data.frame(matrix(ncol = 4,nrow = numberSpecies))
colnames(RSC_mean) <- c("species","mean_RSC","sd","n")
RSC_mean[,1] <- speciesNames
RSC_mean[2] <- round(tapply(RSC$RSC,RSC$species,mean),2)
RSC_mean[3] <- round(tapply(RSC$RSC,RSC$species,sd),3)
RSC_mean[4] <- table(RSC$species)

length(RSC_mean$species) # 89

### Add species from Swain et al., 2016
speciesFromSwain2016 <- c("Diploria labyrinthiformis","Turbinaria reniformis","Favia favus","Montipora foliosa","Montipora digitata")
RSC_Swain <- c(3.39,3.48,3.98,4.03,3.92)
speciesFromSwain2016  <- data.frame(
  species = speciesFromSwain2016,
  mean_RSC = RSC_Swain,
  sd = NA,
  n = 1,
  stringsAsFactors = F
)

RSC_mean <- rbind(RSC_mean,speciesFromSwain2016)
length(RSC_mean$species) # 94

### check for names:
tableNamesChecked <- RSC_mean
speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets = wd_Datasets)
speciesToCheck 

### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp) # 

### correct names in species_colonialityFinal_cut:
species <- tableNamesChecked$species
numberSpecies <- length(species) # 94
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] # only keep the species in vector species
length(speciesChecked$nameSp) # 94
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}

# Check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 1
speciesDuplicated <- tableNamesChecked$species[duplicated(tableNamesChecked$species)]
speciesDuplicated <- levels(as.factor((speciesDuplicated)))
length(speciesDuplicated) # 1

subset(tableNamesChecked,tableNamesChecked$species %in% speciesDuplicated) # Psammocora contigua

### remove it from tableNamesChecked
tableNamesChecked <- subset(tableNamesChecked,tableNamesChecked$species != speciesDuplicated)
length(tableNamesChecked$species) # 92

### get values from RSC and calculate the mean 
speciesOriginalData <- speciesChecked$nameSp[speciesChecked$nameSp_checked == "Psammocora contigua"]
P.contigua <- RSC[RSC$species %in%  speciesOriginalData,] # there is only one value for each name:
P.contigua <- data.frame(
  species = "Psammocora contigua",
  mean_RSC = round(mean(P.contigua$RSC),2),
  sd = round(sd(P.contigua$RSC),2),
  n = 2,
  stringsAsFactors = F
)

### add P.contigua to tableNamesChecked
tableNamesChecked <- rbind(tableNamesChecked,P.contigua)

tableNamesChecked$species <- gsub(" ","_",tableNamesChecked$species)
tableNamesChecked <- tableNamesChecked[order(tableNamesChecked$species),]
tableNamesChecked$species <- gsub("_"," ",tableNamesChecked$species)
length(tableNamesChecked$species) # 93
sum(duplicated(tableNamesChecked$species))

# write.csv(tableNamesChecked,paste(wd_Datasets,"reduced_scattering_coefficient.csv",sep="/"),row.names=F)



