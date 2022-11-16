# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org
## produce sexual_system_original.csv and sexual_system.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_sexual_system <- paste(wd,"/Traits_extra/sexual_system",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### Final result: sexual system: hermaphrodite vs. gonochore (well preserved phylogenetically according to Baird et al., 2009)

### import original dataset from caraltraits.org: (April 18 2018)
# dataSexualSystem_original <- read.csv("https://coraltraits.org/traits/8.csv", as.is=TRUE)
# write.csv(dataSexualSystem_original,paste(wd_Datasets_original,"sexual_system_original.csv",sep="/"),row.names = F)

### import original dataset from folder
dataSexualSystem_original <- read.csv(paste(wd_Datasets_original,"sexual_system_original.csv",sep="/"),header=T,stringsAsFactors = F)
head(dataSexualSystem_original)
names(dataSexualSystem_original)[5] <- "species_name"

nameColumnSelectedSS <- c("species_name","trait_name","standard_id","value")

dataSSCutColumn <- dataSexualSystem_original[nameColumnSelectedSS]
head(dataSSCutColumn)

# check the different traits:
levels(as.factor(dataSSCutColumn$trait_name))

# number of lines :
length(dataSSCutColumn$species_name)    # 340

# different values
levels(as.factor(dataSSCutColumn$value)) # "gonochore"     "hermaphrodite"

### check if there are duplicate species. If so, check if same value. If so, romve the duplicate. If not, show which species
length(levels(as.factor(dataSSCutColumn$species_name))) # 340 --> no duplicate 
sum(duplicated(dataSSCutColumn$species_name)) # 0

# remvoe potential white space at the beginning and end of species names:
dataSSCutColumn$species_name <- trimws(dataSSCutColumn$species_name)

### check for names:
tableNamesChecked <- dataSSCutColumn
speciesToCheck <- checkSpeciesNames(tableNamesChecked$species_name, wd_Datasets = wd_Datasets) # "All species are present in the list

### remove the duplicated species and update the species names for species whose updated name is not already present
speciestoRemove <- c()
x <- 0
for(i in 1:length(speciesToCheck$nameSp)){
  if(speciesPresentFun(speciesToCheck$nameSp_checked[[i]],dataSSCutColumn$species_name)){ # the correct name is already present in dataSSCutColumn$species_name
    x <- x + 1
    speciestoRemove[x] <- speciesToCheck$nameSp[i]
  }else{ # update name
    dataSSCutColumn$species_name[dataSSCutColumn$species_name == speciesToCheck$nameSp[i]] <- speciesToCheck$nameSp_checked[i]
    print(paste("Species",speciesToCheck$nameSp[i]," changed to ",speciesToCheck$nameSp_checked[i]))
  }
}
length(speciestoRemove) # 6

### remove the speciestoRemove
dataSSCutColumn_cut <- dataSSCutColumn[!dataSSCutColumn$species_name %in% speciestoRemove,]
length(dataSSCutColumn_cut$species_name) # 334
sum(duplicated(dataSSCutColumn_cut$species_name)) # 0

names(dataSSCutColumn_cut)[1] <- "species"
names(dataSSCutColumn_cut)[4] <- "sexual_system"

### alphabetic order:
dataSSCutColumn_cut$species <- gsub(" ","_",dataSSCutColumn_cut$species)
dataSSCutColumn_cut <- dataSSCutColumn_cut[order(dataSSCutColumn_cut$species),]
dataSSCutColumn_cut$species <- gsub("_"," ",dataSSCutColumn_cut$species)

# creation csv file in the wd dataCompilation
# write.csv(dataSSCutColumn_cut[,c(1,4)], file = paste(wd_Datasets,"sexual_system.csv",sep="/"),row.names=F)

