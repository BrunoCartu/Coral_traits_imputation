# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org and from Figueiredo et al., 2013, Harriot 1982 
## produce EggDiameter_original.csv, EggDiameter_original_2.csv and egg_diameter.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_Egg_diameter <- paste(wd,"/Traits_extra/egg_diameter",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### Import "egg size" (volume in mm3) orginal dataset from the coraltraits.org: (April 8 2017)
# EggDiameter_original <- read.csv("https://coraltraits.org/traits/217.csv", as.is=TRUE)
# write.csv(EggDiameter_original,paste(wd_Datasets_original,"EggDiameter_original.csv",sep="/"),row.names = F)

### Import orginal dataset from folder:
EggDiameter_original <- read.csv(paste(wd_Datasets_original,"EggDiameter_original.csv",sep="/"),header = T, stringsAsFactors = F)

nameColumnSelected <- c("specie_name","trait_name","value","value_type","precision","precision_type","standard_unit")

EggDiameter_CutColumn <- EggDiameter_original[nameColumnSelected]

## check trait names:
levels(as.factor(EggDiameter_CutColumn$trait_name))
# "Colony area"  "Egg size"  "Habitat type"  "Ratio of female to male gonads"  "Season"  "Water depth"  
# only select "Egg size":
EggDiameter_CutColumn <- subset(EggDiameter_CutColumn, trait_name == "Egg size")       

## check if different standard unit:
levels(factor(EggDiameter_CutColumn$standard_unit))          # "mm^3"

## check if different value type:
levels(as.factor(EggDiameter_CutColumn$value_type))         # "raw_value"

EggDiameter_CutColumn$value <- as.numeric(EggDiameter_CutColumn$value)

## Conversion from volume (mm3) to diameter (mm)
EggDiameter_CutColumn$value <- 2*((3*EggDiameter_CutColumn$value/(4*pi))^(1/3))
EggDiameter_CutColumn$standard_unit <- "mm"
EggDiameter_CutColumn$trait_name <- "egg diameter"

# calculation mean and SE
eggDiameterMean <- tapply(EggDiameter_CutColumn$value,EggDiameter_CutColumn$specie_name,mean)
sd <- tapply(EggDiameter_CutColumn$value,EggDiameter_CutColumn$specie_name,sd)
n <- tapply(EggDiameter_CutColumn$value,EggDiameter_CutColumn$specie_name,length)

eggDiameter_1 <- data.frame(
  species = names(eggDiameterMean),
  eggDiameter_mm = eggDiameterMean,
  sd = sd,
  n = n,
  stringsAsFactors = F
  )

#*****************************************************************************************************************************************
#### Import dataset from Figueiredo et al., 2013 , 20 species, "egg size" = "egg diameter" ########
EggDiameter_Figueriedo_orginal <- read.csv(paste(wd_Egg_diameter,"Figueiredo_et_al_2013_egg_size.csv",sep="/"), header=T,stringsAsFactors = F)

eggDiameter_2 <- data.frame(
  species = EggDiameter_Figueriedo_orginal$species,
  eggDiameter_mm = EggDiameter_Figueriedo_orginal$egg_size_micro_m /1000 ,
  sd = EggDiameter_Figueriedo_orginal$SE / 1000 * sqrt(EggDiameter_Figueriedo_orginal$N),
  n = EggDiameter_Figueriedo_orginal$N,
  stringsAsFactors = F
)

#*****************************************************************************************************************************************
#### Import "mature egg diameter",from the corltraits.org (Harriot 1982  Coral Reefs 2:9-18), 4 species
# EggDiameter_original_2 <- read.csv("https://coraltraits.org/traits/214.csv", as.is=TRUE)
# write.csv(EggDiameter_original_2,paste(wd_Datasets_original,"EggDiameter_original_2.csv",sep="/"),row.names = F)

### Import orgininal datset from folder:
EggDiameter_original_2 <- read.csv(paste(wd_Datasets_original,"EggDiameter_original_2.csv",sep="/"),header=T,stringsAsFactors = F)

nameColumnSelectedMatureEggDiameter <- c("specie_name","trait_name","standard_unit","value","value_type","precision","precision_type")

eggDiameter_cut2 <- EggDiameter_original_2[nameColumnSelectedMatureEggDiameter]

eggDiameter_cut2 <- subset(eggDiameter_cut2, trait_name == "Mature egg diameter")

eggDiameter_cut2$value <- as.numeric(eggDiameter_cut2$value)

## check if different measures for same species:
levels(as.factor(eggDiameter_cut2$standard_unit)) #  "µm"

levels(as.factor(eggDiameter_cut2$value_type)) #  "raw_value"

# conversion in mm:
eggDiameter_cut2$value <- eggDiameter_cut2$value / 1000
eggDiameter_cut2$standard_unit <- "mm"

## calculation of mean and SE:
meanEggSize3 <- tapply(eggDiameter_cut2$value,eggDiameter_cut2$specie_name,mean)
sd <- tapply(eggDiameter_cut2$value,eggDiameter_cut2$specie_name,sd)
n <- tapply(eggDiameter_cut2$value,eggDiameter_cut2$specie_name,length)

eggDiameter_3 <- data.frame(
  species = names(meanEggSize3),
  eggDiameter_mm =meanEggSize3,
  sd = sd,
  n = n,
  stringsAsFactors = F
)

#*****************************************************************************************************************************************
#### combine the 3 tables:
eggDiameter_combined <- rbind(eggDiameter_1,eggDiameter_2,eggDiameter_3)
eggDiameter_combined$species <- gsub(" ","_",eggDiameter_combined$species)
eggDiameter_combined <- eggDiameter_combined[order(eggDiameter_combined$species),]
eggDiameter_combined$species <- gsub("_"," ",eggDiameter_combined$species)

eggDiameter_combined$eggDiameter_mm <- as.numeric(eggDiameter_combined$eggDiameter_mm)
eggDiameter_combined$sd <- as.numeric(eggDiameter_combined$sd)
eggDiameter_combined$n <- as.numeric(eggDiameter_combined$n)
str(eggDiameter_combined)
length(eggDiameter_combined$species) # 26

### check for duplicates:
sum(duplicated(eggDiameter_combined$species)) # 1
duplicatedSpecies <- eggDiameter_combined[duplicated(eggDiameter_combined$species),1]
eggDiameter_combined[eggDiameter_combined$species == duplicatedSpecies,]
# Porites australiensis          0.271 0.05947941 20
# Porites australiensis          0.150         NA  1

### combine the two values:
dtCombined <- meanSDnFund(eggDiameter_combined[eggDiameter_combined$species == duplicatedSpecies,2:4])
dtCombined$species <- duplicatedSpecies
dtCombined <- dtCombined[c("species","mean","sd","n")]
colnames(dtCombined) <- c("species","eggDiameter_mm","sd","n")
# remove Porites australiensis two rows from eggDiameter_combined and add the combined one:
eggDiameter_combined <- subset(eggDiameter_combined,eggDiameter_combined$species != duplicatedSpecies)
length(eggDiameter_combined$species) # 24
eggDiameter_combined <- rbind(eggDiameter_combined,dtCombined)
length(eggDiameter_combined$species) # 25
row.names(eggDiameter_combined) <- NULL

### Check for species names:
tableNamesChecked <- eggDiameter_combined

speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets = wd_Datasets)  # "All species are present in the list"
speciesToCheck
### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp) # 1767

### correct names in species_colonialityFinal_cut:
species <- tableNamesChecked$species
numberSpecies <- length(species) # 25
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] # only keep the species in vector species
length(speciesChecked$nameSp) # 25
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}

### check if there are duplicated names:
sum(duplicated(tableNamesChecked$species)) # 0
length(tableNamesChecked$species) # 25

tableNamesChecked$species <- gsub(" ","_",tableNamesChecked$species)
tableNamesChecked <- tableNamesChecked[order(tableNamesChecked$species),]
tableNamesChecked$species <- gsub("_"," ",tableNamesChecked$species)

# write.csv(tableNamesChecked,paste(wd_Datasets,"egg_diameter.csv",sep="/"),row.names=F)

