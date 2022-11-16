# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org
## produce SizeAtMaturity_original.csv and size_at_maturity.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_Size_at_maturity <- paste(wd,"/Traits_extra/size_at_maturity",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### final output: size at maturity = planar surface area in cm2

### Import original dataset from coraltraits.org: (April8 2017)
# SizeAtMaturity_original <- read.csv("https://coraltraits.org/traits/46.csv", as.is=TRUE)
# write.csv(SizeAtMaturity_original,paste(wd_Datasets_original,"SizeAtMaturity_original.csv",sep="/"),row.names = F)

### Import original dataset from folder:
SizeAtMaturity_original <- read.csv(paste(wd_Datasets_original,"SizeAtMaturity_original.csv",sep="/"),header = T, stringsAsFactors = F)

nameColumnSelectedSAM <- c("specie_name","trait_name","standard_unit","value","value_type","precision","precision_type","replicates","notes")

SizeAtMaturity_cut <- SizeAtMaturity_original[nameColumnSelectedSAM]

## check if different traits:
levels(as.factor((SizeAtMaturity_cut$trait_name)))
# "Corallite width" "Eggs per area" "Habitat type" "Mature egg diameter" "Polyp fecundity" Polyps per area" "Ratio of female to male gonads" "Season" "Size at maturity" "Water depth"

## only select "Size at maturity":
SizeAtMaturity_cut <- subset(SizeAtMaturity_cut,SizeAtMaturity_cut$trait_name =="Size at maturity")

SizeAtMaturity_cut$value <- as.numeric(SizeAtMaturity_cut$value)

## check if several values for a same species:
length(levels(as.factor(SizeAtMaturity_cut$specie_name)))
# 7
length(SizeAtMaturity_cut$specie_name) # 7

### check precision_type
levels(as.factor(SizeAtMaturity_cut$precision_type)) # ""

### check replicates
levels(as.factor(SizeAtMaturity_cut$replicates)) # character()

## checkif different units of measurement:
levels(as.factor(SizeAtMaturity_cut$standard_unit))
# "cm"   "cm^2"
# cm is in Harriott 1983, measured the colony diameter
# cm^2 in Hall and Hughes 1996, S = pi * r1/2 * r2/2   --> surface of an ellipse

# conversio in cm^2 considering circular planar surface area:
SizeAtMaturity_cut$value[SizeAtMaturity_cut$standard_unit == "cm"] <- pi * (SizeAtMaturity_cut$value[SizeAtMaturity_cut$standard_unit == "cm"]/2)^2
SizeAtMaturity_cut$standard_unit[SizeAtMaturity_cut$standard_unit == "cm"] <- "cm^2"

## construction dataset:
sizeAtMaturityFinal <- data.frame(
  species = SizeAtMaturity_cut$specie_name,
  sizeAtMaturity_cm2 = SizeAtMaturity_cut$value,
  stringsAsFactors = F
)

### check for names:
tableNamesChecked <- sizeAtMaturityFinal

speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets = wd_Datasets) # "All species are present in the list"
speciesToCheck

### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp)

### correct names in species_colonialityFinal_cut:
species <- tableNamesChecked$species
numberSpecies <- length(species) # 7
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

# write.csv(tableNamesChecked,paste(wd_Datasets,"size_at_maturity.csv",sep = "/"),row.names = F)

