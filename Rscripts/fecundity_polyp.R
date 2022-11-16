# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org and from Hall and Hughes 1996, 
## produce FecundityPolyp_original.csv, FecundityMesentery_original.csv and fecundity_polyp.csv

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_Fecundity_polyp <- paste(wd,"/Traits_extra/fecundity_polyp",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### Import dataset from coraltraits.org: (April 5 2017)
# FecundityPolyp_original <- read.csv("https://coraltraits.org/traits/12.csv", as.is=TRUE)
# write.csv(FecundityPolyp_original,paste(wd_Datasets_original,"FecundityPolyp_original.csv",sep="/"),row.names = F)

### Import dataset from folder:
FecundityPolyp_original <- read.csv(paste(wd_Datasets_original,"FecundityPolyp_original.csv",sep="/"),header=T,stringsAsFactors = F)

nameColumnSelectedFecundity <- c("specie_name","trait_name","standard_unit","value","value_type","precision","precision_type","replicates")

dataFecundityCutColumn <- FecundityPolyp_original[nameColumnSelectedFecundity]

# check if different types of traits:
levels(as.factor(dataFecundityCutColumn$trait_name))
dataFecundity <- subset(dataFecundityCutColumn,  trait_name == "Polyp fecundity")

levels(as.factor(dataFecundity$precision_type)) 
#  ""     "standard_deviation" "standard_error" 
levels(as.factor(dataFecundity$replicates))
# character(0)
levels(as.factor(dataFecundity$value_type))
# "mean"      "raw_value"

# check number of different species:
length(levels(as.factor(dataFecundity$specie_name)))      # 10 
length(as.factor(dataFecundity$specie_name))              # 44

dataFecundity$value <- as.numeric(dataFecundity$value)

# for Acropora gemmifera,Acropora hyacinthus,Acropora millepora,Acropora nana,Goniastrea retiformis,Stylophora pistillata --> values from Hall and Hughes 1996, figure 2a: n = 9*40 (in text, n = 15*40 though...)
# For Stylophora pistillata, there are many values. Most of them are single observations coming form Hall and Hughes 1996, figure 3c and corresponds to number of planulae (because brooding species) --> to remove
# The few other values come from Liberman et al., 1995 but I don't consider them.
species_Hall <- c("Acropora gemmifera","Acropora hyacinthus","Acropora millepora","Acropora nana","Goniastrea retiformis","Stylophora pistillata")
dataFecundity$replicates[dataFecundity$specie_name %in% species_Hall]
condition1 <- dataFecundity$specie_name == "Stylophora pistillata" & dataFecundity$value_type == "raw_value"
condition2 <- dataFecundity$specie_name != "Stylophora pistillata"
dataFecundity_cut <- subset(dataFecundity, condition1 | condition2)
# Add sample size:
dataFecundity_cut$replicates[dataFecundity_cut$specie_name %in% species_Hall] <- 9*40
### the other species values come from Harriot 1983 and do not have sample size are measure of variation

fecundityPolypFinal <- data.frame(
  species = dataFecundity_cut$specie_name,
  Fecundity_polyp = dataFecundity_cut$value,
  sd = dataFecundity_cut$precision * sqrt(dataFecundity_cut$replicates),
  n = dataFecundity_cut$replicates,
  stringsAsFactors = F
)

fecundityPolypFinal$n[is.na(fecundityPolypFinal$n)] <- 1

#*************************************************************************************************
# Include MESENTERY FECUNDITY (ID 209) --> contains (1)  Eggs per area and (2) Polyps per area  --> (1) / (2) --> polyp fecundity
# But by doing so I cannot obtain a measure of sd because the ratio of two standard normal random variables ($\mu = 0, \sigma = 1) is a Cauchy distribution; 
# the Cauchy has an undefined variance (and hence undefined standard deviation): 
# http://stats.stackexchange.com/questions/58800/what-is-the-mean-and-standard-deviation-of-the-division-of-two-random-variables

### Import dataset from coraltraits.org:
# FecundityMesentery_original <- read.csv("https://coraltraits.org/traits/209.csv", as.is=TRUE)
# write.csv(FecundityMesentery_original,paste(wd_Datasets_original,"FecundityMesentery_original.csv",sep="/")row.names = T)

### Import dataset from folder:
FecundityMesentery_original <- read.csv(paste(wd_Datasets_original,"FecundityMesentery_original.csv",sep="/"),header = T,stringsAsFactors = F)

nameColumnSelectedFecundityMesentery <- c("specie_name","trait_name","standard_unit","value","value_type","precision","precision_type","replicates","notes")
FecundityMesenteryCutColumn <- FecundityMesentery_original[nameColumnSelectedFecundityMesentery]

# check if different traits:
levels(as.factor(FecundityMesenteryCutColumn$trait_name))
traitToKeep <- c("Eggs per area","Polyps per area")
FecundityMesenteryCutColumn_cut <- subset(FecundityMesenteryCutColumn,trait_name %in% traitToKeep)
FecundityMesenteryCutColumn_cut$value <- as.numeric(FecundityMesenteryCutColumn_cut$value)

species <- levels(as.factor(FecundityMesenteryCutColumn_cut$specie_name))
numberSpecies <- length(species) # 3
fecundityPolyp_mesentery <- data.frame(
  species = rep(NA,3),
  Fecundity_polyp = rep(NA,3),
  sd = rep(NA,3),
  n = rep(NA,3),
  stringsAsFactors = F
)
for(i in 1:numberSpecies){
  DTTempo <- subset(FecundityMesenteryCutColumn_cut,FecundityMesenteryCutColumn_cut$specie_name == species[i])
  fecundityPolyp_mesentery$species[i] <-  species[i]
  numbEggArea <- sum(DTTempo$value[DTTempo$trait_name == "Eggs per area"] * DTTempo$replicates[DTTempo$trait_name == "Eggs per area"]) / sum(DTTempo$replicates[DTTempo$trait_name == "Eggs per area"])
  numbPolypArea <- sum(DTTempo$value[DTTempo$trait_name == "Polyps per area"] * DTTempo$replicates[DTTempo$trait_name == "Polyps per area"]) / sum(DTTempo$replicates[DTTempo$trait_name == "Polyps per area"])
  fecundityPolyp_mesentery$Fecundity_polyp[i] <- numbEggArea / numbPolypArea
  fecundityPolyp_mesentery$n[i] <- sum(DTTempo$replicates)/2
}

#### combination of the 2 tables:
fecundityPolypFinalFinal <- rbind(fecundityPolypFinal,fecundityPolyp_mesentery)

fecundityPolypFinalFinal$species <- gsub(" ","zzz",fecundityPolypFinalFinal$species)
fecundityPolypFinalFinal <- fecundityPolypFinalFinal[order(fecundityPolypFinalFinal$species),]
fecundityPolypFinalFinal$species <- gsub("zzz"," ",fecundityPolypFinalFinal$species)
fecundityPolypFinalFinal
length(fecundityPolypFinalFinal$species) # 13

tableNamesChecked <- fecundityPolypFinalFinal

### check for names:
speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets = wd_Datasets) #  "All species are present in the list"
speciesToCheck

### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp)

### correct names in species_colonialityFinal_cut:
species <- tableNamesChecked$species
numberSpecies <- length(species) # 13
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] # only keep the species in vector species
length(speciesChecked$nameSp) # 13
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}

# Check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 0

tableNamesChecked$species <- gsub(" ","_",tableNamesChecked$species)
tableNamesChecked <- tableNamesChecked[order(tableNamesChecked$species),]
tableNamesChecked$species <- gsub("_"," ",tableNamesChecked$species)

# write.csv(tableNamesChecked,paste(wd_Datasets,"fecundity_polyp.csv",sep="/"),row.names = F)



