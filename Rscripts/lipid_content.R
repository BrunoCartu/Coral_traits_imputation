# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org
## produce LipidContent_original.csv and lipid_content.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_Lipid_content <- paste(wd,"/Traits_extra/lipid_content",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

# Def: The amount of lipid in the coral holobiont tissue per unit of skeletal surface area ("mg cm^-2")

### Import original dataset from coraltraits.org:
# LipidContent_original <- read.csv("https://coraltraits.org/traits/131.csv", as.is=TRUE)
# write.csv(AgeAtMaturity_original,paste(wd_Datasets_original,"LipidContent_original.csv",sep="/"),row.names = F)

### Import original dataset from folder:
LipidContent_original <- read.csv(paste(wd_Datasets_original,"LipidContent_original.csv",sep="/"),header = T,stringsAsFactors = F)

nameColumnSelected <- c("specie_name","trait_name","standard_unit","value","value_type","precision","precision_type","replicates","notes")

LipidContent_cut <- LipidContent_original[nameColumnSelected]

### Traits present:
levels(as.factor(LipidContent_cut$trait_name)) # "Habitat type"      "Irradiance"        "Lipid content"     "Season"            "Water depth"       "Water temperature"
LipidContent_cut <- subset(LipidContent_cut,LipidContent_cut$trait_name %in% c("Lipid content"))

### value_type :
levels(as.factor(LipidContent_cut$value_type)) # "mean"

### standard_unit
levels(as.factor(LipidContent_cut$standard_unit)) # "mg cm^-2"

###
levels(as.factor(LipidContent_cut$replicates)) # ""     "10"   "12"   "3"    "4-10" "4-12" "8" 

LipidContent_cut$value <- as.numeric(LipidContent_cut$value)
LipidContent_cut$precision <- as.numeric(LipidContent_cut$precision)

species <- as.character(levels(as.factor(LipidContent_cut$specie_name)))
numberSpecies <- length(species)

LipidContent_final <- data.frame(
  species = species,
  lipid_content_mg.cm2 = rep(NA,numberSpecies),
  sd =  rep(NA,numberSpecies),
  n =  rep(NA,numberSpecies),
  stringsAsFactors = F
)

subset(LipidContent_original,LipidContent_original$specie_name == "Montipora capitata")

#**********************************************************************************************************************************
# Need to check in the references to get exact sample size, and chose value that correspond to pre-spawning period (when appropriate)
#**********************************************************************************************************************************

# For species in Leuzinger et al., 2003, values of mean and sd were obtained from figure 1 and for n from table 1 (see Leuzinger et al., 2003 conversion fig 1.xlsx)
### Update LipidContent_cut: 
LipidContent_cut$reproduction <- NA

i <- 1
LipidContent_cut[LipidContent_cut$specie_name == species[i] ,] # in Middlebrook et al., 2010, effect of fast vs slow heating rate on lipid content --> no reproduction involve so I combine the 2 values
LipidContent_final$species[i] <- species[i]
DF <- LipidContent_cut[LipidContent_cut$specie_name == species[i] ,][c(4,6,8)]
DF$replicates <- as.numeric(DF$replicates)
# converision SE to sd 
DF$precision <- DF$precision * sqrt(DF$replicates)
LipidContent_final[i,2:4] <- meanSDnFund(DF)

i <- 2
LipidContent_cut[LipidContent_cut$specie_name == species[i] ,] # in Leuzinger et al., 2003 --> value before and after spwening --> onlu consider before spawning, n = 6
LipidContent_final$species[i] <- species[i]
LipidContent_final[i,2:4] <- c(4.6,0.711,6)
LipidContent_cut$replicates[LipidContent_cut$specie_name == species[i]] <- 6

i <- 3
LipidContent_cut[LipidContent_cut$specie_name == species[i] ,] # in Leuzinger et al., 2003 --> value before and after spwening --> onlu consider before spawning, n = 6
LipidContent_final$species[i] <- species[i]
LipidContent_final[i,2:4] <- c(11.0,2.370,6)
LipidContent_cut$replicates[LipidContent_cut$specie_name == species[i]] <- 6

i <- 4
LipidContent_cut[LipidContent_cut$specie_name == species[i] ,] # in Rodrigues and Grottoli 2007. I choose n = 5, in accordnace with text but vsampling disgn very not clear,
LipidContent_final$species[i] <- species[i]
LipidContent_final[i,2:4] <- c(LipidContent_cut[LipidContent_cut$specie_name == species[i] ,]["value"],NA,5)

i <- 5
LipidContent_cut[LipidContent_cut$specie_name == species[i] ,] # in Leuzinger et al., 2003 --> value before and after spwening --> onlu consider before spawning, n = 10
LipidContent_final$species[i] <- species[i]
LipidContent_final[i,2:4] <- c(6.8,1.275,10)

i <- 6
LipidContent_cut[LipidContent_cut$specie_name == species[i] ,]  # combine tghe 2 values:
DF <- LipidContent_cut[LipidContent_cut$specie_name == species[i] ,][c(4,6,8)]
# conversion of se to sd for 2nd row:
DF$replicates <- as.numeric(DF$replicates)
DF$precision[2] <- DF$precision[2] * sqrt(DF$replicates[2])
# converision SE to sd 
LipidContent_final[i,2:4] <- meanSDnFund(DF)

i <- 7
LipidContent_cut[LipidContent_cut$specie_name == species[i] ,]   # not sure if that's the value corresponding to the control because it is not the same unit... I just trust the coraltraits.org 
LipidContent_final$species[i] <- species[i]
LipidContent_final[i,2:4] <- subset(LipidContent_cut,LipidContent_cut$specie_name == species[i])[,c("value","precision","replicates")]
# conversion se to sd:
LipidContent_final$n <- as.numeric(LipidContent_final$n)
LipidContent_final[i,3] <- LipidContent_final[i,3] * sqrt(LipidContent_final[i,4])

i <- 8
LipidContent_cut[LipidContent_cut$specie_name == species[i] ,] # in Rodrigues and Grottoli 2007. I choose n = 5, in accordance with text but sampling disign very not clear, SE estimated from figure
LipidContent_final$species[i] <- species[i]
LipidContent_final[i,2:4] <- c(LipidContent_cut[LipidContent_cut$specie_name == species[i] ,]["value"],NA,5)

### check for names:
tableNamesChecked <- LipidContent_final

speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets = wd_Datasets) # "All species are present in the list"
speciesToCheck

### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp)

### correct names in species_colonialityFinal_cut:
species <- tableNamesChecked$species
numberSpecies <- length(species) # 8
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] # only keep the species in vector species
length(speciesChecked$nameSp) # 8
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}
# Check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 0

### alphabetic order:
tableNamesChecked$species <- gsub(" ","_",tableNamesChecked$species)
tableNamesChecked <- tableNamesChecked[order(tableNamesChecked$species),]
tableNamesChecked$species <- gsub("_"," ",tableNamesChecked$species)

# creation csv file in the wd dataCompilation
# write.csv(tableNamesChecked, file = paste(wd_Datasets,"lipid_content.csv",sep="/"),row.names=F)
