# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# NOT USED: because not enough species.
# Goals:
## import and curate species trait data from coraltraits.org and from (Alvarez-Noriega et al., 2016)
## produce fecundity_colony_original.csv

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### Impost dataset from coraltraits.org (April 18 2017)
# fecundity_colony_original <- read.csv("https://coraltraits.org/traits/216.csv", as.is = TRUE)
# View(colonyFecundity_original)
# write.csv(fecundity_colony_original,paste(wd_Datasets_original,"fecundity_colony_original.csv",sep="/"),row.names = F)

### Import original dataset from folder:
fecundity_colony_original <- read.csv(paste(wd_Datasets_original,"fecundity_colony_original.csv",sep="/"),header = T, stringsAsFactors = F)

levels(as.factor(fecundity_colony_original$trait_name))
# "Colony area"      "Colony fecundity" "Habitat type"     "Season"           "Water depth" 

length(levels(as.factor(fecundity_colony_original$specie_name))) # 6
levels(as.factor(fecundity_colony_original$specie_name))
# "Acropora gemmifera"    "Acropora hyacinthus"   "Acropora millepora"    "Acropora nana"         "Goniastrea retiformis" "Stylophora pistillata"

species_coralTrait <- levels(as.factor(fecundity_colony_original$specie_name))

### Species from Alvarez-Noriega et al., 2016:

# tabular (Acropora cytherea and A. hyacinthus), 
# corymbose (A. nasuta and A. spathulata), 
# digitate (A. humilis and A. cf. digitifera) 
# massive (Goniastrea pectinata and G. retiformis;)

species_Alvarez <- c("Acropora cytherea","Acropora hyacinthus","Acropora nasuta","Acropora spathulata","Acropora humilis","Acropora digitifera","Goniastrea pectinata","Goniastrea retiformis")

totalSpecies <- c(species_coralTrait,species_Alvarez)
totalSpecies <- levels(as.factor(totalSpecies))
length(totalSpecies) # 12
totalSpecies
