# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org
## produce life_history_strategy_original.csv and life_history_strategy.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### import original dataset from caraltraits.org: (August 16 2018)
# lhs_original <- read.csv("https://coraltraits.org/traits/233.csv", as.is=TRUE)
# write.csv(lhs_original,paste(wd_Datasets_original,"/life_history_strategy_original.csv",sep=""),row.names = F)

### import original dataset from folder
lhs_original <- read.csv(paste(wd_Datasets_original,"life_history_strategy_original.csv",sep="/"),header=T,stringsAsFactors = F)
head(lhs_original)
names(lhs_original)[5] <- "species_name"

nameColumnSelected <- c("species_name","trait_name","value")

lhs_cut <- lhs_original[nameColumnSelected]
head(lhs_cut)

# check the different traits:
levels(as.factor(lhs_cut$trait_name))

# number of lines :
length(lhs_cut$species_name)    # 143

# different values
levels(as.factor(lhs_cut$value)) # "competitive"     "generalist"      "stress-tolerant" "weedy"          

# remvoe potential white space at the beginning and end of species names:
lhs_cut$species_name <- trimws(lhs_cut$species_name)

### check if there are duplicate species. If so, check if same value. If so, remove the duplicate. If not, show which species
length(levels(as.factor(lhs_cut$species_name))) # 143 --> no duplicate 
sum(duplicated(lhs_cut$species_name)) # 0

### check for names:
tableNamesChecked <- lhs_cut
speciesToCheck <- checkSpeciesNames(tableNamesChecked$species_name, wd_Datasets = wd_Datasets) # "All species are present in the list

### remove the duplicated species and update the species names for species whose updated name is not already present
speciestoRemove <- c()
x <- 0
for(i in 1:length(speciesToCheck$nameSp)){
  if(speciesPresentFun(speciesToCheck$nameSp_checked[[i]],lhs_cut$species_name)){ # TRUE: the correct name is already present in lhs_cut$species_name --> remove the name
    x <- x + 1
    speciestoRemove[x] <- speciesToCheck$nameSp[i]
  }else{ # update name
    lhs_cut$species_name[lhs_cut$species_name == speciesToCheck$nameSp[i]] <- speciesToCheck$nameSp_checked[i]
    print(paste("Species",speciesToCheck$nameSp[i]," changed to ",speciesToCheck$nameSp_checked[i]))
  }
}
length(speciestoRemove) # 0

### remove the speciestoRemove
lhs_cut <- lhs_cut[!lhs_cut$species_name %in% speciestoRemove,]
length(lhs_cut$species_name) # 143
sum(duplicated(lhs_cut$species_name)) # 0

names(lhs_cut)[1] <- "species"
names(lhs_cut)[3] <- "life_history_strategy"

### alphabetic order:
lhs_cut$species <- gsub(" ","_",lhs_cut$species)
lhs_cut <- lhs_cut[order(lhs_cut$species),]
lhs_cut$species <- gsub("_"," ",lhs_cut$species)

# creation csv file in the wd dataCompilation
# write.csv(lhs_cut[,c("species","life_history_strategy")],paste(wd_Datasets,"life_history_strategy.csv",sep="/"),row.names=F)


