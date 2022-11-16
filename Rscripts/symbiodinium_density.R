# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org
## produce SymbiontDensity_original.csv and Symbiodinium_density.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_Symbiont_density <- paste(wd,"/Traits_extra/symbiodinium_density",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### import orginial dataset from coraltraits.org: (April 3 2017)
# SymbiontDensity_original <- read.csv("https://coraltraits.org/traits/135.csv", as.is=TRUE)
# write.csv(SymbiontDensity_original,paste(wd_Datasets_original,"SymbiontDensity_original.csv",sep="/"),row.names = F)

### import original dataset from folder
SymbiontDensity_original <- read.csv(paste(wd_Datasets_original,"SymbiontDensity_original.csv",sep="/"),header=T,stringsAsFactors = F)
# View(SymbiontDensity_original)

length(levels(as.factor(SymbiontDensity_original$specie_name)))   # 35

SymbiontDensity_bleaching <- subset(SymbiontDensity_original,SymbiontDensity_original$trait_name=="Bleaching event")
levels(as.factor(SymbiontDensity_bleaching$value))
# "no"  "yes"

### remove unecessary columns:
colnames(SymbiontDensity_original)
nameColumnSelected <- c("specie_name","trait_name","standard_unit","value","value_type","precision","precision_type","replicates","notes","observation_id")
SymbiontDensity_cut <- SymbiontDensity_original[,nameColumnSelected]

### only keep traits "Symbiodinium density" and "Bleaching event" 
traitsToRemove <- levels(as.factor(SymbiontDensity_original$trait_name))
traitsToRemove <- traitsToRemove[-c(1,6)]
numberTraits <- length(traitsToRemove) # 7
for(i in 1: numberTraits){
  SymbiontDensity_cut <- subset(SymbiontDensity_cut,SymbiontDensity_cut$trait_name != traitsToRemove[i])
}

### check if all the observation of symbiont density are associated to bleaching event value: --> some of them are not
length(subset(SymbiontDensity_cut$specie_name,SymbiontDensity_cut$trait_name=="Bleaching event")) # 3987
length(subset(SymbiontDensity_cut$specie_name,SymbiontDensity_cut$trait_name=="Symbiodinium density")) # 4062
length(levels(as.factor(SymbiontDensity_cut$observation_id)))   #  4062          --> some observation do not have values for bleaching:
levels(as.factor(subset(SymbiontDensity_cut$value,SymbiontDensity_cut$trait_name=="Bleaching event")))

### so I am only going to select the densities not realted to a bleaching event:
# create a column stating of the row should be removed or not
SymbiontDensity_cut$Remove <- "NO"
# create a function that receive a part of the table (onl same observation_od number), check if trait "Bleaching event" is present
# if present then check if value == "yes"
# in this case, put "YES" in column Remove

table <- subset(SymbiontDensity_cut,SymbiontDensity_cut$observation_id==124033)

removeFun <- function(table=data.frame()){
  tableTemporal <- table
  lengthTable <- length(tableTemporal$specie_name)
  if(lengthTable > 1){
    for(i in 1:lengthTable){
      if(tableTemporal$trait_name[i] == "Bleaching event"){
        if(tableTemporal$value[i] == "yes"){
          for(j in 1:lengthTable){
            tableTemporal$Remove[j] <- "YES"
          }
          break
        }
      }
    }
  }
  tableTemporal
}

ObservationID <- as.numeric(levels(as.factor(SymbiontDensity_cut$observation_id)))
numberObservationID <- length(ObservationID) # 4062
for(i in 1:numberObservationID){
  tempTable <- subset(SymbiontDensity_cut,SymbiontDensity_cut$observation_id == ObservationID[i])
  tempTable <- removeFun(tempTable)
  if(tempTable$Remove[1] == "YES"){
    SymbiontDensity_cut$Remove[SymbiontDensity_cut$observation_id ==  ObservationID[i]] <- "YES"
  }
}

subset(SymbiontDensity_cut,SymbiontDensity_cut$Remove=="YES")

SymbiontDensity_cut[,SymbiontDensity_cut$value=="no" && SymbiontDensity_cut$Remove =="YES"]

### checking 
x <- 0
for(i in 1:length(SymbiontDensity_cut$specie_name)){
  if(SymbiontDensity_cut$value[i] =="no" && SymbiontDensity_cut$Remove[i] =="YES"){
    print(SymbiontDensity_cut[i,])
    x <- x +1
  }
}
print(x) # 0 , all good

### remove all the rows with Remove == "YES" ----
SymbiontDensity_cut_noB <- subset(SymbiontDensity_cut,SymbiontDensity_cut$Remove == "NO")
length(levels(as.factor(SymbiontDensity_cut_noB$specie_name)))  # 35    so we did not lose species
# View(SymbiontDensity_cut_noB)

### remove all the row with trait == "Bleaching event" -----
SymbiontDensity_cut_noB <- subset(SymbiontDensity_cut_noB,SymbiontDensity_cut_noB$trait_name == "Symbiodinium density")

levels(as.factor(SymbiontDensity_cut_noB$value_type))  # "raw_value"  --> so I can just I don't need column "value_type","
levels(as.factor(SymbiontDensity_cut_noB$standard_unit))  # "units cm^-2"  --> so I can just I don't need column "standard_unit","
levels(as.factor(SymbiontDensity_cut_noB$precision))  # character(0)  --> so I don't need precision
levels(as.factor(SymbiontDensity_cut_noB$precision_type)) # ""       --> so I don't need precision_type
levels(as.factor(SymbiontDensity_cut_noB$replicates)) # character(0)
levels(as.factor(SymbiontDensity_cut_noB$notes))
TraitsToKeep <- c("specie_name","value","notes","observation_id")

SymbiontDensity_cut_noB <- SymbiontDensity_cut_noB[,TraitsToKeep]

SymbiontDensity_cut_noB$value <- as.numeric(SymbiontDensity_cut_noB$value)

meanSD <- tapply(SymbiontDensity_cut_noB$value,SymbiontDensity_cut_noB$specie_name,mean)
sdSD <- tapply(SymbiontDensity_cut_noB$value,SymbiontDensity_cut_noB$specie_name,sd)
n <- tapply(SymbiontDensity_cut_noB$value,SymbiontDensity_cut_noB$specie_name,function(x){length(x)})
SE.SD <- sdSD/sqrt(n)

Symbiodimium.Density <- cbind(round(meanSD,2),round(sdSD,2),n)
colnames(Symbiodimium.Density) <- c("mean_Symb_density (units/cm2)","sd","n")
Symbiodimium.Density <- as.data.frame((Symbiodimium.Density))
Symbiodimium.Density$species <- rownames(Symbiodimium.Density)

Symbiodimium.Density <- Symbiodimium.Density[,c(4,1,2,3)]
head(Symbiodimium.Density)

tableNamesChecked <- Symbiodimium.Density

### check for names:
speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets = wd_Datasets)
speciesToCheck

### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp)

### correct names in species_colonialityFinal_cut:
species <- tableNamesChecked$species
numberSpecies <- length(species) # 35
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] # only keep the species in vector species
length(speciesChecked$nameSp) # 35
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}

### check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 0
# write.csv(tableNamesChecked,paste(wd_Datasets,"Symbiodinium_density.csv",sep="/"),row.names = F)

