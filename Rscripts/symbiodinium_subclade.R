# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# NOT USED
# Goals:
## import and curate species trait data from coraltraits.org
## produce Symbiodinium_subclass_diversity_original.csv, symbiodinium_subclade.csv and symbiodinium_subclade_noType.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_symbiodinium_subclade <- paste(wd,"Traits_extra/symbiodinium_subclade",sep="/")

source(paste(wd_Rscripts,"functions.R",sep="/"))

# http://www.symbiogbr.org/nomenclature

### ### Import original dataset from coraltraits.org: (March 30 2017)
# SymbiontDiversity_original <- read.csv("https://coraltraits.org/traits/129.csv", as.is=TRUE)
# The genetic identity of Symbiodinium found in coral tissue at the level below clade, but usually above species. 
# This is typically identified using the nuclear ribosomal DNA Internal Transcribed Spacer region (ITS2), but other markers are also used.
# write.csv(SymbiontDiversity_original,paste(wd_Datasets_original,"Symbiodinium_subclass_original.csv",sep="/"),row.names = F)

SymbiontDiversity_original <- read.csv(paste(wd_Datasets_original,"Symbiodinium_subclade_original.csv",sep="/"),header = T, stringsAsFactors = F)
#View(SymbiontDiversity_original)

### remove unecessary columns:
colnames(SymbiontDiversity_original)
nameColumnSelected <- c("specie_name","trait_name","standard_unit","value","value_type","precision","precision_type","replicates","notes","observation_id")
SymbiontDiversity_cut <- SymbiontDiversity_original[,nameColumnSelected]

### traits selection: 
levels(as.factor(SymbiontDiversity_cut$trait_name)) # "Symbiodinium subclade" "Water depth" 
SymbiontDiversity_cut <- subset(SymbiontDiversity_cut,SymbiontDiversity_cut$trait_name=="Symbiodinium subclade")
SymbiontDiversity_cut <- SymbiontDiversity_cut[c("specie_name","value","notes","observation_id")]  # also removed standard_unit, value_type, precision, precision_type, replicates

### number coral species:
length(levels(as.factor(SymbiontDiversity_cut$specie_name))) # 311

### number different symbiodinium clades: 
length(levels(as.factor(SymbiontDiversity_cut$value))) # 248
# View(SymbiontDiversity_cut)

### check for names on WoRMS
# speciesToSendWoRMS <- data.frame(do.call("rbind",strsplit(levels(as.factor(SymbiontDiversity_cut$specie_name))," ",fixed = T)))
# colnames(speciesToSendWoRMS) <- c("Genus","Species")
# head(speciesToSendWoRMS)
# write.csv(speciesToSendWoRMS,paste(wd_symbiodinium_subclade,"speciesToSendWoRMS.csv",sep="/"),row.names=F)
speciesToSendWoRMS <- read.csv(paste(wd_symbiodinium_subclade,"speciesToSendWoRMS.csv",sep="/"),header = T, stringsAsFactors = F)

speciesCoral <- levels(as.factor((SymbiontDiversity_cut$specie_name)))
  
speciesToCheck <- checkSpeciesNames(speciesCoral, wd_Datasets = wd_Datasets) #  "All species are present in the list"

for(i in 1:length(speciesToCheck$nameSp)){
  SymbiontDiversity_cut$specie_name[SymbiontDiversity_cut$specie_name == speciesToCheck$nameSp[i]] <- speciesToCheck$nameSp_checked[i]
}

length(levels(as.factor(SymbiontDiversity_cut$specie_name))) # 309

SymbiontDiversity_cut_correctedNames <- SymbiontDiversity_cut

### creation table with rownames <- coral species and colnames = symbiont clade names:
coralSpecies <- levels(as.factor(SymbiontDiversity_cut_correctedNames$specie_name))
symbiontClade <- levels(as.factor(SymbiontDiversity_cut_correctedNames$value))
numberCoralSpecies <- length(coralSpecies)        # 309
numberSymbiotClade <- length(symbiontClade)       # 248

symbiodinium_subclade_final <- matrix(0,ncol = numberSymbiotClade,nrow=numberCoralSpecies)
colnames(symbiodinium_subclade_final) <- symbiontClade
rownames(symbiodinium_subclade_final) <- coralSpecies
head(symbiodinium_subclade_final)

for(i in 1:numberCoralSpecies){
  temporalTable <- subset(SymbiontDiversity_cut_correctedNames,SymbiontDiversity_cut_correctedNames$specie_name==coralSpecies[i])
  temporalSymbiontClade <- levels(as.factor(temporalTable$value))
  numberTemporalSymbiontClade <- length(temporalSymbiontClade)
  for(j in 1:numberTemporalSymbiontClade){
    for(k in 1:numberSymbiotClade){
      if(temporalSymbiontClade[j] == symbiontClade[k]){
        symbiodinium_subclade_final[i,k] <- 1
        break
      }
    }
  }
}

symbiodinium_subclade_final <- as.data.frame(symbiodinium_subclade_final)
symbiodinium_subclade_final$Number_Clades <- NA
for(i in 1:length(symbiodinium_subclade_final$A1)){
  symbiodinium_subclade_final$Number_Clades[i] <- sum(symbiodinium_subclade_final[i,1:(length(symbiodinium_subclade_final)-1)])
}
# View(symbiodinium_subclade_final)
min(symbiodinium_subclade_final$Number_Clades)

species <- data.frame(species = row.names(symbiodinium_subclade_final))

symbiodinium_subclade_final <- cbind(species,symbiodinium_subclade_final)
# colnames(symbiodinium_subclade_final) <- gsub(" ","_",colnames(symbiodinium_subclade_final))
colnames(symbiodinium_subclade_final)

# write.csv(symbiodinium_subclade_final,paste(wd_Datasets,"symbiodinium_subclade.csv",sep="/"),row.names = F)
# symbiodinium_subclade_final <- read.csv(paste(wd_Datasets,"symbiodinium_subclade.csv",sep="/"),header = T, stringsAsFactors = F)
colnames(symbiodinium_subclade_final)

### some clades had different "type", so I removed the different types and and only consider 1 type of clade, so if a coral species was observed with "C3 (type 1)" and "C3 (type 2)", there is only 1 for "C3" ----
# "C35"           "C35 (type 3)"  
# "C35a"          "C35a (type 1)"   "C35a (type 2)"
# "C3k"           "C3k (type 1)" 
# "C30 (type 2)"
# "C42"           "C42 (type 1)"  "C42 (type 2)" 
# "F2 (type 1)"
# "C33"           "C33 (type 1)"  "C33 (type 2)"
# "F2 (type 1)"

# symbiodinium_subclade_final_noType <- symbiodinium_subclade_final
# symbiodinium_subclade_final_noType$C35[symbiodinium_subclade_final_noType$C35_.type_3. == 1 ] <- 1
# symbiodinium_subclade_final_noType$C35a[symbiodinium_subclade_final_noType$C35a_.type_1. == 1 | symbiodinium_subclade_final_noType$C35a_.type_2. == 1] <- 1
# symbiodinium_subclade_final_noType$C3k[symbiodinium_subclade_final_noType$C3k_.type_1. == 1 ] <- 1
# colnames(symbiodinium_subclade_final_noType)[colnames(symbiodinium_subclade_final_noType) == "C30_.type_2."] <- "C30"
# symbiodinium_subclade_final_noType$C42[symbiodinium_subclade_final_noType$C42_.type_1. == 1 | symbiodinium_subclade_final_noType$C42_.type_2. == 1] <- 1
# colnames(symbiodinium_subclade_final_noType)[colnames(symbiodinium_subclade_final_noType) == "F2 (type 2)"] <- "F2"
# symbiodinium_subclade_final_noType$C33[symbiodinium_subclade_final_noType$C33_.type_1. == 1 | symbiodinium_subclade_final_noType$C33_.type_2. == 1] <- 1
# colnames(symbiodinium_subclade_final_noType)[colnames(symbiodinium_subclade_final_noType) == "F2_.type_1."] <- "F2"
# colnames(symbiodinium_subclade_final_noType)

symbiodinium_subclade_final_noType <- symbiodinium_subclade_final
symbiodinium_subclade_final_noType$C35[symbiodinium_subclade_final_noType$`C35 (type 3)` == 1 ] <- 1
symbiodinium_subclade_final_noType$C35a[symbiodinium_subclade_final_noType$`C35a (type 1)` == 1 | symbiodinium_subclade_final_noType$`C35a (type 2)` == 1] <- 1
symbiodinium_subclade_final_noType$C3k[symbiodinium_subclade_final_noType$`C3k (type 1)` == 1 ] <- 1
colnames(symbiodinium_subclade_final_noType)[colnames(symbiodinium_subclade_final_noType) == "C30 (type 2)"] <- "C30"
symbiodinium_subclade_final_noType$C42[symbiodinium_subclade_final_noType$`C42 (type 1)` == 1 | symbiodinium_subclade_final_noType$`C42 (type 2)` == 1] <- 1
colnames(symbiodinium_subclade_final_noType)[colnames(symbiodinium_subclade_final_noType) == "F2 (type 2)"] <- "F2"
symbiodinium_subclade_final_noType$C33[symbiodinium_subclade_final_noType$`C33 (type 1)` == 1 | symbiodinium_subclade_final_noType$`C33 (type 2)`== 1] <- 1
colnames(symbiodinium_subclade_final_noType)[colnames(symbiodinium_subclade_final_noType) == "F2 (type 1)"] <- "F2"

### remove the columns with (type i)
columnsToRemove <- c("C35 (type 3)","C35a (type 1)","C35a (type 2)","C3k (type 1)","C42 (type 1)","C42 (type 2)","C33 (type 1)","C33 (type 2)","F2 (type 1)")
# columnsToRemove <- gsub(" ","_",columnsToRemove)
# columnsToRemove <- gsub("\\(",".",columnsToRemove)
# columnsToRemove <- gsub("\\)",".",columnsToRemove)

columnsToKeep <- colnames(symbiodinium_subclade_final_noType)[!colnames(symbiodinium_subclade_final_noType) %in% columnsToRemove]
symbiodinium_subclade_final_noType <- symbiodinium_subclade_final_noType[,columnsToKeep]
colnames(symbiodinium_subclade_final_noType)

symbiodinium_subclade_final_noType$Number_Clades <- NA

numberCoralSpecies <- length(symbiodinium_subclade_final_noType$species)
for(i in 1:numberCoralSpecies){
  symbiodinium_subclade_final_noType$Number_Clades[i] <- sum(symbiodinium_subclade_final_noType[i,-c(1,length(symbiodinium_subclade_final_noType))])
}

colnames(symbiodinium_subclade_final_noType)
# write.csv(symbiodinium_subclade_final_noType,paste(wd_Datasets,"symbiodinium_subclade_noType.csv",sep="/"),row.names = F)

