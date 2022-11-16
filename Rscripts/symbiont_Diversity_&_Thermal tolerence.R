# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# NOT USED
# Goals:
## combine thew follwing datasets:
## produce Symbiodinium_subclass_diversity_original.csv, symbiodinium_subclade.csv and symbiodinium_subclade_noType.csv 
## with
## symbiont thermotolerence (from Swain et al. 2016 - Consensus thermotolerance ranking for 110 Symbiodinium phylotypes...)
## with 
## symbiont_List1_List2_Eugenia_Sampayo.csv: list corrected by Dr. Eugenia Sampayo (personal communication, March 20 2017)
## to create:
### symbiont_List1_notList2.csv, symbiont_notList1_List2.csv, and symbiodinium_subcalssDiversity_thermotol.csv

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")

# wd_symbiodinium_subclade <- paste(wd,"Traits_extra/symbiodinium_subclade",sep="/")
# 
# wd <- "~/Dropbox/PhD Bruno/Functional Traits/Traits analysis/Traits"
# wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
# wd_Datasets <- paste(wd,"Datasets",sep="/")
# wd_zooxanthellae <- paste(wd,"Traits_extra/zooxanthellae",sep="/")
# wd_Rscripts <- paste(wd,"Rscripts",sep="/")
wd_symbiont_diversity_thermotolerence <- paste(wd,"Traits_extra/symbiont_diversity_thermotolerence",sep="/")
  
source(paste(wd_Rscripts,"functions.R",sep="/"))

require(plyr)
require(rowr)

# http://www.symbiogbr.org/nomenclature

### get symbiont subclade
symbiodinium_subclade <- read.csv(paste(wd_Datasets,"symbiodinium_subclade_noType.csv",sep="/"),header=T,stringsAsFactors = F,check.names = F)
nrow(symbiodinium_subclade) # 309
ncol(symbiodinium_subclade) # 242
colnames(symbiodinium_subclade)

# pathBis <- "/Users/carturan/Dropbox/PhD Bruno/Functional Traits/Traits analysis/Traits/Traits_extra/symbiont_diversity"
# symbiodinium_subcladeBIS <- read.csv(paste(pathBis,"symbiont_diversity_final_noType.csv",sep="/"),header=T,stringsAsFactors = F,check.names = F)
# nrow(symbiodinium_subcladeBIS) # 309
# ncol(symbiodinium_subcladeBIS) # 243
# colnames(symbiodinium_subcladeBIS)

# colnames(symbiodinium_subcladeBIS)[! colnames(symbiodinium_subcladeBIS) %in% colnames(symbiodinium_subclade)] # "C3-o" "C3-p"
# colnames(symbiodinium_subclade)[! colnames(symbiodinium_subclade) %in% colnames(symbiodinium_subcladeBIS)] # "C3o" "C3p"
# colnames(symbiodinium_subcladeBIS)[(duplicated(colnames(symbiodinium_subcladeBIS)))]

### get symbiont thermotolerence (from Swain et al. 2016 - Consensus thermotolerance ranking for 110 Symbiodinium phylotypes...)
symbiont_thermoTol <- read.csv(paste(wd,"Traits_extra/symbiont_thermotolerance","Swain et al. - 2016 - Symbionte thermotolerance.csv",sep="/"),header=T,stringsAsFactors = F,check.names = F)

symbiontClade_diversity <- colnames(symbiodinium_subclade)[-c(1,243)]
number_symbiontClade_diversity <- length(symbiontClade_diversity) # 241

symbiontClade_thermoTol <- symbiont_thermoTol$Phylotype
number_symbiontClade_thermoTol <- length(symbiontClade_thermoTol) # 110

### clades present in both list
cladeInCommon <- speciesPresentFunBis(symbiontClade_diversity,symbiontClade_thermoTol) 
length(cladeInCommon)

### species present in symbiontClade_diversity but not in symbiontClade_thermoTol
symbiont_List1_notList2 <- speciesNotPresentFun(symbiontClade_diversity,symbiontClade_thermoTol)
symbiont_List1_notList2 <- data.frame(symbiont_List1_notList2, stringsAsFactors = F)
# write.csv(symbiont_List1_notList2,paste(wd_symbiont_diversity_thermotolerence,"symbiont_List1_notList2.csv",sep="/"),row.names = F)

### species present in symbiontClade_thermoTol but not in symbiontClade_diversity
symbiont_notList1_List2 <- speciesNotPresentFun(symbiontClade_thermoTol,symbiontClade_diversity)[order(speciesNotPresentFun(symbiontClade_thermoTol,symbiontClade_diversity))]
symbiont_notList1_List2 <- data.frame(symbiont_notList1_List2, stringsAsFactors = F)
# write.csv(symbiont_notList1_List2,paste(wd_symbiont_diversity_thermotolerence,"symbiont_notList1_List2.csv",sep="/"),row.names = F)

symbiont_unmatchingLists <- cbind(symbiont_List1_notList2$symbiont_List1_notList2,c(symbiont_notList1_List2$symbiont_notList1_List2,rep(NA,nrow(symbiont_List1_notList2) - nrow(symbiont_notList1_List2))))

### import symbiont_List1_List2_Eugenia_Sampayo.csv 
symbiont_List1_List2_Eugenia_Sampayo <- read.csv(paste(wd,"Traits_extra/symbiodinium_subclade/symbiont_List1_List2_Eugenia_Sampayo.csv",sep="/"),header = T, stringsAsFactors = F)
# View(symbiont_List1_List2_Eugenia_Sampayo)

subcladesToRemoveSubcladeDiversity <- symbiont_List1_List2_Eugenia_Sampayo[,3:4][symbiont_List1_List2_Eugenia_Sampayo[,3:4]$comments.eugenia %in% c("never seen this name before",""),]
subcladesToRemoveSubcladeDiversity <- subcladesToRemoveSubcladeDiversity$symbiont_notList1_List2
subcladesToRemoveSubcladeDiversity <- subcladesToRemoveSubcladeDiversity[subcladesToRemoveSubcladeDiversity != ""]
subcladesToRemoveSubcladeDiversity

subcladesToRemoveSubcladeThermoTol <- symbiont_List1_List2_Eugenia_Sampayo[,1:2][symbiont_List1_List2_Eugenia_Sampayo[,1:2]$comments_eugenia != "keep",]
subcladesToRemoveSubcladeThermoTol <- subcladesToRemoveSubcladeThermoTol$symbiont_List1_notList2
subcladesToRemoveSubcladeThermoTol <- subcladesToRemoveSubcladeThermoTol[subcladesToRemoveSubcladeThermoTol != ""]
subcladesToRemoveSubcladeThermoTol

### remove or correct names in symbiont_thermoTol
symbiont_thermoTol_cut <- symbiont_thermoTol[!symbiont_thermoTol$Phylotype %in% subcladesToRemoveSubcladeDiversity,]
length(symbiont_thermoTol_cut$Phylotype) # 98
length(symbiont_thermoTol$Phylotype) # 110
sum(duplicated(symbiont_thermoTol_cut$Phylotype)) # 0
symbiont_thermoTol_cut$Phylotype[symbiont_thermoTol_cut$Phylotype == "C3-ff"] <- "C3ff"
symbiont_thermoTol_cut$Phylotype[symbiont_thermoTol_cut$Phylotype == "D1-4"] <- "D1a"

### remove these clades that Eugenia Sampayo did not ID in symbiodinium_subclade
colnmaes_symbiodinium_subclade <- colnames(symbiodinium_subclade)
length(colnmaes_symbiodinium_subclade) # 243
colnamesToKeep <- colnmaes_symbiodinium_subclade[! colnmaes_symbiodinium_subclade %in% subcladesToRemoveSubcladeThermoTol]
length(colnamesToKeep) # 209

symbiodinium_subclade_cut <- symbiodinium_subclade[,colnamesToKeep]
# View(symbiodinium_subclade_cut)

min(symbiodinium_subclade_cut$Number_Clades) # 1

### determine the symbiont in commons
symbiontCommon <- speciesPresentFunBis(symbiont_thermoTol_cut$Phylotype,colnames(symbiodinium_subclade_cut[,-c(1,length(symbiodinium_subclade_cut))]))
length(symbiontCommon) # 74

### only keep coral species with symbionts present in both list
colnamesToKeep <- c("species",symbiontCommon)
symbiodinium_subclade_cut_cut <- symbiodinium_subclade_cut[,colnamesToKeep]
symbiodinium_subclade_cut_cut$numberClades <- NA
for(i in 1:length(symbiodinium_subclade_cut_cut$species)){
  symbiodinium_subclade_cut_cut$numberClades[i] <- sum(symbiodinium_subclade_cut_cut[i,2:(length(symbiodinium_subclade_cut_cut)-1)])
}
# View(symbiodinium_subclade_cut_cut)
 
symbiodinium_subclade_cut_cut <- symbiodinium_subclade_cut_cut[symbiodinium_subclade_cut_cut$numberClades > 0, ]
length(symbiodinium_subclade_cut_cut$species) # 291
length(symbiodinium_subclade_cut_cut) - 2 # 74
colnames(symbiodinium_subclade_cut_cut)

### change from presence abscence to the thermotolerance ranking
symbiodinium_subclade_cut_cut_thermotTo <- symbiodinium_subclade_cut_cut
for(i in 1:length(symbiodinium_subclade_cut_cut_thermotTo$species)){
  for(j in 2:(length(symbiodinium_subclade_cut_cut_thermotTo)-1)){
    if(symbiodinium_subclade_cut_cut_thermotTo[i,j] == 1 ){
      print(paste(i,j))
      # find the corresponding ranking of symbiont 
      symbiontHere <- colnames(symbiodinium_subclade_cut_cut_thermotTo)[j]
      rank <- symbiont_thermoTol_cut$`Score (R k)`[symbiont_thermoTol_cut$Phylotype == symbiontHere]
      symbiodinium_subclade_cut_cut_thermotTo[i,j] <- rank
    }
  }
}
# View(symbiodinium_subclade_cut_cut)
# View(symbiodinium_subclade_cut_cut_thermotTo)

symbiodinium_subclade_cut_cut_thermotTo$meanThermoRank <- NA
symbiodinium_subclade_cut_cut_thermotTo$maxThermoRank <- NA
for(i in 1:length(symbiodinium_subclade_cut_cut_thermotTo$species)){
  preenceAbsc <- as.numeric(symbiodinium_subclade_cut_cut_thermotTo[i,2:(length(symbiodinium_subclade_cut_cut_thermotTo)-3)])
  symbiodinium_subclade_cut_cut_thermotTo$meanThermoRank[i] <- sum(preenceAbsc)/length(preenceAbsc[preenceAbsc != 0])
  symbiodinium_subclade_cut_cut_thermotTo$maxThermoRank[i] <- max(symbiodinium_subclade_cut_cut_thermotTo[i,2:(length(symbiodinium_subclade_cut_cut_thermotTo)-3)])
}

# write.csv(symbiodinium_subclade_cut_cut_thermotTo[,c("species","numberClades","meanThermoRank","maxThermoRank")],paste(wd_Datasets,"symbiodinium_subcalssDiversity_thermotol.csv",sep="/"),row.names = F)

plot(symbiodinium_subclade_cut_cut_thermotTo$meanThermoRank~symbiodinium_subclade_cut_cut_thermotTo$maxThermoRank)
plot(symbiodinium_subclade_cut_cut_thermotTo$meanThermoRank~symbiodinium_subclade_cut_cut_thermotTo$numberClades)
plot(symbiodinium_subclade_cut_cut_thermotTo$maxThermoRank~symbiodinium_subclade_cut_cut_thermotTo$numberClades)


