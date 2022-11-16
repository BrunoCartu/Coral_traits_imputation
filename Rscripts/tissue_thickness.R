# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org, and from (Loya et al., 2001)
## produce age_at_maturity_original.csv and age_at_maturity.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_Tissue_thickness <- paste(wd,"/Traits_extra/tissue_thickness",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### import original dataset from coraltraits.org: (March 30 2017)
# TissThisk_original <- read.csv("https://coraltraits.org/traits/132.csv", as.is = TRUE)
# View(TissThisk_original)
# write.csv(TissThisk_original,paste(wd_Datasets_original,"Tissue_thickness_original.csv",sep="/"),row.names = F)

### Import original dataset form folder:
TissThisk_original <- read.csv(paste(wd_Datasets_original,"Tissue_thickness_original.csv",sep="/"),header=T,stringsAsFactors = F)

levels(as.factor(TissThisk_original$trait_name))

TissThisk_1 <- subset(TissThisk_original,TissThisk_original$trait_name =="Tissue thickness")

TissThisk_1$value <- as.numeric(TissThisk_1$value)

colnames(TissThisk_1)

colnamesToKeep <- c("specie_name","methodology_name","value","value_type","precision","precision_type","replicates")

TissThisk_2 <- TissThisk_1[,colnamesToKeep]

length(levels(as.factor(TissThisk_2$specie_name))) # 16

levels(as.factor(TissThisk_2$value_type))
# "mean"      "raw_value"
levels(as.factor(TissThisk_2$precision_type))
# ""                   "standard_deviation" "standard_error" 
levels(as.factor(TissThisk_2$replicates))
# "5" "7" "8" "9"
levels(as.factor(TissThisk_1$standard_unit))
"mm"

# creation of final table to determine the mean value for each species:
# mean values are multiplied by sample size when available
# mean values without sample size, medians and raw values are considered as a single observation
species <- levels(as.factor(TissThisk_2$specie_name))
numberSpecies <- length(species) # 16
numberLines <- length(TissThisk_2$specie_name) # 59
TissThisk_3 <- data.frame(species = as.character(rep(NA,numberSpecies)),meanTissueThick = as.numeric(rep(NA,numberSpecies)), sd = as.numeric(rep(NA,numberSpecies)),n = as.numeric(rep(NA,numberSpecies)),stringsAsFactors=FALSE)
for(i in 1:numberSpecies){
  subTable <- subset(TissThisk_2,TissThisk_2$specie_name == species[i])
  numberRaw <- length(subTable$specie_name)
  meanSDn_s <- data.frame(mean=rep(NA,numberRaw),sd=rep(NA,numberRaw),n=rep(NA,numberRaw))
  for(j in 1:numberRaw){
    bool1 <- subTable$value_type[j] == "raw_value"
    bool2 <- subTable$value_type[j] == "mean" && is.na(subTable$replicates[j]) # some values are means but there is not sample size,so they are considered as raw observations
    # bool3 <- TissThisk_2$value_type[j] == "mean" && TissThisk_2$precision_type[j] == "" # none of the raw have means with "" for precision_type
    # TissThisk_2[TissThisk_2$value_type == "mean" & TissThisk_2$precision_type == "",]
    meanSDn_s$mean[j] <- subTable$value[j]
    if(bool1 || bool2){
      meanSDn_s$sd[j] <- 0
      meanSDn_s$n[j] <- 1
    }else {
      meanSDn_s$n[j] <- subTable$replicates[j]
      meanSDn_s$sd[j] <- subTable$precision[j]
      if(subTable$precision_type[j] == "standard_error"){
        meanSDn_s$sd[j] <- meanSDn_s$sd[j] * sqrt(meanSDn_s$n[j])
      }
    }
  }
  TissThisk_3[i,1] <- species[i]
  TissThisk_3[i,c(2:4)] <-  meanSDnFund(meanSDn_s)
}

### species to add from the literature:
# Loya et al., 2001: (n=5)
# Platygyra ryukyuensis 4.0 (SD = 0.096)
# Favites chinensis   9.4 (SD = 0.174)
# Dipsastraea pallida (Favia pallida)  3.1 (SD = 0.011)
# Stylophora pistillata 1.5 (SD = 0.279)

TissThisk_3 <- rbind(TissThisk_3,c("Platygyra ryukyuensis",4.0,0.096/sqrt(5),5))
TissThisk_3 <- rbind(TissThisk_3,c("Favites chinensis",9.4,0.174/sqrt(5),5))
TissThisk_3 <- rbind(TissThisk_3,c("Dipsastraea pallida",3.1,0.011/sqrt(5),5))
TissThisk_3 <- rbind(TissThisk_3,c("Stylophora pistillata",1.5,0.279/sqrt(5),5))

### put in alphabetic order:
TissThisk_3$species <- gsub(" ","zzz",TissThisk_3$species)
TissThisk_3 <- TissThisk_3[order(TissThisk_3$species),]
TissThisk_3$species <- gsub("zzz"," ",TissThisk_3$species)

length(TissThisk_3$species) # 20

tableNamesChecked <- TissThisk_3

### check for names:
speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets =  wd_Datasets)
speciesToCheck

### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp)

### correct names in species_colonialityFinal_cut:
species <- tableNamesChecked$species
numberSpecies <- length(species) # 20
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] # only keep the species in vector species
length(speciesChecked$nameSp) # 20
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}

### check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 0

# write.csv(tableNamesChecked,paste(wd_Datasets,"tissue_thickness.csv",sep="/"),row.names=F)





