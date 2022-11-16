# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals: 
## import and curate taxon-RBI values from Swain et al. 2016 (TableS4) to produce bleaching_response_index_313sp.csv
## combined the latter dataset with values from Marcelino et al., 2013 to produce bleaching_response_index_328sp.csv (NOT USED)

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_taxon_BRI <- paste(wd,"/Traits_extra/bleaching_response_index",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

require(dplyr)
require(stringr)

### Import taxon-RBI values from (Swain et al. 2016_TableS4)
table <- read.csv(paste(wd_Datasets_original,"Swain et al. 2016_TableS4.csv",sep="/"), header= TRUE,stringsAsFactors = F)
head(table)
CoralBleachingIndex <- table[1:374,1:4]
tail(CoralBleachingIndex)
colnames(CoralBleachingIndex) <- c("species","BRI","n",'sd')
head(CoralBleachingIndex)
CoralBleachingIndex <- CoralBleachingIndex[,c(1,2,4,3)]
length(CoralBleachingIndex$species) # 374

### there are raws with only the Genus name --> the one with not " " --> to remove
# library(tidyverse)
CoralBleachingIndex_cut <- CoralBleachingIndex %>% filter(str_detect(species, ' '))
length(CoralBleachingIndex_cut$species) # 316

tableNamesChecked <- CoralBleachingIndex_cut
speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets = wd_Datasets) 
speciesToCheck 

# arange names in alphabetic order
tableNamesChecked$species <- gsub(" ","_",tableNamesChecked$species)
tableNamesChecked <- tableNamesChecked[order(tableNamesChecked$species),]
tableNamesChecked$species <- gsub("_"," ",tableNamesChecked$species)

### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp) # 1767

### correct names in species_colonialityFinal_cut:
sp <- tableNamesChecked$species
numberSpecies <- length(sp) # 316
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% sp,] # only keep the species in vector species
length(speciesChecked$nameSp) # 316
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}

# Check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 3
speciesDuplicated <- tableNamesChecked$species[duplicated(tableNamesChecked$species)]
speciesDuplicated <- levels(as.factor((speciesDuplicated)))
length(speciesDuplicated) # 3

duplicates <- subset(tableNamesChecked,tableNamesChecked$species %in% speciesDuplicated)
duplicates
#                   species      BRI       sd  n
# 73    Colpophyllia natans  2.75    NA  1
# 74    Colpophyllia natans 17.97 16.93 22
# 88   Dichocoenia stokesii  2.75    NA  1
# 89   Dichocoenia stokesii 10.25 15.90  5
# 260 Pocillopora verrucosa 39.84    NA  1
# 265 Pocillopora verrucosa 25.82 20.93 28

### combine the values:
# tableNamesChecked_bis <- tableNamesChecked
# tableNamesChecked <- tableNamesChecked_bis
sp <- as.character(levels(as.factor(duplicates$species)))
numberSpecies <- length(sp) # 3
for(i in 1:numberSpecies){
  tableNamesChecked <- tableNamesChecked[tableNamesChecked$species != sp[i],] # remove the species from tableNamesChecked
  subsetTable <- subset(duplicates,duplicates$species == sp[i])
  combinedValues <- meanSDnFund(subsetTable[,c(2,3,4)])
  combinedValues$species <- sp[i]
  combinedValues <- combinedValues[,c(4,1,2,3)]
  colnames(combinedValues)[2] <- "BRI"
  tableNamesChecked <- rbind(tableNamesChecked,combinedValues)
}

tableNamesChecked$species <- gsub(" ","_",tableNamesChecked$species)
tableNamesChecked <- tableNamesChecked[order(tableNamesChecked$species),]
tableNamesChecked$species <- gsub("_"," ",tableNamesChecked$species)

length(tableNamesChecked$species) # 313
sum(duplicated(tableNamesChecked$species)) # 0

# View(tableNamesChecked)

# write.csv(tableNamesChecked,paste(wd_Datasets,"bleaching_response_index_313sp.csv",sep="/"),row.names=F)

#***********************************************************
#************* with values from Marcelino et al., 2013 *****
#*********** TO NOT USE BECAUSE MANY VALUES ARE AT THE GENUS LEVEL **********
#***********************************************************

### Import list from Marcelino et al., 2013 about RSC
RSC <- read.csv(paste(wd_Datasets_original,"Marcelino_et_al_2013_Table_S1.csv",sep="/"),header=T,stringsAsFactors = F)
head(RSC)
RSC_cut <- RSC[,c(1,4,5,6)]
RSC_cut <- RSC_cut[,c(1,2,4,3)]
head(RSC_cut)
colnames(RSC_cut) <- c("species","BRI","SE","n")

### correct for certain names under column species: 
# Acropora intermedia/nobilis <- Acropora intermedia  (nobilis is not accepted and should be robusta)
# Acropora sp. -> to remove
# Euphyllia sp. -> to remove
# Montipora sp. -> to remove
# Mycedium elephantotus/mancaoi <- Mycedium elephantotus (both are correct names but I have more information for the 1st name)
# Porites sp.  -> to remove
# Siderastrea sp. -> to remove

RSC_cut$species[RSC_cut$species == "Acropora intermedia/nobilis"] <- "Acropora intermedia"
RSC_cut$species[RSC_cut$species == "Mycedium elephantotus/mancaoi"] <- "Mycedium elephantotus"

length(RSC_cut$species)  # 150

RSC_cut <- subset(RSC_cut,species != "Acropora sp.")
RSC_cut <- subset(RSC_cut,species != "Euphyllia sp.")
RSC_cut <- subset(RSC_cut,species != "Montipora sp.")
RSC_cut <- subset(RSC_cut,species != "Porites sp.")
RSC_cut <- subset(RSC_cut,species != "Siderastrea sp.")

length(RSC_cut$species)  # 143

# in this dataset, several rows are replicated but does not correspond to different measurements for the same species --> combine them doing the mean (which is an easy way to get rid of the replicates)
mean_BRI <- tapply(RSC_cut$BRI,RSC_cut$species,mean)
mean_SE <- tapply(RSC_cut$SE,RSC_cut$species,mean)
mean_n <- tapply(RSC_cut$n,RSC_cut$species,mean)
mean_sd <- mean_SE * sqrt(mean_n)

CoralBleachingIndex_2 <- data.frame(
  species = as.character(levels(as.factor(RSC_cut$species))),
  BRI = mean_BRI,
  sd = mean_sd,
  n = mean_n,
  stringsAsFactors = F
)
rownames(CoralBleachingIndex_2) <- NULL
length(CoralBleachingIndex_2$species) # 89
head(CoralBleachingIndex_2)

### Add species from Swain et al., 2016 about RSC that are not in Figueiredo et al., 2013 and Swain et al., 2016 about BRI ----
speciesFromSwain2016 <- c("Turbinaria reniformis","Favia favus")   # Montipora digitata values differe a little between Swain et al., 2016 in BRI and RSC... I keep the one in BRI
BRI_Swain <- c(33.96,27.85)
speciesFromSwain2016  <- data.frame(
  species = speciesFromSwain2016,
  BRI = BRI_Swain,
  sd = NA,
  n = c(19,6),
  stringsAsFactors = F
)

### comnine the 3 datasets:
ds_1 <- CoralBleachingIndex_cut
length(ds_1$species) # 316
ds_2 <- CoralBleachingIndex_2
length(ds_2$species) # 89
ds_3 <- speciesFromSwain2016
length(ds_3$species) # 2

totalDataset_BRI <- merge(ds_1,ds_2,all.x = T, all.y = T,by="species")
totalDataset_BRI <- merge(totalDataset_BRI,ds_3,all.x = T, all.y = T,by="species")
head(totalDataset_BRI)
totalDataset_BRI <- totalDataset_BRI[,c("species","BRI.x","BRI.y","BRI","sd.x","sd.y","sd","n.x","n.y","n")]
head(totalDataset_BRI)
length(totalDataset_BRI$species) # 338

### certain species have different values of BRI depending from which dataset they are from --> average values accounting for sample size:
sum(duplicated(totalDataset_BRI$species)) # 0
colnames <- c("species","BRI","sd","n")
sp <- totalDataset_BRI$species
numberSpecies <- length(sp) # 338
BRI_final <- matrix(ncol=length(colnames),nrow=numberSpecies)
BRI_final <- as.data.frame(BRI_final)
colnames(BRI_final) <- colnames
BRI_final$species <- sp
for(i in 1:numberSpecies){
  subsetBRI <- subset(totalDataset_BRI,totalDataset_BRI$species == sp[i])
  BRI_final$species[i] <- sp[i]
  bool1 <- subsetBRI$BRI.x != subsetBRI$BRI.y & !is.na(subsetBRI$BRI.x) & !is.na(subsetBRI$BRI.y)
  bool2 <- subsetBRI$BRI.x != subsetBRI$BRI & !is.na(subsetBRI$BRI.x) & !is.na(subsetBRI$BRI)
  bool3 <- subsetBRI$BRI.y != subsetBRI$BRI & !is.na(subsetBRI$BRI.y) & !is.na(subsetBRI$BRI)
  bool4 <- (!is.na(subsetBRI$BRI.x) & !is.na(subsetBRI$BRI.y) & !is.na(subsetBRI$BRI)) & ( subsetBRI$BRI.x != subsetBRI$BRI.y | subsetBRI$BRI.x != subsetBRI$BRI | subsetBRI$BRI.y != subsetBRI$BRI ) # not supposed to happen
  if(bool1){
    df_sd <- data.frame(mean=c(subsetBRI$BRI.x,subsetBRI$BRI.y),sd = c(subsetBRI$sd.x,subsetBRI$sd.y),n=c(subsetBRI$n.x,subsetBRI$n.y))
  }else if(bool2){
    df_sd <- data.frame(mean=c(subsetBRI$BRI.x,subsetBRI$BRI),sd = c(subsetBRI$sd.x,subsetBRI$sd),n=c(subsetBRI$n.x,subsetBRI$n))
  }else if(bool3){
    df_sd <- data.frame(mean=c(subsetBRI$BRI.y,subsetBRI$BRI),sd = c(subsetBRI$sd.y,subsetBRI$sd),n=c(subsetBRI$n.y,subsetBRI$n))
  }else if(bool4){
    print(paste("there case where all 3 values are different, for species: ",species[i]))
  }
  if(bool1 | bool2 | bool3){
    # calculate [,1] the mean, [,2] the sd and [,3] the n of multiple samples (see functions.R)
    df_sd_combined <- meanSDnFund(df_sd)
    BRI_final$BRI[i] <- df_sd_combined$mean
    BRI_final$sd[i] <- df_sd_combined$sd
    BRI_final$n[i] <- df_sd_combined$n
  }else if (!bool4){  # cases where only one value is available or it is the same value for different datasets
    if(!is.na(subsetBRI$BRI.x)){
      BRI <- subsetBRI$BRI.x
      sd <- subsetBRI$sd.x
      n <- subsetBRI$n.x
    }else if (!is.na(subsetBRI$BRI.y)){
      BRI <- subsetBRI$BRI.y
      sd <- subsetBRI$sd.y
      n <- subsetBRI$n.y
    }else if (!is.na(subsetBRI$BRI)){
      BRI <- subsetBRI$BRI
      sd <- subsetBRI$sd
      n <- subsetBRI$n
    }
    BRI_final$BRI[i] <- BRI
    BRI_final$sd[i] <- sd
    BRI_final$n[i] <- n
  }
}

length(BRI_final$species) # 338

### check for names:
tableNamesChecked <- BRI_final
speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets = wd_Datasets) 
speciesToCheck 
# arange names in alphabetic order
tableNamesChecked$species <- gsub(" ","_",tableNamesChecked$species)
tableNamesChecked <- tableNamesChecked[order(tableNamesChecked$species),]
tableNamesChecked$species <- gsub("_"," ",tableNamesChecked$species)

### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp) # 1767

### correct names in species_colonialityFinal_cut:
species <- tableNamesChecked$species
numberSpecies <- length(species) # 338
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] # only keep the species in vector species
length(speciesChecked$nameSp) # 338
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}

# Check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 10
speciesDuplicated <- tableNamesChecked$species[duplicated(tableNamesChecked$species)]
speciesDuplicated <- levels(as.factor((speciesDuplicated)))
length(speciesDuplicated) # 10

duplicates <- subset(tableNamesChecked,tableNamesChecked$species %in% speciesDuplicated)
duplicates$species <- gsub(" ","_",duplicates$species)
duplicates <- duplicates[order(duplicates$species),]
duplicates$species <- gsub("_"," ",duplicates$species)
duplicates
#                   species      BRI       sd  n
# 77      Coelastrea aspera 14.27000 16.87000  4
# 142     Coelastrea aspera 14.27000 16.88000  4
# 79    Colpophyllia natans  2.75000       NA  1
# 80    Colpophyllia natans 17.97000 16.93000 22
# 94   Dichocoenia stokesii  2.75000       NA  1
# 95   Dichocoenia stokesii 10.25000 15.90000  5
# 101     Dipsastraea favus 32.05000 28.24000  6
# 122     Dipsastraea favus 27.85000       NA  6
# 107   Dipsastraea pallida 55.04000 20.23000 12
# 124   Dipsastraea pallida 53.56000 20.08292 13
# 42       Isopora palifera 40.85000 22.82541  8
# 166      Isopora palifera 40.85000 22.82000  8
# 202 Montastraea cavernosa 14.42000 12.88000 15
# 204 Montastraea cavernosa 14.42000 12.89703 15
# 203   Orbicella annularis 21.96000 14.73308 19
# 238   Orbicella annularis 18.89000 17.29000 36
# 276 Pocillopora verrucosa 39.84000       NA  1
# 281 Pocillopora verrucosa 25.82000 20.93000 28
# 312   Psammocora contigua 34.95362 21.98485 58
# 314   Psammocora contigua 35.31000 21.87779 55

### combine the values:
# tableNamesChecked_bis <- tableNamesChecked
# tableNamesChecked <- tableNamesChecked_bis
sp <- as.character(levels(as.factor(duplicates$species)))
numberSpecies <- length(sp) # 10
for(i in 1:numberSpecies){
  tableNamesChecked <- tableNamesChecked[tableNamesChecked$species != sp[i],] # remove the species from tableNamesChecked
  subsetTable <- subset(duplicates,duplicates$species == sp[i])
  combinedValues <- meanSDnFund(subsetTable[,c(2,3,4)])
  combinedValues$species <- sp[i]
  combinedValues <- combinedValues[,c(4,1,2,3)]
  colnames(combinedValues)[2] <- "BRI"
  tableNamesChecked <- rbind(tableNamesChecked,combinedValues)
}

tableNamesChecked$species <- gsub(" ","_",tableNamesChecked$species)
tableNamesChecked <- tableNamesChecked[order(tableNamesChecked$species),]
tableNamesChecked$species <- gsub("_"," ",tableNamesChecked$species)

length(tableNamesChecked$species) # 328
sum(duplicated(tableNamesChecked$species)) # 0

# View(tableNamesChecked)

# write.csv(tableNamesChecked,paste(wd_Datasets,"bleaching_response_index_328sp.csv",sep="/"),row.names=F)

#***********************************************************
#************* TO SHOW DISCREPANCY BETWEEN SWAIN ET AL 2016 AND MARCELINO ET AL 2013 BRI DATASETS *****
#***********************************************************------
combineTest <-  merge(CoralBleachingIndex_cut,CoralBleachingIndex_2,by="species",all.x=T,all.y = T)
plot(combineTest$BRI.x,combineTest$BRI.y)
#***********************************************************

BRI_328 <- read.csv(paste(wd_Datasets,"bleaching_response_index_328sp.csv",sep="/"),header = T, stringsAsFactors = F)
BRI_313 <- read.csv(paste(wd_Datasets,"bleaching_response_index_313sp.csv",sep="/"),header = T, stringsAsFactors = F)
RSC <- read.csv(paste(wd_Datasets,"reduced_scattering_coefficient.csv",sep="/"),header = T, stringsAsFactors = F)

speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header = T, stringsAsFactors = F)
speciesNameChecked

combinedDS <- merge(BRI_328,BRI_313,by="species",all.x=T,all.y = T)
combinedDS <- merge(combinedDS,RSC,by="species",all.x=T,all.y = T)
colnames(combinedDS)
# "species"     "BRI.x"       "sd.x"        "n.x"         "BRI.y"       "sd.y"        "n.y"         "mean_RSC"    "sd"          "n"           "BRI.x_logit"

head(combinedDS)

plot(combinedDS$BRI.x ~ combinedDS$BRI.y,ylab="328 species DS",xlab="313 species DS")

points(x=11.99,y=9.098462,col="red")

segments(x0= 10,x1=10,y0=0,y1=40)
segments(x0= 0,x1=40,y0=10,y1=10)

combinedDS[combinedDS$BRI.x > 10 & combinedDS$BRI.y < 10,]

#                   species           BRI.x     sd.x n.x BRI.y  sd.y n.y mean_RSC sd  n
# 61           Agaricia fragilis 28.93172 26.86973  58  7.33 5.72   3     8.52 0.863  2
# 78         Coscinaraea columna 24.80700 26.97543  10  6.51   NA   1     9.80    NA  1
# 84      Cyphastrea chalcidicum 16.38929 25.46170  14  0.00   NA   1     9.64    NA  1
# 86    Cyphastrea microphthalma 16.85643 25.21417  14  6.54   NA   1     6.37 2.800  2
# 186        Merulina scabricula 39.45250 42.23078   8  9.09   NA   1     6.10 1.485  2
# 219        Montipora verrucosa 34.22516 25.85818  95  2.75   NA   1     8.27    NA  1
# 243           Pavona decussata 36.83450 26.27304 111  5.43   NA   1     3.02    NA  1

