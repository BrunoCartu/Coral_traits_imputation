# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org
## produce SkeletalDensity_original.csv and Skelelal_density.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_Skeletal_density <- paste(wd,"/Traits_extra/skeletal_density",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### Import original dataset from coraltraits.org: (April 8 2017)
# skeletalDensity_original <- read.csv("https://coraltraits.org/traits/61.csv", as.is=TRUE)
# write.csv(skeletalDensity_original,paste(wd_Datasets_original,"SkeletalDensity_original.csv",sep="/"),row.names = F)

### Import original dataset from folder:
skeletalDensity_original <- read.csv(paste(wd_Datasets_original,"SkeletalDensity_original.csv",sep="/"),header = T, stringsAsFactors = F)

nameColumnSelectedSkeletalDensity <- c("specie_name","trait_name","value","value_type","precision","precision_type","notes","replicates","standard_unit")

skeletalDensity <- skeletalDensity_original[nameColumnSelectedSkeletalDensity]

# check if other traits:
levels(as.factor(skeletalDensity$trait_name))
skeletalDensity <- subset(skeletalDensity,skeletalDensity$trait_name =="Skeletal density")

# check f different value_type:
levels(as.factor(skeletalDensity$value_type))
#  "mean"      "median"    "raw_value"
levels(as.factor(skeletalDensity$standard_unit))
# "g cm-3"  "g mL^-1"  --> same dimensions
levels(as.factor(skeletalDensity$precision_type))
# ""                   "range"              "standard_deviation" "standard_error"  

# number of species:
length(levels(as.factor(skeletalDensity$specie_name)))
# 54
length(skeletalDensity$specie_name)
# 422

skeletalDensity$value <- as.numeric(skeletalDensity$value)

combinedDFs <- skeletalDensity

# creation of final table to determine the mean value for each species:
# mean or median values without sample size or precison type indicated, as well as raw_values are considered as a single observation
# if "the"precision_type" == range and mean value given because I supposed normal distribution and sd <- range / 6 
# median values with sample size are conidered as mean because I supposed normal distribution
species <- levels(as.factor(combinedDFs$specie_name))
numberSpecies <- length(species) # 54
numberLines <- length(combinedDFs$specie_name) # 422
meanSDnTable <- data.frame(species= rep(NA,numberSpecies),mean_SD=rep(NA,numberSpecies),sd=rep(NA,numberSpecies),n=rep(NA,numberSpecies))
for(i in 1:numberSpecies){
  # make a df for a particular species, each row corresponding to a value (either mean, median or raw values)
  temporalDF <- subset(combinedDFs,combinedDFs$specie_name==species[i])
  numberLines <- length(temporalDF$specie_name)
  meanSDnTable_temporal <- data.frame(mean_SD=rep(NA,numberLines),sd=rep(NA,numberLines),n=rep(NA,numberLines))
  for(j in 1:numberLines){
    vt_mean <- temporalDF$value_type[j] == "mean"
    vt_median <- temporalDF$value_type[j] == "median"
    pt_range <- temporalDF$precision_type[j] == "range"
    pt_sd <- temporalDF$precision_type[j] == "standard_deviation"
    pt_se <- temporalDF$precision_type[j] == "standard_error"
    pt_NA <- temporalDF$precision_type[j] == ""
    replicates_n <- !is.na(temporalDF$replicates[j])
    
    meanSDnTable_temporal$mean_SD[j] <- temporalDF$value[j]
    
    if(replicates_n & !pt_NA & (vt_mean | vt_median)){
      meanSDnTable_temporal$n[j] <- temporalDF$replicates[j]
      # if median but n and precision_type given, then median is considered as the mean, because I assumed a normal distribution
      if(pt_range){
        meanSDnTable_temporal$sd[j] <- temporalDF$precision[j]/6       # sd = range/6
      }else if(pt_sd){
        meanSDnTable_temporal$sd[j] <- temporalDF$precision[j]     
      }else if(pt_se){
        meanSDnTable_temporal$sd[j] <- temporalDF$precision[j]*sqrt(temporalDF$replicates[j])   # sd = sqrt(n)*SE
      }
    }else{
      meanSDnTable_temporal$sd[j] <- NA
      meanSDnTable_temporal$n[j] <- 1
    }
  }
  meanSDnTable$species[i] <- species[i]
  meanSDnTable[i,c(2:4)] <- meanSDnFund(meanSDnTable_temporal)
}
sum(duplicated(meanSDnTable$species)) # 0

length(meanSDnTable$species) # 54

### check for names:
tableNamesChecked <- meanSDnTable

speciesToCheck <- checkSpeciesNames(tableNamesChecked$species,wd_Datasets = wd_Datasets) # "All species are present in the list"
speciesToCheck

### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep='/'),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp) 

### correct names in species_colonialityFinal_cut:
species <- tableNamesChecked$species
numberSpecies <- length(species) # 54
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] # only keep the species in vector species
length(speciesChecked$nameSp) # 54
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}

### check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 0

### alphabetic order:
tableNamesChecked$species <- gsub(" ","zzz",tableNamesChecked$species)
tableNamesChecked <- tableNamesChecked[order(tableNamesChecked$species),]
tableNamesChecked$species <- gsub("zzz"," ",tableNamesChecked$species)
length(tableNamesChecked$species) # 54
sum(duplicated(tableNamesChecked$species)) # 0

colnames(tableNamesChecked) <- c("species","skeletalDensity_gcm3","sd","n")

# write.csv(tableNamesChecked,paste(wd_Datasets,"Skelelal_density.csv",sep="/"),row.names = F)



