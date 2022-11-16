# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals: 
## import and curate species trait data from coraltraits.org
## produce chlorophyll_a_concentration_original.csv and chlorophyll_a_concentration.csv

require(here)     
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_Chlorophyll_a_concentratio <- paste(wd,"/Traits_extra/chlorophyll_a_concentration",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### Import original dataset from coraltraits.org: (April 10 2017)
# Chlorophyll_a_concentration_original <- read.csv("https://coraltraits.org/traits/146.csv", as.is=TRUE)
# Def: The concentration of chlorophyll from the Symbiodinium (and potentially other endolithic algae) per unit of skeletal surface area.
# write.csv(Chlorophyll_a_concentration_original,paste(wd_Datasets_original,"chlorophyll_a_concentration_original.csv",sep="/"),row.names = F)

### Import orginal dataset from folder:
Chlorophyll_a_concentration_original <- read.csv(paste(wd_Datasets_original,"chlorophyll_a_concentration_original.csv",sep="/"),header=T,stringsAsFactors = F)

nameColumnSelected <- c("specie_name","trait_name","standard_unit","value","value_type","precision","precision_type","replicates","notes")

Chlorophyll_a_concentration_cut <- Chlorophyll_a_concentration_original[nameColumnSelected]

###
levels(as.factor(Chlorophyll_a_concentration_cut$trait_name))
Chlorophyll_a_concentration_cut <- subset(Chlorophyll_a_concentration_cut,Chlorophyll_a_concentration_cut$trait_name == "Chlorophyll a")

###
levels(as.factor(Chlorophyll_a_concentration_cut$value_type)) # "mean"      "raw_value"

### 
levels(as.factor(Chlorophyll_a_concentration_cut$precision_type)) # ""                   "95_ci"              "standard_deviation"

### 
levels(as.factor(Chlorophyll_a_concentration_cut$replicates)) #   "3"   "4"   "5"   "6"   "12"  "100"

###
levels(as.factor(Chlorophyll_a_concentration_cut$standard_unit)) # "mg m^-2"  "µg cm^-2"

Chlorophyll_a_concentration_cut$value <- as.numeric(Chlorophyll_a_concentration_cut$value)
Chlorophyll_a_concentration_cut$precision <- as.numeric(Chlorophyll_a_concentration_cut$precision)
Chlorophyll_a_concentration_cut$replicates <- as.numeric(Chlorophyll_a_concentration_cut$replicates)

### conversion from "mg m^-2" to "µg cm^-2"
Chlorophyll_a_concentration_cut[Chlorophyll_a_concentration_cut$standard_unit == "mg m^-2" ,c("value","precision")] <- Chlorophyll_a_concentration_cut[Chlorophyll_a_concentration_cut$standard_unit == "mg m^-2" ,c("value","precision")] / 10
Chlorophyll_a_concentration_cut[Chlorophyll_a_concentration_cut$standard_unit == "mg m^-2" ,"standard_unit"] <- "µg cm^-2"

# creation of final table to determine the mean value for each species:
# mean without sample size or precison type indicated, as well as raw_values are considered as a single observation
combinedDFs <- Chlorophyll_a_concentration_cut
species <- levels(as.factor(combinedDFs$specie_name))
numberSpecies <- length(species) # 22
numberLines <- length(combinedDFs$specie_name) # 110
meanSDnTable <- data.frame(species= rep(NA,numberSpecies),mean_chloro = rep(NA,numberSpecies),sd=rep(NA,numberSpecies),n=rep(NA,numberSpecies))
for(i in 1:numberSpecies){
  # make a df for a particular species, each row corresponding to a value (either mean, median or raw values)
  temporalDF <- subset(combinedDFs,combinedDFs$specie_name == species[i])
  numberLines <- length(temporalDF$specie_name)
  meanSDnTable_temporal <- data.frame(mean_chloro=rep(NA,numberLines),sd=rep(NA,numberLines),n=rep(NA,numberLines))
  for(j in 1:numberLines){
    vt_mean <- temporalDF$value_type[j] == "mean"
    vt_median <- temporalDF$value_type[j] == "median"
    pt_95CI <- temporalDF$precision_type[j] == "95_ci"
    pt_range <- temporalDF$precision_type[j] == "range"
    pt_sd <- temporalDF$precision_type[j] == "standard_deviation"
    pt_se <- temporalDF$precision_type[j] == "standard_error"
    pt_NA <- temporalDF$precision_type[j] == "not_given" | temporalDF$precision_type[j] == ""
    replicates_n <- !is.na(temporalDF$replicates[j])
    
    meanSDnTable_temporal$mean_chloro[j] <- temporalDF$value[j]
    
    if(replicates_n & !pt_NA & (vt_mean | vt_median)){
      meanSDnTable_temporal$n[j] <- temporalDF$replicates[j]
      # if median but n and precision_type given, then median is considered as the mean, because I assumed a normal distribution
      if(pt_95CI){  # “Rule of thumb” 95% confidence interval using these SE values: --> 95%CI = 4SE <--> 95%CI/4 = SE <--> 95%CI/4*sqrt(n) = sd
        meanSDnTable_temporal$sd[j] <- temporalDF$precision[j]/4*sqrt(temporalDF$replicates[j])
      }else if(pt_range){
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

### check for names:
tableNamesChecked <- meanSDnTable
speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets = wd_Datasets) # "All species are present in the list"
speciesToCheck
### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp)

### correct names in species_colonialityFinal_cut:
sp <- tableNamesChecked$species
numberSpecies <- length(sp) # 22
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% sp,] # only keep the species in vector species
length(speciesChecked$nameSp) # 22
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}
# Check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 3

### alphabetic order:
tableNamesChecked$species <- gsub(" ","_",tableNamesChecked$species)
tableNamesChecked <- tableNamesChecked[order(tableNamesChecked$species),]
tableNamesChecked$species <- gsub("_"," ",tableNamesChecked$species)

# creation csv file in the wd dataCompilation
# write.csv(tableNamesChecked, file = paste(wd_Datasets,"chlorophyll_a_concentration.csv",sep="/"),row.names=F)



