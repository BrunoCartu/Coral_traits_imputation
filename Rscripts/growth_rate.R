## Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org
## produce Growth_Rate_original.csv and growth_rate.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_Growth_rate <- paste(wd,"/Traits_extra/growth_rate",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### Final result: mean linear (= diameter) growth rate in mm yr^-1

### import original dataset from caraltraits.org: (April 3 2017)
# dataGrowthRate_original <- read.csv("https://coraltraits.org/traits/60.csv", as.is=TRUE)
# write.csv(dataGrowthRate_original,"paste(wd_Datasets_original,Growth_Rate_original.csv",sep="/"),row.names = F)

### import original dataset from folder
dataGrowthRate_original <- read.csv(paste(wd_Datasets_original,"Growth_Rate_original.csv",sep="/"),header=T,stringsAsFactors = F)

nameColumnSelectedGR <- c("specie_name","trait_name","standard_id","standard_unit","value","value_type","precision","precision_type","replicates","notes")

dataGRCutColumn <- dataGrowthRate_original[nameColumnSelectedGR]

# check the different traits:
levels(as.factor(dataGRCutColumn$trait_name))
dataGRCutColumn <- subset(dataGRCutColumn, trait_name =="Growth rate")

# number of lines :
length(dataGRCutColumn$specie_name)    # 1552

# check the different value_type:
levels(as.factor(dataGRCutColumn$value_type))     #    ""          "mean"      "median"    "raw_value"

# check the different precision_type
levels(as.factor(dataGRCutColumn$precision_type))
# ""                   "95_ci"              "not_given"          "range"              "standard_deviation" "standard_error"  

# check the different replicates
levels(as.factor(dataGRCutColumn$replicates))

levels(as.factor(dataGRCutColumn$standard_unit))
# "mm d^-1", "mm month^-1", "mm per 3 months^-1","mm per 4 months^-1", "mm per 6 months^-1", "mm per 7 months^-1" "mm per 8 months^-1", "mm yr^-1" 

levels(as.factor(dataGRCutColumn$standard_id))
# "38" "49" "60" "62" "63" "80" "81" "83" "84" "85" "86" --> either radial or linear --> covert all of them in linear 

### conversion of radial measures into linear (diameters) <- x2
numberRows <- length(dataGRCutColumn$specie_name) # 1552
dataGRCutColumn$value <- as.numeric(dataGRCutColumn$value)
for(i in 1:numberRows){
  linear <- dataGRCutColumn$standard_id[i] %in% c("38","49","62","63","80","83","84","85","86")
  radial <- dataGRCutColumn$standard_id[i] %in% c("60","81")
  if(radial){
    dataGRCutColumn$value[i] <- dataGRCutColumn$value[i] * 2
  }else if(!linear){
    print(paste("problem for row ",i))
  }
}

### Conversions in the mm.year-1:
## mm d^-1 :
dataGR_Day_DF <- subset(dataGRCutColumn, standard_unit =="mm d^-1")
dataGR_Day_DF$value <- as.numeric(dataGR_Day_DF$value) * 365
dataGR_Day_DF$standard_unit <- "mm yr^-1"

# mm month^-1
dataGR_1_Month_DF <- subset(dataGRCutColumn, standard_unit =="mm month^-1")
dataGR_1_Month_DF$value <- as.numeric(dataGR_1_Month_DF$value) * 12
dataGR_1_Month_DF$standard_unit <- "mm yr^-1"
# !!!! correction needs to be made for Acropora palmata: wrongly converted from cm to mm (from Bak et al., 2009). The 4 values need to 
# be / 10 !!!!!!!! TO CHECK IF THEY HAVE DONE THE CORRECTIONS !!!!!
# --> wrong values : 972 900 744 828 "mm yr^-1"
dataGR_1_Month_DF$value[dataGR_1_Month_DF$specie_name =="Acropora palmata"] <-  dataGR_1_Month_DF$value[dataGR_1_Month_DF$specie_name =="Acropora palmata"] /10

# mm per 3 months^-1
dataGR_3_Month_DF <- subset(dataGRCutColumn, standard_unit =="mm per 3 months^-1")
dataGR_3_Month_DF$value <- as.numeric(dataGR_3_Month_DF$value) * 4
dataGR_3_Month_DF$standard_unit <- "mm yr^-1"

# mm per 4 months^-1
dataGR_4_Month_DF <- subset(dataGRCutColumn, standard_unit =="mm per 4 months^-1")
dataGR_4_Month_DF$value <- as.numeric(dataGR_4_Month_DF$value) * 3
dataGR_4_Month_DF$standard_unit <- "mm yr^-1"

# mm per 6 months^-1
dataGR_6_Month_DF <- subset(dataGRCutColumn, standard_unit =="mm per 6 months^-1")
dataGR_6_Month_DF$value <- as.numeric(dataGR_6_Month_DF$value) * 2
dataGR_6_Month_DF$standard_unit <- "mm yr^-1"

# mm per 7 months^-1
dataGR_7_Month_DF <- subset(dataGRCutColumn, standard_unit =="mm per 7 months^-1")
dataGR_7_Month_DF$value <- as.numeric(dataGR_7_Month_DF$value) * 12 / 7
dataGR_7_Month_DF$standard_unit <- "mm yr^-1"

# mm per 8 months^-1
dataGR_8_Month_DF <- subset(dataGRCutColumn, standard_unit =="mm per 8 months^-1")
dataGR_8_Month_DF$value <- as.numeric(dataGR_8_Month_DF$value) * 3 / 2
dataGR_8_Month_DF$standard_unit <- "mm yr^-1"

# mm yr^-1
dataGR_Year_DF <- subset(dataGRCutColumn, standard_unit =="mm yr^-1")

## binding of all these data.frames and put names in alphabetoc order:
combinedDFs <- rbind(dataGR_Day_DF, dataGR_1_Month_DF, dataGR_3_Month_DF, dataGR_4_Month_DF,dataGR_6_Month_DF,
                     dataGR_7_Month_DF, dataGR_8_Month_DF, dataGR_Year_DF)

combinedDFs$specie_name <- gsub(" ","zzz",combinedDFs$specie_name)
combinedDFs <- combinedDFs[order(combinedDFs$specie_name),]
combinedDFs$specie_name <- gsub("zzz"," ",combinedDFs$specie_name)

# check if number of lines preserved:
length(combinedDFs$specie_name)    #    1552

# check different value_types:
levels(as.factor(combinedDFs$value_type))     #  ""          "mean"      "median"    "raw_value"

combinedDFs$value <- as.numeric(combinedDFs$value)

# creation of final table to determine the mean value for each species:
# mean or median values without sample size or precison type indicated, as well as raw_values are considered as a single observation
# if "the"precision_type" == range and mean value given --> suppose normal distribution and sd <- range / 6 
# median values with sample size are conidered as mean because I supposed normal distribution
species <- levels(as.factor(combinedDFs$specie_name))
numberSpecies <- length(species) # 130
numberLines <- length(combinedDFs$specie_name) # 152
meanSDnTable <- data.frame(species= rep(NA,numberSpecies),mean_GR_mm.yr=rep(NA,numberSpecies),sd=rep(NA,numberSpecies),n=rep(NA,numberSpecies))
for(i in 1:numberSpecies){
  # make a df for a particular species, each row corresponding to a value (either mean, median or raw values)
  temporalDF <- subset(combinedDFs,combinedDFs$specie_name==species[i])
  numberLines <- length(temporalDF$specie_name)
  meanSDnTable_temporal <- data.frame(mean_GR_mm.yr=rep(NA,numberLines),sd=rep(NA,numberLines),n=rep(NA,numberLines))
  for(j in 1:numberLines){
    vt_mean <- temporalDF$value_type[j] == "mean"
    vt_median <- temporalDF$value_type[j] == "median"
    pt_95CI <- temporalDF$precision_type[j] == "95_ci"
    pt_range <- temporalDF$precision_type[j] == "range"
    pt_sd <- temporalDF$precision_type[j] == "standard_deviation"
    pt_se <- temporalDF$precision_type[j] == "standard_error"
    pt_NA <- temporalDF$precision_type[j] == "not_given" | temporalDF$precision_type[j] == ""
    replicates_n <- !is.na(temporalDF$replicates[j])
    
    meanSDnTable_temporal$mean_GR_mm.yr[j] <- temporalDF$value[j]
    
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
species <- tableNamesChecked$species
numberSpecies <- length(species) # 130
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] # only keep the species in vector species
length(speciesChecked$nameSp) # 130
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}
# Check for duplicates:
sum(duplicated(tableNamesChecked$species)) # 3
speciesDuplicated <- tableNamesChecked$species[duplicated(tableNamesChecked$species)]
subset(tableNamesChecked,tableNamesChecked$species %in% speciesDuplicated)

# combine the values:
tableNamesChecked_backup <- tableNamesChecked
numberSpeciesToCombine <- length(speciesDuplicated) # 3
for(i in 1:numberSpeciesToCombine){
  tableTempo <- subset(tableNamesChecked_backup,tableNamesChecked_backup$species == speciesDuplicated[i])
  combination <- meanSDnFund(tableTempo[,2:4])
  combination$species <- speciesDuplicated[i]
  combination <- combination[,c(4,1,2,3)]
  colnames(combination) <- c("species","mean_GR_mm.yr","sd","n")
  # remove the duplicated species
  tableNamesChecked_backup <- subset(tableNamesChecked_backup,tableNamesChecked_backup$species != speciesDuplicated[i])
  # add the species name with updated stats:
  tableNamesChecked_backup <- rbind(tableNamesChecked_backup,combination)
}
sum(duplicated(tableNamesChecked_backup$species))
length(tableNamesChecked_backup$species) # 127

### alphabetic order:
tableNamesChecked_backup$species <- gsub(" ","_",tableNamesChecked_backup$species)
tableNamesChecked_backup <- tableNamesChecked_backup[order(tableNamesChecked_backup$species),]
tableNamesChecked_backup$species <- gsub("_"," ",tableNamesChecked_backup$species)

# creation csv file in the wd dataCompilation
# write.csv(tableNamesChecked_backup, file = paste(wd_Datasets,"growth_rate.csv",sep="/"),row.names=F)




