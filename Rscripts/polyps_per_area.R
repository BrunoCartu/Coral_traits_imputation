# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## import and curate species trait data from coraltraits.org and check (Hall and Hughes 1996) for sample sizes
## produce PolypPerarea_original.csv and polyps_per_area.csv 

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### Import dataset from coraltraits.org: (27 November)
# polypPerarea_original <- read.csv("https://coraltraits.org/traits/148.csv", as.is = TRUE)
# write.csv(polypPerarea_original,paste(wd_Datasets_original,"PolypPerarea_original.csv",sep="/"),row.names = F)

### Import dataset from folder:
polypPerarea_original <- read.csv(paste(wd_Datasets_original,"polypPerarea_original.csv",sep="/"),header=T,stringsAsFactors = F)
head(polypPerarea_original)
nameColumn <- c("specie_name","trait_name","value","standard_unit","value_type","precision","precision_type","replicates")
polypPerarea_cut <- polypPerarea_original[,nameColumn]
polypPerarea_cut <- polypPerarea_cut[polypPerarea_cut$trait_name == "Polyps per area",]

levels(as.factor(polypPerarea_cut$standard_unit)) # "units cm^-2"
levels(as.factor(polypPerarea_cut$value_type)) #  ""          "mean"      "raw_value"
levels(as.factor(polypPerarea_cut$precision_type)) #  ""  "standard_deviation" "standard_error"  
levels(as.factor(polypPerarea_cut$replicates)) #  ""  "standard_deviation" "standard_error"  

tableNamesChecked <- polypPerarea_cut

### check for names:
speciesToCheck <- checkSpeciesNames(tableNamesChecked$specie_name, wd_Datasets = wd_Datasets) #  "All good"

### conversion in numerical values 
polypPerarea_cut$value <- as.numeric(polypPerarea_cut$value)
polypPerarea_cut$precision <- as.numeric(polypPerarea_cut$precision)

# Creation of final table to determine the mean value for each species:
# mean or median values without sample size or precison type indicated, as well as raw_values are considered as a single observation
# if "the"precision_type" == range and mean value given --> suppose normal distribution and sd <- range / 6 
# median values with sample size are conidered as mean because I supposed normal distribution
species <- levels(as.factor(polypPerarea_cut$specie_name))
numberSpecies <- length(species) # 26
numberLines <- length(polypPerarea_cut$specie_name) # 56
meanSDnTable <- data.frame(species= rep(NA,numberSpecies),mean_NumbPolyp_per_cm2=rep(NA,numberSpecies),sd=rep(NA,numberSpecies),n=rep(NA,numberSpecies))
for(i in 1:numberSpecies){
  # make a df for a particular species, each row corresponding to a value (either mean, median or raw values)
  temporalDF <- subset(polypPerarea_cut,polypPerarea_cut$specie_name==species[i])
  numberLines <- length(temporalDF$specie_name)
  meanSDnTable_temporal <- data.frame(mean_GR_mm.yr=rep(NA,numberLines),sd=rep(NA,numberLines),n=rep(NA,numberLines))
  for(j in 1:numberLines){
    vt_mean <- temporalDF$value_type[j] == "mean"
    vt_rv <- temporalDF$value_type[j] == "raw_value" || temporalDF$value_type[j] == ""
    pt_sd <- temporalDF$precision_type[j] == "standard_deviation"
    pt_se <- temporalDF$precision_type[j] == "standard_error"
    pt_NA <- temporalDF$precision_type[j] == ""
    replicates_n <- !is.na(temporalDF$replicates[j])
    
    meanSDnTable_temporal$mean_GR_mm.yr[j] <- temporalDF$value[j]
    
    if(replicates_n & !pt_NA & vt_mean){
      meanSDnTable_temporal$n[j] <- temporalDF$replicates[j]
      # if median but n and precision_type given, then median is considered as the mean, because I assumed a normal distribution
      if(pt_sd){
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

### alphabetic order:
meanSDnTable$species <- gsub(" ","_",meanSDnTable$species)
meanSDnTable <- meanSDnTable[order(meanSDnTable$species),]
meanSDnTable$species <- gsub("_"," ",meanSDnTable$species)

# creation csv file in the wd dataCompilation
# write.csv(meanSDnTable, file = paste(wd_Datasets,"polyps_per_area.csv",sep="/"),row.names=F)
