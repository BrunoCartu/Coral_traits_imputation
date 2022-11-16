# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## Combine all the trait values in one dataset
## produce:
### traits_compilation.csv : 1499 coral species, 22 traits
##  traits_compilation_zooxanthellate.csv : same as above but only "zooxanthellate" or "both" coral species (n = 828)

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_Datasets_trait_data_compilation <- paste(wd,"Traits_extra/trait_data_compilation",sep="/")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### Creation of traitsdataset:
traitData <- data.frame(species = NA)

### Import datasets:
traitData <- read.csv(paste(wd_Datasets,"age_at_maturity.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2]
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"aggressiveness.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"bleaching_response_index_313sp.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"chlorophyll_a_concentration.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"coloniality.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"colony_max_diameter.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"corallite_area.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"dark_respiration.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"egg_diameter.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"fecundity_polyp.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"growth_form.csv",sep="/"),header = T, stringsAsFactors = F)[,1:3],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"growth_rate.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"lipid_content.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"mode_larval_development.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"reduced_scattering_coefficient.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"size_at_maturity.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"skelelal_density.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"symbiodimium_density.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"tissue_thickness.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"sexual_system.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)
traitData <- merge(traitData,read.csv(paste(wd_Datasets,"zooxanthellate.csv",sep="/"),header = T, stringsAsFactors = F)[,1:2],all.x=TRUE,all.y = TRUE)

colnames(traitData)<-c("species","age_at_maturity","aggressiveness","bleaching_response_index","chlorophyll_a_concentration","coloniality","colony_max_diameter",
                        "corallite_area","dark_respiration","egg_diameter","fecundity_polyp","growth_form","growthFormClasses","growth_rate",
                        "lipid_content","mode_larval_development","reduced_scattering_coefficient","size_at_maturity","skeletal_density","symbiodinium_density","tissue_thickness",
                        "sexual_system","zooxanthellae")

sum(duplicated(traitData$species)) # 0
length(traitData$species) # 1499

# write.csv(traitData,paste(wd_Datasets,"traits_compilation.csv",sep="/"),row.names = F)
traits_compilation <- read.csv(paste(wd_Datasets,"traits_compilation.csv",sep="/"),header = T, stringsAsFactors = F)

### Creation dataset with only zooxanthellate species:
levels(as.factor(traits_compilation$zooxanthellae))
sum(is.na(traits_compilation$zooxanthellae)) # 22  --> check their family and order with WoRMS

# check the species with NA for zooxanthellae --> some might not be coral species
traitData_zoo_NA <- traits_compilation[is.na(traits_compilation$zooxanthellae),]

### send species to WoRMS
# library(stringr)
# speciesTSend <- str_split_fixed(traitData_zoo_NA$species," ",2)
# speciesToSendToWoRMS <- data.frame(
#   Genus = speciesTSend[,1],
#   Species = speciesTSend[,2]
# )
# 
# write.csv(speciesToSendToWoRMS,paste(wd_Datasets_trait_data_compilation,"speciesToSendToWoRMS.csv",sep="/"),row.names = F)

speciestosendtoworms_matched <- read.csv(paste(wd_Datasets_trait_data_compilation,"speciestosendtoworms_matched.csv",sep="/"),header=T,stringsAsFactors = F)
levels(as.factor(speciestosendtoworms_matched$Order))
# "Alcyonacea"    "Anthoathecata" "Helioporacea"  "Scleractinia" 

# species whose original name was not correct:
speciestosendtoworms_matched$ScientificName[speciestosendtoworms_matched$ScientificName != speciestosendtoworms_matched$ScientificName_accepted]
#  "Turbinaria immersa"  --> nomen dubium

# scleractinian species 
speciesScleractian <- speciestosendtoworms_matched$ScientificName[speciestosendtoworms_matched$Order == "Scleractinia"]
speciesScleractian

traits_compilation[traits_compilation$species %in% speciesScleractian,]
# only 2 species with only aggressiveness as trait value --> we don't care.

### selection zooxanthellae species
traitData_zoo <- subset(traits_compilation,traits_compilation$zooxanthellae %in% c("both","zooxanthellate"))
length(traitData_zoo$species) # 828

# write.csv(traitData_zoo,paste(wd_Datasets,"traits_compilation_zooxanthellate.csv",sep="/"),row.names = F)

# 
traitData_azoo <- subset(traits_compilation,! traits_compilation$zooxanthellae %in% c("both","zooxanthellate"))
length(traitData_azoo$species) # 671

traitData_azoo_aggress <- traitData_azoo[,c("species","aggressiveness")]
traitData_azoo_aggress <- traitData_azoo_aggress[!is.na(traitData_azoo_aggress$aggressiveness),]

# CHECK FOR NUMBER OF SPECIES FOR EACH TRAIT:

### Number of species available for each trait in traits_compilation_zooxanthellate.csv:
traits_compilation_zooxanthellate <- read.csv(paste(wd_Datasets,"traits_compilation_zooxanthellate.csv",sep="/"),header = T, stringsAsFactors = F)
length(traits_compilation_zooxanthellate$species) # 828
str(traits_compilation_zooxanthellate)

traits_compilation_zooxanthellateBis <- traits_compilation_zooxanthellate[!is.na(traits_compilation_zooxanthellate$aggressiveness),]

traits_compilation_zooxanthellateBis[,c("species","aggressiveness")][na.omit(traits_compilation_zooxanthellateBis$aggressiveness) == max(traits_compilation_zooxanthellateBis$aggressiveness),]

length(traits_compilation_zooxanthellate)
numberTraits <- length(colnames(traits_compilation_zooxanthellate))
for(i in 1:numberTraits){
  print(paste(colnames(traits_compilation_zooxanthellate)[i],"",length(na.omit(traits_compilation_zooxanthellate[,i]))))
}









