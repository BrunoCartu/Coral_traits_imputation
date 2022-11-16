# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## To create 
### functionalTraitDF_model.csv: the dataset for the 798 coral species and six functional groups of algae used in the model to input trait values; exported in Simulations/Datasets_for_model and pasted in coralreef2/coralreef2/data

## From 
### functionTrait_DS_Imputed.csv                       : the coral trait dataset of imputed traits for 798 species and 11 traits; created in imputation_traits_missForest.R
### bleaching_probability.csv                          : the species-specific index of bleaching susceptibility (Appendix S4: 1); created in BRI_responseTraits_associations_Maveraging.R
### Initial_benthic_composition_Martinique_FB_PB_IR.csv: the initial coral community composition in the thre Caribbean sites; created in Empirical_datasets_Martinique_calibration > Rscripts > prep_for_model.R. NOTE: it is the same composition as the one in the three 1st communitis in coralreef2 > coralreef2 > data > Initial_benthic_composition.csv
### growthRate_randomRadiusConversion.csv              : converted from growthRate_randomRadiusConversion.xls which I manually created to convert real grouth rates into radius size (Appendix S2: 6.1)

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")

wd_Empirical_datasets_Martinique_calibration <- paste(gsub("/Traits_and_imputation","",wd),"Empirical_datasets_Martinique_calibration",sep="/")
wd_Empirical_datasets_Martinique_calibration_Datasets_for_model <- paste(wd_Empirical_datasets_Martinique_calibration,"Datasets_for_model",sep="/")

wd_coralreef2 <-  paste(gsub("/Traits_and_imputation","",wd),"coralreef2/coralreef2",sep="/")
wd_coralreef2_data <- paste(wd_coralreef2,"data",sep="/")

wd_Simulations <-  paste(gsub("/Traits_and_imputation","",wd),"Simluations",sep="/")
wd_Simulations_Datasets_for_model <- paste(wd_Simulations,"Datasets_for_model",sep="/")

source(paste(wd_Rscripts,"functions.R",sep="/"))

# Import the imputed coral traits database: -----
traitsCoralImputed <- read.csv(paste(wd_Datasets,"functionTrait_DS_Imputed.csv",sep="/"),header=T,stringsAsFactors = F)
bleaching_probability <- read.csv(paste(wd_Datasets,"bleaching_probability.csv",sep="/"),header = T, stringsAsFactors = F)

functionTraitDF <- traitsCoralImputed
functionTraitDF$bleaching_probability <- bleaching_probability$bleaching_probability

length(functionTraitDF$species) # 798

# Create functionTraitDF_model.csv --------
## = the dataframe that is going to contain functional characteristics of all the different agents 
functionTraitDF_model <- functionTraitDF

functionTraitDF_model$size_maturity <- 35  # not used anymore 
functionTraitDF_model$age_maturity <- 3 # value in Guest et al., 2014 for Acropora millepora --> smaller than the 3 values found in the coraltrait.org data base
functionTraitDF_model$red <- NA
functionTraitDF_model$green <- NA
functionTraitDF_model$blue <- NA

colnames(functionTraitDF_model)
str(functionTraitDF_model)
head(functionTraitDF_model)
colnamesInRightOrderForModel <- c("species","aggressiveness","coloniality","colony_max_diameter","corallite_area","egg_diameter","fecundity_polyp",
                                  "growthFormClasses","growth_rate","mode_larval_development","reduced_scattering_coefficient","size_maturity",
                                  "age_maturity","red","green","blue","bleaching_probability","sexual_system")

functionTraitDF_model <- functionTraitDF_model[,colnamesInRightOrderForModel]
head(functionTraitDF_model)

## add column for correction coefficient for probability (Cc) of a polyp to be fecund (in function of colony size) --> only for small coral species (Appendix S2: 3.1.1.b)
# for most species is it 0
# only for the smallest species it is different from 0 
functionTraitDF_model$correction_coeff_polypFecundity <- 0
# the maximum max_planar_surface_area of the small species
Smax <- exp((logitTrans(dataProportion = 0.9,epsilon = 0,dataInPercentage = F) - 8.626)/1.682) * 10000
Smax # 218.8176 cm2
# determine biggest maximum diameter:
Dmax <- sqrt(4*Smax/pi)
Dmax #  16.69153 cm
smallSp <- functionTraitDF_model$species[functionTraitDF_model$colony_max_diameter < Dmax] 
length(smallSp) # 50
# for these species the Cc is expressed as:
for(i in 1 : length(smallSp)){
  Smax_ln <- log(pi/4*(functionTraitDF_model$colony_max_diameter[functionTraitDF_model$species == smallSp[i]]/100)^2)
  functionTraitDF_model$correction_coeff_polypFecundity[functionTraitDF_model$species == smallSp[i]] <- (logitTrans(0.9,dataInPercentage = F)-8.626)/1.682 - Smax_ln
}

namesColumns <- colnames(functionTraitDF_model)
functionTraitDF_model$substrate_category <- "Coral"
functionTraitDF_model$substrate_subcategory <- "LiveCoral"

functionTraitDF_model <- functionTraitDF_model[c("substrate_category","substrate_subcategory",namesColumns)]
functionTraitDF_model$species <- gsub(" ","_",functionTraitDF_model$species)

head(functionTraitDF_model)

colourBorrenGround <- c(160,160,160)
colourMA <- c(20,160,20)
colourAMA <- c(14,125,14)
colourHalimeda <- c(35,243,32)
colourTurf <- c(178,255,102)
colourCCA <- c(229,192,69)
colourACA <- c(243,78,144)
colourBleachedCoral <- c(245,245,245)
colourDeadCoral <- c(255,255,255)

listColourTaken <- list(colourBorrenGround,colourMA,colourAMA,colourHalimeda,colourTurf,colourCCA,colourACA,colourBleachedCoral,colourDeadCoral)
sum(duplicated(listColourTaken)) # 0

## Give species Martinique specific color for display
## combine all the species in the 3 sites considered
martiniqueSpecies <- read.csv(paste(wd_Empirical_datasets_Martinique_calibration_Datasets_for_model,"Initial_benthic_composition_Martinique_FB_PB_IR.csv",sep="/"),header = T, stringsAsFactors = F)
martiniqueSpecies <- martiniqueSpecies[c(1,5,9),3:14]
martiniqueSpecies <- levels(as.factor(c(as.character(martiniqueSpecies[1,]),as.character(martiniqueSpecies[2,]),as.character(martiniqueSpecies[3,]))))

length(martiniqueSpecies) # 19
martiniqueSpecies <- data.frame(
  species = martiniqueSpecies,
  R = c(10,10,10,10,10,102,76,178,153,255,153,255,255,204,255,204,100,100,180),
  G = c(102,153,51,102,10,102,10,102,10,102,10,102,102,102,178,204,100,150,105),
  B = c(102,153,102,204,204,255,153,255,10,153,255,76,78,10,102,10,160,165,36),
  stringsAsFactors = F
)

## to check if color to add is already present in the table
# colSp <- c(180,105,36)
# for(i in 1:length(martiniqueSpecies$R)){
#   if(sum(colSp == martiniqueSpecies[i,]) == 3){
#     print("already present")
#   }
# }

j <- length(listColourTaken) # 9
for(i in 1:length(martiniqueSpecies$species)){
  j <- j + 1
  listColourTaken[[j]] <- as.numeric(martiniqueSpecies[i,2:4])
}
length(listColourTaken) # 28
sum(duplicated(listColourTaken)) # 0
length(unique(listColourTaken)) # 28

# attribute the R, G, B colors values to the caribben species in functionTraitDF_model and te non-caribbean species in a randomly
numberSpecies <- length(functionTraitDF_model$species) # 798
for(i in 1:numberSpecies){
  if(functionTraitDF_model$species[i] %in% martiniqueSpecies$species){
    functionTraitDF_model$red[i] <- martiniqueSpecies$R[martiniqueSpecies$species == functionTraitDF_model$species[i]]
    functionTraitDF_model$green[i] <- martiniqueSpecies$G[martiniqueSpecies$species == functionTraitDF_model$species[i]]
    functionTraitDF_model$blue[i] <- martiniqueSpecies$B[martiniqueSpecies$species == functionTraitDF_model$species[i]]
  }else{
    validColour <- F     # for no Caribbean coral sp, determine color randomly
    while(!validColour){
      colourCoral <- c(sample(seq(10,255),1),sample(seq(10,255),1),sample(seq(10,255),1))
      allGood <- T
      for(k in 1: length(listColourTaken)){
        if(identical(colourCoral,listColourTaken[[k]])){
          allGood <- F
          break
        }
      }
      if(allGood){
        functionTraitDF_model[c("red","green","blue")][i,] <- colourCoral
        validColour <- T
        j <- j + 1
        listColourTaken[[j]] <- colourCoral
      }
    }
  }
  print(i)
}

# View(functionTraitDF_model)
functionTraitDF_model[functionTraitDF_model$species == martiniqueSpecies$species[1] ,]

max(functionTraitDF_model$growth_rate) # 179.3333 mm/y-1 (diameter/linear GR)
min(functionTraitDF_model$growth_rate) # 0.8 mm/y-1

### add barrenground and algae and bleach and dead coral
# the growth rates are entered in mm.yr-1
coralDead <- c("Coral","DeadCoral",NA,0,NA,0,0,0,0,NA,0,NA,0,0,0,colourDeadCoral,0,NA,0)
names(coralDead) <- colnames(functionTraitDF_model)
functionTraitDF_model <- rbind(coralDead,functionTraitDF_model)
coralBleach <- c("Coral","BleachedCoral",NA,0,NA,0,0,0,0,NA,0,NA,0,0,0,colourBleachedCoral,0,NA,0)
functionTraitDF_model <- rbind(coralBleach,functionTraitDF_model)
CCA <- c("Algae","CCA","CCA",0,"None",100000,0,0,0,"None",24,"None",0,0,0,colourCCA,0,"None",0)    # 0.01 to 0.07 mm.day-1 --> 0.05mm.d-1 --> 18.25mm.y-1 marginal --> 36.5mm.y-1 diameter (?)  "Colonization and growth of crustose coralline algae (Corallinales, Rhodophyta) on the Rocas Atoll
# (Adey and Vassar 1975): 0.9 to 2.3 mm/month  --> 10.8 to 27.6 mm.yr-1 
# (Boas et al., 2005): 0.01 to 0.05 mm/day     --> 3.65 to 18.25 mm.yr-1
# (Steneck 1986) crustise coralline algae are, In any given habitat, corallines are usually the least productive group of algae (10, 69, 73, 77, 80); they grow slowly (24, 49).
# (Martone 2010) Marginal extension rates for encrusting corallines may exceed 2– 3mm Æ month)1 (Gardiner 1931, Adey and Vassar 1975, Steneck 1985) but are typically <1 mm Æ month)1 (Littler 1972, Steneck and Adey 1976, Johansen 1981, Matsuda 1989).--> < 12 mm.yr-1
# ***** final call: 12  mm.yr-1 radial extension --> 24 mm.yr-1 diameter extension *******
functionTraitDF_model <- rbind(CCA,functionTraitDF_model)
ACA <- c("Algae","ACA","ACA",0,"None",100000,0,0,0,"None",42,"None",0,0,0,colourACA,0,"None",0)    # 
# Bossiella (?)
# Corallina (?)
# Calliartlzron tuberculosum (Johansen and Austin 1970)(not a coral reef): 1.7 mm.month-1 --> 20.4 mm.yr-1
# (Martone 2010) Erect articulated corallines may grow much faster, up to 5 mm Æ month)1 (Haas et al. 1935, Smith 1972), but on average grow 1.5–2 mm Æ month)1 (Johansen and Austin 1970, Pearse 1972, Andrake and Johansen 1980, Blake and Maggs 2003) --> 18 to 24 mm.yr-1
# ***** final call: 21 mm.yr-1 ***** --> 42 mm.yr-1 diameter extention ***********
functionTraitDF_model <- rbind(ACA,functionTraitDF_model)
Turf <- c("Algae","Turf","Turf",0,"None",100000,0,0,0,"None",500,"None",0,0,0,colourTurf,0,"None",0)               
# Turf (mainly Polysiphonia and Ceramium) (Golber et al., 2006) (values taken with cage --> no grazing, measure of hight): 0.0492 to 0.7692 mm.d-1 --> 17.958 to 280.758 mm.yr-1
# Nostococaceae (?)
# Lyngbya (Lyngbya majuscula) (Ahern et al., 2007): 0.0178 mm.d-1 (control) to 0.5758 mm.d-1 (environment rich) --> 6.497 mm.yr-1 to 210.167 mm.yr-1
# Oscillatoria (?)
# Cladophora (?)
# ***** final call 250 mm.yr-1 ***** --> 500 mm.yr-1 diameter extention ***********
functionTraitDF_model <- rbind(Turf,functionTraitDF_model)
Halimeda <- c("Algae","Halimeda","Halimeda",0,"None",100000,0,0,0,"None",300,"None",0,0,0,colourHalimeda,0,"None",0)    
# (Drew 1983), gr = sqrt(0.75cm2)*0.16.d-1 = 0.1386 cm.d-1 = 505.89 mm.yr-1
# (Vroom et al., 2003): segment lenght: 1.95 to 10.02 mm --> 5.985 mm ; number segments produced.d-1: 0.2925 to 2.7515 -->  638.9734 to 6010.721 mm.yr-1
# ***** final call:  300 mm.yr-1 diameter extention (same value for MA, AMA, Halimeda)  ***********
functionTraitDF_model <- rbind(Halimeda,functionTraitDF_model)
AMA <- c("Algae","AMA","AMA",0,"None",100000,0,0,0,"None",300,"None",0,0,0,colourAMA,0,"None",0)
# Chlorodesmis (?)
# Lobophora (?)
# ***** final call:  300 mm.yr-1 diameter extention (same value for MA, AMA, Halimeda)  ***********
functionTraitDF_model <- rbind(AMA,functionTraitDF_model)
MA <- c("Algae","Macroalgae","Macroalgae",0,"None",100000,0,0,0,"None",300,"None",0,0,0,colourMA,0,"None",0)  
# Sargassum polycystum (May-Lin and Ching-Lee 2013) (seasonal variation --> only take mean max value):  -3.75 to 2.54 mm.d-1 -->  927.1  mm.yr-1
# Sargassum binderi (May-Lin and Ching-Lee 2013) (seasonal variation --> only take mean max value):     -3.23 to 1.89 mm.d-1 -->  689.85 mm.yr-1
# Sargassum siliquosum (May-Lin and Ching-Lee 2013) (seasonal variation --> only take mean max value):  -5.22 to 4.08 mm.d-1 --> 1489.2  mm.yr-1
# Sargassum pteropleuron (Prince and O'Neal 1979) (seasonaly variation, only take mean max value):      - 10...to 7.1 mm.d-1 --> 2591.5 mm.yr-1
# Hydroclathrus (?) 
# Dictyota dichotoma (Kuhlenkamp et al., 2001) (colder water 12.5 to 17 C): 0.0286 to 0.1571 mm.d-1 --> 10.43 to 57.34 mm.yr-1
# ***** final call:  300 mm.yr-1 diameter extention (same value for MA, AMA, Halimeda)  ***********
functionTraitDF_model <- rbind(MA,functionTraitDF_model)
barrenGround <- c("BarrenGround","BarrenGround","BarrenGround",0,"None",0,0,0,0,"None",0,"None",0,0,0,colourBorrenGround,0,"None",0)
functionTraitDF_model <- rbind(barrenGround,functionTraitDF_model)

head(functionTraitDF_model,12)
str(functionTraitDF_model)
functionTraitDF_model$growth_rate <- as.numeric(functionTraitDF_model$growth_rate)
### Necessity to convert the growth rate (linear/diamter) mm.y-1 in random radius length (because the model is discrete in space and time):
# 1st convert diameter GR in radial GR:
functionTraitDF_model$growth_rate <- functionTraitDF_model$growth_rate / 2
max(functionTraitDF_model$growth_rate) # 250 mm/y-1

# 2nd: conversion in cm/y-1:
functionTraitDF_model$growth_rate <- functionTraitDF_model$growth_rate / 10
max(functionTraitDF_model$growth_rate) # 25 cm/y-1

# 3rd, import the csv conversion file: growthRate_randomRadiusConversion.csv
converionGR <- read.csv(paste(wd_Datasets,"growthRate_randomRadiusConversion.csv",sep="/"),header=T)
max(converionGR$GR) # 27.45 cm/y
#View(converionGR)
convertedGR <- rep(NA,numberSpecies)

# Conversion:
numberSpecies <- length(functionTraitDF_model$species) # 807
numberGR <- length(converionGR$RF) # 550
for(i in 1:numberSpecies){     # need to find the smallest difference between the value in converionGR$GR and the actual value 
  RF <- converionGR$RF[1]
  ABS <- abs(converionGR$GR[1] - functionTraitDF_model$growth_rate[i])
  for(j in 1:numberGR){
    ABS_j <- abs(converionGR$GR[j] - functionTraitDF_model$growth_rate[i])
    if(ABS_j < ABS){
      ABS <- ABS_j
      RF <- converionGR$RF[j]
    }else{  # in this case ABS_j > ABS and we are going through the grouth rate values that are > functionTraitDF_model$growth_rate[i] --> no need to go further
      if(converionGR$GR[j] - functionTraitDF_model$growth_rate[i] > 0){
        break()
      }
    }
  }
  convertedGR[i] <- RF
}

plot(functionTraitDF_model$growth_rate~convertedGR) # all good

functionTraitDF_model$growth_rate <- convertedGR

colnames(functionTraitDF_model)[colnames(functionTraitDF_model)=="growthFormClasses"] <- "growth_form"

head(functionTraitDF_model)

### check for size at maturity (cm2) vs. maximum colony diameter --> not relevent anymore
max_planarArea <- pi * as.numeric(functionTraitDF_model$colony_max_diameter)^2
functionTraitDF_model$size_maturity <- as.numeric(functionTraitDF_model$size_maturity)
functionTraitDF_model[functionTraitDF_model$size_maturity > max_planarArea,]

hist(log10(max_planarArea),xlim=c(0,12),breaks = seq(0,12,0.5))
sizeMaturity_log10 <- log10(functionTraitDF_model$size_maturity[20])
segments(x0=sizeMaturity_log10,y0=0,x1=sizeMaturity_log10,y1=200,col="red",lwd=3)

# write.csv(functionTraitDF_model,paste(wd_Simulations_Datasets_for_model,"functionalTraitDF_model.csv",sep="/"),row.names = F)
# write.csv(functionTraitDF_model,paste(wd_coralreef2_data,"functionalTraitDF_model.csv",sep="/"),row.names = F) # in the folder where data is inputed in the model during simulations

DF <- read.csv(paste(wd_Datasets_inputs,"functionalTraitDF_model.csv",sep="/"),header=T)
