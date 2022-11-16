# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## to import and curate species trait data from coraltraits.org and from the following references: Reynaud-Vaganay et al. 2001, Hennige et al., 2010, Cooper et al., 2011, Ulstrup et al., 2011 
## To create:
### produce Dark_respiration_original.csv
### dark_respiration.csv

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_Dark_respiration <- paste(wd,"/Traits_extra/dark_respiration",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

### import dataset from coraltrait.org: (April 3 2017)
# Dark_respiration_original <- read.csv("https://coraltraits.org/traits/133.csv", as.is = TRUE)
# write.csv(Dark_respiration_original,paste(wd_Datasets_original,"Dark_respiration_original.csv",sep="/"),row.names = F)

### import original dataset: ----
Dark_respiration_original <- read.csv(paste(wd_Datasets_original,"Dark_respiration_original.csv",sep="/"),header=T,stringsAsFactors = F)

levels(as.factor(Dark_respiration_original$trait_name))

Dark_respiration_1 <- subset(Dark_respiration_original,Dark_respiration_original$trait_name =="Dark respiration")

colnames(Dark_respiration_1)

colnamesToKeep <- c("specie_name","trait_name","standard_unit","methodology_name","value","value_type","precision","precision_type","replicates","notes")

Dark_respiration_2 <- Dark_respiration_1[,colnamesToKeep]

length(levels(as.factor(Dark_respiration_2$specie_name))) # 11

levels(as.factor(Dark_respiration_2$value_type))
#  "raw_value"
levels(as.factor(Dark_respiration_2$precision_type))
# NA
levels(as.factor(Dark_respiration_2$replicates))
# NA
levels(as.factor(Dark_respiration_2$standard_unit))
# "µmol O2 cm^-2 h^-1"
levels(as.factor(Dark_respiration_2$notes))

# creation of final table to determine the mean value for each species (there are only raw values here:
Dark_respiration_2$value <- as.numeric(Dark_respiration_2$value)

meanDarkRespiration <- tapply(Dark_respiration_2$value,Dark_respiration_2$specie_name,mean)
sd <- tapply(Dark_respiration_2$value,Dark_respiration_2$specie_name,sd)
lengthDR <- tapply(Dark_respiration_2$value,Dark_respiration_2$specie_name,function(x){length(x)})
SE <- sd/sqrt(lengthDR)

Dark_respiration_3 <- data.frame(
  species = rownames(meanDarkRespiration),
  dark_respiration = meanDarkRespiration,
  sd = sd,
  n = lengthDR,
  stringsAsFactors = F
)

### additional references (check excel sheet "Data from extra references.xlsx":
# Reynaud-Vaganay et al. 2001 -----
# measures are taken from text and figure 2
# "Stylophora pistillata"
# value = 0.15 "µmol O2 cm^-2 h^-1"
mean1 <- 0.14
mean2 <- 0.16
SE1 <- 0.056/2
SE2 <- 0.063/2
n <- 20
data <- data.frame(
  species = "Stylophora pistillata",
  dark_respiration = c(mean1,mean2),
  sd = round(c(SE1*sqrt(n),SE2*sqrt(n)),3),
  n = n
)
finalData <- meanSDnFund(data[,-1])
finalData$species <- data$species[1]
finalData <- finalData[,c(4,1,2,3)]
colnames(finalData)[2] <- "dark_respiration"
Dark_respiration_3 <- rbind(Dark_respiration_3,finalData)

# Hennige et al., 2010 (averaged over the different sites) -----
# measures taken from figure 4
# "Porites lutea"
means <- c(20.87,19.57,23.04)
n <- 3
sd <- c(rep(3.043,3))/2*sqrt(n)
data <- data.frame(
  species = "Porites lutea",
  dark_respiration = means,
  sd = round(sd,3),
  n = n
)
# conversion from "µmol O2 cm^-2 d^-1" to "µmol O2 cm^-2 h^-1":
data$dark_respiration <- data$dark_respiration/24
data$sd <- data$sd/24
finalData <- meanSDnFund(data[,-1])
finalData$species <- data$species[1]
finalData <- finalData[,c(4,1,2,3)]
colnames(finalData)[2] <- "dark_respiration"
Dark_respiration_3 <- rbind(Dark_respiration_3,finalData)

# "Favites abdita"
means <- 40.65
n <- 3
sd <- 3.913/2*sqrt(n)
data <- data.frame(
  species = "Favites abdita",
  dark_respiration = means,
  sd = round(sd,3),
  n = n
)
# conversion from "µmol O2 cm^-2 d^-1" to "µmol O2 cm^-2 h^-1":
data$dark_respiration <- data$dark_respiration/24
data$sd <- data$sd/24
Dark_respiration_3 <- rbind(Dark_respiration_3,data)

# "Porites lobata"
means <- c(15.43,17.39,34.35)
n <- 3
sd <- c(rep(3.043,3))/2*sqrt(n)
data <- data.frame(
  species = "Porites lobata",
  dark_respiration = means,
  sd = round(sd,3),
  n = n
)
# conversion from "µmol O2 cm^-2 d^-1" to "µmol O2 cm^-2 h^-1":
data$dark_respiration <- data$dark_respiration/24
data$sd <- data$sd/24
finalData <- meanSDnFund(data[,-1])
finalData$species <- data$species[1]
finalData <- finalData[,c(4,1,2,3)]
colnames(finalData)[2] <- "dark_respiration"
Dark_respiration_3 <- rbind(Dark_respiration_3,finalData)

# "Goniastrea aspera"
means <- c(63.48,30.87,39.57)
n <- 3
sd <- c(3.478,3.043,2.826)/2*sqrt(n)
data <- data.frame(
  species = "Goniastrea aspera",
  dark_respiration = means,
  sd = round(sd,3),
  n = n
)

# conversion from "µmol O2 cm^-2 d^-1" to "µmol O2 cm^-2 h^-1":
data$dark_respiration <- data$dark_respiration/24
data$sd <- data$sd/24
finalData <- meanSDnFund(data[,-1])
finalData$species <- data$species[1]
finalData <- finalData[,c(4,1,2,3)]
colnames(finalData)[2] <- "dark_respiration"
Dark_respiration_3 <- rbind(Dark_respiration_3,finalData)

# "Porites cylindrica"  # ALREADY IN TABLE SO I DON T ADD IT
# (20+13)/2=16.5
# value = 0.6875 "µmol O2 cm^-2 h^-1"
# Dark_respiration_3 <- rbind(Dark_respiration_3,c("Porites cylindrica",0.6875,NA))

# "Acropora formosa"
means <- c(36.96,22.39)
n <- 3
sd <- c(10,12.609)/2*sqrt(n)
data <- data.frame(
  species = "Acropora formosa",
  dark_respiration = means,
  sd = round(sd,3),
  n = n
)
# conversion from "µmol O2 cm^-2 d^-1" to "µmol O2 cm^-2 h^-1":
data$dark_respiration <- data$dark_respiration/24
data$sd <- data$sd/24
finalData <- meanSDnFund(data[,-1])
finalData$species <- data$species[1]
finalData <- finalData[,c(4,1,2,3)]
colnames(finalData)[2] <- "dark_respiration"
Dark_respiration_3 <- rbind(Dark_respiration_3,finalData)

# Cooper et al., 2011 ----
# I estimated the mean and standard deviation values based from the box plots assuming a normal distribution of the data:
# I only considered the values from 10  to 45 meters, after values decrease a lot
# from figures 3 and 4
# "Pachyseris speciosa"
means <- c(1.06,0.88,1.06)
n <- 5
sd <- c(0.143,0.142,0.198)/2*sqrt(n)
data <- data.frame(
  species = "Pachyseris speciosa",
  dark_respiration = means,
  sd = round(sd,3),
  n = n
)
finalData <- meanSDnFund(data[,-1])
finalData$species <- data$species[1]
finalData <- finalData[,c(4,1,2,3)]
colnames(finalData)[2] <- "dark_respiration"
Dark_respiration_3 <- rbind(Dark_respiration_3,finalData)

# "Seriatopora hystrix"
means <- c(0.76,0.54,0.76,0.54,0.59)
n <- 5
sd <- c(0.10,0.03,0.16,0.06,0.13)/2*sqrt(n)
data <- data.frame(
  species = "Seriatopora hystrix",
  dark_respiration = means,
  sd = round(sd,3),
  n = n
)
finalData <- meanSDnFund(data[,-1])
finalData$species <- data$species[1]
finalData <- finalData[,c(4,1,2,3)]
colnames(finalData)[2] <- "dark_respiration"
Dark_respiration_3 <- rbind(Dark_respiration_3,finalData)

# Ulstrup et al., 2011 -----
# I averaged the values over seasons (winter, summer) and geographic regions (Center vs Northern GBR)
# values from table
# "Pocillopora damicornis" # ALREADY IN THE TABLE SO I DON T ADD IT UP
# (0.36+0.28+0.64+0.43)/4 = 0.4275 "nmol O2 cm–2 s–1"
# value = 1.539 "µmol O2 cm^-2 h^-1"

# "Turbinaria reniformis"
means <- c(0.91,0.51,0.22,0.37)
n <- 4
sd <- c(0.12,0.11,0.05,0.06)/2*sqrt(n)
data <- data.frame(
  species = "Turbinaria reniformis",
  dark_respiration = means,
  sd = round(sd,3),
  n = n
)
# conversion from "nmol O2 cm–2 s–1" to "µmol O2 cm^-2 h^-1":
data$dark_respiration <- data$dark_respiration /1000 * 3600
data$sd <- data$sd /1000 * 3600
finalData <- meanSDnFund(data[,-1])
finalData$species <- data$species[1]
finalData <- finalData[,c(4,1,2,3)]
colnames(finalData)[2] <- "dark_respiration"
Dark_respiration_3 <- rbind(Dark_respiration_3,finalData)

### put species is order -----
Dark_respiration_3$species <- gsub(" ","zzz",Dark_respiration_3$species)
Dark_respiration_3 <- Dark_respiration_3[order(Dark_respiration_3$species),]
Dark_respiration_3$species <- gsub("zzz"," ",Dark_respiration_3$species)

### check for duplicates: 
sum(duplicated(Dark_respiration_3$species)) # 0 all good

tableNamesChecked <- Dark_respiration_3

### check for names:
speciesToCheck <- checkSpeciesNames(tableNamesChecked$species, wd_Datasets = wd_Datasets)
speciesToCheck

### import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp) # 

### correct names in species_colonialityFinal_cut:
species <- tableNamesChecked$species
numberSpecies <- length(species) # 20
speciesChecked <- speciesNameChecked[speciesNameChecked$nameSp %in% species,] # only keep the species in vector species
length(speciesChecked$nameSp) # 20
for(i in 1:numberSpecies){
  tableNamesChecked$species[tableNamesChecked$species == speciesChecked$nameSp[i]] <-  speciesChecked$nameSp_checked[i]
}

### check if there are duplicated names:
sum(duplicated(tableNamesChecked$species)) # 0
length(tableNamesChecked$species) # 20

# write.csv(tableNamesChecked,paste(wd_Datasets,"dark_respiration.csv",sep="/"),row.names=F)


