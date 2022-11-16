# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## Check and eventually correct the species names for each trait and phylogeny considered using WoRMS (http://www.marinespecies.org/) and eventually coraltraits.org
## Produce speciesNameChecked.csv

## I updated (added more species names) speciesNameChecked.csv as I was currating traits datasets.
## In each trait R script (e.g., coloniality.R, growth_form.R), I checked the species name using speciesNameChecked.csv. If species were 
## not present in speciesNameChecked.csv, I would create a speciesToSendWoRMS.csv (I then removed the code) which I then imported in WoRMS.
## The site would return speciestosendtoworms_matched.xlsx, which I would convert to speciestosendtoworms_matched.csv (I would eventually create
## a namesCheck.csv). I would then update speciesNameChecked.csv in this R script and then go back to the trait R script where I was 
## checking species names and update the names.

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")

wd_Symbiont_density <- paste(wd,"/Traits_extra/symbiodinium_density",sep="")
wd_Coloniality <- paste(wd,"/Traits_extra/coloniality",sep="")
wd_CoralliteArea <- paste(wd,"/Traits_extra/corallite_area",sep="")
wd_growthForm <- paste(wd,"/Traits_extra/growth_form",sep="")
wd_zooxanthellae <- paste(wd,"/Traits_extra/zooxanthellae",sep="")
wd_Fecundity_polyp <- paste(wd,"/Traits_extra/fecundity_polyp",sep="")
wd_Aggressiveness <- paste(wd,"/Traits_extra/aggressiveness",sep="")
wd_Egg_diameter <- paste(wd,"/Traits_extra/egg_diameter",sep="")
wd_Skeletal_density <- paste(wd,"/Traits_extra/skeletal_density",sep="")
wd_Lipid_content <- paste(wd,"/Traits_extra/lipid_content",sep="")
wd_Chlorophyll_a_concentratio <- paste(wd,"/Traits_extra/chlorophyll_a_concentration",sep="")
wd_phylogeny_coral <- paste(wd,"/Traits_extra/phylogeny_coral",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

listFinal <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(listFinal$nameSp) # 1786
length(levels(as.factor(listFinal$nameSp_checked))) 

list1 <- read.csv(paste(wd,"Traits_extra/dark_respiration/namesCheck.csv",sep="/"),header = T, stringsAsFactors = F)
length(list1$nameSp) # 20

list2 <- read.csv(paste(wd,"Traits_extra/bleaching_response_index/namesCheck.csv",sep="/"),header = T, stringsAsFactors = F)
length(list2$nameSp) # 316

list3 <- read.csv(paste(wd,"Traits_extra/symbiodinium_density/namesCheck.csv",sep="/"),header = T, stringsAsFactors = F)
length(list3$nameSp) # 35

list4 <- read.csv(paste(wd,"Traits_extra/tissue_thickness/namesCheck.csv",sep="/"),header = T, stringsAsFactors = F)
length(list4$nameSp) # 20

list5 <- read.csv(paste(wd,"/Traits_extra/reduced_scattering_coefficient/namesCheck.csv",sep="/"),header = T, stringsAsFactors = F)
length(list5$nameSp) # 89

list6 <- read.csv(paste(wd,"Traits_extra/coloniality/namesCheck.csv",sep="/"),header = T, stringsAsFactors = F)
length(list6$nameSp) # 823

list7 <- read.csv(paste(wd_CoralliteArea,"speciestosendtoworms_matched.csv",sep="/"),header = T, stringsAsFactors = F)
list7 <- list7[,c(6,8)]
length(list7$ScientificName) # 18  no need to change
colnames(list7) <- colnames(listFinal)

list8 <- read.csv(paste(wd_growthForm,"speciestosendtoworms_matched.csv",sep="/"),header = T, stringsAsFactors = F)
length(list8$X.Genus.) # 56
subset(list8,list8$ScientificName == "")
speciesToAdd <- data.frame(
  ScientificName = c("Dipsastraea mirabilis","Barabattoia mirabilis","Scolymia wellsi"),
  ScientificName_accepted = c("Dipsastraea amicorum","Dipsastraea amicorum","Scolymia wellsii"),
  stringsAsFactors = F
)
# "Dipsastraea mirabilis"   Barabattoia mirabilis Dipsastraea amicorum
#  "Scolymia wellsi"      Scolymia wellsii
list8_cut <- list8[,c("ScientificName","ScientificName_accepted")]
# remove rows with list8$ScientificName == ""
list8_cut <- subset(list8_cut,list8_cut$ScientificName != "")
# all the species with list8$ScientificName != "" but list8$ScientificName_accepted == "" are "nomen nudum" --> I keep them
list8_cut$ScientificName_accepted[list8_cut$ScientificName_accepted == ""] <- list8_cut$ScientificName[list8_cut$ScientificName_accepted == ""]
# add speciesToAdd
list8_cut <- rbind(list8_cut,speciesToAdd)
# correct for one name: "Acropora surculosa var. turbinata" --> "Acropora surculosa"
list8_cut$ScientificName[list8_cut$ScientificName ==  "Acropora surculosa var. turbinata"]  <-  "Acropora surculosa"
colnames(list8_cut) <- c("nameSp","nameSp_checked")

list9 <- read.csv(paste(wd_zooxanthellae,"speciestosendtoworms_matched.csv",sep="/"),header = T, stringsAsFactors = F)
length(list9$X.Genus.) # 675
head(list9)
# combine the 1st 2 columns:
list9$species <- with(list9,paste0(X.Genus.,X.Species.))
list9$species <- gsub('""'," ",list9$specie)
list9$species <- gsub('"',"",list9$specie)
list9_cut <- list9[c("species","ScientificName_accepted")]
head(list9_cut)
# checked species wiht list9$ScientificName_accepted == ""
speciesToCheck <- subset(list9_cut,list9_cut$ScientificName_accepted== "")[,1] 
n <- length(speciesToCheck) # 31
# names checked one by one on WoRMS:
speciesCheked <- c("Anthemiphyllia patera costata","Anthemiphyllia patera patera","Caryophyllia (Caryophyllia) ambrosia","Caryophyllia (Caryophyllia) ambrosia","Caryophyllia balaenacea","Caryophyllia cincticulatus","Caryophyllia scillaemorpha","Coenocyathus parvulus",
                   "Culicia tenella natalensis","Culicia tenella tenella","Deltocyathus inusitiatus","Euphyllia fimbriata","Flabellum (Ulocyathus) apertum","Flabellum (Ulocyathus) apertum","Flabellum (Flabellum) transversale","Flabellum (Flabellum) transversale",
                   "Flabellum (Flabellum) transversale","Fungiacyathus (Fungiacyathus) pusillus","Fungiacyathus (Fungiacyathus) pusillus","Goniopora mauritiensis","Heterocyathus aequicostatus","Oulangia stokesiana","Oulangia stokesiana","Paraconotrochus antarcticus",
                   "Paracyathus indicus","Paracyathus indicus","Phyllangia americana","Phyllangia americana",     
                   "Rhizosmilia multipalifera","Truncatoflabellum zuluense","Tubastraea micranthus")
DFToAdd <- data.frame(
  species = speciesToCheck,
  ScientificName_accepted = speciesCheked,
  stringsAsFactors = F
)
n <- length(DFToAdd$species) # 31
for(i in 1:n){
  list9_cut$ScientificName_accepted[list9_cut$species == DFToAdd$species[i]] <- DFToAdd$ScientificName_accepted[i]
}
# add 2 more species:
list9_cut <- rbind(list9_cut,data.frame(species = c("Heterocyathus japonicus","Tubastrea micrantha"), ScientificName_accepted = c("Heterocyathus aequicostatus","Tubastraea micranthus"),stringsAsFactors = F))
colnames(list9_cut) <- c("nameSp","nameSp_checked")

list10 <- read.csv(paste(wd_Aggressiveness,"speciestosendtoworms_matched.csv",sep="/"),header = T, stringsAsFactors = F)
length(list10$X.Genus.) # 16
head(list10)
list10$species <- with(list10,paste0(X.Genus.,X.Species.))
list10$species <- gsub('""'," ",list10$specie)
list10$species <- gsub('"',"",list10$specie)
### some species names to change :
list10$ScientificName[list10$species == "Dipsastrea paillida"] <- "Dipsastraea pallida"
list10$ScientificName_accepted[list10$species == "Dipsastrea paillida"] <- "Dipsastraea pallida"
list10$ScientificName[list10$species == "Sacrophyton glaucum"] <- "Sarcophyton glaucum"
list10$ScientificName_accepted[list10$species == "Sacrophyton glaucum"] <- "Sarcophyton glaucum"
list10$ScientificName_accepted[list10$species == "Turbinaria immersa"] <- "Turbinaria immersa"    # nomen dubium
list10_cut <- list10[c("species","ScientificName_accepted")]
colnames(list10_cut) <- c("nameSp","nameSp_checked")
### add Favia stelligera var. fanningensis in list10_cut$species:
list10_cut <- rbind(list10_cut,c("Favia stelligera var. fanningensis","Goniastrea stelligera"))

list11 <- read.csv(paste(wd_phylogeny_coral,"speciestosendtoworms_matched.csv",sep="/"),header = T, stringsAsFactors = F)
length(list11$X.Genus.) # 162
head(list11)
list11$species <- with(list11,paste0(X.Genus.,X.Species.))
list11$species <- gsub('""'," ",list11$specie)
list11$species <- gsub('"',"",list11$specie)

speciesNamesToChange <- c("Acropora lovelli","Acropora wallaceae","Acropora wallacea","Isopora cylindrica","Isopora meridiana","Coenocyathus brooksi","Deltocyathus agassizii","Paracyathus monteryensis","Dendrophyllia fotojiku","Heteropsammia cochleata",
                          "Leptopsammia brittanica","Favia albidus","Favia danae","Favia maxima","Favia truncatus","Montastrea colemani","Phymastrea salebrosa","Truncatoflabellum viginifarium","Isophyllastrea rigida","Lobophyllia dentatus",
                          "Lobophyllia serratus","Scolymia australis","Pocillopora setichelli","Pocillopora setichellii","Seriatopora dendritica","Seriatopora guttatus","Astrangia astrata","Psammocora decussata","Astreopora montiporina")

speciesNamesCorrected <- c("Acropora loveli","Acropora samoensis","Acropora samoensis","Isopora togianensis","Isopora brueggemanni","Coenocyathus brooki","Deltocyathus agassizi","Paracyathus montereyensis","Dendrophyllia futojiku","Heteropsammia cochlea",
                           "Leptopsammia britannica","Dipsastraea albida","Dipsastraea danai","Dipsastraea maxima","Dipsastraea truncatus","Favites colemani","Paramontastraea salebrosa","Truncatoflabellum vigintifarium","Isophyllia rigida","Lobophyllia dentata",
                           "Lobophyllia serrata","Homophyllia australis","Pocillopora brevicornis","Pocillopora brevicornis","Seriatopora dentritica","Seriatopora guttata","Astrangia atrata","Psammocora contigua","Astreopora monteporina")

# species that should not be changed, but no validated name was suggested:
nameToNotChange <- c("Montipora vaughani","Anthemiphyllia patera","Madracis fragilis","Caryophyllia ambrosia","Caryophyllia aspera","Caryophyllia concreta","Caryophyllia corona","Caryophyllia crypta","Caryophyllia balanacea","Caryophyllia huinayensis",
                     "Caryophyllia laevigata","Caryophyllia oblonga","Caryophyllia tangaroae","Ceratotrochus laxus","Phyllangia american","Stephanocyathus imperialis","Stephanocyathus isabellae","Trochocyathus laboreli","Trochocyathus wellsi",
                     "Balanophyllia dilatata","Balanophyllia helenae","Balanophyllia japonica","Balanophyllia malouinensis","Balanophyllia spongiosa","Balanophyllia striata","Balanophyllia vanderhorsti","Dendrophyllia granosa","Dendrophyllia minima",
                     "Eguchipsammia strigosa","Heteropsammia moretonensis","Thecopsammia elongata","Flabellum apertum","Flabellum transversale","Placotrochides minuta","Truncatoflabellum cumingi","Fungiacyathus pusillus","Meandrina jacksoni","Acanthastrea maxima",
                     "Culicia tenella","Psammocora ramosa","Sphenotrochus cuneolus","Sphenotrochus lindstroemi","Tropidocyathus lessoni","Phyllangia americana")

numberSpToChange <- length(speciesNamesToChange) # 29
for(i in 1:numberSpToChange){
  list11$ScientificName_accepted[list11$species == speciesNamesToChange[i]] <- speciesNamesCorrected[i]
}
numberSpCorrect <- length(nameToNotChange) # 44
for(i in 1:numberSpCorrect){
  list11$ScientificName_accepted[list11$species == nameToNotChange[i]] <- nameToNotChange[i]
}
speciesNotPresentFun(nameToNotChange,list11$species)

list11_cut <- list11[,c("species","ScientificName_accepted")]
colnames(list11_cut) <- c("nameSp","nameSp_checked")

### complete final liste:
listFinal <- list1
length(list1$nameSp) # 20
listFinal <- merge(listFinal,list2,all.x=TRUE,all.y = TRUE)
sum(duplicated(listFinal$nameSp)) # 0
length(listFinal$nameSp) # 319
listFinal <- merge(listFinal,list3,all.x=TRUE,all.y = TRUE)
sum(duplicated(listFinal$nameSp)) # 0
length(listFinal$nameSp) # 324
listFinal <- merge(listFinal,list4,all.x=TRUE,all.y = TRUE)
sum(duplicated(listFinal$nameSp)) # 0
length(listFinal$nameSp) # 325
listFinal <- merge(listFinal,list5,all.x=TRUE,all.y = TRUE)
sum(duplicated(listFinal$nameSp)) # 0
length(listFinal$nameSp) # 344
listFinal <- merge(listFinal,list6,all.x=TRUE,all.y = TRUE)
sum(duplicated(listFinal$nameSp)) # 0
length(listFinal$nameSp) # 823
### from Growth rate dataset:
listFinal <- rbind(listFinal,c("Endopachys bulbosa","Lophelia pertusa"))
listFinal <- rbind(listFinal,c("Lophelia pertusa","Lophelia pertusa"))
length(listFinal$nameSp) # 825
### from Colony maximum diameter:
listFinal <- rbind(listFinal,c("Colpophyllia breviserialis","Colpophyllia natans"))
listFinal <- rbind(listFinal,rep("Helioseris cucullata",2))
listFinal <- rbind(listFinal,rep("Oculina arbuscula",2))
listFinal <- rbind(listFinal,rep("Phacelocyathus flos",2))
listFinal <- rbind(listFinal,rep("Placotrochides frustum",2))
listFinal <- rbind(listFinal,rep("Porites decasepta",2))
length(listFinal$nameSp) # 831
### from corallite area dataset
listFinal <- merge(listFinal,list7,all.x=TRUE,all.y = TRUE)
length(listFinal$nameSp) # 849
### from growth from :
listFinal <- merge(listFinal,list8_cut,all.x=TRUE,all.y = TRUE)
length(listFinal$nameSp) # 906
### from zooxanthellae
listFinal <- merge(listFinal,list9_cut,all.x=TRUE,all.y = TRUE)
length(listFinal$nameSp) # 1583
### from aggressiveness:
listFinal <- merge(listFinal,list10_cut,all.x=TRUE,all.y = TRUE)
length(listFinal$nameSp) # 1600
### From egg diameter:
listFinal <- rbind(listFinal,c("Fungia horrida","Danafungia horrida"))
listFinal <- rbind(listFinal,c("Leptastrea purpurae","Leptastrea purpurea"))
listFinal <- rbind(listFinal,c("Physogyra lichtenseini","Physogyra lichtensteini"))
### correction for " Australophyllia wilsoni" -> "Australophyllia wilsoni" in listFinal$nameSp_checked
listFinal$nameSp_checked[listFinal$nameSp_checked == " Australophyllia wilsoni"] <- "Australophyllia wilsoni"
listFinal <- rbind(listFinal,c("Australophyllia wilsoni","Australophyllia wilsoni"))
### add Favia favus --> ] "Dipsastraea favus" (from Swainet al. 2016 about RSC vs BRI)
listFinal <- rbind(listFinal,c("Favia favus","Dipsastraea favus"))
### add species not already present from coral super tree (Huang and Roy 2015)
listFinal <- merge(listFinal,list11_cut,all.x=TRUE,all.y = TRUE)
length(listFinal$nameSp) 
length(levels(as.factor(listFinal$nameSp_checked))) 
### in coralTrait,org, "Duncanopsammia peltata" and "Turbinaria peltata" are the same species and in WoRMS in the latter is validated and the former is not referenced --> "Duncanopsammia peltata" --> "Turbinaria peltata"
listFinal[listFinal$nameSp == "Duncanopsammia peltata",]$nameSp_checked <- "Turbinaria peltata"
### modofiction of species in Martinique datasets:
listFinal <- rbind(listFinal,
                   c("Agaricia sp","Agaricia sp"),
                   c("Coral juvenile", "Coral juvenile"),
                   c("Coral_nonID","Coral_nonID"),
                   c("Dichocenia stokesii","Dichocoenia stokesii"),
                   c("Millepora sp","Millepora sp"),
                   c("Millepora squarrosa","Millepora squarrosa"),
                   c("Millipora alcicornis","Millepora alcicornis"),
                   c("Millipora complanata","Millepora complanata"),
                   c("Millipora squarrosa","Millepora squarrosa"),
                   c("Montastrea faveolata","Orbicella faveolata"),
                   c("Montastrea franksi","Orbicella franksi"),
                   c("Scolymia sp","Scolymia sp"),
                   c("Stephanocenia michelinii","Stephanocoenia intersepta"),
                   c("Stephanocoenia michelini","Stephanocoenia intersepta"),
                   c("Stephanocoenia mechelinii","Stephanocoenia intersepta"),
                   c("Stephanocoenia michelinii","Stephanocoenia intersepta"),
                   c("Tubastraea aurea","Tubastraea coccinea"))
### Add species in Chancerelle 2000
listFinal <- rbind(listFinal,
                   c("Synaraea rus","Porites rus"),
                   c("Fungia scutaria","Lobactis scutaria"))

listFinal$nameSp <- gsub(" ","_",listFinal$nameSp)
listFinal <- listFinal[order(listFinal$nameSp),]
listFinal$nameSp <- gsub("_"," ",listFinal$nameSp)

sum(duplicated(listFinal$nameSp)) # 0

listFinal[listFinal$nameSp==listFinal$nameSp[duplicated(listFinal$nameSp)],]

length(listFinal$nameSp) # 1786
length(levels(as.factor(listFinal$nameSp_checked))) # 1588

### check for "(" characters in the speciesNameChecked.csv:
listFinal[grepl("\\(",listFinal$nameSp_checked),]
sum(grepl("\\(",listFinal$nameSp_checked)) # 234
# remove the characters in brakets:
listFinal$nameSp_checked <- gsub("\\s*\\([^\\)]+\\)","",listFinal$nameSp_checked)
sum(grepl("\\(",listFinal$nameSp_checked)) 

# remove white space at begining and end of species names:
trimws(listFinal$nameSp_checked)

# sum(grepl("  ",listFinal$nameSp_checked)) # 0
# sum(grepl("\" ",listFinal$nameSp_checked)) # 0     space at the beginning    --> use trimws insead 
# sum(grepl(" \"",listFinal$nameSp_checked)) # 0     space at the end

length(listFinal$nameSp) 
length(levels(as.factor(listFinal$nameSp_checked))) 

# write.csv(listFinal,paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),row.names = F)
