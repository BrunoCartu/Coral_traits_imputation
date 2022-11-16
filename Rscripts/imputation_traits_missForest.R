# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals:
## imput missing trait data (Appendix S1: 2 and 3)
## To create: 
### FTTable_Phylo_randomAdd_9PC_combined.csv: trait data + 9 principal components from phylogeny
## functionTrait_DS_Imputed.csv             : imputed coral traits dataset for 789 species and 11 traits
## functionTrait_DS_Imputed_OOerror.csv     : normalized root mean squared error (= 0.0985) and proportion of falsely classified (= 0.0735)

## From:
### Traits_compilation_zooxanthellate_cut.csv: 798 species, 22 traits, created in phylogeny_coral_and_functional_traits.R
### supertree_names_validated_1000_798sp.tre : the updated and cut supertree from (Huang and Roy 2015), produced inphylogeny_coral_and_functional_traits.R


require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")

require(missForest)

# source(paste(wd_Rscripts,"functions_coralModel.R",sep="/"))
source(paste(wd_Rscripts,"functions.R",sep="/"))

### Import the coral trait dataset with only zooxanthellae ----
# and without species for which (1) not in coral phylogeny and (2) only had value for categorical traits --> 29 species
# and "Sclerophyllia maxima" because (1) not in the phylogeny and (2) the only available numerical value is for corallite area --> 1 
# - 30 species 
# Functional_traits_zooxan_cut <- read.csv(paste(wd_Datasets_original,"Traits_compilation_zooxanthellate_cut.csv",sep = "/"),header=T,stringsAsFactors = F)
Functional_traits_zooxan_cut <- read.csv(paste(wd_Datasets,"Traits_compilation_zooxanthellate_cut.csv",sep = "/"),header=T,stringsAsFactors = F)
length(Functional_traits_zooxan_cut$species) # 798

### Import updated and cut supertree for coral from (Huang and Roy 2015) --> names have been updated -----
# 1000 trees 
# 11 species from the trait dataset have been added at the genus level (randomly to avoid creation polytomies): see phylogeny_coral_and_functional_traits.R
# 713 species that were not in Traits_compilation_zooxanthellate_cut.csv have been removed
supertree <- read.tree(paste(wd_Datasets,"supertree_names_validated_1000_798sp.tre",sep="/"))
length(supertree[[1]]$tip.label) # 798
head(supertree[[1]]$tip.label)

# remove "_" between genus and species name in trees
numberTrees <- length(supertree)
for(i in 1:numberTrees){
  supertree[[i]]$tip.label <- gsub("_"," ",supertree[[i]]$tip.label)
}
head(supertree[[1]]$tip.label)

# check if same species in FT table and no duplicated names:
length(speciesPresentBothList(Functional_traits_zooxan_cut$species,supertree[[1]]$tip.label)) # 798
sum(duplicated(Functional_traits_zooxan_cut$species)) # 0

### do PCA on average distance matrices ------
## creation of a list of 1000 distance matrices (1 for each tree): (see Swenson 2014, Functional and Phylogenetic ecology in R, p. 152):
listDisMatrices <- list()
numberTrees <- length(supertree)  # 1000
for(i in 1:numberTrees){
  listDisMatrices[[i]] <- cophenetic(supertree[[i]])
}

### creation of distance matrix that is the mean of the 1000 distance matrices -----
numberSpecies <- length(supertree[[1]]$tip.label) # 798
distanceMatrix <- matrix(0,nrow=numberSpecies,ncol=numberSpecies)
for(i in 1:numberTrees){
  distanceMatrix <- distanceMatrix + listDisMatrices[[i]]
}
distanceMatrix <- distanceMatrix / numberTrees
distanceMatrix[1:10,1:10]

### PCA on the distanceMatrix
pr <-  prcomp(distanceMatrix)   # princomp did not work ?!
summary_pr <- summary(pr)
summary_pr$importance[,1:10]

# calcualte the variance, its proportion and cumulative proportion
princ <- princomp(distanceMatrix)
princ$sdev
var <- princ$sdev * princ$sdev
prop_var <- var / sum(var)
sum_prop_var <- cumsum(prop_var)

summary_PCA <- rbind(princ$sdev,var,prop_var,sum_prop_var)
rownames(summary_PCA) <- c("sd","variance","prop_variance","cumulative_proportion")
summary_PCA[,1:10]

### Determination of the number of eigenvector to include using the "broken stick" cumulative percentage BSCP and broken stick percentages BSP:
BSCP <- c()
BSCP[1] <- 0.5
BAP <- c()
BAP[1] <- 0.5
for(i in 2: 15){
  BSCP[i] <- BSCP[i-1] + (1-BSCP[i-1])/2
  BAP[i] <- (1-BSCP[i-1])/2
}
BSCP
summary_PCA[4,1:15]

for(i in 1:15){
  if(summary_PCA[4,i] < BSCP[i]){
    print(i)
    cat(summary_PCA[4,i],BSCP[i])
    break
  }
}
#  10
# 0.998864 0.9990234

# plot the % variance explained by 1st 9 eigenvectors 
graphics.off()
par(mar=c(4.5,3.8,0.5,0.5))
barplot(100*summary_PCA[3,1:9], las=2, xlab='', ylab='',ylim=c(0,100))
mtext('% Variance Explained',side=2,line=2.5,cex=1)
text(x=c(0.8,2,3.2,4.4,5.5,6.7,7.9,9.1,10.3),y=100*summary_PCA[3,1:9]+5,labels=paste(as.character(formatC(round(100*summary_PCA[3,1:9],2),format="f",digits = 2)),"%"),cex=0.8)

par(mar=c(5,4.5,2,2))
barplot(100*summary_PCA[3,1:9], las=2, xlab='', ylab='% Variance Explained',ylim=c(0,100))
text(x=c(0.8,2,3.2,4.4,5.5,6.7,7.9,9.1,10.3),y=100*summary_PCA[3,1:9]+5,labels=paste(as.character(formatC(round(100*summary_PCA[3,1:9],2),format="f",digits = 2))," %"),cex=0.9)

table9PC <- data.frame(
  species = rownames(princ$scores) ,
  comp1 = round(as.numeric(princ$scores[,1]),2) ,
  comp2 = round(as.numeric(princ$scores[,2]),2) ,
  comp3 = round(as.numeric(princ$scores[,3]),2) ,
  comp4 = round(as.numeric(princ$scores[,4]),2) ,
  comp5 = round(as.numeric(princ$scores[,5]),2) ,
  comp6 = round(as.numeric(princ$scores[,6]),2) ,
  comp7 = round(as.numeric(princ$scores[,7]),2) ,
  comp8 = round(as.numeric(princ$scores[,8]),2) ,
  comp9 = round(as.numeric(princ$scores[,9]),2) ,
  stringsAsFactors = FALSE
)

### combination of Functional_traits_zooxan_cut and table9PC: ----
## combination of the 2 tables:
Functional_traits_zooxan_cut_phylo9PC <- merge(Functional_traits_zooxan_cut,table9PC,all.x=TRUE,all.y = TRUE)

### write csv file:
# write.csv(Functional_traits_zooxan_cut_phylo9PC,paste(wd_Datasets_original,"FTTable_Phylo_randomAdd_9PC_combined.csv",sep="/"),row.names = FALSE)

### filling up gaps in the trait data base:
### remove traits not in use in the model:
# columnsToKeep <- c("species","aggressiveness","bleaching_response_index","coloniality","colony_max_diameter","corallite_area","egg_diameter",
#                    "fecundity_polyp","growthFormClasses","growth_rate","mode_larval_development","reduced_scattering_coefficient","skeletal_density",
#                    "comp1","comp2","comp3","comp4","comp5","comp6","comp7","comp8","comp9")

columnsToKeep <- c("species","aggressiveness","coloniality","colony_max_diameter","corallite_area","egg_diameter",
                   "fecundity_polyp","growthFormClasses","growth_rate","mode_larval_development","reduced_scattering_coefficient","sexual_system",
                   "comp1","comp2","comp3","comp4","comp5","comp6","comp7","comp8","comp9")

Functional_traits_zooxan_cut_phylo9PC_cut <- Functional_traits_zooxan_cut_phylo9PC[,columnsToKeep]

# need to convert categorical variable into factor:
Functional_traits_zooxan_cut_phylo9PC_cut$coloniality <- as.factor(Functional_traits_zooxan_cut_phylo9PC_cut$coloniality)
Functional_traits_zooxan_cut_phylo9PC_cut$growthFormClasses <- as.factor(Functional_traits_zooxan_cut_phylo9PC_cut$growthFormClasses)
Functional_traits_zooxan_cut_phylo9PC_cut$mode_larval_development <- as.factor(Functional_traits_zooxan_cut_phylo9PC_cut$mode_larval_development)
Functional_traits_zooxan_cut_phylo9PC_cut$sexual_system <- as.factor(Functional_traits_zooxan_cut_phylo9PC_cut$sexual_system)

# in-filling values w:-----
Functional_traits_zooxan_cut_phylo9PC_cut_Imp <- missForest(Functional_traits_zooxan_cut_phylo9PC_cut[,-1], verbose = TRUE)

Functional_traits_zooxan_cut_phylo9PC_cut_Imp$OOBerror
#    NRMSE        PFC 
# 0.09424169 0.07742333 

functionTrait_DS_Imputed <- cbind(Functional_traits_zooxan_cut_phylo9PC$species,Functional_traits_zooxan_cut_phylo9PC_cut_Imp$ximp[,c("aggressiveness","coloniality","colony_max_diameter","corallite_area","egg_diameter",
                                                                                                                                      "fecundity_polyp","growthFormClasses","growth_rate","mode_larval_development",
                                                                                                                                      "reduced_scattering_coefficient","sexual_system")])
colnames(functionTrait_DS_Imputed)[1] <- "species"

### write file:
# write.csv(functionTrait_DS_Imputed,paste(wd_Datasets,"functionTrait_DS_Imputed.csv",sep="/"),row.names = F)
# write.csv(Functional_traits_zooxan_cut_phylo9PC_cut_Imp$OOBerror,paste(wd_Datasets,"functionTrait_DS_Imputed_OOerror.csv",sep="/"),row.names = T)
