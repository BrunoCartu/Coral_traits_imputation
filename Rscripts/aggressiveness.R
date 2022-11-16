# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goals: 
## create the aggressivenes traits from six sources: Lang 1973, Sheppard 1979, Logan 1984, Dai 1990, Abelson and Loya 1999, Connell et al., 1999
## create 
### aggressiveness.csv: the aggressive ranking score for 132 species, values are in [0,100], higer values for more aggressive species

## From: 
### aggressiveness_0_to_1.csv: the normalized ranking I did manually from the six different study; I copy pasted the values from Literature review.xls

require(here)              
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_Aggressiveness <- paste(wd,"/Traits_extra/aggressiveness",sep="")

source(paste(wd_Rscripts,"functions.R",sep="/"))

require(Kendall) 

#*********************************** NOTE about aggressiveness ranking ***********************************
# Species are ranked accordingly to there aggressivity. Ranking data comes from six studies that used different metric to determine the aggressivity of coral species:
# (1) number of subordinates, (2) Coral Index: (#wins - #losses)/(#wins + #losses), (3) the % of winning interactions.
# I first manually normalized the ranking in each study from 0 to 1 (1 being the most aggressive species): I created in Traits_extra > aggressiveness > Literature review.xls then aggressiveness_0_to_1.csv then aggressiveness.csv
# I then computed the iterative partial rank aggregation pivoting algorithm (IPRAPA) developed and explained in details by (Swain et al., 2016. Functional Ecology, Vol. 31).
# The final ranking values values are comprised between 0 (not aggressive) and 100 (most aggressive)
# Data soure in wd_Aggressiveness:
#   Literature review.xlsx
#   aggressiveness_0_to_1.csv
#*******************************************************************************************************

### Import dataset from folder: -------
speciesHierarchies <- read.csv(paste(wd_Aggressiveness,"aggressiveness_0_to_1.csv",sep="/"),header = TRUE,stringsAsFactors = F)
speciesHierarchies[speciesHierarchies$species == "",]
speciesHierarchies <- speciesHierarchies[speciesHierarchies$species != "",]

### update species names
speciesToCheck <- checkSpeciesNames(speciesHierarchies$species, wd_Datasets = wd_Datasets)
speciesToCheck
for(i in 1:length(speciesToCheck$nameSp)){
  speciesHierarchies$species[speciesHierarchies$species == speciesToCheck$nameSp[i]] <- speciesToCheck$nameSp_checked[i]
}

### check for duplicated species 
TableDuplicates <- speciesHierarchies[speciesHierarchies$species %in% speciesHierarchies[duplicated(speciesHierarchies$species),],]
speciesHierarchies[speciesHierarchies$species %in% speciesHierarchies[duplicated(speciesHierarchies$species),],]$Abelson_and_Loya_1999[1] <- TableDuplicates$Abelson_and_Loya_1999[2]
speciesHierarchies[speciesHierarchies$species %in% speciesHierarchies[duplicated(speciesHierarchies$species),],]$ Connell_et_al_.2004[2] <- TableDuplicates$ Connell_et_al_.2004[1]
speciesHierarchies[speciesHierarchies$species %in% speciesHierarchies[duplicated(speciesHierarchies$species),],]

speciesHierarchies <- speciesHierarchies[!duplicated(speciesHierarchies$species),]
sum(duplicated(speciesHierarchies$species))

### iterative partial-rank aggregation algorithm (IPRARA from Swain et al., 2016):
## Prior to the first iteration, all lists are sorted by the number of cohorts (i.e., levels of aggressiveness) and then phylotypes, and are assessed in order of decreasing information content
list_Lang1973 <- speciesHierarchies[,c("species","Lang_1973")]
list_Lang1973 <- subset(list_Lang1973,!is.na(list_Lang1973$Lang_1973))
# number cohorts (hierarchy levels):
length(levels(as.factor(list_Lang1973$Lang_1973)))           # 16
list_Lang1973$species <- as.factor(list_Lang1973$species)
# number species:
length(list_Lang1973$species)                               # 27

list_Sheppard_1979 <- speciesHierarchies[,c("species","Sheppard_1979")]
list_Sheppard_1979 <- subset(list_Sheppard_1979,!is.na(list_Sheppard_1979$Sheppard_1979))
# number cohorts (hierarchy levels):
length(levels(as.factor(list_Sheppard_1979$Sheppard_1979)))    # 21
list_Sheppard_1979$species <- as.factor(list_Sheppard_1979$species)
# number species:
length(list_Sheppard_1979$species)                             # 26

list_Logan_1984 <- speciesHierarchies[,c("species","Logan_1984")]
list_Logan_1984 <- subset(list_Logan_1984,!is.na(list_Logan_1984$Logan_1984))
# number cohorts (hierarchy levels):
length(levels(as.factor(list_Logan_1984$Logan_1984)))       # 12
list_Logan_1984$species <- as.factor(list_Logan_1984$species)
# number species:
length(list_Logan_1984$species)                             # 16

list_Dai_1990 <- speciesHierarchies[,c("species","Dai_1990")]
list_Dai_1990 <- subset(list_Dai_1990,!is.na(list_Dai_1990$Dai_1990))
# number cohorts (hierarchy levels):
length(levels(as.factor(list_Dai_1990$Dai_1990)))         # 5
list_Dai_1990$species <- as.factor(list_Dai_1990$species)
# number species:
length(list_Dai_1990$species)                             # 76

# 
# levels(as.factor(list_Dai_1990$Dai_1990)) #"0.0975" "0.2975" "0.4975" "0.6975" "0.9" 
# # list_Dai_1990 need to be changed to go from 0 to 1 (and not from 0.0975 to 0.9)
# max <- max(list_Dai_1990$Dai_1990)
# min <- min(list_Dai_1990$Dai_1990)
# numberSpecies <- length(list_Dai_1990$species)
# for(i in 1:numberSpecies){
#   list_Dai_1990$Dai_1990[i] <- (list_Dai_1990$Dai_1990[i] - min)/(max-min)
# }
# levels(as.factor(list_Dai_1990$Dai_1990))

list_Abelson_and_Loya_1999 <- speciesHierarchies[,c("species","Abelson_and_Loya_1999")]
list_Abelson_and_Loya_1999 <- subset(list_Abelson_and_Loya_1999,!is.na(list_Abelson_and_Loya_1999$Abelson_and_Loya_1999))
# number cohorts (hierarchy levels):
length(levels(as.factor(list_Abelson_and_Loya_1999$Abelson_and_Loya_1999)))         # 17
list_Abelson_and_Loya_1999$species <- factor(list_Abelson_and_Loya_1999$species)
# number species:
length(list_Abelson_and_Loya_1999$species)                                          # 33

list_Connell_et_al_.2004 <- speciesHierarchies[,c("species","Connell_et_al_.2004")]
list_Connell_et_al_.2004 <- subset(list_Connell_et_al_.2004,!is.na(list_Connell_et_al_.2004$Connell_et_al_.2004))
# number cohorts (hierarchy levels):
list_Connell_et_al_.2004$species <- factor(list_Connell_et_al_.2004$species)
length(levels(as.factor(list_Connell_et_al_.2004$Connell_et_al_.2004)))           # 13
# number species:
length(list_Connell_et_al_.2004$species)                                          # 13

### Function that gives the number or the name(s) of species common in 2 list of names: -----
speciesCommonFun <- function(liste1 = character(),liste2 = character(),Number_or_species = c("number","species")){
  size1 <- length(liste1)
  size2 <- length(liste2)
  numbercommonSpecies <- 0
  speciesList <- c()
  for(i in 1: size1){
    for(j in 1 : size2){
      if(liste1[i] == liste2[j]){
        numbercommonSpecies <- numbercommonSpecies + 1
        speciesList[numbercommonSpecies] <- liste1[i]
        #print(liste1[i])
        break 
      }
    }
  }
  if(Number_or_species == "number"){
    numbercommonSpecies
  }else{
    speciesList
  }
}

## Function that creates a nxn matrices displaying the number of common species between 2 lists of species names:
# The lists come in a list ------
speciesCommonMatrixFun <- function(listLists = list()){
  numberLists <- length(listLists)
  listSizeVect <- c()
  matriceList <- matrix(ncol = numberLists, nrow = numberLists) 
  for(i in 1:numberLists){
    listSizeVect[i] <- length(listLists[[i]])
  }
  for(i in 1:numberLists){
    for(j in 1 : numberLists){
      matriceList[i,j] <- speciesCommonFun(listLists[[i]],listLists[[j]],"number")
      if(i == j){
        matriceList[i,j] <- 0
      }
    }
  }
  matriceList
}

## Function that transforms a table so that ranking score is place instead of a value from 0 (= not aggressive --> low ranking score) to 1 (= aggressive --> high ranking score)
## 1 --> high ranking score (most aggressive), and 0 --> lowest ranking score 
## the table input has 2 columns, the 1st with species names and the 2nd wiht the ranking score
rankingFunction <- function(listScores = data.frame()){
  numberRanks <- length(levels(as.factor(listScores[,2])))
  numberSpecies <- length(listScores[,1])
  dataFrame <- listScores
  levelsNum <- as.numeric(levels(as.factor(listScores[,2])))  # this orders the values from the lowest to the biggest
  for(i in 1:numberSpecies){
    for(j in 1:numberRanks){
      if(round(dataFrame[i,2],4) == round(levelsNum[j],4)){
        dataFrame[i,2] <- j
        break
      }
    }
  }
  dataFrame
}

## Function that returns true if a species name is present in a list of names:
speciesPresent <- function(nameSpecies = character(1),listNames = character()){
  numberSpecies <- length(listNames)
  bool <- FALSE
  for(i in 1:numberSpecies){
    if(nameSpecies == listNames[i]){
      bool <- TRUE
      break
    }
  }
  bool
}

## Function that takes a list of ranked-species-tables and creates a new list of tables with new ranking made from R scores (the initial tables are extended)
## the new Rscores and related rankings are added to the tables to keep track of the change
# THe "pivot" is a species present in 2 consecutive lists that is used to compare the speices between the 2 lists. If more than one species, then the species with the least uncertainty is chosen. if still several pivot candidates, then the 1st species of the list is chosen (alternatively it could be random but not imoplemented)
RscoreRanking <- function(listRankList = list(),iterationNumber = numeric()){
  ## add a empty colun in each data frame for the R scores:
  numberLists <- length(listRankList)
  if(iterationNumber > 1){
    numberLists <- numberLists - 1         # to not consider the meanNromR table places at tge end of the list of the lists
  }
  for(i in 1: numberLists){ # add a NA column to each list 
    numberSpeciesList <- length(listRankList[[i]][1])
    listRankList[[i]] <- cbind(listRankList[[i]],rep(NA,numberSpeciesList))
    colnames(listRankList[[i]])[2*iterationNumber + 1] <- paste("R",iterationNumber,sep="_")
  }
  #
  for(i in 1: numberLists){
    if(i == 1){             # attribution of the R scores for the 1st list, which is different because there is no pivot species involved yet
      numberRanks <- max(listRankList[[i]][,2*iterationNumber])
      numberSpecies <- length(listRankList[[i]][,1])
      for(j in 1: numberSpecies){
        listRankList[[i]][j,2*iterationNumber+1] <- 100 * listRankList[[i]][j,2*iterationNumber] / numberRanks
      }
      # cat("Transition list ",i-1," to list ",i,"\n")
    }else{
      # selection of the pivot species: it is a species common in both consecutive lists (i and i-1) with the highest cohort value in the current list i
      # if more than one species, then the species with the least uncertainty is chosen
      # if still several pivot candidates, then the 1st species of the list is chosen 
      numberCommonSpecies <- speciesCommonFun(listRankList[[i]][,1],listRankList[[i-1]][,1],"number")
      pivot <- speciesCommonFun(listRankList[[i]][,1],listRankList[[i-1]][,1],"species")
      numberPreviousList <- i - 1
      cat("Transition list ",i-1," to list ",i,"\n")
      cat("numberCommonSpecies: ",numberCommonSpecies,"\n")
      cat("pivot:",pivot,"\n")
      if(numberCommonSpecies > 1){    # if more than 1 common species
        # chose species with highest cohort ranking in list i:
        pivotsCohortScores <- subset(listRankList[[i]],listRankList[[i]][,1] %in% pivot)
        pivot <- subset(pivotsCohortScores[,1],pivotsCohortScores[,2*iterationNumber] == max(pivotsCohortScores[,2*iterationNumber]))
        if(length(pivot) > 1){ # still more than one species -> calculate uncertainty
          # calculation of uncertainty:
          UncertaintySpecies <- data.frame(species=pivot,  U = numeric(length(pivot)))        # data frame containing the value of the uncertainty for each species 
          for(k in 1:length(pivot)){
            squaresRatios <- rep(NA,numberPreviousList)             
            for(l in 1: numberPreviousList){
              if(speciesPresent(pivot[k],listRankList[[i-l]][,1])){
                squaresRatios[l] <- (100/max(listRankList[[i-l]][,2*iterationNumber]))^2
              }
            }
            UncertaintySpecies[k,2] <- sqrt(sum(squaresRatios,na.rm = TRUE)) / length(na.omit(squaresRatios))
          }
          pivot <- as.character(subset(UncertaintySpecies[,1],UncertaintySpecies[,2] == min(UncertaintySpecies[,2],na.rm = TRUE)))
          if(length(pivot) > 1){ # if still more than one species, chose the 1st in the list
            pivot <- pivot[1]
          }
          print(UncertaintySpecies)
        }
        cat("pivot:",pivot,"\n")
      }
      # calculation of pivot score A for the pivot
      Ascore <- numeric(1)
      Rscores <- rep(NA,numberPreviousList)
      for(l in 1: numberPreviousList){
        if(speciesPresent(pivot,listRankList[[i-l]][,1])){
          Rscores[l] <- subset(listRankList[[i-l]][,2*iterationNumber+1],listRankList[[i-l]][,1] == pivot)
          print(subset(listRankList[[i-l]],listRankList[[i-l]][,1] == pivot))
        }
      }
      Ascore <- sum(Rscores,na.rm = TRUE) / length(na.omit(Rscores))
      cat("Ascore: ",Ascore,"\n")
      # calculation of the R for other species in list i based on A of pivot:
      numberSpeciesList <- length(listRankList[[i]][,1])
      for(m in 1: numberSpeciesList){
        listRankList[[i]][m,2*iterationNumber+1] <- Ascore * listRankList[[i]][m,2*iterationNumber] / subset(listRankList[[i]][,2*iterationNumber], listRankList[[i]][,1] == pivot)
      }
    }
    print(listRankList[[i]])
    cat("\n")
  }
  # determination of average iteration score for each species:
  Rcol <- 2*iterationNumber + 1
  combineTables <- listRankList[[1]][,c(1,Rcol)]
  numberLists <- length(listRankList)
  if(iterationNumber > 1){
    numberLists <- numberLists - 1         # to remove the meanNromR table
  }
  for(i in 2: numberLists){
    combineTables <- merge(combineTables,listRankList[[i]][,c(1,Rcol)],all.x = TRUE, all.y = TRUE)
  }
  meanRscores <- tapply(combineTables[,2],combineTables$species,mean)
  MaxRscore <- max(meanRscores)
  Iteration1Rscores <- data.frame(
    species = as.character(rownames(meanRscores)),
    meanNormR = round(as.numeric(meanRscores) * 100/MaxRscore,1)
  )
  
  ### place all ranked tables and meanR table in a list "Iteration1":
  Iteration1 <- listRankList   # change name     "Iteration1" 
  lengthList <- length(Iteration1)
  if(iterationNumber == 1){             # a table that contains the mean R scores for each species is added to list Iteration1
    Iteration1[[lengthList + 1]] <- Iteration1Rscores
    lengthList <- lengthList + 1
  }else{                               # here the table 
    Iteration1[[lengthList]] <- merge(Iteration1[[lengthList]],Iteration1Rscores)
  }
  colnames(Iteration1[[lengthList]])[iterationNumber + 1] <- paste("meanNormR",iterationNumber,sep="_")
  
  ### replace the R values in each individual tables with the meanNormR scores:
  lengthList <- lengthList - 1
  for(i in 1:lengthList){
    numberSpecies <- length(Iteration1[[i]][,1])
    for(j in 1:numberSpecies){
      Iteration1[[i]][j,2*iterationNumber + 1] <- subset(Iteration1[[lengthList+1]][,iterationNumber+1],Iteration1[[lengthList+1]][,1] == Iteration1[[i]][j,1])
    }
    colnames(Iteration1[[i]])[2*iterationNumber + 1] <- paste("meanNormR",iterationNumber,sep="_")
  }
  
  ### add a new ranking column in each table based on the new Rscores:
  numberCol <- length(colnames(Iteration1[[1]]))
  for(i in 1:lengthList){
    Iteration1[[i]][,numberCol+1] <- rankingFunction(Iteration1[[i]][,c(1,numberCol)])[,2]
    colnames(Iteration1[[i]])[numberCol+1] <- paste("rankIt",iterationNumber,sep="_")
  }
  Iteration1
}

## Function that returns the data frame of % reversed pair associations between 2 hierachy matrices (composed of 1 for dominance, 0 for equality and -1 for domination) 
percentageSimilarity <- function(matrix1 = matrix(),matrix2 = matrix()){
  matrixDif <- matrix1 - matrix2
  numberReversed <- 0
  dim <- length(matrixDif[,1])
  for(i in 1:dim){
    for(j in 1:dim){
      if(matrixDif[i,j] == 2 || matrixDif[i,j] == -2){
        numberReversed <- numberReversed + 1 
      }
    }
  }
  numberReversed <- numberReversed / 2
  numberAssociations <- dim * (dim - 1) / 2
  proportionReversed <- round(numberReversed / numberAssociations * 100,2)
  proportionReversed
}

### function that determines what is the proportion of reversed pair-relationship in each dataset (loss or gain of equality is not considered):
proportionReversedFun <- function(rankingIteration = list()){
  numberLists <- length(rankingIteration) - 1
  NumberIteration <- (length(rankingIteration[[1]][1,]) - 2) / 2
  proportionReverseTable <- as.data.frame(matrix(NA,ncol= NumberIteration,nrow = numberLists))  #  table that contains the proportion of reversed pair-associations for numberLists x NumberIteration
  for(j in 1: numberLists){
    rownames(proportionReverseTable)[j] <- paste("table",j,sep="_")
  }
  for(j in 1: NumberIteration){
    colnames(proportionReverseTable)[j] <- paste("iter",j,sep="_")
  }
  for(i in 1 : numberLists){
    listMatrices <- list()         # list that is going to contain the different pair-association matrices
    numberSpecies <- length(rankingIteration[[i]][,1])
    matrixOriginal <- matrix(,ncol = numberSpecies, nrow = numberSpecies)
    # make the pair-association matrix for original ranking:
    for(j in 1:numberSpecies){
      for(k in 1:numberSpecies){
        if(rankingIteration[[i]][j,2] > rankingIteration[[i]][k,2]){
          matrixOriginal[j,k] <- 1
        }else if(rankingIteration[[i]][j,2] == rankingIteration[[i]][k,2]){
          matrixOriginal[j,k] <- 0
        }else {
          matrixOriginal[j,k] <- -1
        }
      }
    }
    # make the pair-association matrices for the iterated ranking:
    listMatrices[[1]] <- matrixOriginal       # the 1st matrix in the list is the one form the original data set against the other matrices are going to be compared 
    for(l in 1:NumberIteration){
      matrixIt <- matrix(,ncol = numberSpecies, nrow = numberSpecies)
      for(j in 1:numberSpecies){
        for(k in 1:numberSpecies){
          if(rankingIteration[[i]][j,(l+1)*2] > rankingIteration[[i]][k,(l+1)*2]){
            matrixIt[j,k] <- 1
          }else if(rankingIteration[[i]][j,(l+1)*2] == rankingIteration[[i]][k,(l+1)*2]){
            matrixIt[j,k] <- 0
          }else {
            matrixIt[j,k] <- -1
          }
        }
      }
      listMatrices[[1 + l]] <- matrixIt
    }
    # calculate the proportion of reversed pair-associations between orginal martix and the iterated ones:
    listProportionReversed <- c()
    for(l in 1:NumberIteration){
      listProportionReversed[l] <- percentageSimilarity(listMatrices[[1]],listMatrices[[1+NumberIteration]])
    }
    proportionReverseTable[i,] <- listProportionReversed
  }
  # add another row with cumulative % of reversed relationships:
  numberTables <- length(proportionReverseTable[,1])
  numberIterations <- length(proportionReverseTable[1,])
  sumPercentReversed <- numeric()
  for(i in 1:numberIterations){
    sum <- 0
    for(j in 1:numberTables){
      sum <- sum + proportionReverseTable[j,i]
    }
    sumPercentReversed[i] <- sum
  }
  proportionReverseTable <- rbind(proportionReverseTable,sumPercentReversed)
  rownames(proportionReverseTable)[numberTables+1] <- "sum"
  proportionReverseTable
}

### Funciton that determines the iterative partial ranking on a list of tables and there iterated R scores, for s given number of iterations
iterativePartialRankingFun <- function(listRankList = list(), numberIterations = numeric()){
  listRankListIt <- listRankList
  for(i in 1:numberIterations){
    listRankListIt_I <- RscoreRanking(listRankListIt,i)
    listRankListIt <- listRankListIt_I
  }
  listRankListIt
}

### Function that computes the Kendall's Tau instead of proportion of reversed relationships
KendallTauFun <- function(rankingIteration = list()){
  numberLists <- length(rankingIteration) - 1
  NumberIteration <- (length(rankingIteration[[1]][1,]) - 2) / 2
  KendallTauTable <- as.data.frame(matrix(NA,ncol= NumberIteration,nrow = numberLists))  # table that contains the proportion of reversed pair-associations for numberLists x NumberIteration
  for(i in 1:numberLists){
    for(j in 1:NumberIteration){
      sp <- rankingIteration[[i]][,1]
      numberSpecies <- length(sp) # 
      initialRanking <- rankingIteration[[i]][,2]
      rankingIteration[[7]][,1] <- as.character(rankingIteration[[7]][,1])
      iteartionRanking_cut <- rankingIteration[[7]][,c(1,j+1)][rankingIteration[[7]]$species %in% sp,]
      finalRanking <- rankingFunction(iteartionRanking_cut)[,2]
      Kendall <- Kendall(initialRanking,finalRanking)
      KendallTauTable[i,j] <- round(Kendall$tau,2)[1]
      colnames(KendallTauTable)[j] <- paste("iter",j,sep="_")
      rownames(KendallTauTable)[i] <- paste("table",i,sep="_")
    }
  }
  # add another row with sum of Tau by iteration through all the tables:
  numberTables <- length(KendallTauTable[,1])
  numberIterations <- length(KendallTauTable[1,])
  sumKendallTau <- numeric()
  for(i in 1:numberIterations){
    sum <- 0
    for(j in 1:numberTables){
      sum <- sum + KendallTauTable[j,i]
    }
    sumKendallTau[i] <- sum
  }
  KendallTauTable <- rbind(KendallTauTable,sumKendallTau)
  rownames(KendallTauTable)[numberTables+1] <- "sum"
  KendallTauTable
}

### Function that compute the "coherence" value (cf. Fang et al., 2010. Frontiers in Algorithmics. 6213 LNCS)
# the coherence considers both the normalized Kendall-τ distance and the size of overlap between two partial rankings
# the coherence is calculated between the initial ranking of each list and the total ranking with all species, for each iterations
coherenceFun <- function(rankingIteration = list()){
  numberLists <- length(rankingIteration) - 1
  NumberIteration <- (length(rankingIteration[[1]][1,]) - 2) / 2
  KendallTauTable <- as.data.frame(matrix(NA,ncol= NumberIteration,nrow = numberLists))  # table that contains the proportion of reversed pair-associations for numberLists x NumberIteration
  coherenceTable <- as.data.frame(matrix(NA,ncol= NumberIteration,nrow = numberLists))
  for(i in 1:numberLists){
    for(j in 1:NumberIteration){
      sp <- rankingIteration[[i]][,1]
      numberSpecies <- length(sp) # 26
      initialRanking <- rankingIteration[[i]][,2]
      rankingIteration[[7]][,1] <- as.character(rankingIteration[[7]][,1])
      iteartionRanking_cut <- rankingIteration[[7]][,c(1,j+1)][rankingIteration[[7]]$species %in% sp,]
      finalRanking <- rankingFunction(iteartionRanking_cut)[,2]
      Kendall <- Kendall(initialRanking,finalRanking)
      coherence <- numberSpecies *  (1 - Kendall$tau[1] / (factorial(numberSpecies)/(2*factorial(numberSpecies - 2))))
      coherenceTable[i,j] <- coherence
      #coherenceTable[i,j] <- round(coherence,2)
      colnames(coherenceTable)[j] <- paste("iter",j,sep="_")
      rownames(coherenceTable)[i] <- paste("table",i,sep="_")
    }
  }
  # add another row with sum of coherence by iteration through all the tables:
  numberTables <- length(coherenceTable[,1])
  numberIterations <- length(coherenceTable[1,])
  sumCoherence <- numeric()
  for(i in 1:numberIterations){
    sum <- 0
    for(j in 1:numberTables){
      sum <- sum + coherenceTable[j,i]
    }
    sumCoherence[i] <- sum
  }
  coherenceTable <- rbind(coherenceTable,sumCoherence)
  rownames(coherenceTable)[numberTables+1] <- "sum"
  coherenceTable
}

## Classification of the lists by 1st their number of cohorts and then phylotypes: 
list1 <- list_Sheppard_1979               # 21     26
list2 <- list_Abelson_and_Loya_1999       # 17     33
list3 <- list_Lang1973                    # 16     27
list4 <- list_Connell_et_al_.2004         # 13     13
list5 <- list_Logan_1984                  # 12     16
list6 <- list_Dai_1990                    # 5      76

list1$species <- as.character(list1$species)
list2$species <- as.character(list2$species)
list3$species <- as.character(list3$species)
list4$species <- as.character(list4$species)
list5$species <- as.character(list5$species)
list6$species <- as.character(list6$species)

## check if pivot species present between 2 consecutive lists:
speciesCommonFun(list1$species,list2$species,"number")    # 8           1, 2
speciesCommonFun(list2$species,list3$species,"number")    # 0           1, 2, 3 --> no
#
speciesCommonFun(list2$species,list4$species,"number")    # 4           1, 2, 4 --> ok
speciesCommonFun(list4$species,list3$species,"number")    # 1           1, 2, 4, 3 --> ok but not top
speciesCommonFun(list4$species,list3$species,"species")  # Porites 5sp
speciesCommonFun(list3$species,list5$species,"number")    # 7           1, 2, 4, 3, 5 --> ok
speciesCommonFun(list6$species,list5$species,"number")    # 0           1, 2, 4, 3, 5, 6 --> no, neither is  1, 2, 4, 3, 6, 5 

speciesCommonFun(list4$species,list5$species,"number")    # 0           1, 2, 4, 5 --> no
speciesCommonFun(list4$species,list6$species,"number")    # 5           1, 2, 4, 6 --> possibly
speciesCommonFun(list6$species,list3$species,"number")    # 0           1, 2, 4, 6, 3 --> no
speciesCommonFun(list6$species,list5$species,"number")    # 0           1, 2, 4, 6, 5 --> no

speciesCommonFun(list2$species,list5$species,"number")    # 0          1, 2, 5 --> no

speciesCommonFun(list2$species,list6$species,"number")    # 17          1, 2, 6 --> possibly
speciesCommonFun(list4$species,list6$species,"number")    # 5           1, 2, 6, 4 --> possibly
speciesCommonFun(list4$species,list3$species,"number")    # 1           1, 2, 6, 4, 3 --> ok but not top
speciesCommonFun(list4$species,list5$species,"number")    # 0           1, 2, 6, 4, 5 --> no
speciesCommonFun(list3$species,list5$species,"number")    # 7           1, 2, 6, 4, 3, 5 -->              POSSIBILITY: but not top

listLists <- list(list1$species,list2$species,list3$species,list4$species,list5$species,list6$species)

speciesCommonMatrixFun(listLists)
      [,1] [,2] [,3] [,4] [,5] [,6]
[1,]    0    8    0    4    0   14
[2,]    8    0    0    3    0   17
[3,]    0    0    0    1    7    0
[4,]    4    3    1    0    0    5
[5,]    0    0    7    0    0    0
[6,]   14   17    0    5    0    0

## final order:  that's the only possible order so that at least one species is common between 2 consecutive lists
list1Final <- list1
list2Final <- list2
list3Final <- list6
list4Final <- list4
list5Final <- list3
list6Final <- list5

## transformation of these tables so that ranking score is place instead of a value from 0 to 1
## 1 --> high rank number (most aggressive), and 0 --> lowest: 1
rankL1 <- rankingFunction(list1Final)
rankL2 <- rankingFunction(list2Final)
rankL3 <- rankingFunction(list3Final)
rankL4 <- rankingFunction(list4Final)
rankL5 <- rankingFunction(list5Final)
rankL6 <- rankingFunction(list6Final)

## creation of a list that contains the rankLi:
listRankList <- list(rankL1,rankL2,rankL3,rankL4,rankL5,rankL6)  # this list is going to change

### Determination of the ranking scores for 10 iterations:
listRankList_10 <- iterativePartialRankingFun(listRankList,10)

# % of pair-associations whose hierarchy is reversed after iteration (change in equality is no considered as a change) (bigger is worth):
proportionReversedFun(listRankList_10)
        iter_1 iter_2 iter_3 iter_4 iter_5 iter_6 iter_7 iter_8 iter_9 iter_10
table_1   5.54   5.54   5.54   5.54   5.54   5.54   5.54   5.54   5.54    5.54
table_2   7.20   7.20   7.20   7.20   7.20   7.20   7.20   7.20   7.20    7.20
table_3   4.81   4.81   4.81   4.81   4.81   4.81   4.81   4.81   4.81    4.81
table_4  19.23  19.23  19.23  19.23  19.23  19.23  19.23  19.23  19.23   19.23
table_5   3.13   3.13   3.13   3.13   3.13   3.13   3.13   3.13   3.13    3.13
table_6   2.50   2.50   2.50   2.50   2.50   2.50   2.50   2.50   2.50    2.50
sum      42.41  42.41  42.41  42.41  42.41  42.41  42.41  42.41  42.41   42.41

# Kendall's tau between the original ranking and each ranking obtained after an iteration (bigger is better):
KendallTauFun(listRankList_10)
        iter_1 iter_2 iter_3 iter_4 iter_5 iter_6 iter_7 iter_8 iter_9 iter_10
table_1   0.44   0.44   0.44   0.44   0.44   0.44   0.44   0.44   0.44    0.44
table_2   0.85   0.84   0.84   0.83   0.83   0.83   0.83   0.83   0.83    0.83
table_3   0.34   0.34   0.33   0.33   0.33   0.33   0.33   0.33   0.33    0.33
table_4   0.62   0.62   0.62   0.62   0.62   0.62   0.62   0.62   0.62    0.62
table_5   0.95   0.94   0.93   0.93   0.93   0.93   0.93   0.93   0.93    0.93
table_6   0.93   0.93   0.93   0.93   0.93   0.93   0.93   0.93   0.93    0.93
sum       4.13   4.11   4.09   4.08   4.08   4.08   4.08   4.08   4.08    4.08

# coherence between the original rankings and each ranking obtained after an iteration (bigger is better):
coherenceFun(listRankList_10)
          iter_1    iter_2    iter_3    iter_4    iter_5    iter_6    iter_7    iter_8    iter_9   iter_10
table_1  25.96516  25.96466  25.96466  25.96466  25.96466  25.96466  25.96466  25.96466  25.96466  25.96466
table_2  32.94715  32.94739  32.94763  32.94788  32.94788  32.94788  32.94788  32.94788  32.94788  32.94788
table_3  75.99102  75.99092  75.99126  75.99126  75.99126  75.99126  75.99126  75.99126  75.99126  75.99126
table_4  12.89744  12.89744  12.89744  12.89744  12.89744  12.89744  12.89744  12.89744  12.89744  12.89744
table_5  26.92703  26.92781  26.92871  26.92871  26.92871  26.92871  26.92871  26.92871  26.92871  26.92871
table_6  15.87577  15.87577  15.87577  15.87577  15.87577  15.87577  15.87577  15.87577  15.87577  15.87577
sum     190.60356 190.60399 190.60547 190.60571 190.60571 190.60571 190.60571 190.60571 190.60571 190.60571

#**************************************************************************************************
# There is very little difference between iterations.
# I take iteration nb 5
#**************************************************************************************************

### creation of figures:
par(mar=c(4,3.5,0.5,0.5))
cex <- 1.2
iteration <- 5
positionLetter <- "topright"
positionStats <- "bottomleft"
layout(matrix(seq(1:6),nrow=2,ncol=3,byrow = T), widths=c(1.17,1,1),heights=c(1,1.17))
par(mfrow = c(2,3))
numberTables <- length(listRankList_10) - 1
kandall <- KendallTauFun(listRankList_10)
kandall <- kandall[,iteration]
# Ppa <- proportionReversedFun(listRankList_10)
# Ppa <- Ppa[,iteration]
nameList <- c("A","B","C","D","E","F")
for(i in 1:numberTables){
  y <- listRankList_10[[i]][,2*(iteration + 1)]
  x <- listRankList_10[[i]][,2]
  n <- length(listRankList_10[[i]][,1])
  kandallList <- kandall[i]
  # PpaList <- Ppa[i]
  plot(jitter(y)~ x,las = 1,xlab="",ylab="",col = "darkgrey",cex = 1.5,lwd = 1.5)
  mtext("Original ranking",side=1,line=2.1,cex=0.9)
  mtext("Iterated ranking",side=2,line=2.2,cex=0.9)
  ymax <- max(y)
  xmax <- max(x)
  lines(x=c(0,xmax),y=c(0,ymax))
  lm <- lm(y~x)
  abline(lm,lty=2)
  R2 <- as.character(round(as.numeric(summary(lm)[8]),2))
  legend("bottomright",c(paste("R2:",R2),paste("tau:",kandallList),paste("n:",n)),cex=1.2,adj = 0,bty = "n")
  legend("topleft",nameList[i],cex=1.5,bty = "n")
}

### comupted aggressivness:
aggressivenessComputed <- listRankList_10[[7]][,c(1,iteration + 1)]
aggressivenessComputed$species <- as.character(aggressivenessComputed$species)

### several rows only have genera name --> to remove:
speciesToCheck <- checkSpeciesNames(aggressivenessComputed$species, wd_Datasets = wd_Datasets)

speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)

speciesGeneraToRemove <- c("Agaricia spp.","Astreopora sp.","Fungia 3sp.","Gonjopora spp.","Lopbophyllia sp.","Madracis spp.","Millepora sp.","Siderastrea spp.","Montipora spp.","Oculina spp.","Porites 5sp.","Stephanocoenia spp.","Turbinaria sp.")
length(speciesGeneraToRemove) # 13

### remove these species 
aggressivenessComputed <- aggressivenessComputed[!aggressivenessComputed$species %in% speciesGeneraToRemove,]
length(aggressivenessComputed$species) # 132
sum(duplicated(aggressivenessComputed$species)) # 0
colnames(aggressivenessComputed)[2] <- "aggressiveness"

# write.csv(aggressivenessComputed,paste(wd_Datasets,"aggressiveness.csv",sep="/"),row.names = F)

### figure Final -----
### need to update names in each original list 
## import list of names to change:
speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
length(speciesNameChecked$nameSp) # 1786
# need to create listRankList
aggressivness <- read.csv(paste(wd_Datasets,"aggressiveness.csv",sep="/"),stringsAsFactors = T)
length(aggressivness$species) # 132

listRankList_updated <- listRankList

layout(matrix(seq(1:6),nrow=2,ncol=3,byrow = T), widths=c(1.13,1,1),heights=c(1,1.15))
letter <- 0
for(i in 1:length(listRankList)){
  # prep DF original ----
  print(names(listRankList_updated[[i]][2]))
  ### correct names in species_colonialityFinal_cut:
  speciesNameChecked_cut <- speciesNameChecked[speciesNameChecked$nameSp %in% listRankList[[i]]$species,]
  speciesToCheck <- checkSpeciesNames(listRankList[[i]]$species, wd_Datasets = wd_Datasets)
  # remove species names with only genera given --> there are not present in speciesNameChecked$nameSp
  speciesToRemoveInOriginalList <- speciesNotPresentFun(listRankList[[i]]$species,speciesNameChecked_cut$nameSp)
  listRankList_updated[[i]] <- listRankList_updated[[i]][listRankList_updated[[i]]$species %in% speciesNameChecked_cut$nameSp,]
  # species names if updated name is not already present, otherwise remove it
  if(length(speciesToCheck$nameSp) > 0){
    for(j in 1:length(speciesToCheck$nameSp)){
      if(speciesPresentFun(speciesToCheck$nameSp_checked[j],listRankList_updated[[i]]$species)){
        print(paste(speciesToCheck$nameSp[j],"is removed because updated name is allready present"))
        listRankList_updated[[i]] <- listRankList_updated[[i]][listRankList_updated[[i]]$species != speciesToCheck$nameSp[j],]
      }else{
        listRankList_updated[[i]]$species[listRankList_updated[[i]]$species == speciesToCheck$nameSp[j]] <- speciesToCheck$nameSp_checked[j]
      }
    }
  }
  if(sum(duplicated(listRankList_updated[[i]]$species)) > 0){
    print("PROBLEMMMMMM: there are duplicaed species after updating for species name")
  }
  # update ranking order 
  listRankList_updated[[i]] <- listRankList_updated[[i]][order(listRankList_updated[[i]][,2]),]
  listRankList_updated[[i]]$Original_updated <- NA
  listRankList_updated[[i]]$Original_updated[1] <- 1
  for(j in 2:length(listRankList_updated[[i]]$species)){
    if(listRankList_updated[[i]][j,2] == listRankList_updated[[i]][j-1,2]){
      listRankList_updated[[i]]$Original_updated[j] <- listRankList_updated[[i]]$Original_updated[j-1]
    }else{
      listRankList_updated[[i]]$Original_updated[j] <- listRankList_updated[[i]]$Original_updated[j-1] + 1
    }
  }
  # prep inputed DF -----
  ### import the inputed ranking score
  aggressivness_cut <- aggressivness[aggressivness$species %in% listRankList_updated[[i]]$species,]
  if(length(aggressivness_cut$species) != length(listRankList_updated[[i]]$species)){
    print("PROBLEMMMMMM: different number of species in the originla updated list and the computed list")
  }
  # # update ranking order
  aggressivness_cut <- aggressivness_cut[order(aggressivness_cut$aggressiveness),]
  aggressivness_cut$aggressiveness_updated <- NA
  aggressivness_cut$aggressiveness_updated[1] <- 1
  for(j in 2:length(aggressivness_cut$species)){
    if(aggressivness_cut$aggressiveness[j] == aggressivness_cut$aggressiveness[j-1]){
      aggressivness_cut$aggressiveness_updated[j] <- aggressivness_cut$aggressiveness_updated[j-1]
    }else{
      aggressivness_cut$aggressiveness_updated[j] <- aggressivness_cut$aggressiveness_updated[j-1] + 1
    }
  }
  ### add the imputed ranking in listRankList_updated
  listRankList_updated[[i]] <- merge(listRankList_updated[[i]],aggressivness_cut,all.x = T, all.y = T)
  
  # plot ----
  x <- listRankList_updated[[i]]$Original_updated
  y <- listRankList_updated[[i]]$aggressiveness_updated
  bottom <- 2
  left <- 2
  if(i > 3){
    bottom <- 3.7
  }
  if(i == 1 | i == 4){
   left <- 3.5
  }
  par(mar=c(bottom,left,0.5,0.5))
  plot(jitter(y)~jitter(x), xlim = c(0,round(max(x)+max(x)/10,0)), ylim=c(0,max(y)),las=1,col="darkgrey",pch=1,cex=1.5,lwd=2,ylab="",xlab="")
  abline(a = 0, b = 1, lty = 2)
  if(i > 3){
    mtext("Original ranking",side=1,line=2.4,cex=1)
  }
  if(i == 1 | i == 4){
    mtext("Computed ranking",side=2,line=2.2,cex=1)
  }
  # Spearman ranking coefficient
  p <- correlationValFun(y,x,method="spearman",output="p.value",numberDecimals=2,p="stars") # 
  r <- correlationValFun(y,x,method="spearman",output="r",numberDecimals=2,p="stars") # 
  
  pfull <- correlationValFun(y,x,method="spearman",output="p.value",numberDecimals=2,p = "num") # 

  
  n <- length(x)
  legend("topleft",LETTERS[letter <- letter + 1],cex=1.5,bty = "n")
  print(paste(LETTERS[letter <- letter + 1],"Pvalue = ",pfull))
  
  legend("bottomright",c(eval(bquote(expression('r'[s]*' = ' ~ .(r) ~ .(p)))),paste("n = ",n)), cex=1.5,bty = "n")
# ----
}

for(i in 1:length(listRankList)){
  print(paste(names(listRankList[[i]][2]),length(listRankList[[i]]$species),"vs",length(listRankList_updated[[i]]$species)))
}

