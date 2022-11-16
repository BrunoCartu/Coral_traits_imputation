# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# Goal:  Producing the functions used to manipulate the trait and phylogenetic data.

require(ape)
require(apTreeshape)
require(ape)
require(vegan)
require(picante)
require(phangorn)
require(PHYLOGR)

## methods that compares 3 vectors of names and create 3 vectors names species names places in a list: listCoral[[1]]: those that are reefspecies, [[2]]: non-reef species and [[3]] those that are not in the both of the 2 1st lists
nameListCompareasonFun <- function(listReefSpecies = character(),listNonReefSpecies = character(),listSpeciesFT = character()){
  
  FT_table_all <- listSpeciesFT
  reefSpecies <- listReefSpecies
  NonReefSpecies <- listNonReefSpecies
  
  numberSpeciesInFDTable <- length(FT_table_all)
  numberSpeciesInNoneReefSpecies <- length(NonReefSpecies)
  numberSpeciesInReefSpecies <- length(reefSpecies)
  listCoral <- list()
  reefSpeciesVec <- character()
  nonreefSpeciesVec <- character()
  notInTreeSpeciesVec <- FT_table_all
  reefSpeciesNum <- 0
  nonreefSpeciesNum <- 0
  for(i in 1:numberSpeciesInFDTable){
    for(j in 1:numberSpeciesInReefSpecies){
      if(FT_table_all[i] == reefSpecies[j]){
        reefSpeciesNum <- reefSpeciesNum + 1
        reefSpeciesVec[reefSpeciesNum] <- FT_table_all[i]
        notInTreeSpeciesVec[i] <- NA
        break
      }
    }
    for(k in 1: numberSpeciesInNoneReefSpecies){
      if (FT_table_all[i] == NonReefSpecies[k]){
        nonreefSpeciesNum <- nonreefSpeciesNum + 1
        nonreefSpeciesVec[nonreefSpeciesNum] <- FT_table_all[i]
        notInTreeSpeciesVec[i] <- NA
        break
      }
    }
  }
  notInTreeSpeciesVec <- subset(notInTreeSpeciesVec,!is.na(notInTreeSpeciesVec))
  listCoral[[1]] <- reefSpeciesVec
  listCoral[[2]] <- nonreefSpeciesVec
  listCoral[[3]] <- notInTreeSpeciesVec
  listCoral
}

## methods that return FALSE if the name of a species entered as a 1st paramter is not in the list of names entered as a 2nd parameter:
speciesPresentFun <- function(nameSpecies = character(1),vectNames = character()){
  numberSpeciesVectNames <- length(vectNames)
  found <- FALSE
  for(i in 1:numberSpeciesVectNames){
    if(vectNames[i] == nameSpecies){
      found <- TRUE
      break
    }
  }
  found
}

## method that takes 2 vectors of names and returns a character() of names that are present in both vectors (opposit of speciesNotPresentFun)
speciesPresentFunBis <- function(listnames1 = character(), listnames2 = character()){
  numberSpeicesLN1 <- length(listnames1)
  numberSpeicesLN2 <- length(listnames2)
  namesToBeReturned <- character()
  numberNonCommonSpecies <- 0
  for(i in 1:numberSpeicesLN1){
    present <- 0
    for(j in 1:numberSpeicesLN2){
      if(listnames1[i] == listnames2[j]){
        present <- 1
        break
      }
    }
    if(present == 1){
      numberNonCommonSpecies <- numberNonCommonSpecies + 1
      namesToBeReturned[numberNonCommonSpecies] <- listnames1[i]
    }
  }
  namesToBeReturned
}


## method that takes 2 vectors of names and returns a character() of names that are present in the 1st list but not in the second:
speciesNotPresentFun <- function(listnames1 = character(), listnames2 = character()){
  numberSpeicesLN1 <- length(listnames1)
  numberSpeicesLN2 <- length(listnames2)
  namesToBeReturned <- character()
  numberNonCommonSpecies <- 0
  for(i in 1:numberSpeicesLN1){
    present <- 0
    for(j in 1:numberSpeicesLN2){
      if(listnames1[i] == listnames2[j]){
        present <- 1
        break
      }
    }
    if(present == 0){
      numberNonCommonSpecies <- numberNonCommonSpecies + 1
      namesToBeReturned[numberNonCommonSpecies] <- listnames1[i]
    }
  }
  namesToBeReturned
}

## method that takes 2 vectors of names and returns a character() of names that are present in both vectors (opposit of speciesNotPresentFun)
speciesPresentBothList <- function(listnames1 = character(), listnames2 = character()){
  numberSpeicesLN1 <- length(listnames1)
  numberSpeicesLN2 <- length(listnames2)
  namesToBeReturned <- character()
  numberNonCommonSpecies <- 0
  for(i in 1:numberSpeicesLN1){
    present <- 0
    for(j in 1:numberSpeicesLN2){
      if(listnames1[i] == listnames2[j]){
        present <- 1
        break
      }
    }
    if(present == 1){
      numberNonCommonSpecies <- numberNonCommonSpecies + 1
      namesToBeReturned[numberNonCommonSpecies] <- listnames1[i]
    }
  }
  namesToBeReturned
}

## method that checks if strings in 2 vectors are indentical for a same position in the vector (the vector are of same length)
# and return a list of names present in the 1st vectors that differ in the same position in the second vector
# in case of a difference, the name in the 2nd vector is "valide" as opposed to the one in the 1st vect "not valide" 
speciesSameLineIdentical <- function(vecString1 = character(),vecString2 = character()){
  numLines1 <- length(vecString1)
  numLines2 <- length(vecString2)
  if(numLines1 != numLines1){
    print("Vectors not same length")
  }
  nameCorrected <- as.data.frame(matrix(NA,nrow=numLines1,ncol=2))
  colnames(nameCorrected) <- c("not_valide","valide")
  x <- 0
  for(i in 1:numLines1){
    if(vecString1[i] != vecString2[i]){
      x <- x + 1
      nameCorrected[x,1] <- vecString1[i]
      nameCorrected[x,2] <- vecString2[i]
    }
  }
  if(x ==0){
    print("no difference")
    nameCorrected <- NA
  }else{
    print(paste("there are ", x," differences"))
    nameCorrected <- na.omit(nameCorrected)
  }
  nameCorrected
}

### function that returns a 1 row dataframe with the mean, the sd and the n from another dataframe with same column (IN THIS ORDER) but each raw corresponding to different samples:
# treats single observations as well
meanSDnFund <- function(meanSDn = data.frame()){
  
  output <- data.frame(mean = NA,sd = NA,n = NA)
  ### table for when n > 1:
  nSup1_sub <- subset(meanSDn,meanSDn[,3] > 1)
  colnames(nSup1_sub) <- c("mean","sd","n")
  ### table for when n = 1:
  nEqu1_sub <- subset(meanSDn,meanSDn[,3] == 1)
  colnames(nEqu1_sub) <- c("mean","sd","n")
  
  ### if there are rows with sd values and n > 1:
  if(length(row.names(nSup1_sub)) > 0){
    numberSamples <- length(nSup1_sub[,1])
    nTot <- sum(nSup1_sub[,3])
    output$n <- nTot
    ### calculation of the mean of the means:
    weigthedMean <- c()
    for(i in 1:numberSamples){
      weigthedMean[i] <- nSup1_sub[i,1]*nSup1_sub[i,3]
    }
    output$mean <- sum(weigthedMean)/output$n
    ### cacluation sd of the sds:
    # see http://stats.stackexchange.com/questions/55999/is-it-possible-to-find-the-combined-standard-deviation for formula
    operations <- c()
    for(i in 1:numberSamples){
      operations[i] <- nSup1_sub[i,3]*nSup1_sub[i,2]*nSup1_sub[i,2] + nSup1_sub[i,3]*(nSup1_sub[i,1] - output$mean)*(nSup1_sub[i,1] - output$mean)
    }
    sumOperations <- sum(na.omit(operations))
    output$sd <- sqrt(sumOperations/output$n)
    
    ### if there ALSO also rows with single observations (n == 1), then add these observation incrementally (http://stats.stackexchange.com/questions/105773/how-to-calculate-sd-of-sample-for-one-new-observation)
    if(length(row.names(nEqu1_sub)) > 0){
      numberSingleObs <- length(row.names(nEqu1_sub))
      for(i in 1: numberSingleObs){
        output$n <- output$n + 1
        previousMean <- output$mean
        output$mean <- (nEqu1_sub$mean[i] + (output$n-1)*output$mean) / output$n
        output$sd <- sqrt(((output$n-1)*output$sd^2 + output$n*(previousMean - output$mean)^2 + (nEqu1_sub$mean[i] - output$mean)^2)/output$n)  # http://math.stackexchange.com/questions/102978/incremental-computation-of-standard-deviation/103025#103025
      }
    }
  }else{  # the dataset only has unic observations:
    output$mean <- mean(meanSDn[,1])
    output$sd <- sd(meanSDn[,1])
    output$n <- length(meanSDn[,1])
  }
  output
}

### Method that checks species names based on my personal list of species names cheked on WoRMS and coraltraits.org; returns a list of names that are not already present and need to be checked
checkSpeciesNames <- function(spNames = character(),wd_Datasets){
  # get dataset with species names checked:
  speciesNameChecked <- read.csv(paste(wd_Datasets,"speciesNameChecked.csv",sep="/"),header=T,stringsAsFactors = F)
  # species present in the 1st list and not in the 2nd:
  speciesNotPresentin_nameSp <- speciesNotPresentFun(spNames,speciesNameChecked$nameSp) # species not present in speciesNameChecked$nameSp
  speciesPresentin_nameSp_checked <- speciesPresentFunBis(spNames,speciesNameChecked$nameSp_checked) # species present in speciesNameChecked$nameSp_checked
  speciesToCheck <- speciesNotPresentin_nameSp[!speciesNotPresentin_nameSp %in% speciesPresentin_nameSp_checked] # species not present in spNames but present in nameSp_checked need to be remove because name is correct
  numberspeices <- length(speciesToCheck)
  if(numberspeices == 0){
    print("All species are present in the list")
  }else{
    print(paste("there are: ",numberspeices," species not present in the list:"))
    for(i in 1:length(speciesToCheck)){
      print(speciesToCheck[i])
    }
    
  }
  # print the species whose name need to be updated:
  speciesToUpdate <- speciesNotPresentFun(spNames,speciesNameChecked$nameSp_checked)
  correctNames <- speciesNameChecked[speciesNameChecked$nameSp %in% spNames,]
  correctNames <- correctNames[! correctNames$nameSp %in% correctNames$nameSp_checked,]
  if(length(speciesToUpdate) == 0){
    print("All names are good")
  }else{
    print("Here are the species that needs an updated name:")
    print(correctNames)
  }
  correctNames
}

### function that takes a vector of species names with Genus and species and return a dataframe with the genus and species names in different column
# it allows to chose the type of sepator used in the initial table
speciesGenusSplitFun <- function(vectNames = charactcer(),separator = " "){
  numberSpecies <- length(vectNames)
  DF <- data.frame(do.call("rbind",strsplit(vectNames,separator,fixed = TRUE)))
  colnames(DF) <- c("Genus","Species")
  DF
}

### Modified version of Jason Pither BIOL202 stats courses' stripchart.CI.fun
stripchart.CI.fun <- function(yvar, group.var, bars=c("confidence","SE"),confidence.lev = 0.95, groupnames=levels(group.var),
                              las.graph = 1,ylab.char="",xlab.char="",vertical=T,ylabVertiF=10,Yaxt="s",Xaxt="s",ylimit=NULL,cex.axis=1){
  options(warn=-1);
  # "groupnames" is a vector of labels for the groups
  # "bars" tells it which type of error bars to use, confidence or +/- one standard error (SE)
  
  # calcualte group means
  group.means <- tapply(yvar,group.var,mean,na.rm=T)
  # calculate group sizes
  group.size <- tapply(yvar,group.var,function(x){length(na.omit(x))})
  group.size[is.na(group.size)] <- 0
  # calculate confidence limits for group means
  lower.cls <- numeric()
  upper.cls <- numeric()
  for(i in 1:length(group.size)){
    if(group.size[i] == 0 | group.size[i] == 1){
      lower.cls[i] <- NA
      upper.cls[i] <- NA
    }else{
      x <- yvar[group.var == groupnames[i]]
      lower.cls[i] <-  t.test(x,conf.level=confidence.lev)$conf[1]
      upper.cls[i] <-  t.test(x,conf.level=confidence.lev)$conf[2]
    }
  }
  
  # lower.cls <-  tapply(yvar,group.var,function(x){t.test(x,conf.level=confidence.lev)$conf[1]})
  # upper.cls <-  tapply(yvar,group.var,function(x){t.test(x,conf.level=confidence.lev)$conf[2]})
  
  # calculate standard error of the mean:
  se.vals <- tapply(yvar,group.var,function(x){sd(x)/sqrt(length(x))})
  
  # calculate an appropriate offset distance for the confidence bands
  numberGroupVar <- length(levels(group.var))
  offset.by <- .4/numberGroupVar
  
  # determine offset of ylimit in case not null
  if(!is.null(ylimit)){
    offset.y <- (max(ylimit) - min(ylimit))/10
    ylimit <- c(min(ylimit)-offset.y,max(ylimit)+offset.y)
  }
  
  # now create main figure
  if(vertical == T){
    stripchart(yvar~group.var,method="jitter",group.names=groupnames,
               pch=1,cex=1.2,lwd=1.5,col="darkgrey",
               xlim=c(0.5,length(levels(group.var))+0.5), 
               ylim=ylimit,cex.axis = cex.axis,
               las=las.graph,vertical=vertical,
               main="",yaxt=Yaxt,xaxt=Xaxt,
               ylab=ylab.char, xlab=xlab.char);
    if(Yaxt == "n"){
      axis(side=2,labels=F) 
    }
    if(Xaxt == "n"){
      axis(side=1,labels=F,seq(1,numberGroupVar,1))   # TO CHECK IF NOT REVERSED WITH Y  
    }
  } else{
    stripchart(yvar~group.var,method="jitter",group.names=groupnames,
               pch=1,cex=1.2,lwd=1.5,col="darkgrey",
               ylim=c(0.5,length(levels(group.var))+0.5), 
               xlim=ylimit,cex.axis = cex.axis,
               las=las.graph,vertical=vertical,
               main="",yaxt=Yaxt,xaxt=Xaxt,
               ylab="", xlab= ylab.char);
    mtext(xlab.char, 2, line = ylabVertiF, cex = 1.2)
    if(Yaxt == "n"){
      axis(side=2,labels=F,seq(1,numberGroupVar,1)) 
    }
    if(Xaxt == "n"){
      axis(side=1,labels=F)  
    }
  }
  
  # loop through each group and add associated group means and confidence intervals
  for (i in 1:numberGroupVar){
    
    if(bars=="confidence"){
      if(vertical == T){
        X0 <- i+offset.by
        Y0 <- lower.cls[i]
        X1 <- i+offset.by
        Y1 <- upper.cls[i]
      } else{
        X0 <- lower.cls[i] 
        Y0 <- i+0.2
        X1 <- upper.cls[i] 
        Y1 <- i+0.2
      }
      if(group.size[i] > 2){
        arrows(X0,Y0,X1,Y1,code=3, length=0.1,angle=90,lwd=1.5)
      }
    }else {
      if(vertical == T){
        X0 <- i+offset.by
        Y0 <- group.means[i]-se.vals[i]
        X1 <- i+offset.by
        Y1 <- group.means[i]+se.vals[i]
      }else{
        X0 <- group.means[i]-se.vals[i]
        Y0 <- i+0.2
        X1 <- group.means[i]+se.vals[i]
        Y1 <- i+0.2
      }
      if(group.size[i] > 2){
        arrows(X0,Y0,X1,Y1,code=3,length=0.1,angle=90,lwd=1.5)
      }
    }
    if(vertical == T){
      X <- i+offset.by
      Y <- group.means[i]
    }else{
      X <- group.means[i] 
      Y <- i+0.2
    }
    if(group.size[i] > 1){
      points(X, Y,col="black",pch=16,cex=1.4)
    }
  }
  options(warn=1);
}

## method that changes
changeCharacterTipLabelPhyloTreefun <- function(phyloTree = multiPhylo(),charToRemove = character(),charToPlace = character()){
  numberTree <-  length(phyloTree)
  numberSpecies <-  length(phyloTree[[1]]$tip.label)
  for(i in 1:numberTree){
    for(j in 1 : numberSpecies){
      phyloTree[[i]]$tip.label[j] <- gsub(charToRemove,charToPlace,phyloTree[[i]]$tip.label[j])
    }
  }
}

### Function that go get the trait csv files and creates a list of these. IN PROGRESS
barFun <- function(xvar = c("age_at_maturity","aggressiveness","coral_bleaching_index","chlorophyll_a_concentration","coloniality","colony_max_diameter","corallite_area","dark_respiration",
                            "egg_diameter","fecundity_polyp","growth_form","growth_rate","lipid_content","mode_larval_development","reduced_scattering_coefficient","size_at_maturity","skelelal_density","symbiodinium_density",
                            "tissue_thickness","zooxanthellate"), 
                   yvar = c("age_at_maturity","aggressiveness","coral_bleaching_index","chlorophyll_a_concentration","coloniality","colony_max_diameter","corallite_area","dark_respiration",
                            "egg_diameter","fecundity_polyp","growth_form","growth_rate","lipid_content","mode_larval_development","reduced_scattering_coefficient","size_at_maturity","skelelal_density","symbiodinium_density",
                            "tissue_thickness","zooxanthellate"),
                   species = character(), x_log = F, y_log = F){
  
  listTraits <- traitsListFun(wd)
  listTraits_cut <- listTraits[c(xvar,yvar)]
  dfx <-  listTraits_cut[[1]][listTraits_cut[[1]]$species %in% species,]
  dfy <-  listTraits_cut[[2]][listTraits_cut[[2]]$species %in% species,]
  
  numberSpecies <- length(species)
  segmentCoordinates <- data.frame(
    x0 = rep(NA,numberSpecies),
    x1 = rep(NA,numberSpecies),
    y0 = rep(NA,numberSpecies),
    y1 = rep(NA,numberSpecies))
  
  # xvar corresponds to a numerical trait that has a sd and n:
  x_numeric_var_bool <- xvar %in% c("chlorophyll_a_concentration","dark_respiration","egg_diameter","fecundity_polyp","growth_rate","lipid_content","reduced_scattering_coefficient",
                                    "skelelal_density","symbiodinium_density","tissue_thickness")
  
  x_coralliteArea_bool <- xvar == "corallite_area"
  y_categorical_bool <- yvar %in% c("coloniality","mode_larval_development","zooxanthellate")
  y_numeric_var_bool <- yvar %in% c("chlorophyll_a_concentration","dark_respiration","egg_diameter","fecundity_polyp","growth_rate","lipid_content","reduced_scattering_coefficient",
                                    "skelelal_density","symbiodinium_density","tissue_thickness")
  y_growth_form_bool <- yvar == "growth_form"
  
  ### x variable:
  if(x_numeric_var_bool){
    if(x_log){
      dfx[,2:3] <- log(dfx[,2:3])
    }
    # x0:
    segmentCoordinates[,1] <- dfx[,2] - dfx[,3] / sqrt(dfx[,4])
    segmentCoordinates[,1][is.na(segmentCoordinates[,1])] <- dfx[,2][is.na(dfx[,3])]
    # x1:
    segmentCoordinates[,2] <- dfx[,2] + dfx[,3] / sqrt(dfx[,4])
    segmentCoordinates[,2][is.na(segmentCoordinates[,2])] <- dfx[,2][is.na(dfx[,3])]
  }else if (x_coralliteArea_bool){
    if(x_log){
      dfx[,2:4] <- log(dfx[,2:4])
    }
    # x0:
    segmentCoordinates[,1] <- dfx[,3]
    segmentCoordinates[,2] <- dfx[,4]
    segmentCoordinates[,1][is.na(segmentCoordinates[,1])] <- dfx[,2][is.na(dfx[,3])]
    segmentCoordinates[,2][is.na(segmentCoordinates[,2])] <- dfx[,2][is.na(dfx[,4])]
  }
  ### y variable:
  if(y_growth_form_bool){
    dfy$growthFormClasses <- ordered(dfy$growthFormClasses,levels=rev(c("branching","tables_or_plates","digitate","corymbose","laminar","columnar","massive","encrusting_long_uprights","encrusting")))
    numberCategories <- length(levels(as.factor(dfy[,3]))) # 9
    #y0
    segmentCoordinates[,3] <- as.integer(dfy$growthFormClasses)
    segmentCoordinates[,4] <- segmentCoordinates[,3]
  }else if(y_numeric_var_bool){
    if(y_log){
      dfy[,2:3] <- log(dfy[,2:3])
    }
    # x0:
    segmentCoordinates[,1] <- dfy[,2] - dfy[,3] / sqrt(dfy[,4])
    segmentCoordinates[,1][is.na(segmentCoordinates[,1])] <- dfy[,2][is.na(dfy[,3])]
    # y1:
    segmentCoordinates[,2] <- dfy[,2] + dfy[,3] / sqrt(dfy[,4])
    segmentCoordinates[,2][is.na(segmentCoordinates[,2])] <- dfy[,2][is.na(dfy[,3])]
  }else if(Y_colony_max_diameter){
    if(y_log){
      dfy[,2] <- log(dfy[,2])
    }
    segmentCoordinates[,3] <- dfy[,2]
    segmentCoordinates[,4] <- dfy[,2]
  }
  
  for(i in 1:numberSpecies){
    segments(x0 = segmentCoordinates[i,1],y0 = segmentCoordinates[i,3],x1 = segmentCoordinates[i,2],y1 = segmentCoordinates[i,4],col="darkgrey",lwd=1.5,lty = 1)
  }
  
}

## method that read.nexus a tree file and retun a multiphylo object with tip.lables that can be changed
read.nexus.fun <- function(nameFileSuperTree = character()){
  trees <- read.nexus(file  = nameFileSuperTree)
  trees_readable <-  rtreeshape(n=length(trees),tip.number=length(trees[[1]]$tip.label),model="yule")
  numberTrees <- length(trees)
  for (i in 1:numberTrees){
    trees_readable[[i]] <- trees[[i]]
  }
  class(trees_readable) <- "multiPhylo"
  trees_readable
}

# function that calculate the mean, sd and SE of a distance matrix. 
distMatrixSumStats <- function(myMatrix = matrix){
  distances <- c()
  ndistances <- 0
  sd.wRF_dis <- 0
  mean.wRF_dist.2 <- 0
  for(i in 1:(numberTrees-1)){
    jstart <- i + 1
    for(j in jstart:numberTrees){
      ndistances <- ndistances + 1
      distances[ndistances] <- myMatrix[i,j]
    }
    sd.wRF_dis <- sd(distances)    #  323.4839
    mean.wRF_dist.2 <- mean(distances)  # 5627.14
    SE.wRF_dis <- sd.wRF_dis / sqrt(numberTrees*(numberTrees-1)/2) # 0.4577042
  }
  distMatrixSummary <- data.frame(
    mean= mean.wRF_dist.2,
    sd= sd.wRF_dis,
    SE=SE.wRF_dis
  )
}

# This function tests for phylogenetic signal with categorical traits. CHECK R SCRIPT JSI.R. 
phylo.signal.disc <- function(trait,phy,rep = 999,cost=NULL)  {
  lev <- attributes(factor(trait))$levels
  if (length(lev) == length(trait)) 
    stop("Are you sure this variable is categorical?")
  if(is.null(cost)){ 
    cost1 <- 1-diag(length(lev))
  }
  else {
    if (length(lev) != dim(cost)[1]) 
      stop("Dimensions of the character state transition matrix do not agree with the number of levels")
    cost1<- t(cost)
  }
  dimnames(cost1) <- list(lev,lev)		
  trait <- as.numeric(trait)
  attributes(trait)$names <- phy$tip.label
  NULL.MODEL <- matrix(NA,rep,1)
  obs <- t(data.frame(trait))
  obs <- phyDat(t(obs),type="USER",levels=attributes(factor(obs))$levels)
  OBS <- parsimony(phy,obs,method="sankoff",cost=cost1)	
  for (i in 1:rep){
    null <- sample(as.numeric(trait))
    attributes(null)$names <- attributes(trait)$names  
    null <- t(data.frame(null))
    null <- phyDat(t(null),type="USER",levels=attributes(factor(null))$levels)
    NULL.MODEL[i,]<-parsimony(phy,null,method="sankoff",cost=cost1)
    P.value <- sum(OBS >= NULL.MODEL)/(rep + 1)
  }
  
  #hist(NULL.MODEL,xlab="Transitions.in.Randomizations",xlim=c(min(c(min(NULL.MODEL,OBS-1))),max(NULL.MODEL)+1), main="",
  #     lwd=2)
  #mtext("(b)", side=3, font=4, line=1, adj=0, cex=1.5)
  #arrows(OBS,rep/10,OBS,0,angle=20,col="red",lwd=4)
  phy$tip.label <- rep(".",length(trait))
  
  OUTPUT1 <- t(data.frame(Number.of.Levels = length(attributes(factor(trait))$levels), Evolutionary.Transitions.Observed=OBS,Evolutionary.Transitions.Randomization.Median=median(NULL.MODEL),Evolutionary.Transitions.Randomization.Min=min(NULL.MODEL),Evolutionary.Transitions.Randomization.Max=max(NULL.MODEL),P.value))
  
  if(is.null(cost)){ 
    list(.Randomization.Results=OUTPUT1,.Levels= lev,.Costs.of.character.state.transition.UNORDERED.PARSIMONY = t(cost1))
  }
  else {
    list(.Randomization.Results=OUTPUT1,.Levels= lev,.Costs.of.character.state.transition.FROM.ROW.TO.COL = t(cost1))		}
}

### function that take the fucntion trait table of categorical traits and a multiphylo object and run the phylo.signal.disc() command on each tree 
# and then output a dataframe with n p values (n being the number of trees) for each trait. 
# It deals with NAs
cat_trait_signal.fun <- function(multiphylo,catTraitTable){
  numberTrees <- length(multiphylo) # 1000
  numberTraits <- length(FFT.categ) # 3
  cat_trat_signal <- matrix(nrow=numberTrees,ncol=numberTraits)
  colnames(cat_trat_signal) <- colnames(catTraitTable)
  for(i in 1 : numberTraits){
    FFT.temp <- subset(FFT.categ,!is.na(catTraitTable[,i]))
    speciesToKeep <- rownames(FFT.temp)
    speciesToRemove <- rownames(subset(FFT.categ,is.na(catTraitTable[,i])))
    
    phylo <- multiphylo
    for(j in 1 : numberTrees){
      phylo[[j]] <- drop.tip(multiphylo[[j]],speciesToRemove)
      
      phylo.signal <- phylo.signal.disc(FFT.temp[phylo[[j]]$tip.label,i],phylo[[j]])
      cat_trat_signal[j,i] <- phylo.signal$.Randomization.Results[6]
    }
  }
  cat_trat_signal
}

### perform a logit transformation, and potentially implement a epsilon value in case of negative value inside log()
logitTrans <- function(dataProportion = numeric(),epsilon = 0,dataInPercentage = T){
  
  # conversion from % to proportion
  if(dataInPercentage){
    dataProportion <- dataProportion / 100
  }
  
  # cases where epsilon need to be different from 0: proportion = 0 or 1
  number_zero <- sum(dataProportion == 0)
  number_one <- sum(dataProportion == 1)
  
  if(number_zero > 0 || number_one > 0){
    epsilon <- sort(dataProportion)[1]
    if(epsilon == 0){
      j <- 1
      while(epsilon == 0){
        j <- j + 1
        epsilon <- sort(dataProportion)[j]
      }
    }
    print(paste("Epsilon = ",epsilon))
  }
  
  logitTransResult <- log((dataProportion + epsilon)/(1 - dataProportion + epsilon))
  logitTransResult
}

logit.inverseTrans <- function(x){
  y <- exp(x)/(1+exp(x))
  y
}

### fucnton that performs a regression anaylis
regressionFun <- function(yvar,xvar,colorPoint="darkgrey",colorSeg="black",ylabel=character(),xlabel=character(),
                          displayRegression=T,displayAssumptionTest=F,categoryvar=NA,main="",Yaxt="s",Xaxt="s"){
  
  lm_temporal <- lm(yvar~xvar)
  
  if(displayRegression==T){
    
    layout(matrix(seq(1:1),nrow=1,ncol=1,byrow = T), widths=c(1),heights=c(1))
    
    if(!is.na(categoryvar)){
      cat("Result regression for ****** ",categoryvar," *******")
    }else{
      print("Result regression")
    }
    print(summary(lm_temporal))
    
    print(confint(lm_temporal))     # 95% confidence interval of the slope
    
    maxX <- max(xvar)
    minX <- min(xvar)
    maxY <- max(yvar)
    minY <- min(yvar)
    offsetX <- (maxX - minX)/10
    offsetY <- (maxY - minY)/10
    plot(yvar~xvar,col=colorPoint,las=1,xlab=xlabel,ylab=ylabel,cex=1.5,pch=1,lwd=1.5,xlim=c(minX-offsetX,maxX+offsetX),ylim=c(minY-offsetY,maxY+offsetY),main=main,yaxt=Yaxt,xaxt=Xaxt)
    if(Yaxt == "n"){
      axis(side=2,labels=F) 
    }
    if(Xaxt == "n"){
      axis(side=1,labels=F)   # TO CHECK IF NOT REVERSED WITH Y  
    }
    
    ### for segements:
    min <- min(xvar)-offsetX
    max <- max(xvar)+offsetX
    extrems <- data.frame(xvar=c(min,max))
    extrems$predicted_values <-  predict(lm_temporal,newdata = extrems)
    
    segments(x0=extrems[1,1], y0=extrems[1,2], 
             x1=extrems[2,1], y1=extrems[2,2],
             col = colorSeg, lty = 1, lwd = 4)
    
    ### for condidence bands:
    maxValues <- data.frame(xvar=seq(from=min,to=max,length=10))
    confi.lim <- as.data.frame(predict(lm_temporal,maxValues,level=0.95,interval="confidence"))
    lines(cbind(maxValues,confi.lim$lwr),col=colorSeg,lty="dashed",lwd=2)
    lines(cbind(maxValues,confi.lim$upr),col=colorSeg,lty="dashed",lwd=2)
  }
  
  if(displayAssumptionTest==T){
    #### Assumptions testing
    print("Assumptions testing")
    lm_temporal.resids <- residuals(lm_temporal)
    print(shapiro.test(lm_temporal.resids))
    
    layout(matrix(seq(1:2),nrow=1,ncol=2,byrow = T), widths=c(1,1),heights=c(1))
    plot(lm_temporal.resids~xvar,las=1,cex=1.2,lwd=2, xlab=xlabel, ylab="Residuals",main=main,yaxt=Yaxt,xaxt=Xaxt)
    abline(0,0,lty=2)
    qqnorm(lm_temporal.resids,las=1,main = "")
    qqline(lm_temporal.resids)
    
    if(Yaxt == "n"){
      axis(side=2,labels=F) 
    }
    if(Xaxt == "n"){
      axis(side=1,labels=F)   # TO CHECK IF NOT REVERSED WITH Y  
    }
  }
  sum <- summary(lm_temporal)
  sum
}

### function that performs a linear regressions between 2 numerical variables plotted on 1 only graph and test assumptions
regression2Fun <- function(yvar=numeric(),xvar=numeric(),catevar=character(),color=c("deepskyblue4","firebrick"),ylabel=character(),xlabel=character(),
                           displayRegression=T,displayAssumptionTest=T,legendPosition="topleft"){
  
  tableTempo <- data.frame(
    yvar = yvar,
    xvar = xvar,
    catevar = as.character(catevar),
    stringsAsFactors = F)
  
  categories <- levels(as.factor(as.character(tableTempo$catevar)))
  numberCategories <- length(categories)
  
  if(displayRegression==T){
    maxX <- max(xvar)
    minX <- min(xvar)
    maxY <- max(yvar)
    minY <- min(yvar)
    offsetX <- (maxX - minX)/20
    offsetY <- (maxY - minY)/20
    plot(0,0,xlim=c(minX-offsetX,maxX+offsetX),ylim=c(minY-offsetY,maxY+offsetY),col="white",las=1,xlab=xlabel,ylab=ylabel,cex=1.5,pch=1,lwd=1.5)
    # diplsay one regression after another
    for(i in 1:numberCategories){
      tableTempo_subset <- subset(tableTempo,tableTempo$catevar==categories[i])
      lm_temporal <- lm(yvar~xvar,data=tableTempo_subset)
      cat("Result regression for: ****** ",categories[i]," ********")
      print(summary(lm_temporal))
      print(confint(lm_temporal))     # 95% confidence interval of the slope
      numberSpeciesTemporal <- length(tableTempo_subset$yvar)
      # Plot the points
      for(j in 1:numberSpeciesTemporal){
        y <- tableTempo_subset$yvar[j]
        x <- tableTempo_subset$xvar[j]
        points(y=y,x=x,pch=1,cex=1.5,lwd=2,col=color[i])
      }
      ### for segements:
      min <- min(tableTempo_subset$xvar)-offsetX
      max <- max(tableTempo_subset$xvar)+offsetX
      extrems <- data.frame(xvar=c(min,max))
      extrems$predicted_values <-  predict(lm_temporal,newdata = extrems)
      segments(x0=extrems[1,1], y0=extrems[1,2], 
               x1=extrems[2,1], y1=extrems[2,2],
               col = color[i], lty = 1, lwd = 4)
      ### for condidence bands:
      maxValues <- data.frame(xvar=seq(from=min,to=max,length=10))
      confi.lim <- as.data.frame(predict(lm_temporal,maxValues,level=0.95,interval="confidence"))
      lines(cbind(maxValues,confi.lim$lwr),col=color[i],lty="dashed",lwd=2)
      lines(cbind(maxValues,confi.lim$upr),col=color[i],lty="dashed",lwd=2)
    } 
    legend(legendPosition, legend = categories, bty = "n",lwd = 2, cex = 1, col = color, lty = NA, pch = 1)
  }
  if(displayAssumptionTest==T){
    # diplsay one regression after another
    for(i in 1:numberCategories){
      tableTempo_subset <- subset(tableTempo,tableTempo$catevar==categories[i])
      lm_temporal <- lm(yvar~xvar,data=tableTempo_subset)
      #### Assumptions testing
      cat("Assumptions testing for ",categories[i]," *******")
      lm_temporal.resids <- residuals(lm_temporal)
      print(shapiro.test(lm_temporal.resids))
      # QQplot
      plot(lm_temporal.resids~tableTempo_subset$xvar,las=1,cex=1.2,lwd=2, xlab=X.label, ylab="Residuals")
      abline(0,0,lty=2)
    }
  }
}

### function that perform many linear regressions and display the output on a same splited screen:
multiple.regression.fun <- function(yvar=numeric(),xvar=numeric(),catevar=character(),outputRegression = T,assumptionTesting=F,
                                    color=charatcer(),Y.label=character(),X.label=character(),legendPosition="topleft"){
  
  dataFrameTempo <- data.frame(
    yvar=yvar,
    xvar=xvar,
    catevar=as.character(catevar),
    stringsAsFactors = F
  )
  categories <- levels(as.factor(dataFrameTempo$catevar))
  numberCategories <- length(categories)
  
  if(outputRegression){
    nRow <- 1
    nCol <- 2
    if(numberCategories == 3 || numberCategories == 4){
      nRow <- 2
    }else if(numberCategories == 5 || numberCategories == 6){
      nRow <- 3
    }else if(numberCategories == 7 || numberCategories == 8 || numberCategories == 9){
      nRow <- 3
      nCol <- 3
    }else if(numberCategories == 10 || numberCategories == 11 || numberCategories == 12){
      nRow <- 4
      nCol <- 3
    }
    par(mfrow=c(nRow,nCol))
    for(i in 1:numberCategories){
      dataFrameTempo_subset <- subset(dataFrameTempo,dataFrameTempo$catevar==categories[i])
      regressionFun(yvar=dataFrameTempo_subset$yvar,xvar=dataFrameTempo_subset$xvar,categoryvar=categories[i],
                    color=color[1],ylabel=Y.label,xlabel=X.label,main=categories[i],
                    displayRegression=T,displayAssumptionTest=F)
      #legend(legendPosition, legend = categories[i], bty = "n",lwd = 2, cex = 1, lty = NA)
      print(cat("sample size for ",categories[i]," :",length(dataFrameTempo_subset$yvar)," "))
    }
  }
  
  if(assumptionTesting){
    nRow <- 1
    nCol <- 2
    if(numberCategories == 3 || numberCategories == 4){
      nRow <- 2
    }else if(numberCategories == 5 || numberCategories == 6){
      nRow <- 3
    }else if(numberCategories == 7 || numberCategories == 8 || numberCategories == 9){
      nRow <- 3
      nCol <- 3
    }else if(numberCategories == 10 || numberCategories == 11 || numberCategories == 12){
      nRow <- 4
      nCol <- 3
    }
    par(mfrow=c(nRow,nCol))
    for(i in 1:numberCategories){
      print(cat("***** ASSUMPTION FOR ",categories[i]," *******"))
      dataFrameTempo_subset <- subset(dataFrameTempo,dataFrameTempo$catevar==categories[i])
      print(cat("sample size:",length(dataFrameTempo_subset$yvar)))
      regressionFun(yvar=dataFrameTempo_subset$yvar,xvar=dataFrameTempo_subset$xvar,categoryvar=categories[i],
                    color=color[1],ylabel=Y.label,xlabel=X.label,main=categories[i],
                    displayRegression=F,displayAssumptionTest=T)
      #legend(legendPosition, legend = categories[i], bty = "n",lwd = 2, cex = 1, lty = NA)
    } 
  }
}

### fucnton that performs a Pearson / Spearman correlation anaylisis
correlationFun <- function(yvar,xvar,colorPoint="darkgrey",colorSeg="black",ylabel=character(),xlabel=character(),
                          displayRegression=T,main="",Yaxt="s",Xaxt="s",pch = 1,xlimit=NULL,ylimit=NULL){
  
  if(is.null(xlimit)){
    maxX <- max(xvar)
    minX <- min(xvar)
  }else{
    maxX <- max(xlimit)
    minX <- min(xlimit)
  }
  if(is.null(ylimit)){
    maxY <- max(yvar)
    minY <- min(yvar)
  }else{
    maxY <- max(ylimit)
    minY <- min(ylimit)
  }

  offsetX <- (maxX - minX)/10
  offsetY <- (maxY - minY)/10
  
  if(displayRegression==T){
    lm_temporal <- lm(yvar~xvar)

    plot(yvar~xvar,col=colorPoint,las=1,xlab=xlabel,ylab=ylabel,cex=1.5,pch=pch,lwd=1.5,xlim=c(minX-offsetX,maxX+offsetX),ylim=c(minY-offsetY,maxY+offsetY),main=main,yaxt=Yaxt,xaxt=Xaxt)
    if(Yaxt == "n"){
      axis(side=2,labels=F) 
    }
    if(Xaxt == "n"){
      axis(side=1,labels=F)   # TO CHECK IF NOT REVERSED WITH Y  
    }
    
    ### for segements:
    min <- min(xvar)-offsetX
    max <- max(xvar)+offsetX
    extrems <- data.frame(xvar=c(min,max))
    extrems$predicted_values <-  predict(lm_temporal,newdata = extrems)
    
    segments(x0=extrems[1,1], y0=extrems[1,2], 
             x1=extrems[2,1], y1=extrems[2,2],
             col = colorSeg, lty = 1, lwd = 4)
    
  }else{
    plot(yvar~xvar,col=colorPoint,las=1,xlab=xlabel,ylab=ylabel,cex=1.5,pch=pch,lwd=1.5,xlim=c(minX-offsetX,maxX+offsetX),ylim=c(minY-offsetY,maxY+offsetY),main=main,yaxt=Yaxt,xaxt=Xaxt)
    if(Yaxt == "n"){
      axis(side=2,labels=F) 
    }
    if(Xaxt == "n"){
      axis(side=1,labels=F)   # TO CHECK IF NOT REVERSED WITH Y  
    }
  }
  print(cor.test(yvar, xvar, method="spearman"))
  print(cor.test(yvar, xvar, method="pearson"))
}

### function that return different values of correlation:
correlationValFun <- function(yvar,xvar,method="pearson",output=c("t_statistic","df","p.value","r","conf.int"),numberDecimals=2,p=c("num","stars")){
  correlation <- cor.test(yvar, xvar, method=method)
  value <- 0
  if(output == "t_statistic"){
    value <- correlation$statistic
  }else if(output == "df"){
    value <- correlation$parameter
  }else if(output == "p.value"){
    value <- correlation$p.value
  }else if(output == "r"){
    value <- correlation$estimate
  }else if(output == "conf.int"){
    value <- correlation$conf.int
  }
  print(method)
  print(output)
  print(value)
 
  if(output == "p.value" &&  p == "stars"){
    stars <- ""
    if(value < 0.05){
      stars <- "*"
    }
    if(value < 0.01){
      stars <- "**"
    }
    if(value < 0.001){
      stars <- "***"
    }
    stars
  }else{
    round(value,numberDecimals)
  }
  
}

### function that calculate the z score values of a numerical vector:
zscoreFun <- function(num=numeric()){
  z_score <- (num - mean(na.omit(num)))/sd(na.omit(num))
  z_score
}

### Function that creates a list containing all the traits csv data datasets, with the measure of varibility they have:
traitsListFun <- function(wd = character){
  
  wd_toSetBack <- getwd()
  setwd(paste(wd,"/Datasets",sep=""))
  
  listTraiTDatasets <- list()
  
  ### Import datasets:
  i <- 1
  listTraiTDatasets[[i]] <- read.csv("age_at_maturity.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("aggressiveness.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("chlorophyll_a_concentration.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("coloniality.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("colony_max_diameter.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("corallite_area.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("dark_respiration.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("egg_diameter.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("fecundity_polyp.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("growth_form.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("growth_rate.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("lipid_content.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("mode_larval_development.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("reduced_scattering_coefficient.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("response_bleaching_index.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("size_at_maturity.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("skelelal_density.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("symbiodimium_density.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("tissue_thickness.csv",header = T, stringsAsFactors = F)
  i <- i + 1
  listTraiTDatasets[[i]] <-read.csv("zooxanthellate.csv",header = T, stringsAsFactors = F)
  i # 20
  
  namesList <- c("age_at_maturity","aggressiveness","chlorophyll_a_concentration","coloniality","colony_max_diameter","corallite_area","dark_respiration",
                 "egg_diameter","fecundity_polyp","growth_form","growth_rate","lipid_content","mode_larval_development","reduced_scattering_coefficient","coral_bleaching_index","size_at_maturity","skelelal_density","symbiodinium_density",
                 "tissue_thickness","zooxanthellate")
  
  length(namesList) # 20
  
  names(listTraiTDatasets) <- namesList
  
  setwd(wd_toSetBack)
  
  listTraiTDatasets
}

### function that display the least square regression line on a pre-existing scatter plot:
linearRegressionLineFun <- function(yvar,xvar,colorSeg="black",lwd = 4,origin = F){
  # for margines
  maxX <- max(xvar)
  minX <- min(xvar)
  maxY <- max(yvar)
  minY <- min(yvar)
  lm_temporal <- lm(yvar~xvar)
  if(origin == T){
    lm_temporal <- lm(yvar~xvar - 1)
  }
  offsetX <- (maxX - minX)/10
  offsetY <- (maxY - minY)/10
  ### for segements:
  min <- min(xvar)-offsetX
  max <- max(xvar)+offsetX
  extrems <- data.frame(xvar=c(min,max))
  extrems$predicted_values <-  predict(lm_temporal,newdata = extrems)
  segments(x0=extrems[1,1], y0=extrems[1,2], 
           x1=extrems[2,1], y1=extrems[2,2],
           col = colorSeg, lty = 1, lwd = lwd)
}

### function that displat error bars on an existing graph, need to give the min and max values, center value and the center value for the opposit axis
barFun <- function(center = numeric(),upperL = numeric(),lowerL = numeric(),axis = c("x","y"),otherCoordValues = numerical(),colorSeg = "darkgrey",lwd = 2){
  numberSpecies <- length(center)
  for(i in 1:numberSpecies){
    max <- upperL[i]
    min <- lowerL[i]
    otherCoordValue <- otherCoordValues[i]
    #    if(is.na(max)){
    #     max <- center[i]
    #    }
    #    if(is.na(min)){
    #      min <- center[i]
    #    }
    if(axis == "x"){
      x <- center[i]
      y <- otherCoordValues[i]
      x0 <- min
      y0 <- y
      x1 <- max
      y1 <- y
    }else{
      x <- otherCoordValues[i]
      y <- center[i]
      x0 <- x
      y0 <- min
      x1 <- x
      y1 <- max
    }
    if(!is.na(max) & !is.na(min)){
      segments(x0=x0,y0=y0,x1=x1,y1=y1,col = colorSeg, lty = 1, lwd = lwd)
    }
    points(x=x,y=y,cex=1.5,lwd=1.5,lty = 1,col = colorSeg)
  }
}

### methods that take an vector of number as an interval and transforms the values so they are comprised between [0,1]
zeroTooneInterval.fun <- function(interval = numeric()){
  min <- min(interval)
  unitInteval <- interval
  unitInteval <- unitInteval - min 
  max <- max(unitInteval)
  unitInteval <- unitInteval / max
}


# Function that perform a PIC on 2 numerical functional traits. The yvar and xvar are numerical vectors NAMED (with species), and no NA, and same species for all arguments
# It calculates the mean R2 squares value and the SE and determine proporion of p > 0.05, <0.05, < 0.01, < 0.001 
# fail to reject H0: proportion(p > 0.05) > 0.8,
# Ha *** : proportion(p <= 0.001) > 0.8
# Ha **  : proportion(p <= 0.001) +  proportion(0.001 < p <= 0.01) > 0.8 AND !HA ***
# Ha *   : proportion(p <= 0.001) +  proportion(0.001 < p <= 0.01) + proportion(0.01 < p <= 0.05) > 0.8 AND !Ha **


### Function that calculates the Pearson or spearman r or the linear least sqaures regression R2 and p value for the number of trees in multiphyloTree
# then calculate the mean r or R2
# returns an approximate p value corresponding to the range when the proportion of p values > threshold; if none, return "sig. ?"

PIC.multiplylo.fun <- function(multiphyloTree,yvar,xvar,species,histo = F,model=c("lm","Pearson","Spearman"),threshold = 0.7){
  numberTree <- length(multiphyloTree)
  r_p_table <- as.data.frame(matrix(nrow = numberTree,ncol = 2))
  colnames(r_p_table) <- c("r_coeff","p")
  names(yvar) <- species
  names(xvar) <- species
  for(i in 1:numberTree){
    pic.y <- pic(yvar[multiphyloTree[[i]]$tip.label],multiphyloTree[[i]])
    pic.x <- pic(xvar[multiphyloTree[[i]]$tip.label],multiphyloTree[[i]])
    if(model == "lm"){
      pic.output <- lm(pic.y ~ pic.x - 1)
      p <- summary(pic.output)$coefficients[,4]
      r_coeff <- summary(pic.output)$r.squared
    }else if(model == "Pearson"){
      p <- correlationValFun(pic.y,pic.x,method="pearson",output="p.value",numberDecimals=2,p="num")
      r_coeff <- correlationValFun(pic.y,pic.x,method="pearson",output="r",numberDecimals=2,p="num")
    }else if(model == "Spearman"){
      p <- correlationValFun(pic.y,pic.x,method="spearman",output="p.value",numberDecimals=2,p="num")
      r_coeff <- correlationValFun(pic.y,pic.x,method="spearman",output="r",numberDecimals=2,p="num")
    }else{
      print("Wromg model, lm chosen")
      pic.output <- lm(pic.y~pic.x - 1)
      p <- summary(pic.output)$coefficients[,4]
      r_coeff <- summary(pic.output)$r.squared
    }
    r_p_table[i,] <- c(r_coeff,p)
  }
  if(histo){  # for interpretation of p values distributions: http://varianceexplained.org/statistics/interpreting-pvalue-histogram/
    hist(r_p_table$p)
  }
  r_p_table_summary <- data.frame(
    mean_r = round(mean(r_p_table$r_coeff),2),
    SE_r = round(sd(r_p_table$r_coeff)/sqrt(numberTree),3),
    p_prop_nonSig = round(sum(with(r_p_table,p > 0.05))/numberTree,3),
    p_prop_oneStar = round(sum(with(r_p_table,p <= 0.05 & p > 0.01))/numberTree,3),
    p_prop_twoStar = round(sum(with(r_p_table,p <= 0.01 & p > 0.001))/numberTree,3),
    p_prop_threeStar = round(sum(with(r_p_table,p <= 0.001))/numberTree,3),
    p_final = NA
  )
  # p_majority <- colnames(r_p_table_summary[,3:6])[apply(r_p_table_summary[,3:6],1,which.max)]
  H0 <- r_p_table_summary$p_prop_nonSig > threshold
  Ha_threeStars <- r_p_table_summary$p_prop_threeStar > threshold
  Ha_twoStars <- (r_p_table_summary$p_prop_threeStar + r_p_table_summary$p_prop_twoStar) > threshold && !Ha_threeStars
  Ha_oneStar <- (r_p_table_summary$p_prop_threeStar + r_p_table_summary$p_prop_twoStar + r_p_table_summary$p_prop_oneStar) > threshold && !Ha_twoStars && !Ha_threeStars
  
  if(H0){
    r_p_table_summary$p_final <- ""
  }else if(Ha_oneStar){
    r_p_table_summary$p_final <- "*"
  }else if(Ha_twoStars){
    r_p_table_summary$p_final <- "**"
  }else if(Ha_threeStars){
    r_p_table_summary$p_final <- "***"
  }else{
    r_p_table_summary$p_final <- "(sig.?)"
  }
  r_p_table_summary
}

### function that performs a phylogenetic-informed general least square (PGLS) ANOVA tests over a multiphylo object
# then calculate the mean F
# returns an approximate p value corresponding to the range when the proportion of p values > threshold; if none, return "sig. ?"
PGLS.multiplylo.fun <- function(multiphyloTree,yvar,group.var,species,histo = F,threshold = 0.7){
  numberTree <- length(multiphyloTree)
  F_p_table <- as.data.frame(matrix(nrow = numberTree,ncol = 2))
  colnames(F_p_table) <- c("F.stats","p")
  names(yvar) <- species
  names(group.var) <- species
  for(i in 1:numberTree){
    supertree_1_cut <- multiphyloTree[[i]]
    trait.y <- yvar[supertree_1_cut$tip.label]
    trait.group <- group.var[supertree_1_cut$tip.label]
    
    head(supertree_1_cut$tip.label)
    head(trait.y)
    
    cor.bm <- corBrownian(,phy = supertree_1_cut)
    pgls.output <- gls(trait.y ~ trait.group, correlation = cor.bm)
    # summary(pgls.output)
    anova.gls <- anova(pgls.output)
    F_p_table$F.stats[i] <- anova.gls$`F-value`[2]
    F_p_table$p[i] <- anova.gls$`p-value`[2]
    print(i)
  }
  F_p_table_summary <- data.frame(
    mean_F = round(mean(F_p_table$F.stats),2),
    SE_F = round(sd(F_p_table$F.stats)/sqrt(numberTree),3),
    p_prop_nonSig = round(sum(with(F_p_table,p > 0.05))/numberTree,3),
    p_prop_oneStar = round(sum(with(F_p_table,p <= 0.05 & p > 0.01))/numberTree,3),
    p_prop_twoStar = round(sum(with(F_p_table,p <= 0.01 & p > 0.001))/numberTree,3),
    p_prop_threeStar = round(sum(with(F_p_table,p <= 0.001))/numberTree,3),
    p_final = NA
  )
  # p_majority <- colnames(F_p_table_summary[,3:6])[apply(F_p_table_summary[,3:6],1,which.max)]
  H0 <- F_p_table_summary$p_prop_nonSig > threshold
  Ha_threeStars <- F_p_table_summary$p_prop_threeStar > threshold
  Ha_twoStars <- (F_p_table_summary$p_prop_threeStar + F_p_table_summary$p_prop_twoStar) > threshold && !Ha_threeStars
  Ha_oneStar <- (F_p_table_summary$p_prop_threeStar + F_p_table_summary$p_prop_twoStar + F_p_table_summary$p_prop_oneStar) > threshold && !Ha_twoStars && !Ha_threeStars
  
  if(H0){
    F_p_table_summary$p_final <- ""
  }else if(Ha_oneStar){
    F_p_table_summary$p_final <- "*"
  }else if(Ha_twoStars){
    F_p_table_summary$p_final <- "**"
  }else if(Ha_threeStars){
    F_p_table_summary$p_final <- "***"
  }else{
    F_p_table_summary$p_final <- "(sig.?)"
  }
  F_p_table_summary
}

### function that recieve a vector of "Genus species" characters and returns a dataframe with a column "Genus" and  column "Species", ready to be sent to WoRMS
# in WoRMS: Linefeed (LF) and comma (,)
genus.species.DF.fun <- function(genusSpecies = character()){
  reeftree_tosendWoRMS <- data.frame(do.call("rbind",strsplit(genusSpecies," ",fixed = TRUE)))
  colnames(reeftree_tosendWoRMS) <- c("Genus","Species")
  reeftree_tosendWoRMS$Genus <- as.character(reeftree_tosendWoRMS$Genus)
  reeftree_tosendWoRMS$Species <- as.character(reeftree_tosendWoRMS$Species)
  reeftree_tosendWoRMS
}

# if "star", function that returns "", "*", "**" or "***" depending on the value of the p value given as argument
# if "range", "> 0.01" for instance
pStart.fun <- function(pvalue=numeric(),notation = c("star","range","exact"),numDecimals=3){
  pNotation <- data.frame(star="",range="P > 0.05")
  output <- character()
  if(pvalue < 0.05 & pvalue >= 0.01){
    pNotation$star <- "*"
    pNotation$range <- "P < 0.05" 
  }else if(pvalue < 0.01 & pvalue >= 0.001){
    pNotation$star <- "**"
    pNotation$range <- "P < 0.01" 
  }else if(pvalue < 0.001){
    pNotation$star <- "***"
    pNotation$range <- "P < 0.001" 
  }
  if(pvalue > 1 | pvalue < 0){
    print("wrong p value")
    start <- NA
  }
  if(notation == "range"){
    output <- pNotation$range
  }else if(notation == "star"){
    output <- pNotation$star
  }else if (notation == "exact"){
    output <- paste("P = ",round(pvalue,numDecimals))
  }
  as.character(output)
}

pStart.fun(0.741213,"exact")
# function that performs a permutation test and returns the p value
permutationTest.fun <- function(yvar,group.var,numBPermu=10000){
  obsermeans <- tapply(yvar,group.var,mean)
  obsMeanDiff <- obsermeans[1] - obsermeans[2]
  permResults <- vector()
  for(i in 1:numBPermu){
    permSample <- sample(yvar,replace = F)
    permMeans <- tapply(permSample,group.var,mean)
    permResults[i] <- permMeans[1] - permMeans[2]
  }
  if(obsMeanDiff < 0){
    p <- 2 * (sum(permResults <= obsMeanDiff) / numBPermu) 
  }else{
    p <- 2 * (sum(permResults >= obsMeanDiff) / numBPermu) 
  }
  p
}

# function that remove species to remove for each tree of a multiphylo object:
drop.tip.multiphylo.fun <- function(multiphylo,speciesToRemove=character()){
  supertree_cut <- multiphylo
  numberTrees <- length(supertree_cut)
  for(j in 1: numberTrees){
    supertree_cut[[j]] <- drop.tip(multiphylo[[j]],speciesToRemove)
    print(j)
  }
  supertree_cut
}

### transform values so that the extrem values are excluded from interval, e.g.: [0,1] --> ]0,1[, from (Smithson and Verkuilen 2006)
# minPossible = the smallest possible value in the range, vice versa for maxPossible
y.transf.betareg <- function(var = numeric(), minPossible = numeric(), maxPossible = numeric()){
  # 1st reduce the interval to [0,1]
  y1 <- (var - minPossible)/(maxPossible - minPossible)
  # 2nd exclude 0 and 1 --> ]0,1[
  n <- length(var)
  y2 <- (y1 * (n - 1) + 0.5 ) / n
  y2
}

# ### function that gives the number of days corresponding to the biginning of the month, starting January (e.g., january --> 0, Feb --> 31, etc.) It consider leap years
# months have to be entered as numbers
numberDayFromMonths <- function(month,year){
  nbDays_months <- 0
  feb <- 28
  leapYear = year %in% seq(1900,2100,4) 
  if(leapYear){
    feb = 29
  }
  if(month == 1){
    nbDays_months <- 0                         # useless
  }else if(month == 2){       # note that both 2005 and 2006 have Feb with feb days
    nbDays_months <- 31
  }else if(month == 3){
    nbDays_months <- 31 + feb
  }else if(month == 4){
    nbDays_months <- 31 + feb + 31
  }else if(month == 5){
    nbDays_months <- 31 + feb + 31 + 30
  }else if(month == 6){
    nbDays_months <- 31 + feb + 31 + 30 + 31
  }else if(month == 7){
    nbDays_months <- 31 + feb + 31 + 30 + 31 + 30
  }else if(month == 8){
    nbDays_months <- 31 + feb + 31 + 30 + 31 + 30 + 31
  }else if(month == 9){
    nbDays_months <- 31 + feb + 31 + 30 + 31 + 30 + 31 + 31
  }else if(month == 10){
    nbDays_months <- 31 + feb + 31 + 30 + 31 + 30 + 31 + 31 + 30
  }else if(month == 11){
    nbDays_months <- 31 + feb + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31
  }else if(month == 12){
    nbDays_months <- 31 + feb + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30
  }
  nbDays_months
}

### function that calculates the number of days between 2 years (the upper year being not completed ! -> there are 366 days between 2000 and 2001), considers leap years. If year0 > year --> value is negative
# numberDaysFromYears(year0 = 2008, year= 2000)
numberDaysFromYears <- function(year0,year){
  nbDays_years <- 0
  leapYears <- seq(1900,2100,4)
  numberYears <- if(year >= year0){ year - year0 } else { year0 - year }
  leapYearsPresent <- if(year >= year0){ leapYears[leapYears < year] }else{ leapYears[leapYears >= year] }    # year is not included !
  leapYearsPresent <- if(year >= year0){ leapYearsPresent[leapYearsPresent >= year0] } else { leapYearsPresent[leapYearsPresent < year0]  }
  numberLeapYears <- length(leapYearsPresent)
  numberNoneLeapYears <- numberYears - numberLeapYears
  nbDays_years <- 365 * numberNoneLeapYears + 366 * numberLeapYears
  if(year < year0){ nbDays_years <- -1 * nbDays_years }
  nbDays_years
}

### function that gives the total number of days between two dates, day0 is not being included
# uses numberDayFromMonths()
# uses numberDaysFromYears()
# calculateNumberDaysFun(day=25,month=02,year=2000,day0=03,month0=03,year0=2000)
calculateNumberDaysFun <- function(day,month,year,day0,month0,year0){
  count <- 0
  numberDaysTotal <- numeric()
  for(i in 1:length(day)){
    ndDays_from_days <- day[i] - day0  # to include day0 
    nbDays_from_months <- numberDayFromMonths(month=month[i],year=year[i]) - numberDayFromMonths(month=month0,year=year0) # give the number of days from january 1st (included) of the year, considering possibiily of leap year 
    nbDays_from_years <- numberDaysFromYears(year0 = year0, year = year[i]) # give the number of days between years0 and year, year being not completed !
    count <- count + 1
    numberDaysTotal[count] <- ndDays_from_days + nbDays_from_months + nbDays_from_years
  }
  numberDaysTotal
}

### same as calculateNumberDaysFun but date are giveng as the following characters: "24-04-2010"
# uses calculateNumberDaysFun()
calculateNumberDays.CharacterDates.Fun <- function(date,date0){
  
  dateCut <- data.frame(do.call("rbind",strsplit(date,"-",fixed = TRUE)))
  colnames(dateCut) <- c("day","month","year")
  dateCut$day <- as.numeric(as.character(dateCut$day))
  dateCut$month <- as.numeric(as.character(dateCut$month))
  dateCut$year <- as.numeric(as.character(dateCut$year))
  
  date0Cut <- data.frame(do.call("rbind",strsplit(date0,"-",fixed = TRUE)))
  colnames(date0Cut) <- c("day0","month0","year0")
  date0Cut$day0 <- as.numeric(as.character(date0Cut$day0))
  date0Cut$month0 <- as.numeric(as.character(date0Cut$month0))
  date0Cut$year0 <- as.numeric(as.character(date0Cut$year0))
  
  numberDaysTotal <- calculateNumberDaysFun(dateCut$day,dateCut$month,dateCut$year,date0Cut$day0,date0Cut$month0,date0Cut$year0)
  numberDaysTotal
}

### function that receive a numerical vector of months and return the corresponding 3-letter vector character of months (e.g., "Nov"), or vice versa
monthNumToCharaFun <- function(month){
  month_out <- rep(0,length(month))
  monthYear <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  for(j in 1: length(month)){
    if(class(month[j]) == "numeric"){
      if(month[j] > 12 || month[j] < 1){ print("Yo, this month does not exist")}
      for(i in 1:12){
        if(month[j] == i){
          month_out[j] <- monthYear[i]
          break
        }
      }
    }else if(class(month[j]) == "character"){
      if(! month[j] %in% monthYear){print("Yo, this month does not exist")}
      for(i in 1:12){
        if(month[j] == monthYear[i]){
          month_out[j] <- i
        }
      }
    }else{
      Print("yo, wrong class for month")
    }
  }
  month_out
}

### function that take a vector string month name and returns the corresponding numerical month name
month.char.to.numb.fun <- function(month=character()){
  monthNumb <- c()
  count <- 0
  for(i in 1:length(month)){
    count <- count + 1
    if(month[i] == "Jan" || month[i] == "January"){
      monthNumb[count] <- 01
    }else if(month[i] == "Feb" || month[i] == "February"){
      monthNumb[count] <- 02
    }else if(month[i] == "Mar" || month[i] == "March"){
      monthNumb[count] <- 03
    }else if(month[i] == "Apr" || month[i] == "Aprile"){
      monthNumb[count] <- 04
    }else if(month[i] == "May"){
      monthNumb[count] <- 05
    }else if(month[i] == "Jun" || month[i] == "June"){
      monthNumb[count] <- 06
    }else if(month[i] == "Jul" || month[i] == "July"){
      monthNumb[count] <- 07
    }else if(month[i] == "Aug" || month[i] == "August"){
      monthNumb[count] <- 08
    }else if(month[i] == "Sep" || month[i] == "September"){
      monthNumb[count] <- 09
    }else if(month[i] == "Oct" || month[i] == "October"){
      monthNumb[count] <- 10
    }else if(month[i] == "Nov" || month[i] == "November"){
      monthNumb[count] <- 11
    }else if(month[i] == "Dec" || month[i] == "December"){
      monthNumb[count] <- 12
    }else{
      print(paste("Wrong spelling for: ",month[i],sep=""))
      monthNumb[count] <- NA
    }
  }
  monthNumb
}

### Function that gives the number or the name(s) of species common in 2 list of names
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