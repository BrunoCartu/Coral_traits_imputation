# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# This is the same analysis done in Appendix 4 - Implementation bleaching response.R
# Goals:
## create the intrinsic bleaching probability from coral bleaching resistance traits (see Appendix S4)
## To produce:
### BRI_ResistanceTraits_NineModels_averaging.csv : AIC and AIC weight of the the beta regression models generated
### BRI_ResistanceTraits_10out65_averaging.csv    : same as above but only the 10 1st models = Table S2 in Appendix S4
### bleaching_probability.csv                     : the intrinsic bleaching probability predicted for the 798 species using averaged betagression models 
### functionTrait_DS_Imputed_probaBleaching.csv   : traits_compilation_zooxanthellate_cut.csv + bleaching_probability.csv

## From:
### traits_compilation_zooxanthellate_cut.csv : non-imputed traits dataset with 798 species and 22 traits; created in trait_data_compilation.R
### functionTrait_DS_Imputed.csv              : imputed traits dataset with 798 species and 11 traits ; created in Imputation_traits_missForest.R
 
require(here)
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")

require(betareg)
require(MuMIn)

source(paste(wd_Rscripts,"functions.R",sep="/"))

# Import dataset ------
### Import non-imputed trait dataset:
traitDS_nonImputed <- read.csv(paste(wd_Datasets,"traits_compilation_zooxanthellate_cut.csv",sep="/"),header=T,stringsAsFactors = T)
head(traitDS_nonImputed)

traitDS_Imputed <- read.csv(paste(wd_Datasets,"functionTrait_DS_Imputed.csv",sep="/"),header=T,stringsAsFactors = T)
head(traitDS_Imputed)

# Data preparation -----
### coral species with RBI available:
speciesWithBRI <- traitDS_nonImputed$species[!is.na(traitDS_nonImputed$bleaching_response_index)]
length(speciesWithBRI)  # 304

### subset traitDS_Imputed with only these species:
traitDS_Imputed_cut <- traitDS_Imputed[traitDS_Imputed$species %in% speciesWithBRI,]
length(traitDS_Imputed_cut$species) # 304

traitDS_nonImputed_cut <- traitDS_nonImputed[traitDS_nonImputed$species %in% speciesWithBRI,]
length(traitDS_nonImputed_cut$species) # 304

traitDS_Imputed_cut$bleaching_response_index <- traitDS_nonImputed_cut$bleaching_response_index

### select the response traits:
colnames(traitDS_Imputed_cut)

colnamesToKeep <- c("species","bleaching_response_index","colony_max_diameter","corallite_area","growthFormClasses","growth_rate","reduced_scattering_coefficient")
traitDS_Imputed_cut <- traitDS_Imputed_cut[,colnamesToKeep]
traitDS_Imputed_cut$growthFormClasses <- as.factor(traitDS_Imputed_cut$growthFormClasses)

# log transfomations
traitDS_Imputed_cut$colony_max_diameter_log <- log(traitDS_Imputed_cut$colony_max_diameter)
traitDS_Imputed_cut$corallite_area_log <- log(traitDS_Imputed_cut$corallite_area)
traitDS_Imputed_cut$growth_rate_log <- log(traitDS_Imputed_cut$growth_rate)
traitDS_Imputed_cut$reduced_scattering_coefficient_log <- log(traitDS_Imputed_cut$reduced_scattering_coefficient)
# convert taxon-BRI [0,100] to ]0,1[ for betareg
traitDS_Imputed_cut$bleaching_response_index_ratio_betareg <- y.transf.betareg(var=traitDS_Imputed_cut$bleaching_response_index,minPossible=0,maxPossible=100)
# also # convert taxon-BRI with logit (not used for betareg)
traitDS_Imputed_cut$bleaching_response_index_logit <- logitTrans(traitDS_Imputed_cut$bleaching_response_index) # not used here

## order GF from the most to least complex 
GF_ordered <- c("branching","tables_or_plates","corymbose","digitate","laminar","columnar","massive","encrusting_long_uprights","encrusting")
traitDS_Imputed_cut$growthFormClasses <- ordered(traitDS_Imputed_cut$growthFormClasses,levels=GF_ordered)

# The global model (Beta regression) -----
## assumption testing on the complete model
# http://www.statisticssolutions.com/assumptions-of-multiple-linear-regression/
# http://support.sas.com/kb/57/480.html

formula.0 <- as.formula("bleaching_response_index_ratio_betareg ~ 
                        colony_max_diameter_log +
                        corallite_area_log + 
                        growth_rate_log + 
                        reduced_scattering_coefficient_log + 
                        growthFormClasses +
                        I(colony_max_diameter_log^2) +
                        I(corallite_area_log^2) +
                        I(growth_rate_log^2) +
                        I(reduced_scattering_coefficient_log^2) +
                        colony_max_diameter_log * corallite_area_log +
                        colony_max_diameter_log * growthFormClasses +
                        colony_max_diameter_log * growth_rate_log
                        colony_max_diameter_log * reduced_scattering_coefficient_log +
                        corallite_area_log * growthFormClasses +
                        corallite_area_log * growth_rate_log +
                        growthFormClasses * growth_rate_log")

beta.0 <- betareg(formula.0, data= traitDS_Imputed_cut, link = "logit")
summary(beta.0)

# Cook's distance (to ID influencial observation/species) ------
# species that influence the model too much https://www.r-bloggers.com/outlier-detection-and-treatment-with-r/
# Cook distance (graph top right) computes the influence exerted by each data point (row) on the predicted outcome.
# data points with Cook distance value > 4*mean(Cook dist) are influencial
cook.beta.0 <- cooks.distance(beta.0)

# setwd(wd_Figures)
# pdf(file = "CookDistance_beta.0.pdf",height = 3, width = 4)
par(mar=c(4,4.5,0.8,0.5))
plot(cook.beta.0~as.numeric(names(cook.beta.0)),pch=1,cex=1,main="",ylab="Leverage (Cook's distance)",xlab="",ylim=c(0,max(cook.beta.0)+0.2),las=1,col="darkgrey",lwd=1.5)
mtext("Species reference number",side=1,line=2.5)
# abline(h=4*mean(cook.beta.0),col="red")
abline(h=1,col="red")
#dev.off()
length(names(cook.beta.0[ cook.beta.0 > 4*mean(cook.beta.0, na.rm=T)])) # 18
length(names(cook.beta.0[ cook.beta.0 > 1])) # 0

# Selection of link function for mean ------
link.fun <- c("logit","loglog","cauchit","probit","cloglog","log")
listModels <- list()
for(i in 1:length(link.fun)){
  listModels[[i]] <- betareg(formula.0, data=traitDS_Imputed_cut, link=link.fun[i])
  names(listModels)[i] <- paste("beta.0.",link.fun[i],sep="")
}

# plot predicted vs. observed values for each models
# setwd(wd_Figures)
# pdf(file = "BetaReg.linkFun.mu_Predict.VS.Observed.pdf",height = 3.5, width = 5.5)
layout(matrix(seq(1:6),nrow=2,ncol=3,byrow = T), widths=c(1.3,1,1),heights=c(1,1.3))
obsers <- traitDS_Imputed_cut$bleaching_response_index_ratio_betareg
for(i in 1:length(listModels)){
  model <- listModels[[i]]
  predicted <- predict(model,type="response") 
  Xaxt <- "n"
  Yaxt <- "n"
  side1 <- 0.5
  side2 <- 0.5
  if(i %in% c(1,4)){
    Yaxt <- "s"
    side2 <- 4.2
  }
  if(i > 3){
    Xaxt <- "s"
    side1 <- 4
  }
  par(mar=c(side1,side2,0.5,0.5))
  plot(predicted~obsers,col="darkgrey",lwd=1.5,las=1,xlim=c(0,0.8),ylim=c(0,0.8),cex=1.5,cex.lab=1,xaxt=Xaxt,yaxt=Yaxt,xlab="",ylab="")
  #plot(predicted~obsers,col="darkgrey",lwd=1.5,las=1,cex=1.5,cex.lab=1,xaxt=Xaxt,yaxt=Yaxt,xlab="",ylab="")
  abline(0,1,lty=2)
  if(i %in% c(1,4)){
    mtext("Predicted values",side = 2, line = 2.5)
  }
  if(i > 3){
    mtext("Observed values",side = 1, line = 2.5)
  }
  pseudoR2 <- round(model$pseudo.r.squared,2)
  aic <- round(AIC(model),2)
  legend("topleft",c(paste("link:",link.fun[i]),paste("AIC =",aic),eval(bquote(expression("pseudo R"^2*" = " ~ .(pseudoR2))))), cex=1,bty = "n")
}
# dev.off()

### plot standardized/Pearson residuals vs predicted values
# setwd(wd_Figures)
# pdf(file = "BetaReg.linkFun.mu_PearsonResiduals.vs.predicted.pdf",height = 3.5, width = 5.5)
layout(matrix(seq(1:6),nrow=2,ncol=3,byrow = T), widths=c(1.3,1,1),heights=c(1,1.3))
obsers <- traitDS_Imputed_cut$bleaching_response_index_ratio_betareg
for(i in 1:length(listModels)){
  model <-   listModels[[i]]
  residual.std <- residuals(model,"pearson")
  predicted <- predict(model,type="response")
  # predicted <-  predict(model,type="link")   # linear predictors scale
  ## std.residuals vs Obs number --> not needed because observation are indenpendent from one another (not taken along a temporal or spatial gradiant)
  Xaxt <- "n"
  Yaxt <- "n"
  side1 <- 0.5
  side2 <- 0.5
  if(i %in% c(1,4)){
    Yaxt <- "s"
    side2 <- 4.2
  }
  if(i > 3){
    Xaxt <- "s"
    side1 <- 4
  }
  ## std.residuals vs fitted values --> for overdispersion (Zuur et al., 2009, p. 230)
  par(mar=c(side1,side2,0.5,0.5))
  #plot(residual.std ~ predicted,col="darkgrey",lwd=1.5,las=1,,cex=1,cex.lab=1,xlab="",ylab="",xaxt=Xaxt,yaxt=Yaxt,xlim=c(0,0.6),ylim=c(-2,4))
  plot(residual.std ~ predicted,col="darkgrey",lwd=1.5,las=1,,cex=1,cex.lab=1,xlab="",ylab="",xaxt=Xaxt,yaxt=Yaxt)
  if(i %in% c(1,4)){
    mtext("Pearson residuals",side = 2, line = 2.5)
  }
  if(i > 3){
    mtext("Predicted values",side = 1, line = 2.5)
  }
  abline(0,0,lty=2)
  legend("topright",paste("link:",link.fun[i]), cex=1,bty = "n")
}
# dev.off()

### plot standardized/Pearson residuals vs observed values
# setwd(wd_Figures)
# pdf(file = "BetaReg.linkFun.mu_PearsonResiduals.vs.Observed.pdf",height = 3.5, width = 5.5)
layout(matrix(seq(1:6),nrow=2,ncol=3,byrow = T), widths=c(1.3,1,1),heights=c(1,1.3))
obsers <- traitDS_Imputed_cut$bleaching_response_index_ratio_betareg
for(i in 1:length(listModels)){
  model <-   listModels[[i]]
  residual.std <- residuals(model,"pearson")
  predicted <- predict(model,type="response")
  ## std.residuals vs Obs number --> not needed because observation are indenpendent from one another (not taken along a temporal or spatial gradiant)
  Xaxt <- "n"
  Yaxt <- "n"
  side1 <- 0.5
  side2 <- 0.5
  if(i %in% c(1,4)){
    Yaxt <- "s"
    side2 <- 4.2
  }
  if(i > 3){
    Xaxt <- "s"
    side1 <- 4
  }
  ## std.residuals vs predicted/fitted values --> for overdispersion (Zuur et al., 2009, p. 230)
  par(mar=c(side1,side2,0.5,0.5))
  plot(residual.std ~ predicted,col="darkgrey",lwd=1.5,las=1,,cex=1,cex.lab=1,xlab="",ylab="",xaxt=Xaxt,yaxt=Yaxt,xlim=c(0,0.6),ylim=c(-2,4))
  if(i %in% c(1,4)){
    mtext("Pearson residuals",side = 2, line = 2.5)
  }
  if(i > 3){
    mtext("Observed values",side = 1, line = 2.5)
  }
  abline(0,0,lty=2)
  legend("topright",paste("link:",link.fun[i]), cex=1,bty = "n")
}
# dev.off()

## plot half-normal plot: sdt.residuals vs normal quantiles
# setwd(wd_Figures)
# pdf(file = "BetaReg.linkFun.mu_PearsonResiduals.vs.NormalQuantiles.pdf",height = 6.5, width = 8)
layout(matrix(seq(1:6),nrow=2,ncol=3,byrow = T), widths=c(1,1,1),heights=c(1,1.15))
obsers <- traitDS_Imputed_cut$bleaching_response_index_ratio_betareg
for(i in 1:length(listModels)){
  model <-   listModels[[i]]
  residual.std <- residuals(model,"pearson")
  predicted <- predict(model)
  Xaxt <- "n"
  Yaxt <- "s"
  side1 <- 0.5
  side2 <- 4.4
  if(i > 3){
    Xaxt <- "s"
    side1 <- 4.2
  }
  ## std.residuals vs fitted values --> for overdispersion (Zuur et al., 2009, p. 230)
  par(mar=c(side1,side2,2,0.5))
  # plot(residual.std ~ predicted,col="darkgrey",lwd=1.5,las=1,,cex=1,cex.lab=1,xlab="",ylab="",xaxt=Xaxt,yaxt=Yaxt,xlim=c(0,0.6),ylim=c(-2,4))
  plot(model,which = 5, type = "pearson",cex=1.5,col="darkgrey",lwd=1.5,las=1,cex.lab=1.5,xaxt=Xaxt,yaxt=Yaxt)
  abline(0,0,lty=2)
  legend("topleft",paste("link:",link.fun[i]), cex=1.5,bty = "n")
}
# dev.off()

# Dispersion parameter phi: which factor to implement with which link function -----
link.fun.mu <- "cloglog"
beta.1 <-  betareg(formula.0, data=traitDS_Imputed_cut, link = link.fun.mu)
linkFun.phi <- c("identity","log","sqrt")
## chose variables
independVar <- c("colony_max_diameter_log","corallite_area_log","growthFormClasses","growth_rate_log","reduced_scattering_coefficient_log")
dependentVar <- "bleaching_response_index_ratio_betareg" # bleaching_response_index_ratio_betareg

# ## build up the different models
regMat <- expand.grid(c(TRUE,FALSE), c(TRUE,FALSE),c(TRUE,FALSE),c(TRUE,FALSE),c(TRUE,FALSE))
names(regMat) <- independVar

### Try the different combination of factors and link functions for phi and only retain the one that are significant using likelihood ratio-test (Cribari and Zeileis 2009)
# 1st try with each and only one factor
# retain the models that are significant
# then for each model that are significant, add another factor and retain the significant ones
# do that 5 times in total
# do that for each link function (3)
# note: only identity can be used when growthform is a predictor
listModelSelected <- list()
index <- 1
formula.0.char <- paste(formula.0[2],formula.0[1],formula.0[3])
link.mu <- link.fun.mu
library(lmtest) # for lieklihood test ratio
for(i in 1:length(linkFun.phi)){
  regMat.here <- regMat
  link.fun.phi.here <- linkFun.phi[i]
  print(paste("***** link function:",link.fun.phi.here))
  for(j in 1:length(independVar)){
    factor.here.1 <- independVar[j]
    rest.factors.here.1 <- independVar[independVar !=  factor.here.1]
    # regMat.here <- regMat.here[regMat.here[,factor.here.1],]
    if(!(link.fun.phi.here != "identity" && factor.here.1 == "growthFormClasses")){ # only the identity function can be used with categorical vriable (?) 
      formula.here.1.char <- paste(formula.0.char,"|",factor.here.1,sep=" ")
      formula.here.1 <- as.formula(formula.here.1.char)
      betaREg.here.1 <-  betareg(formula.here.1, data= traitDS_Imputed_cut, link = link.mu,link.phi = link.fun.phi.here)
      # likeluhood ration test
      lhrt <- lrtest(beta.1,betaREg.here.1)
      print(lhrt$`Pr(>Chisq)`[2])
      if(lhrt$`Pr(>Chisq)`[2] < 0.05){
        factors <- paste(factor.here.1)
        print(paste("*** significant factor:",factors))
        listModelSelected[[index <- index + 1]] <- betaREg.here.1
        names(listModelSelected)[index-1] <- paste(link.fun.phi.here,"-",factors)
        # remove the corresponding row in regMat
        rowToRemove <- rep(T,5)
        names(rowToRemove) <- independVar
        rowToRemove[names(rowToRemove) %in% rest.factors.here.1] <- F
        regMat.here <- regMat.here[-rowToRemove,]
        # add another factor
        for(l in 1:length(rest.factors.here.1)){
          factor.here.2 <- rest.factors.here.1[l]
          rest.factors.here.2 <- rest.factors.here.1[rest.factors.here.1 !=  factor.here.2]
          boolean1 <- !(link.fun.phi.here != "identity" && factor.here.2 == "growthFormClasses")  # only the identity function can be used with categorical vriable (?) 
          # check if this combination of factor has not been tested yet --> if not, it should be in regMat.here
          rowTest <- rep(T,5)
          names(rowTest) <- independVar
          rowTest[rest.factors.here.2] <- F
          boolean2 <- sum(duplicated(rbind(regMat.here,rowTest))) > 0
          if(boolean1 && boolean2){ 
            formula.here.2.char <- paste(formula.0.char,"|",factor.here.1,"+",factor.here.2,sep=" ")
            formula.here.2 <- as.formula(formula.here.2.char)
            betaREg.here.2 <-  betareg(formula.here.2, data= traitDS_Imputed_cut, link = link.mu,link.phi = link.fun.phi.here)
            # likelihood ratio test
            lhrt <- lrtest(betaREg.here.1,betaREg.here.2)
            print(lhrt$`Pr(>Chisq)`[2])
            if(lhrt$`Pr(>Chisq)`[2] < 0.05){
              factors <- paste(factor.here.1,"+",factor.here.2)
              print(paste("*** significant factor:",factors))
              listModelSelected[[index <- index + 1]] <- betaREg.here.2
              names(listModelSelected)[index-1] <- paste(link.fun.phi.here,"-",factors)
              # remove the corresponding row in regMat
              rowToRemove <- rep(T,5)
              names(rowToRemove) <- independVar
              rowToRemove[names(rowToRemove) %in% rest.factors.here.2] <- F
              regMat.here <- regMat.here[-rowToRemove,]
              # add another factor
              for(m in 1:length(rest.factors.here.2)){
                factor.here.3 <- rest.factors.here.2[m]
                rest.factors.here.3 <- rest.factors.here.2[rest.factors.here.2 !=  factor.here.3]
                boolean1 <- !(link.fun.phi.here != "identity" && factor.here.3 == "growthFormClasses")
                # check if this combination of factor has not been tested yet --> if not, it should be in regMat.here
                rowTest <- rep(T,5)
                names(rowTest) <- independVar
                rowTest[rest.factors.here.3] <- F
                boolean2 <- sum(duplicated(rbind(regMat.here,rowTest))) > 0
                if(boolean1 && boolean2){
                  formula.here.3.char <- paste(formula.0.char,"|",factor.here.1,"+",factor.here.2,"+",factor.here.3,sep=" ")
                  formula.here.3 <- as.formula(formula.here.3.char)
                  betaREg.here.3 <-  betareg(formula.here.3, data= traitDS_Imputed_cut, link = link.mu,link.phi = link.fun.phi.here)
                  # likeluhood ration test
                  lhrt <- lrtest(betaREg.here.2,betaREg.here.3)
                  print(lhrt$`Pr(>Chisq)`[2])
                  if(lhrt$`Pr(>Chisq)`[2] < 0.05){
                    factors <- paste(factor.here.1,"+",factor.here.2,"+",factor.here.3)
                    print(paste("*** significant factor:",factors))
                    listModelSelected[[index <- index + 1]] <- betaREg.here.3
                    names(listModelSelected)[index-1] <- paste(link.fun.phi.here,"-",factors)
                    # remove the corresponding row in regMat
                    rowToRemove <- rep(T,5)
                    names(rowToRemove) <- independVar
                    rowToRemove[names(rowToRemove) %in% rest.factors.here.3] <- F
                    regMat.here <- regMat.here[-rowToRemove,]
                    # add another factor
                    for(n in 1:length(rest.factors.here.3)){
                      factor.here.4 <- rest.factors.here.3[n]
                      rest.factors.here.4 <- rest.factors.here.3[rest.factors.here.3 !=  factor.here.4]
                      boolean1 <- !(link.fun.phi.here != "identity" && factor.here.4 == "growthFormClasses")
                      # check if this combination of factor has not been tested yet --> if not, it should be in regMat.here
                      rowTest <- rep(T,5)
                      names(rowTest) <- independVar
                      rowTest[rest.factors.here.4] <- F
                      boolean2 <- sum(duplicated(rbind(regMat.here,rowTest))) > 0
                      if(boolean1 && boolean2){
                        formula.here.4.char <- paste(formula.0.char,"|",factor.here.1,"+",factor.here.2,"+",factor.here.3,"+",factor.here.4,sep=" ")
                        formula.here.4 <- as.formula(formula.here.4.char)
                        betaREg.here.4 <-  betareg(formula.here.4, data= traitDS_Imputed_cut, link = link.mu,link.phi = link.fun.phi.here)
                        # likeluhood ration test
                        lhrt <- lrtest(betaREg.here.3,betaREg.here.4)
                        print(lhrt$`Pr(>Chisq)`[2])
                        if(lhrt$`Pr(>Chisq)`[2] < 0.05){
                          factors <- paste(factor.here.1,"+",factor.here.2,"+",factor.here.3,"+",factor.here.4)
                          print(paste("*** significant factor:",factors))
                          listModelSelected[[index <- index + 1]] <- betaREg.here.4
                          names(listModelSelected)[index-1] <- paste(link.fun.phi.here,"-",factors)
                          # remove the corresponding row in regMat
                          rowToRemove <- rep(T,5)
                          names(rowToRemove) <- independVar
                          rowToRemove[names(rowToRemove) %in% rest.factors.here.4] <- F
                          regMat.here <- regMat.here[-rowToRemove,]
                          # add another factor (last potential one)
                          factor.here.5 <- rest.factors.here.4[1]
                          boolean1 <- !(link.fun.phi.here != "identity" && factor.here.5 == "growthFormClasses")
                          # check if this combination of factor has not been tested yet --> if not, it should be in regMat.here
                          rowTest <- rep(T,5)
                          names(rowTest) <- independVar
                          boolean2 <- sum(duplicated(rbind(regMat.here,rowTest))) > 0
                          if(boolean1 && boolean2){
                            formula.here.5.char <- paste(formula.0.char,"|",factor.here.1,"+",factor.here.2,"+",factor.here.3,"+",factor.here.4,"+",factor.here.5,sep=" ")
                            formula.here.5 <- as.formula(formula.here.5.char)
                            betaREg.here.5 <-  betareg(formula.here.5, data= traitDS_Imputed_cut, link = link.mu,link.phi = link.fun.phi.here)
                            # likeluhood ration test
                            lhrt <- lrtest(betaREg.here.4,betaREg.here.5)
                            print(lhrt$`Pr(>Chisq)`[2])
                            if(lhrt$`Pr(>Chisq)`[2] < 0.05){
                              factors <- paste(factor.here.1,"+",factor.here.2,"+",factor.here.3,"+",factor.here.4,"+",factor.here.5)
                              print(paste("*** significant factor:",factors))
                              listModelSelected[[index <- index + 1]] <- betaREg.here.5
                              names(listModelSelected)[index-1] <- paste(link.fun.phi.here,"-",factors)
                              # remove the corresponding row in regMat
                              rowToRemove <- rep(T,5)
                              names(rowToRemove) <- independVar
                              regMat.here <- regMat.here[-rowToRemove,]
                            }
                          }
                        }
                      }
                    }
                  }
                } 
              }
            }
          }
        }
      }
    }
  }
}

### model selection based on the 95% condidence set --------
link.fun.mu <- "cloglog"

beta.1 <- betareg(formula.0, data= traitDS_Imputed_cut, link = link.fun.mu)

### data dredging ()
options(na.action = "na.fail")

dd <- dredge(global.model = beta.1 , rank="AIC",
             subset = dc(colony_max_diameter_log,{I(colony_max_diameter_log^2)}) && 
                      dc(corallite_area_log,{I(corallite_area_log^2)}) &&
                      dc(growth_rate_log,{I(growth_rate_log^2)}) &&
                      dc(reduced_scattering_coefficient_log,{I(reduced_scattering_coefficient_log^2)}))

options(na.action = "na.omit")
head(dd,10)
length(dd$`(Intercept)`) # 504

# chedck if there are cases where there is the Xi^2 but not Xi --> there are many
for(i in 1:length(dd$`(Intercept)`)){
  boolean_MCD <- is.na(dd$colony_max_diameter_log[i]) && !is.na(dd$`I(colony_max_diameter_log^2)`[i])
  boolean_CA <- is.na(dd$corallite_area_log[i]) && !is.na(dd$`I(corallite_area_log^2)`[i])
  boolean_GR <- is.na(dd$growth_rate_log[i]) && !is.na(dd$`I(growth_rate_log^2)`[i])
  boolean_mRSC <- is.na(dd$reduced_scattering_coefficient_log[i]) && !is.na(dd$`I(reduced_scattering_coefficient_log^2)`[i])
  if(boolean_MCD || boolean_CA || boolean_GR || boolean_mRSC){
    print(i)
  }
  if(boolean_MCD){
    print("CMD TOO !!!")
  }
  if(boolean_CA){
    print("CA TOO !!!")
  }
  if(boolean_GR){
    print("GR TOO !!!")
  }
  if(boolean_mRSC){
    print("mRSC TOO !!!")
  }
}
# all good

### model selection
modelSelected_total <- model.sel(object = dd)

### model selection based on Delta AIC < 4
modelSelectedAIC <- get.models(modelSelected_total,subset = delta <= 4) 
summary(modelSelectedAIC)
length(modelSelectedAIC) # 37

### select the models that are in the 95% confidence set
modelSelected95CI <- get.models(modelSelected_total,subset = cumsum(weight) <= .95)
summary(modelSelected95CI)
length(modelSelected95CI) # 65

# Model averaging based on the 95% confidence set -------
modelSelected95CI
modelAveraged <- model.avg(modelSelected95CI,fit = T)

modelAveraged$importance 
modelAveraged$sw # Table S1 in Appendix S4
modelAveraged$msTable[1:10,] # Table S2 in Appendix S4

### write csv for the stats relative to the nine models
# write.csv(modelAveraged$msTable,paste(wd_Datasets,"BRI_ResistanceTraits_NineModels_averaging.csv",sep="/"),row.names = T)
# write.csv(modelAveraged$msTable[1:10,],paste(wd_Datasets,"BRI_ResistanceTraits_10out65_averaging.csv",sep="/"),row.names = T)

summary(modelAveraged)
coef(modelAveraged)
confint(modelAveraged)
formula(modelAveraged)
vcov(modelAveraged)

### Note: the function provide coefficients the "full average"/"full" and the "conditional avergae"/"substet"
# https://stats.stackexchange.com/questions/208724/interpreting-model-averaging-results-in-r
# in Grueber et al., 2011: the "natural average" (= "conditional average" or "subset"): the parameter estimate is averaged only from the parameter values of the model in which it is present, --> it is the "conditional average" or the "substet" above
# at the opposit, "zero method" (= "full average"/"full"): 0 are replaced as value of parameter in the models where the variable is not present 
# I chose the "natural average"/"conditional average"/"subset" parameter values in order to avoid shrinking the factors only present in certain of the nine models
summary(modelAveraged) # Table S3 in Appendix S4

formula.averaged <- as.formula("bleaching_response_index_ratio_betareg ~ 
                                colony_max_diameter_log +
                                corallite_area_log + 
                                growth_rate_log + 
                                reduced_scattering_coefficient_log + 
                                I(colony_max_diameter_log^2)+
                                I(corallite_area_log^2) +
                                I(reduced_scattering_coefficient_log^2) +
                                I(growth_rate_log^2) +
                                colony_max_diameter_log * corallite_area_log +
                                colony_max_diameter_log * growth_rate_log")

beta.average.toChange <-  betareg(formula.averaged, data= traitDS_Imputed_cut, link = link.fun.mu)
beta.average.sub <-beta.average.toChange
beta.average.full <- beta.average.toChange

coefficient <- as.data.frame(modelAveraged$coefficients)
coefficientSub <- coefficient[row.names(coefficient)=="subset",]
coefficientFull <- coefficient[row.names(coefficient)=="full",]

for(i in 1 :length(coefficientSub)){
  variable <-  colnames(coefficientSub)[i]
  val.original <-  beta.average.sub$coefficients$mean[names(beta.average.sub$coefficients$mean) == variable]
  val.updatedSub <- coefficientSub[1,i]
  val.updatedFull <- coefficientFull[1,i]
  print(paste(variable,":",round(val.original,2),"subset:",round(val.updatedSub,2),"full:",round(val.updatedFull,2),sep=" "))
  beta.average.sub$coefficients$mean[names(beta.average.sub$coefficients$mean) == variable] <- val.updatedSub
  beta.average.full$coefficients$mean[names(beta.average.full$coefficients$mean) == variable] <- val.updatedFull
}
beta.average.sub$coefficients$precision <- coefficientSub[1,"(phi)"]
beta.average.full$coefficients$precision <- coefficientFull[1,"(phi)"]

### prediocted values
predictedVal_sub <- predict(beta.average.sub,new = traitDS_Imputed_cut)
predictedVal_full <- predict(beta.average.full,new = traitDS_Imputed_cut)

graphics.off()
plot(predictedVal_full ~ predictedVal_sub) # good

# phi 
phi <- modelAveraged$coefficients["subset","(phi)"]

# residuals 
residuals.sub <- traitDS_Imputed_cut$bleaching_response_index_ratio_betareg - predictedVal_sub
residuals.full <- traitDS_Imputed_cut$bleaching_response_index_ratio_betareg - predictedVal_full
plot(residuals.sub ~ residuals.full)

# std residuals (Ferrai and Cribari-neto 2004)
varY.sub <- predictedVal_sub * (1 - predictedVal_sub) / (1 + phi)     ### from (Ferrari Cribari-Neto 2004)
residuals.Std.sub <- residuals.sub / sqrt(varY.sub)
varY.full <- predictedVal_full * (1 - predictedVal_full) / (1 + phi)     ### from (Ferrari Cribari-Neto 2004)
residuals.Std.full <- residuals.full / sqrt(varY.full)
plot(residuals.Std.sub ~ residuals.Std.full)

### plot std residuasl against fitted values and each of the factors (including growth form) Sub model ----
# setwd(wd_Figures)
# pdf(file = "BRI.validation.average.model.sub.pdf",height = 5, width = 7)
layout(matrix(seq(1:6),nrow=2,ncol=3,byrow = T), widths=c(1.2,1,1),heights=c(1,1))
let <- 0
side2_l <- 4 
side2_s <- 0.5
side1 <- 4
# plot std residuals vs fitted values 
par(mar=c(side1,side2_l,0.5,0.5))
plot(residuals.Std.sub ~ predictedVal_sub,col="darkgrey",lwd=1.5,las=1,cex=1.5,cex.lab=1,xaxt="s",yaxt="s",xlab="",ylab="")
abline(0,0,lty=2)
mtext("Predicted values", side = 1, line = 2.5)
mtext("Pearson residuals", side = 2, line = 2.5)
legend("topright",LETTERS[let <- let + 1], cex=1.2,bty = "n")
# std residuals vs factor1 
par(mar=c(side1,side2_s,0.5,0.5))
xvar <- traitDS_Imputed_cut$colony_max_diameter_log
plot(residuals.Std.sub ~ xvar,col="darkgrey",lwd=1.5,las=1,cex=1.5,cex.lab=1,xaxt="s",yaxt="s",xlab="",ylab="",yaxt = "n")
axis(side=2,labels=F) 
abline(0,0,lty=2)
mtext(expression('Colony maximum diameter (cm)(ln)'),side=1,line=2.5)
# mtext("Pearson residuals", side = 2, line = 2.5)
legend("topright",LETTERS[let <- let + 1], cex=1.2,bty = "n")
# std residuals vs factor2
par(mar=c(side1,side2_s,0.5,0.5))
xvar <- traitDS_Imputed_cut$corallite_area_log
plot(residuals.Std.sub ~ xvar,col="darkgrey",lwd=1.5,las=1,cex=1.5,cex.lab=1,xaxt="s",yaxt="s",xlab="",ylab="",yaxt = "n")
axis(side=2,labels=F) 
abline(0,0,lty=2)
mtext(expression('Corallite area (cm'^2*')(ln)'),side=1,line=2.5)
# mtext("Pearson residuals", side = 2, line = 2.5)
legend("topright",LETTERS[let <- let + 1], cex=1.2,bty = "n")
# std residuals vs factor3
par(mar=c(side1,side2_l,0.5,0.5))
xvar <- traitDS_Imputed_cut$growth_rate_log
plot(residuals.Std.sub ~ xvar,col="darkgrey",lwd=1.5,las=1,cex=1.5,cex.lab=1,xaxt="s",yaxt="s",xlab="",ylab="",yaxt = "s")
axis(side=2,labels=F) 
abline(0,0,lty=2)
mtext(expression('Growth rate (mm.yr'^-1*')(ln)'),side=1,line=2.5)
mtext("Pearson residuals", side = 2, line = 2.5)
legend("topright",LETTERS[let <- let + 1], cex=1.2,bty = "n")
# std residuals vs factor4
par(mar=c(side1,side2_s,0.5,0.5))
xvar <- traitDS_Imputed_cut$reduced_scattering_coefficient_log
plot(residuals.Std.sub ~ xvar,col="darkgrey",lwd=1.5,las=1,cex=1.5,cex.lab=1,xaxt="s",yaxt="s",xlab="",ylab="",yaxt = "n")
axis(side=2,labels=F) 
abline(0,0,lty=2)
mtext(expression('mRSC (mm'^-1*')(ln)'),side=1,line=2.5)
# mtext("Pearson residuals", side = 2, line = 2.5)
legend("topright",LETTERS[let <- let + 1], cex=1.2,bty = "n")
# std residuals vs factor4
par(mar=c(side1,side2_s,0.5,0.5))
xvar <- traitDS_Imputed_cut$growthFormClasses
stripchart.CI.fun(yvar = residuals.Std.sub,group.var = xvar ,confidence.lev = 0.95, groupnames=levels(xvar),las.graph = 1,bars="SE",ylab.char="",xlab.char="",vertical=T,ylabVertiF=8,Yaxt="n",Xaxt="n")
abline(0,0,lty=2)
mtext("Growth forms",side=1,line=2.5,cex=1)
mtext("Simple     to      Complex",side=1,line=1,cex=0.8)
# mtext("Pearson residuals", side = 2, line = 2.5)
legend("topright",LETTERS[let <- let + 1], cex=1.2,bty = "n")

### plot std residuasl against fitted values and each of the factors (including growth form) FULL model ----
# setwd(wd_Figures)
# pdf(file = "BRI.validation.average.model.full.pdf",height = 5, width = 7)
layout(matrix(seq(1:6),nrow=2,ncol=3,byrow = T), widths=c(1.2,1,1),heights=c(1,1))
let <- 0
side2_l <- 4 
side2_s <- 0.5
side1 <- 4
# plot std residuals vs fitted values 
par(mar=c(side1,side2_l,0.5,0.5))
plot(residuals.Std.full ~ predictedVal_full,col="darkgrey",lwd=1.5,las=1,cex=1.5,cex.lab=1,xaxt="s",yaxt="s",xlab="",ylab="")
abline(0,0,lty=2)
mtext("Predicted values", side = 1, line = 2.5)
mtext("Pearson residuals", side = 2, line = 2.5)
legend("topright",LETTERS[let <- let + 1], cex=1.2,bty = "n")
# std residuals vs factor1 
par(mar=c(side1,side2_s,0.5,0.5))
xvar <- traitDS_Imputed_cut$colony_max_diameter_log
plot(residuals.Std.full ~ xvar,col="darkgrey",lwd=1.5,las=1,cex=1.5,cex.lab=1,xaxt="s",yaxt="s",xlab="",ylab="",yaxt = "n")
axis(side=2,labels=F) 
abline(0,0,lty=2)
mtext(expression('Colony maximum diameter (cm)(ln)'),side=1,line=2.5)
# mtext("Pearson residuals", side = 2, line = 2.5)
legend("topright",LETTERS[let <- let + 1], cex=1.2,bty = "n")
# std residuals vs factor2
par(mar=c(side1,side2_s,0.5,0.5))
xvar <- traitDS_Imputed_cut$corallite_area_log
plot(residuals.Std.full ~ xvar,col="darkgrey",lwd=1.5,las=1,cex=1.5,cex.lab=1,xaxt="s",yaxt="s",xlab="",ylab="",yaxt = "n")
axis(side=2,labels=F) 
abline(0,0,lty=2)
mtext(expression('Corallite area (cm'^2*')(ln)'),side=1,line=2.5)
# mtext("Pearson residuals", side = 2, line = 2.5)
legend("topright",LETTERS[let <- let + 1], cex=1.2,bty = "n")
# std residuals vs factor3
par(mar=c(side1,side2_l,0.5,0.5))
xvar <- traitDS_Imputed_cut$growth_rate_log
plot(residuals.Std.full ~ xvar,col="darkgrey",lwd=1.5,las=1,cex=1.5,cex.lab=1,xaxt="s",yaxt="s",xlab="",ylab="",yaxt = "s")
axis(side=2,labels=F) 
abline(0,0,lty=2)
mtext(expression('Growth rate (mm.yr'^-1*')(ln)'),side=1,line=2.5)
mtext("Pearson residuals", side = 2, line = 2.5)
legend("topright",LETTERS[let <- let + 1], cex=1.2,bty = "n")
# std residuals vs factor4
par(mar=c(side1,side2_s,0.5,0.5))
xvar <- traitDS_Imputed_cut$reduced_scattering_coefficient_log
plot(residuals.Std.full ~ xvar,col="darkgrey",lwd=1.5,las=1,cex=1.5,cex.lab=1,xaxt="s",yaxt="s",xlab="",ylab="",yaxt = "n")
axis(side=2,labels=F) 
abline(0,0,lty=2)
mtext(expression('mRSC (mm'^-1*')(ln)'),side=1,line=2.5)
# mtext("Pearson residuals", side = 2, line = 2.5)
legend("topright",LETTERS[let <- let + 1], cex=1.2,bty = "n")
# std residuals vs factor4
par(mar=c(side1,side2_s,0.5,0.5))
xvar <- traitDS_Imputed_cut$growthFormClasses
stripchart.CI.fun(yvar = residuals.Std.full,group.var = xvar ,confidence.lev = 0.95, groupnames=levels(xvar),las.graph = 1,bars="SE",ylab.char="",xlab.char="",vertical=T,ylabVertiF=8,Yaxt="n",Xaxt="n")
abline(0,0,lty=2)
mtext("Growth forms",side=1,line=2.5,cex=1)
mtext("Simple     to      Complex",side=1,line=1,cex=0.8)
# mtext("Pearson residuals", side = 2, line = 2.5)
legend("topright",LETTERS[let <- let + 1], cex=1.2,bty = "n")

#### plot QQplot submodel ------
graphics.off()
# setwd(wd_Figures)
# pdf(file = "BRI.validation.average.model.QQplot.sub.full.pdf",height = 3, width = 5)
layout(matrix(seq(1:2),nrow=1,ncol=2,byrow = T), widths=c(1,1),heights=c(1,1))
par(mar=c(4,4,0.5,0.5))
qqnorm(residuals.Std.sub,col="darkgrey",lwd=1.5,las=1,cex=1.5,main="")
qqline(residuals.Std.sub,lty=2)
legend("topleft","A", cex=1.2,bty = "n")
# dev.off()
shapiro.test(residuals.Std.sub) # W = 0.94235, p-value = 1.62e-09
#### plot QQplot Full model ------
# graphics.off()
# setwd(wd_Figures)
# pdf(file = "BRI.validation.average.model.QQplot.full.pdf",height = 4, width = 4)
par(mar=c(4,4,0.5,0.5))
qqnorm(residuals.Std.full,col="darkgrey",lwd=1.5,las=1,cex=1.5,main="")
qqline(residuals.Std.full,lty=2)
legend("topleft","B", cex=1.2,bty = "n")
shapiro.test(residuals.Std.full) # W = 0.94235, p-value = 1.62e-09

# Predict the a species-specific coral bleaching response for all the species with  -------
dat.pred <- data.frame(species = traitDS_Imputed$species,
                       colony_max_diameter_log = log(traitDS_Imputed$colony_max_diameter),
                       growth_rate_log = log(traitDS_Imputed$growth_rate),
                       corallite_area_log = log(traitDS_Imputed$corallite_area),
                       reduced_scattering_coefficient_log = log(traitDS_Imputed$reduced_scattering_coefficient),
                       growthFormClasses = traitDS_Imputed$growthFormClasses,
                       stringsAsFactors = T)

predictedVal.798.sub <- predict(beta.average.sub,new=dat.pred)
names(predictedVal.798.sub) <- dat.pred$species
predictedVal.304.sub <- predictedVal.798.sub[names(predictedVal.798.sub) %in% traitDS_Imputed_cut$species]

predictedVal.798.full <- predict(beta.average.full,new=dat.pred)
names(predictedVal.798.full) <- dat.pred$species
predictedVal.304.full <- predictedVal.798.full[names(predictedVal.798.full) %in% traitDS_Imputed_cut$species]

plot(predictedVal.798.full ~ predictedVal.798.sub)

# Calculate pseudo R2  ----- 
# "the pseudo R2 is defined as the square of the sample correlation coefficient between "linears predictors" and g(y)" (Ferrari and Cribari-Neto 2004).
# "squared correlation of linear predictor and link-transformed response" Zeileis et al., 2016, Package "betareg"
cloglog.fun <- function(x){log(-log(1-x))}
cloglog.inverse.fun <- function(x){1 - exp(-exp(x))}

# x <- seq(0.01,0.99,0.01)
# plot(cloglog.inverse.fun(cloglog.fun(x)) ~ x)
# plot(cloglog.fun(x)~x)
# plot(logitTrans(x)~x)

linear.predictors.sub <- cloglog.fun(predictedVal.304.sub)
linear.predictors.full <- cloglog.fun(predictedVal.304.full)

linear.predictors.sub <- logitTrans(predictedVal.304.sub)
linear.predictors.full <- logitTrans(predictedVal.304.full)

pseudo.R2.dataFrame <- data.frame(
  linear.predictors.sub = linear.predictors.sub,
  linear.predictors.full = linear.predictors.full,
  reponseVar.cloglog = cloglog.fun(traitDS_Imputed_cut$bleaching_response_index_ratio_betareg)
)

# pseudo.R2.dataFrame <- data.frame(
#   linear.predictors.sub = linear.predictors.sub,
#   linear.predictors.full = linear.predictors.full,
#   reponseVar.cloglog = logitTrans(traitDS_Imputed_cut$bleaching_response_index_ratio_betareg)
# )

# plot(pseudo.R2.dataFrame$linear.predictors.sub ~  pseudo.R2.dataFrame$reponseVar.cloglog)
# plot(pseudo.R2.dataFrame$linear.predictors.full ~ pseudo.R2.dataFrame$reponseVar.cloglog)

pseudo.R2.sub <- (cor.test(pseudo.R2.dataFrame$linear.predictors.sub , pseudo.R2.dataFrame$reponseVar.cloglog , method = "pearson")$estimate)^2
pseudo.R2.sub #  0.1485445 
pseudo.R2.full<- (cor.test(pseudo.R2.dataFrame$linear.predictors.full , pseudo.R2.dataFrame$reponseVar.cloglog , method = "pearson")$estimate)^2
pseudo.R2.full #  0.1847146

# Plot validation figure for subset model -------
graphics.off()
# setwd(wd_Figures)
# pdf(file = "BRI.vs.IBP.sub.pdf",height = 3.5, width = 7)
zones <- matrix(c(1,5,2, 
                  3,4,2), ncol = 3, byrow = TRUE)
layout(zones, widths=c(4,1,5), heights = c(1,4))
# layout.show(n=5)
nClass <- 20
let <- 0
xvar <- traitDS_Imputed_cut$bleaching_response_index_ratio_betareg
yvar <- predictedVal.304.sub 
hist.BRI <- hist(xvar,breaks = seq(0,1,0.05),plot=F)
hist.IBP <- hist(yvar,breaks = seq(0,1,0.05),plot=F)
top.count <- max(hist.BRI$counts,hist.IBP$counts)
# 1 - plot hist.BRI
par(mar=c(0,4,0.5,0.5))
barplot(hist.BRI$counts, axes = FALSE,space = 0, horiz = F,ylim=c(0,top.count))
# 2 -plot hist.IBP.tot 
par(mar=c(4,6,2,2))
hist(predictedVal.798.sub ,breaks = seq(0,1,0.05),main="",xlim=c(0,1),las=1,nclass = nClass,cex.lab=1.2,xlab="",xaxt="s",yaxt="s",plot=T,col = "grey",ylab="")
mtext(expression("Intrinsic probability of bleaching"),side=1,line=2.8,cex = 1)
mtext(expression("Frequency"),side=2,line=2.5,cex = 1)
legend("top",paste("n = ",length(traitDS_Imputed$species)), cex=1.2,bty = "n")
legend("topright","B", cex=1.2,bty = "n")
# 3 - plot scatter plot
par(mar=c(4,4,0.5,0.5))
plot(yvar ~ xvar,col="darkgrey",lwd=1.5,las=1,cex=1.5,cex.lab=1,xaxt="s",yaxt="s",xlab="",ylab="",ylim=c(0,1),xlim=c(0,1))
abline(0,1,lty=2)
mtext(expression("Intrinsic probability of bleaching"),side=2,line=2.5,cex = 1)
mtext(expression("taxon-BRI"["]0,1["]),side=1,line=2.8,cex = 1)
legend("topleft",c(eval(bquote(expression("pseudo R"^2*" = " ~ .(round(pseudo.R2.sub,2))))),paste("n = ",length(traitDS_Imputed_cut$species))), cex=1.2,bty = "n")
# 4 - plot hist.IBP
par(mar=c(4,0,0.5,0))
barplot(hist.IBP$counts, axes = F,space = 0, horiz = T,xlim=c(0,top.count))
# 5 - plot legend A
par(mar=c(0,0,2,0))
plot(0,1,col="white",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
legend("top",c("A"), cex=1.2,bty = "n")
dev.off()

# Plot validation figure for fullset model -------
graphics.off()
# setwd(wd_Figures)
# pdf(file = "BRI.vs.IBP.full.pdf",height = 3.5, width = 7)
zones <- matrix(c(1,5,2, 
                  3,4,2), ncol = 3, byrow = TRUE)
layout(zones, widths=c(4,1,5), heights = c(1,4))
# layout.show(n=5)
nClass <- 20
let <- 0
xvar <- traitDS_Imputed_cut$bleaching_response_index_ratio_betareg
yvar <- predictedVal.304.full 
hist.BRI <- hist(xvar,breaks = seq(0,1,0.05),plot=F)
hist.IBP <- hist(yvar,breaks = seq(0,1,0.05),plot=F)
top.count <- max(hist.BRI$counts,hist.IBP$counts)
# 1 - plot hist.BRI
par(mar=c(0,4,0.5,0.5))
barplot(hist.BRI$counts, axes = FALSE,space = 0, horiz = F,ylim=c(0,top.count))
# 2 -plot hist.IBP.tot 
par(mar=c(4,6,2,2))
hist(predictedVal.798.full ,breaks = seq(0,1,0.05),main="",xlim=c(0,1),las=1,nclass = nClass,cex.lab=1.2,xlab="",xaxt="s",yaxt="s",plot=T,col = "grey",ylab="",ylim=c(0,300))
mtext(expression("Intrinsic probability of bleaching"),side=1,line=2.8,cex = 1)
mtext(expression("Frequency"),side=2,line=2.5,cex = 1)
legend("top",paste("n = ",length(traitDS_Imputed$species)), cex=1.2,bty = "n")
legend("topright","B", cex=1.2,bty = "n")
# 3 - plot scatter plot
par(mar=c(4,4,0.5,0.5))
plot(yvar ~ xvar,col="darkgrey",lwd=1.5,las=1,cex=1.5,cex.lab=1,xaxt="s",yaxt="s",xlab="",ylab="",ylim=c(0,1),xlim=c(0,1))
abline(0,1,lty=2)
mtext(expression("Intrinsic probability of bleaching"),side=2,line=2.5,cex = 1)
mtext(expression("taxon-BRI"["]0,1["]),side=1,line=2.8,cex = 1)
legend("topleft",c(eval(bquote(expression("pseudo R"^2*" = " ~ .(round(pseudo.R2.full,2))))),paste("n = ",length(traitDS_Imputed_cut$species))), cex=1.2,bty = "n")
# 4 - plot hist.IBP
par(mar=c(4,0,0.5,0))
barplot(hist.IBP$counts, axes = F,space = 0, horiz = T,xlim=c(0,top.count))
# 5 - plot legend A
par(mar=c(0,0,2,0))
plot(0,1,col="white",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
legend("top",c("A"), cex=1.2,bty = "n")

# Write csv ----

traitDS_Imputed$bleaching_probability <- predictedVal.798.sub

bleachingProbability <- traitDS_Imputed[,c("species","bleaching_probability")]

# write.csv(traitDS_Imputed,paste(wd_Datasets,"functionTrait_DS_Imputed_probaBleaching.csv",sep="/"),row.names = F)
# write.csv(bleachingProbability,paste(wd_Datasets,"bleaching_probability.csv",sep="/"),row.names = F)

