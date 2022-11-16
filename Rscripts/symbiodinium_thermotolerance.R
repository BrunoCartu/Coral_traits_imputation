# Carturan_et_al_2020_Combining_agent-based_trait-based_and_demographic_approaches_to_model_coral_community_dynamics
# NOT USED
# Goals:
## check dataset from Swain et al. 2016 - Consensus thermotolerance ranking for 110 Symbiodinium phylotypes...
## I did not use the dataset.

require(here)               
# dr_here(show_reason = TRUE)
wd_Rscripts <- getwd()
wd <- gsub("/Rscripts","",wd_Rscripts)
wd_Datasets <- paste(wd,"Datasets",sep="/")
wd_Datasets_original <- paste(wd,"Datasets_original",sep="/")
wd_symbiont_thermotolerance <- paste(wd,"Traits_extra/symbiont_thermotolerance",sep="/")

# http://www.symbiogbr.org/nomenclature

symbThermoTol_original <- read.csv(paste(wd_symbiont_thermotolerance,"Swain et al. - 2016 - Symbionte thermotolerance.csv",sep="/"),header=T,stringsAsFactors = F)
head(symbThermoTol_original)


