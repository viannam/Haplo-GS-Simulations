#==================================================================================
#Running simulation SNP and Haplo GS over 30 years
# Pipeline adapted from from Dr. Marco A. Peixoto

###>>>----------------------------------------------------###
###> 1. Setting the base population and initial variables ###
###>>>----------------------------------------------------###
rm(list=ls())

setwd("/blue/mresende/share/viannam/Simu_v2")

require('AlphaSimR')
require('BGLR')
require('AGHmatrix')

options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
rep <- as.numeric(args[1])

#---- 2. Setting the global parameters and creating files for results
# Load global parameters
source("GlobalParameters.R")

{
  # Creating the files for record the results
  MeanG_pop = matrix(NA, 35)
  MeanP_pop = matrix(NA, 35)
  Accuracy  = matrix(NA, 35)
  MeanA_pop = matrix(NA, 35)
  VarA_pop  = matrix(NA, 35)
  VarG_pop  = matrix(NA, 35)
  GenicVA_pop  = matrix(NA, 35)
  GenicVG_pop = matrix(NA, 35)
  LDB_pop  = matrix(NA, 35)
  LDW_pop  = matrix(NA, 35)
  LDT_pop  = matrix(NA, 35)
  covG_L_pop = matrix(NA, 35)
  LD_pop = matrix(NA, 35)
  covG_HW_pop = matrix(NA, 35)
  inbreeding = matrix(NA, 35)
  inbreedingQTL = matrix(NA, 35)
  relIbdMean = matrix(NA, 35)
  meanHyb = matrix(NA, 35)
}

#---- 3. Creating parents and filling the initial pipeline

# Create initial parents and set testers and hybrid parents
source("CreateParents.R")

# Fill breeding pipeline with unique individuals from initial parents
source("FillPipeline.R")

# p-values for GxY effects
P = runif(burninYears+futureYears)


###>>>-------------------------------------------------
###> 2. Burn-in period - do not change it
###>>>-------------------------------------------------

# 2.1 Cycle years
for(year in 1:burninYears){ #Change to any number of desired years
  cat("Working on year:",year,"\n")
  p = P[year]
  source("UpdateParents.R")  # Pick new parents based on last year's data
  source("UpdateTesters.R")  # Pick new testers and hybrid parents
  source("AdvancePheno.R")   # Advances yield trials by a year
  source("WriteRecordsGS.R") # Write records for GS predictions
  source("UpdateResults.R")  # Track summary data

}

# 2.2 Save burn-in to load later use
save.image(paste0("BURNIN_",rep,".RData"))

###>>>-------------------------------------------------
###> 3. Scenario 1 - Genomic selection program based on SNPs
###>>>-------------------------------------------------

# 3.0 Loading the scenarios
load(paste0("BURNIN_",rep,".RData"))


# 3.1 Looping through the years
source("GenomicModel.R")

cat("Working on Scenario 1\n")
for(year in (burninYears+1):(burninYears+20)){

  cat("Working on year:",year,"\n")

  p = P[year]
  source("UpdateParents_GS.R") #Pick new parents based on last year's data
  source("UpdateTesters.R") #Pick new testers and hybrid parents
  source("AdvanceDHGS.R") #Advances yield trials by a year
  source('WriteRecordsGS.R') # Fill data for next cycle GS models
  source("UpdateResultsGS.R") #Track summary data
}

# 3.2 Recording results
output1 = data.frame(rep=rep(rep, 35),
                     scenario=rep("SNP_300", 35),
                     MeanG_pop,
                     MeanP_pop,
                     Accuracy,
                     MeanA_pop,
                     VarA_pop,
                     VarG_pop,
                     GenicVA_pop,
                     GenicVG_pop,
                     LDB_pop,
                     LDW_pop,
                     LDT_pop,
                     covG_L_pop,
                     LD_pop,
                     covG_HW_pop,
                     inbreeding,
                     inbreedingQTL,
                     relIbdMean,
                     meanHyb,
                     stringsAsFactors=FALSE)

# 3.3 Saving the results as RDS
saveRDS(output1,paste0("Results_300_SNP_",rep,".rds"))


###>>>-------------------------------------------------
###> 4. Scenario 2 - Haplotypes
###>>>-------------------------------------------------
# 4.0 Loading the scenarios
load(paste0("BURNIN_",rep,".RData"))

# 4.1 Looping through the years
source("GenomicModel.R")

cat("Working on Scenario 2\n")
for(year in (burninYears+1):(burninYears+20)){

  cat("Working on year:",year,"\n")

  p = P[year]
  source("UpdateParents_Haplo.R") #Pick new parents based on last year's data
  source("UpdateTesters.R") #Pick new testers and hybrid parents
  source("AdvanceHaplo.R") #Advances yield trials by a year
  source('WriteRecordsGS.R') # Fill data for next cycle GS models
  source("UpdateResultsGS.R") #Track summary data
}


# 4.2 Recording results
output2 = data.frame(rep=rep(rep, 35),
                     scenario=rep("Haplo_300", 35),
                     MeanG_pop,
                     MeanP_pop,
                     Accuracy,
                     MeanA_pop,
                     VarA_pop,
                     VarG_pop,
                     GenicVA_pop,
                     GenicVG_pop,
                     LDB_pop,
                     LDW_pop,
                     LDT_pop,
                     covG_L_pop,
                     LD_pop,
                     covG_HW_pop,
                     inbreeding,
                     inbreedingQTL,
                     relIbdMean,
                     meanHyb,
                     stringsAsFactors=FALSE)

# 4.3 Saving the results as RDS
saveRDS(output2, paste0("Results_300_Haplo_",rep,".rds"))

###>>>-------------------------------------------------
###> 3. Scenario 3 - Genomic selection program based on QTLs- SNPs
###>>>-------------------------------------------------

# 3.0 Loading the scenarios
load(paste0("BURNIN_",rep,".RData"))


# 3.1 Looping through the years
source("QTL_GenomicModel.R")

cat("Working on Scenario 1\n")
for(year in (burninYears+1):(burninYears+20)){

  cat("Working on year:",year,"\n")

  p = P[year]
  source("UpdateParents_GS.R") #Pick new parents based on last year's data
  source("UpdateTesters.R") #Pick new testers and hybrid parents
  source("AdvanceDHGS.R") #Advances yield trials by a year
  source('WriteRecordsGS.R') # Fill data for next cycle GS models
  source("UpdateResultsGS.R") #Track summary data
}

# 3.2 Recording results
output1 = data.frame(rep=rep(rep, 35),
                     scenario=rep("SNP_QTL_300", 35),
                     MeanG_pop,
                     MeanP_pop,
                     Accuracy,
                     MeanA_pop,
                     VarA_pop,
                     VarG_pop,
                     GenicVA_pop,
                     GenicVG_pop,
                     LDB_pop,
                     LDW_pop,
                     LDT_pop,
                     covG_L_pop,
                     LD_pop,
                     covG_HW_pop,
                     inbreeding,
                     inbreedingQTL,
                     relIbdMean,
                     meanHyb,
                     stringsAsFactors=FALSE)

# 3.3 Saving the results as RDS
saveRDS(output1,paste0("Results_300_SNP_QTLs_",rep,".rds"))


###>>>-------------------------------------------------
###> 4. Scenario 4 - GS Haplotypes QTLs
###>>>-------------------------------------------------
# 4.0 Loading the scenarios
load(paste0("BURNIN_",rep,".RData"))

# 4.1 Looping through the years
source("QTL_GenomicModel.R")

cat("Working on Scenario 2\n")
for(year in (burninYears+1):(burninYears+20)){

  cat("Working on year:",year,"\n")

  p = P[year]
  source("UpdateParents_Haplo.R") #Pick new parents based on last year's data
  source("UpdateTesters.R") #Pick new testers and hybrid parents
  source("AdvanceHaplo.R") #Advances yield trials by a year
  source('WriteRecordsGS.R') # Fill data for next cycle GS models
  source("UpdateResultsGS.R") #Track summary data
}


# 4.2 Recording results
output2 = data.frame(rep=rep(rep, 35),
                     scenario=rep("Haplo_300", 35),
                     MeanG_pop,
                     MeanP_pop,
                     Accuracy,
                     MeanA_pop,
                     VarA_pop,
                     VarG_pop,
                     GenicVA_pop,
                     GenicVG_pop,
                     LDB_pop,
                     LDW_pop,
                     LDT_pop,
                     covG_L_pop,
                     LD_pop,
                     covG_HW_pop,
                     inbreeding,
                     inbreedingQTL,
                     relIbdMean,
                     meanHyb,
                     stringsAsFactors=FALSE)

# 4.3 Saving the results as RDS
saveRDS(output2, paste0("Results_300_Haplo_QTL_",rep,".rds"))


###>>>-------------------------------------------------
###> 6. Removing the temporary files
###>>>-------------------------------------------------

# # Delete tmp file
# file.remove(paste0("BURNIN_",rep,"_S.RData"))
