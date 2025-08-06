# Simulation Pipeline for Double Haploid Sweet Corn Population
rm(list = ls())
setwd("/blue/mresende/share/viannam/Simulations/SweetCorn/")

# Load required libraries
library('AlphaSimR')
library(BGLR)
library('AGHmatrix')

# Arguments
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)

# Default parameters

rept = as.numeric(args[1])
##rept = 1

### Loading
load(paste0("/blue/mresende/share/viannam/Simulations/SweetCorn/Base_pop_5K_", rept,".RData"))

#G matrix SNP
GMat_SNP <- as.matrix(read.table(paste0("/blue/mresende/share/viannam/Simulations/output/5K_SNPs/GMat_SNP_5K_", rept, ".txt")))
ETA_SNP <- list(A = list(K = GMat_SNP, model="RKHS"))

# Gamtrix haplotype-blocks 
GMat_Haplo <- as.matrix(read.table(paste0("/blue/mresende/share/viannam/Simulations/output/5K_SNPs/GMat_Haplo_5K_", rept, ".txt")))
ETA_Haplo <- list(A = list(K = GMat_Haplo, model="RKHS"))

#Gmatrix haplotype LD based
GMat_LD <- as.matrix(read.table(paste0("/blue/mresende/share/viannam/Simulations/output/5K_SNPs/GMat_LD_5K_", rept, ".txt")))
ETA_LD <- list(A = list(K = GMat_LD, model = "RKHS"))

#G_matrix haplotype_Slide window
GMat_SW <- as.matrix(read.table(paste0("/blue/mresende/share/viannam/Simulations/output/5K_SNPs/GMat_SW_5K_", rept, ".txt")))
ETA_SW <- list(A = list(K = GMat_SW, model = "RKHS"))

#Gmatrix fix window
GMat_FW <- as.matrix(read.table(paste0("/blue/mresende/share/viannam/Simulations/output/5K_SNPs/GMat_FW_5K_", rept, ".txt")))
ETA_FW <- list(A = list(K = GMat_FW, model = "RKHS"))

# Cross-validation parameters
nIter = 10000 
burnIn = 600
nReps = 10
nFolds = 5

# Run GBLUP with cross-validation
run_GBLUP_CV <- function(ETA, model_name) {
  pred_results <- data.frame(Repetition = numeric(), Correlation = numeric(), Model = character(0))
  set.seed(1234)
  for (Rep in 1:nReps) {
    cat("Running Rep:", Rep, "\n")
    yHatCV = rep(NA, length(Pop@pheno))
    tstFold = sample(1:nFolds, size=Pop@nInd, replace=TRUE)
    for (fold in 1:nFolds) {
      yNA = Pop@pheno
      yNA[which(tstFold == fold)] <- NA
      fm = BGLR(y = as.matrix(yNA), ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE)
      yHatCV[which(tstFold == fold)] = fm$yHat[which(tstFold == fold)]
    }
    cor_val = cor(yHatCV, Pop@pheno, use="pairwise.complete.obs")
    pred_results = rbind(pred_results, data.frame(Repetition = Rep, Correlation = cor_val, Model = model_name))
  }
  write.table(pred_results, paste0("/blue/mresende/share/viannam/Simulations/output/5K_SNPs/Pred_", model_name,"_5K_", rept, ".txt"), row.names = FALSE, quote = FALSE)
}

run_GBLUP_CV(ETA_SNP, "SNP")
run_GBLUP_CV(ETA_Haplo, "Haplo_blocks")
run_GBLUP_CV(ETA_LD,"LD")
run_GBLUP_CV(ETA_SW, "SW")
run_GBLUP_CV(ETA_FW, "FW")
