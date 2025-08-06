args <- commandArgs(trailingOnly = TRUE)
rept <- as.numeric(args[1])
nSnp_input <- as.numeric(args[2])  # SNPs per chromosome (user input)

library(AlphaSimR)
cat("=æ Running make_base_pop_MV.R for rep:", rept, " | nSnp per chr (input):", nSnp_input, "\n")

source("GlobalParameters.R")
# Override any default nSnp with the user specified value:
assign("nSnp", nSnp_input, envir = .GlobalEnv)
cat(" Overridden nSnp in GlobalEnv with input:", nSnp_input, "per chromosome\n")

library(BGLR)
library(AlphaSimR)
library(AGHmatrix)

# Run base pipeline scripts.
source("CreatParents.R")
cat(" Loaded CreatParents.R\n")

source("FillPipeline.R")
cat(" Loaded FillPipeline.R\n")

source("UpdateParents_private.R")
source("Advance_cycle_private.R")
source("UpdateTesters.R")
source("UpdateParents_DH.R")
source("Advance_cycle.DH.R")
cat(" All pipeline steps completed\n")

# Check for the expected population object.
if (!exists("PublicInbredYT2")) {
  stop("L ERROR: PublicInbredYT2 does not exist. Pipeline may have failed.")
}
cat(" Found PublicInbredYT1. Setting Pop...\n")
Pop <- PublicInbredYT1 
Pop <- setPheno(PublicInbredYT1, H2 = 0.5)

# Check if SP exists. SP is the simulation parameters object used by functions like pullSnpGeno.
if (exists("SP")) {
  cat(" SP object found.\n")
  mySave <- list(Pop = Pop, SP = SP)
} else {
  warning("  SP object not found. Some downstream functions may fail!")
  mySave <- list(Pop = Pop)
}

cat("=Ð Pop dimensions: ", dim(pullSnpGeno(Pop)), "\n")
cat(">ì Pop nInd:", Pop@nInd, " | Total SNPs:", nSnp_input * 10, "\n")

# Save the base pop (and SP) into an RDS file.
save_file <- paste0("SweetCorn/Base_pop_", nSnp_input, "_SNPs_", rept, ".RData")
cat("=¾ Saving Pop (and SP) to: ", save_file, "\n")
saveRDS(mySave, file = save_file)
cat("< Base population successfully saved!\n")
