## Generating G_matrixes
rm(list = ls())
setwd("/blue/mresende/share/viannam/Simulations/SweetCorn/")

# Load required libraries
library('AlphaSimR')
library(BGLR)
library('AGHmatrix')
source("/blue/mresende/share/viannam/Simulations/Haplo/functionsHaplo.R")


# Arguments
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)

# Default parameters

rept = as.numeric(args[1])
##rept = 2

### Loading
load(paste0("/blue/mresende/share/viannam/Simulations/SweetCorn/Base_pop_5K_", rept,".RData"))
source("/blue/mresende/share/viannam/Simulations/SweetCorn/haplo_functions_MV.R")

## Extract SNP-based genomic relationship matrix
SnpMat = pullSnpGeno(Pop)
GMat_SNP = AGHmatrix::Gmatrix(SnpMat, method="VanRaden", ploidy=2, maf=0.05)

## Extract haplotype-blocks based genomic relationship matrix
haploData = getHaploBlock(popSim = Pop, nBlock = 6, overlapSeg = 4, outPattern = TRUE)
GMat_Haplo = (calc.K(matrix=t(haploData), haplotypes=TRUE, ploidy=2))[[1]]

# ##haplotype LD based
haploData_LD <- create_GMat_LD(popSim = Pop, ld_threshold = 0.2)
GMat_LD <- (calc.K(matrix = (haploData_LD), haplotypes = TRUE, ploidy = 2))[[1]]

##haplotype_Slide windonw
haploData_SW <- create_GMat_SlidingWindow(popSim = Pop, window_size = 10, overlapSeg = 4)
GMat_SW <- calcHaploRelMat(haploData_SW)

source("/blue/mresende/share/viannam/Simulations/SweetCorn/Haplo_Block_fix_window.R")

# Step 1: Create block results from the SNP matrix.
blocks <- haploblock_fixed_window(PopSim = Pop, window_size = 31, min_block_frequency = 0.05)

# Step 2: Convert block results into a haplotype matrix.
haplo_matrix <- convert_block_results_to_matrix(blocks)

# Step 3: Compute the relationship matrix.
GMat_FW <- compute_relationship_matrix(haplo_matrix)

##View a portion of the relationship matrix.
print(GMat_FW[1:10, 1:10])

# ---- Print Range and Diagonal Summaries for all G_matrices ----
matrices_list <- list(
  GMat_SNP = GMat_SNP,
  GMat_Haplo = GMat_Haplo,
  GMat_LD = GMat_LD,
  GMat_SW = GMat_SW,
  GMat_FW = GMat_FW
)

for(name in names(matrices_list)) {
  current_mat <- matrices_list[[name]]
  cat("Summary for", name, ":\n")
  cat("  Range: ", min(current_mat), "to", max(current_mat), "\n")
  cat("  Dimensions: ", paste(dim(current_mat), collapse = " x "), "\n")
  cat("  Mean(diag): ", round(mean(diag(current_mat)), 3), "\n\n")
  # Saving each matrix as a text file
  filename <- paste0("/blue/mresende/share/viannam/Simulations/output/", name, "_5K_", rept, ".txt")
  write.table(current_mat, file=filename, sep="\t", quote=FALSE, col.names=NA)

}


save.image(paste0("/blue/mresende/share/viannam/Simulations/output/G_matrixes_5K_", rept, ".RData"))
