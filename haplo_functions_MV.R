##Install Slection tools package:
# if (!requireNamespace("SelectionTools", quietly = TRUE)) {
#   +     install.packages("SelectionTools", repos = "https://population-genetics.uni-giessen.de/~software/SelectionTools_23.1_R_x86_64-pc-linux-gnu.tar.gz")
#   + }
#M: SelectionTools 23.1
library(SelectionTools)
library(SNPRelate)
library(AGHmatrix)
library(AlphaSimR)

##Reconstucting the genetic map
# Total number of markers (assuming no filtering)
total_markers <- sum(Pop@nLoci)

# Create a chromosome vector based on the number of loci per chromosome
chr_vector_full <- rep(1:Pop@nChr, times = Pop@nLoci)

# Create a position vector; assuming markers are equally spaced, use sequential numbers
pos_vector_full <- unlist(lapply(Pop@nLoci, function(n) seq(1, n)))

# Combine into a data.frame representing the genetic map
genetic_map <- data.frame(
  SNP = paste0("SNP", 1:total_markers),
  Chromosome = chr_vector_full,
  Position = pos_vector_full
)

head(genetic_map)

create_GMat_LD <- function(popSim, ld_threshold = 0.2) {
  # Extract SNP matrix and assign names
  SnpMat <- pullSnpGeno(popSim)
  nInd <- nrow(SnpMat)  # Number of individuals (e.g., 1250)
  rownames(SnpMat) <- paste0("Ind", 1:nInd)
  colnames(SnpMat) <- paste0("SNP", 1:ncol(SnpMat))
  
  # Reconstruct genetic map information
  chr_vector <- rep(1:popSim@nChr, times = popSim@nLoci)
  pos_vector <- unlist(lapply(popSim@nLoci, function(n) seq(1, n)))
  selected_indices <- as.numeric(gsub("SNP", "", colnames(SnpMat)))
  chr_vector <- chr_vector[selected_indices]
  pos_vector <- pos_vector[selected_indices]
  
  # Create a temporary GDS file and perform LD pruning
  gds_file <- tempfile(fileext = ".gds")
  snpgdsCreateGeno(gds.fn = gds_file, 
                   genmat = as.matrix(SnpMat),
                   sample.id = rownames(SnpMat), 
                   snp.id = colnames(SnpMat),
                   snp.chromosome = chr_vector,
                   snp.position = pos_vector,
                   snpfirstdim = FALSE)
  
  gds <- snpgdsOpen(gds_file)
  snpset <- snpgdsLDpruning(gds, ld.threshold = ld_threshold)
  pruned_snps <- unlist(snpset)
  snpgdsClose(gds)
  file.remove(gds_file)
  
  # Subset SNP matrix to the pruned markers
  pruned_SnpMat <- SnpMat[, pruned_snps, drop = FALSE]
  
  # Instead of one row per individual, allocate a matrix with 2 rows per individual.
  new_nInd <- 2 * nInd  # Total haplotypes (e.g., 2500 if nInd is 1250)
  haplo_matrix_encoded <- matrix(NA_character_, nrow = new_nInd, ncol = ncol(pruned_SnpMat))
  
  # For each SNP marker, convert genotype to two haplotype values:
  #   If genotype == 0 -> ("H0", "H0")
  #   If genotype == 1 -> ("H0", "H1")
  #   If genotype == 2 -> ("H1", "H1")
  for (j in seq_len(ncol(pruned_SnpMat))) {
    geno_vector <- pruned_SnpMat[, j]
    hap_A <- ifelse(geno_vector == 2, "H1", "H0")
    hap_B <- ifelse(geno_vector == 0, "H0", "H1")
    
    # Assign haplotype A to the odd rows and haplotype B to the even rows.
    haplo_matrix_encoded[seq(1, new_nInd, by = 2), j] <- hap_A
    haplo_matrix_encoded[seq(2, new_nInd, by = 2), j] <- hap_B
  }
  
  # Optionally, set column names for clarity
  snp_names <- colnames(pruned_SnpMat)
  colnames(haplo_matrix_encoded) <- snp_names
  
  return(haplo_matrix_encoded)
}


library(ff)

create_GMat_SlidingWindow <- function(popSim, window_size, overlapSeg, chr_test = NULL) {
  SnpMat <- pullSnpGeno(popSim)
  nInd <- nrow(SnpMat)
  chr_vector_full <- rep(1:popSim@nChr, times = popSim@nLoci)
  selected_indices <- 1:ncol(SnpMat)
  chr_vector <- chr_vector_full[selected_indices]
  
  block_list <- list()
  
  chromosomes <- if (is.null(chr_test)) sort(unique(chr_vector)) else chr_test
  
  for (ch in chromosomes) {
    idx_chr <- which(chr_vector == ch)
    SnpMat_chr <- SnpMat[, idx_chr, drop = FALSE]
    nSNPs <- ncol(SnpMat_chr)
    
    if (nSNPs < window_size) {
      cat("Chromosome", ch, ": Not enough SNPs for window size. Skipping.\n")
      next
    }
    
    starts <- seq(1, nSNPs - window_size + 1, by = window_size - overlapSeg)
    total_windows <- length(starts)
    cat("Chromosome", ch, "has", total_windows, "blocks.\n")
    
    haplo_chr <- matrix(NA_integer_, nrow = nInd, ncol = total_windows)
    
    for (w in seq_along(starts)) {
      window_data <- SnpMat_chr[, starts[w]:(starts[w] + window_size - 1), drop = FALSE]
      
      # Factorization (numerical encoding) to reduce memory
      haplo_encoded <- as.integer(factor(do.call(paste, c(as.data.frame(window_data), sep = ""))))
      haplo_chr[, w] <- haplo_encoded
    }
    
    block_list[[paste0("chr", ch)]] <- haplo_chr
  }
  
  # Combine numerically encoded blocks
  haplo_matrix_encoded <- do.call(cbind, block_list)
  
  cat("\nSliding Window Haplotype G-Matrix Summary:\n")
  cat("Total Individuals:", nInd, "\n")
  cat("Total Blocks:", ncol(haplo_matrix_encoded), "\n")
  
  return(haplo_matrix_encoded)
}

calcHaploRelMat <- function(haplo_matrix) {
  nInd <- nrow(haplo_matrix)
  nMarkers <- ncol(haplo_matrix)
  
  # Center the haplotype matrix by columns (markers)
  Z <- scale(haplo_matrix, center = TRUE, scale = FALSE)
  
  # Compute the genomic relationship matrix (VanRaden method, normalized)
  K <- (Z %*% t(Z)) / sum(apply(Z, 2, var))
  
  cat("Haplotype Relationship Matrix computed:\n")
  cat("Number of Individuals:", nInd, "\n")
  cat("Number of Haplotype Blocks:", nMarkers, "\n")
  
  return(K)
}
