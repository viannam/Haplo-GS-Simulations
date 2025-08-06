generate_all_gmatrices <- function(Pop, density, rept, output_base) {
  dir.create(file.path(output_base, density), showWarnings = FALSE, recursive = TRUE)
  output_path <- file.path(output_base, density)
  
  SnpMat <- pullSnpGeno(Pop)
  
  # SNP-based G matrix
  GMat_SNP <- AGHmatrix::Gmatrix(SnpMat, method = "VanRaden", ploidy = 2, maf = 0.05)
  
  # Haploblock G matrix
  haploData <- getHaploBlock(popSim = Pop, nBlock = 6, overlapSeg = 4, outPattern = TRUE)
  GMat_Haplo <- (calc.K(matrix = t(haploData), haplotypes = TRUE, ploidy = 2))[[1]]
  
  # LD-based G matrix
  GMat_LD <- (calc.K(matrix = create_GMat_LD(Pop), haplotypes = TRUE, ploidy = 2))[[1]]
  
  # Sliding window G matrix
  GMat_SW <- calcHaploRelMat(create_GMat_SlidingWindow(Pop, window_size = 10, overlapSeg = 4))
  
  # Fixed window G matrix
  blocks <- haploblock_fixed_window(Pop, window_size = 31, min_block_frequency = 0.05)
  GMat_FW <- compute_relationship_matrix(convert_block_results_to_matrix(blocks))
  
  mats <- list(SNP = GMat_SNP, Haplo = GMat_Haplo, LD = GMat_LD, SW = GMat_SW, FW = GMat_FW)
  for (name in names(mats)) {
    write.table(mats[[name]],
                file = file.path(output_path, paste0("GMat_", name, "_", density, "_", rept, ".txt")),
                sep = "\t", quote = FALSE, col.names = NA)
  }
  
  return(mats)
}
