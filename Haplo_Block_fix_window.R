haploblock_fixed_window <- function(PopSim, window_size = 31, min_block_frequency = 0.05) {
  # snp_matrix: a matrix with individuals in rows and SNP markers in columns.
  # window_size: number of SNPs per block.
  # min_block_frequency: minimum fraction of individuals required for a haplotype to be kept.
  SnpMat <- pullSnpGeno(PopSim)
  n_ind <- nrow(SnpMat)
  n_snps <- ncol(SnpMat)
  n_windows <- floor(n_snps / window_size)
  
  block_results <- vector("list", n_windows)
  
  for (w in 1:n_windows) {
    start <- (w - 1) * window_size + 1
    end <- start + window_size - 1
    window_data <- SnpMat[, start:end, drop = FALSE]
    
    # Collapse the SNP values for each individual into one string (the haplotype).
    hap_strings <- apply(window_data, 1, paste0, collapse = "")
    
    # Count the frequency of each unique haplotype.
    hap_table <- table(hap_strings)
    
    # Define the threshold (number of individuals) based on the minimum frequency.
    threshold <- min_block_frequency * n_ind
    frequent_haps <- names(hap_table)[hap_table >= threshold]
    
    # Store the results for this window.
    block_results[[w]] <- list(
      window_range = c(start, end),
      hap_table = hap_table,
      frequent_haplotypes = frequent_haps,
      hap_strings = hap_strings
    )
  }
  
  return(block_results)
}

convert_block_results_to_matrix <- function(block_results) {
  n_windows <- length(block_results)
  n_ind <- length(block_results[[1]]$hap_strings)
  haplo_matrix <- matrix(NA_character_, nrow = n_ind, ncol = n_windows)
  
  for (w in 1:n_windows) {
    haplo_matrix[, w] <- block_results[[w]]$hap_strings
  }
  return(haplo_matrix)
}

compute_relationship_matrix <- function(haplo_matrix) {
  # haplo_matrix: a matrix with individuals in rows and haploblocks in columns.
  # Each element is a haplotype string for that window.
  
  nInd <- nrow(haplo_matrix)
  nBlocks <- ncol(haplo_matrix)
  
  # Initialize an empty matrix for storing the similarity values.
  rel_mat <- matrix(0, nrow = nInd, ncol = nInd)
  
  # Loop over individuals and compare each pair.
  for (i in 1:nInd) {
    for (j in i:nInd) {
      # Count the number of windows where the haplotypes are identical.
      same <- sum(haplo_matrix[i, ] == haplo_matrix[j, ])
      similarity <- same / nBlocks  # Proportion of windows with identical haplotypes.
      
      rel_mat[i, j] <- similarity
      rel_mat[j, i] <- similarity  # Fill symmetric entry.
    }
  }
  
  # Optionally, set row/column names if available.
  rownames(rel_mat) <- rownames(haplo_matrix)
  colnames(rel_mat) <- rownames(haplo_matrix)
  
  return(rel_mat)
}


