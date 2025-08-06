run_GS_CV <- function(Pop, GMat, model_name, density, rept, repID, output_base,
                      nIter = 10000, burnIn = 600, nFolds = 5) {
  yHatCV <- rep(NA, length(Pop@pheno))
  tstFold <- sample(1:nFolds, size = Pop@nInd, replace = TRUE)
  set.seed(1234 + repID)
  
  for (fold in 1:nFolds) {
    yNA <- Pop@pheno
    yNA[tstFold == fold] <- NA
    # Create a unique output prefix for this rep-model-density combo
    save_prefix <- file.path(output_base, density, paste0("GS_", model_name, "_", density, "_", rept, "_"))
    
    # Ensure the output folder exists
    dir.create(dirname(save_prefix), showWarnings = FALSE, recursive = TRUE)
    
    fm <- BGLR(y = as.matrix(yNA),
               ETA = list(A = list(K = GMat, model = "RKHS")),
               nIter = nIter, burnIn = burnIn, verbose = FALSE,
               saveAt = save_prefix)
    yHatCV[tstFold == fold] <- fm$yHat[tstFold == fold]
  }
  
  cor_val <- cor(yHatCV, Pop@pheno, use = "pairwise.complete.obs")
  
  # Use the provided density string directly (e.g., "50_SNPs")
  out_dir <- file.path(output_base, density)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  out_file <- file.path(out_dir, paste0("Pred_", model_name, "_", density, "_", rept, ".txt"))
  
  write.table(data.frame(Repetition = repID, Correlation = cor_val, Model = model_name),
              out_file, row.names = FALSE, quote = FALSE)
}
