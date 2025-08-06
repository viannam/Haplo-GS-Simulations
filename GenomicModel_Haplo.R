

#' `runGBLUP`
#' Function to run a BayesB multivarite and univariate using BGLR in populations from AlphaSimR
#'
#'
#'
#' @param trainPop population where to train the model
#' @param targetPop Target population from AlphaSimR
#' @param nTraits number of traits in the simulation
#' @param Model the model implemented should be additive or additive and dominance model
#' 
#' @return estimated SNP effects for the individuals in the population
#' 
#' @import BGLR
#' 
#' @author Marco A Peixoto 


runGBLUP <- function(trainPop, targetPop, nTraits = 1, Model = NULL){
  
  library('AlphaSimR')
  library(BGLR)
  library('AGHmatrix')
  source("/blue/mresende/share/viannam/Simulations/Haplo/functionsHaplo.R")
  
  # 1.1 Phenotypic data
  y1 = as.matrix(pheno(trainPop))
  rownames(y1) = c(trainPop@id)
  
  y2 = matrix(ncol = nTraits, nrow = length(targetPop@id))
  rownames(y2) = targetPop@id
  
  y = rbind(y1,y2)
  
  # 1.2 SNP panel from base population
  M1 = pullSnpGeno(trainPop) #Access SNP genotype data
  rownames(M1) = c(trainPop@id)
  
  # 1.3 SNP panel from target population
  M2 = pullSnpGeno(targetPop) #Access SNP genotype data
  rownames(M2) = targetPop@id
  
  Markers = rbind(M1,M2)
  
  if(Model == "Haplo"){
    # Additive matrix
    relMat = AGHmatrix::Gmatrix(SNPmatrix = Markers)
    
    ETA = list(list(model='RKHS', K=relMat))
    
  }else if(Model == "Haplo"){
    # Additive matrix
    haploData = getHaploBlock(popSim = Pop, nBlock = 6, overlapSeg = 4, outPattern = TRUE)
    GMat_Haplo = (calc.K(matrix=t(haploData), haplotypes=TRUE, ploidy=2))[[1]]
    
    # Kernel
    ETA = list(A = list(model='RKHS', K=GMat_Haplo))
    
  }
  
  if(nTraits == 1){
    # Model
    gmodel = BGLR::BGLR(y = y,
                        ETA=ETA,
                        nIter = 10000,
                        burnIn = 1000,
                        thin = 5)
    
    
    ###>>>----------- 5. Pushing data back to AlphaSimR code
    indEbv = as.matrix(gmodel$yHat)[((trainPop@nInd+1):(trainPop@nInd + targetPop@nInd))]
    
  }else if(nTraits > 1){
    # Model
    gmodel = BGLR::Multitrait(y = y,
                              ETA=ETA,
                              nIter = 10000,
                              burnIn = 1000,
                              thin = 5)
    
    
    ###>>>----------- 5. Pushing data back to AlphaSimR code
    indEbv = as.matrix(gmodel$ETAHat)[((trainPop@nInd+1):(trainPop@nInd + targetPop@nInd)),]
    
  }
  
  
  
  return(indEbv)
  
}

