#######################################################################################
####### Topic Haplotypes and haplobocks
####### Marco Peixoto
####### 9/18/2024
#######################################################################################

###-----------------------------------------------
### Function getHaploBlock 
### Author Marco Peixoto
###-----------------------------------------------


#' Function to create the haploblocks from haplotypes
#'
#' @description Using p haplotypes per individual (where p is ploidy), 
#' it creates a haploblock matrix from a population simulated in AlphaSim or from a matrix.
#'
#' @param popSim population from alphasimR
#' @param haploMat matrix with the information on haplotypes
#' @param nBlock integer indicating the size of the block
#' @param overlapSeg integer indicating how many snps should overlap in each read of haploblock
#' @param outPattern boulean, if true, the output will be coded in AB, otherwise, 01 is used. Default is false.
#'
#' @return A matrix with the haploblocks
#' 
#' @export
#'

getHaploBlock = function(popSim, haploMat = NULL, nBlock = 6, overlapSeg = 4, outPattern = FALSE){
  
  if(!is.null(popSim)){
    suppressMessages(requireNamespace('AlphaSimR'))
    
    segSHaplo = pullSnpHaplo(popSim, simParam=SP)
  }else{
    
    segSHaplo = haploMat
  }
  
  
  # Function to concatenate every nStep number of columns
  createBlocks <- function(hapMatrix, nBlockSize, overSeg, sep = "") {
    nCols = ncol(hapMatrix)
    
    start_points <- seq(1, (nCols+1) - nBlockSize, by = nBlockSize - overSeg)
    stepsInit <- lapply(start_points, function(start) {
      seq(start, length.out = nBlockSize)
    })
    
    
    outBlock = lapply(stepsInit, FUN = function(.a){
      apply(hapMatrix[, .a], 1, paste, collapse = sep)
    })
    
    finalHaploblock = do.call(rbind, outBlock)
    return(finalHaploblock)
    
  }
  
  # Applying the function
  if(outPattern){
    segSHaplo[segSHaplo == 0] <- "A"
    segSHaplo[segSHaplo == 1] <- "B"
    
  }
  
  hapBlockMatrix = createBlocks(hapMatrix = segSHaplo, nBlockSize=nBlock-1, overSeg = overlapSeg)
  return(data.frame(hapBlockMatrix))
  
}




###-----------------------------------------------
### Function calc.K and dosage.X
### Author Alejandro Navarro
### Paper: Multiallelic models for QTL mapping in diverse polyploid populations
### https://doi.org/10.1186/s12859-022-04607-z
###-----------------------------------------------

#' Calculation of realized distance matrix (K)
#'
#' @description Using dosage or haplotype scores, a distance matrix is calculated
#' such that the average distance of an individual with itself is 1, and
#' the average with an unrelated individual is 0. Based on the "Realized
#' Relationship" model found in \href{https://dl.sciencesocieties.org/publications/tpg/abstracts/9/2/plantgenome2015.08.0073}{Rosyara et al. 2016}
#'
#' @param matrix Numeric matrix with individuals on rows and markers on columns.
#' @param haplotypes logical, whether haplotypes are present in the matrix. If T,
#' a matrix is expected to have p columns per individual, where p is ploidy.
#' @param ploidy integer, the ploidy number p of all individuals. Only used if
#' haplotypes = T.
#'
#' @return A numeric matrix nxn where n is the number of rows.
#' @export
#'
#' @examples
#' #Create 10 tetraploid individuals with 200 markers each
#' inds <- lapply(1:10,function(i) round(runif(200,0,4)) )
#'
#' #Put them in a matrix or data.frame with individuals in rows
#' geno <- do.call(rbind,inds)
#'
#' #K matrix of 10 individuals (each row is one individual)
#' K<-calc.K(geno)
#'
#' #K matrix of 2 individuals (every 5 rows is one individual)
#' K <- calc.K(geno, haplotye = T, ploidy = 5)
#'
calc.K <- function(  matrix,  haplotypes=F,  ploidy=NULL){
  #When haplotypes are given we can still perform K!
  if(haplotypes){
    if(is.null(ploidy))
      stop("For haplotype distance calculation ploidy must be defined.")
    #Create an ANOVA type matrix for all haplotypes
    matrix<-lapply(1:ncol(matrix),function(i){
      dosage.X(matrix[,i],haplotype = T,ploidy=ploidy,normalize = F)
    })
    matrix<-do.call(cbind,matrix)
  } else {
    #impute NAs. Let's leave this here, useful for our imputator
    matrix <- imputeNA(matrix)
  }
  
  #Substract mean dosage for each marker
  M<-apply(matrix,2,function(x) x-mean(x))
  K<-M%*%t(M) #Calculate distance matrix
  K<-K/mean(diag(K)) #average the center
  colnames(K)<-rownames(matrix)
  rownames(K)<-rownames(matrix)
  return(list(K = K,
              M = M,
              matrix))
}


#' X matrix calculator
#'
#' @description Using a vector of dosages per individual, or p haplotypes per
#' individual (where p is ploidy), it creates a design matrix X, with or without
#' normalization.
#'
#' @param genotypes vector of dosages per individual, or p haplotypes per individual
#' (where p is ploidy).
#' @param ploidy integer indicating ploidy.
#' @param haplotype logical, whether the vector contains haplotype or genotype data
#' @param normalize logical, should the X matrix be normalized?
#'
#' @return if haplotypes = F, a matrix is returned whith a single column
#' containing the same values that were provided. If haplotypes = T, a design matrix
#' is created with ncol = number of haplotypes and nrow = number of individuals.
#' @export
#'
#' @examples
#'
#' dos <- round(runif(100,0,4))
#' dosage.X(dos)
#'
#' hap <- round(sample(400,1:20))
#' dosage.X(haplotype,haplotype = T, ploidy = 4)
dosage.X <- function(genotypes, haplotype=F, ploidy=NULL, normalize = F ){
  
  if(!haplotype){
    alcount <- matrix(genotypes,ncol=1)
    if(normalize) alcount <- (alcount-mean(alcount))/sd(alcount)
  }else{
    #we obtain the different alleles present
    unals <- unique(genotypes)
    #we obtain a design matrix indicating the allele of each chromosome
    #atn: changed this function to add a NA column
    match <- sapply(unals, function(x) {
      if (!is.na(x)) {
        m <- x == genotypes
        m[is.na(m)] <- F
      } else{
        m <- is.na(genotypes)
      }
      return(m)
    })
    
    
    #we count the number of each allele for each individual
    #For haplotypes we need to combine ploidy columns into one count
    n <- length(genotypes)/ploidy
    alcount <- t(sapply(1:n, function(x){
      colSums(match[1:ploidy + (x - 1) * ploidy, ,drop=F])
    }))
    
    #this line will not give the correct answer if we have a single individual
    if(nrow(alcount)==1) alcount <- t(alcount)
    
    if(normalize & ncol(alcount)!= 1) alcount <- apply(alcount,2,function(a) (a-mean(a))/sd(a))
    else if(normalize) alcount <- apply(alcount,2,function(a) (a-mean(a)))
    
    # #This part is also not nice
    # inds <- sapply(1:n,function(i){
    #   j <- 1:ploidy + (i - 1) * ploidy
    #   s <- names(genotypes)[j]
    # })
    
    inds <- unique(substr(names(genotypes), 1, nchar(names(genotypes)) - 2))
    rownames(alcount) <- inds
    colnames(alcount) <- unals
  }
  
  return(alcount)
}





# Imputation ----------------------------------

#' NA replacement by mean score per marker
#'
#' @description For some operations (e.g. calculating correlation) the
#' presence of missing values is not allowed. If the level of missingness
#' per marker is not too high (below 25%), missing values can be imputed
#' using the mean score per marker, so that allele frequency is not affected.
#'
#' @param m a numeric matrix, with markers in rows and
#' individuals in columns.
#' @param miss a number used for missing values instead of NA.
#' The function expects missing values are indicated by NA. Alternatively,
#' since the inpute matrix is a numeric matrix, you could indicate
#' missing values by a number.
#'
#' @return a matrix with no missing values.
#' @keywords internal
imputeNA <- function(m,
                     miss = NULL) {
  
  if(!is.null(miss) & is.numeric(miss)) m[m==miss] <- NA
  
  meanDose <- rowMeans(m, na.rm=T)
  na <- is.na(m)
  means <- na * meanDose
  m[na] <- means[na]
  
  # if (any(is.na(m))) stop("There are still NAs!!")
  return(m)
}
