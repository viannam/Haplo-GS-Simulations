#Set initial yield trials with unique individuals

P = runif(6) #p-values for GxY effect
year=1
for(year in 1:6){
  cat("FillPipeline year:",year,"of 6\n")

  #Year 1
  MaleF1 = randCross(MaleParents,nCrosses)
  FemaleF1 = randCross(FemaleParents,nCrosses)

  MaleDH = makeDH(MaleF1,nDH)
  FemaleDH = makeDH(FemaleF1,nDH)

  #Year 2
  if(year<6){
    p = P[6-year]
    
    MaleYT1 = MaleDH
    FemaleYT1 = FemaleDH
    
    MaleHybridYT1 = hybridCross(FemaleTester1, MaleYT1)
    FemaleHybridYT1 = hybridCross(FemaleYT1, MaleTester1)
    
    MaleHybridYT1 = setPheno(MaleHybridYT1, varE = varE, reps = repYT1, p = P[year])
    FemaleHybridYT1 = setPheno(FemaleHybridYT1, varE = varE, reps = repYT1, p = P[year])
    
    MaleYT1@pheno = as.matrix(calcGCA(MaleHybridYT1)$GCAm[,2])
    FemaleYT1@pheno = as.matrix(calcGCA(FemaleHybridYT1)$GCAf[,2])
    
    
  
 }

  #Year 3
  if(year<5){
    p = P[5-year]

    MaleYT2 = selectInd(selectWithinFam(MaleYT1, famMaxYT1, use = "pheno"), 
                        nInbred2, use = "pheno")
    
    FemaleYT2 = selectInd(selectWithinFam(FemaleYT1, famMaxYT1, 
                                          use = "pheno"), nInbred2, use = "pheno")
    
    
    MaleHybridYT2 = hybridCross(FemaleTester2, MaleYT2)
    FemaleHybridYT2 = hybridCross(FemaleYT2, MaleTester2)
    
    MaleHybridYT2 = setPheno(MaleHybridYT2, varE = varE, reps = repYT2, p = P[year])
    FemaleHybridYT2 = setPheno(FemaleHybridYT2, varE = varE, reps = repYT2, p = P[year])
    
    MaleYT2@pheno = as.matrix(calcGCA(MaleHybridYT2)$GCAm[,2])
    FemaleYT2@pheno = as.matrix(calcGCA(FemaleHybridYT2)$GCAf[,2])
    

  }

  #Year 4
  if(year<4){
    p = P[4-year]

    MaleYT3 = selectInd(MaleYT2, nInbred2, use = "pheno")
    FemaleYT3 = selectInd(FemaleYT2, nInbred2, use = "pheno")
    
    MaleHybridYT3 = hybridCross(FemaleTester3, MaleYT3)
    FemaleHybridYT3 = hybridCross(FemaleYT3, MaleTester3)
    
    MaleHybridYT3 = setPheno(MaleHybridYT3, varE = varE, reps = repYT3, p = P[year])
    FemaleHybridYT3 = setPheno(FemaleHybridYT3, varE = varE, reps = repYT3, 
                               p = P[year])
    
    MaleYT3@pheno = as.matrix(calcGCA(MaleHybridYT3)$GCAm[,2])
    FemaleYT3@pheno = as.matrix(calcGCA(FemaleHybridYT3)$GCAf[,2])
    

  }

  #Year 5
  if(year<3){
    p = P[3-year]

    MaleHybridYT4 = selectInd(MaleHybridYT3,nYT4)
    FemaleHybridYT4 = selectInd(FemaleHybridYT3,nYT4)
    
    MaleHybridYT4 = setPheno(MaleHybridYT4,varE=varE,reps=repYT4,p=P[year])
    FemaleHybridYT4 = setPheno(FemaleHybridYT4,varE=varE,reps=repYT4,p=P[year])
    
    MaleYT4 = MaleYT3[
      MaleYT3@id%in%MaleHybridYT4@father
    ]
    
    FemaleYT4 = FemaleYT3[
      FemaleYT3@id%in%FemaleHybridYT4@mother
    ]

  }

  #Year 6
  if(year<2){
    p = P[2-year]

    MaleHybridYT5 = selectInd(MaleHybridYT4,nYT5, use="pheno")
    FemaleHybridYT5 = selectInd(FemaleHybridYT4,nYT5, use="pheno")
    
    MaleHybridYT5 = setPheno(MaleHybridYT5,varE=varE,
                             reps=repYT5,p=P[year])
    FemaleHybridYT5 = setPheno(FemaleHybridYT5,varE=varE,
                               reps=repYT5,p=P[year])
    
    MaleYT5 = MaleYT4[
      MaleYT4@id%in%MaleHybridYT5@father
    ]
    
    FemaleYT5 = FemaleYT4[
      FemaleYT4@id%in%FemaleHybridYT5@mother
    ]

  }

  #Year 7, release
}










