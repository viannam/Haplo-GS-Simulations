#Advance breeding program by 1 year
#Works backwards through pipeline to avoid copying data

####>>>>>-------------- Year 7
# Release top hybrids hybrid

#Year 6
MaleHybridYT5 = selectInd(MaleHybridYT4,nYT5)
FemaleHybridYT5 = selectInd(FemaleHybridYT4,nYT5)

MaleHybridYT5 = setPheno(MaleHybridYT5,reps=repYT5,p=p,varE = varE)
FemaleHybridYT5 = setPheno(FemaleHybridYT5,reps=repYT5,p=p,varE = varE)

MaleYT5 = MaleYT4[
  MaleYT4@id%in%MaleHybridYT5@father
]
FemaleYT5 = FemaleYT4[
  FemaleYT4@id%in%FemaleHybridYT5@mother
]

HybridPop = c(MaleHybridYT5, FemaleHybridYT5) 

#Year 5
MaleHybridYT4 = selectInd(MaleHybridYT3,nYT4)
FemaleHybridYT4 = selectInd(FemaleHybridYT3,nYT4)

MaleHybridYT4 = setPheno(MaleHybridYT4,reps=repYT4,p=p,varE = varE)
FemaleHybridYT4 = setPheno(FemaleHybridYT4,reps=repYT4,p=p,varE = varE)

MaleYT4 = MaleYT3[
  MaleYT3@id%in%MaleHybridYT4@father
]
FemaleYT4 = FemaleYT3[
  FemaleYT3@id%in%FemaleHybridYT4@mother
]

####>>>>>-------------- Year 4
MaleYT3 = selectInd(MaleYT2, nInbred2, use = "pheno")
FemaleYT3 = selectInd(FemaleYT2, nInbred2, use = "pheno")

MaleHybridYT3 = hybridCross(FemaleTester3, MaleYT3)
FemaleHybridYT3 = hybridCross(FemaleYT3, MaleTester3)

MaleHybridYT3 = setPheno(MaleHybridYT3, varE = varE, reps = repYT3, p = P[year])
FemaleHybridYT3 = setPheno(FemaleHybridYT3, varE = varE, reps = repYT3, p = P[year])

MaleYT3@pheno = as.matrix(calcGCA(MaleHybridYT3)$GCAm[,2])
FemaleYT3@pheno = as.matrix(calcGCA(FemaleHybridYT3)$GCAf[,2])


####>>>>>-------------- Year 3
MaleYT2 = selectInd(selectWithinFam(MaleYT1, famMaxYT1, use = "pheno"), nInbred2, use = "pheno")
FemaleYT2 = selectInd(selectWithinFam(FemaleYT1, famMaxYT1, use = "pheno"), nInbred2, use = "pheno")

MaleHybridYT2 = hybridCross(FemaleTester2, MaleYT2)
FemaleHybridYT2 = hybridCross(FemaleYT2, MaleTester2)

MaleHybridYT2 = setPheno(MaleHybridYT2, varE = varE, reps = repYT2, p = P[year])
FemaleHybridYT2 = setPheno(FemaleHybridYT2, varE = varE, reps = repYT2, p = P[year])

MaleYT2@pheno = as.matrix(calcGCA(MaleHybridYT2)$GCAm[,2])
FemaleYT2@pheno = as.matrix(calcGCA(FemaleHybridYT2)$GCAf[,2])


####>>>>>-------------- Year 2
MaleYT1 = MaleDH
FemaleYT1 = FemaleDH

MaleHybridYT1 = hybridCross(FemaleTester1, MaleYT1)
FemaleHybridYT1 = hybridCross(FemaleYT1, MaleTester1)

MaleHybridYT1 = setPheno(MaleHybridYT1, varE = varE, reps = repYT1, p = P[year])
FemaleHybridYT1 = setPheno(FemaleHybridYT1, varE = varE, reps = repYT1, p = P[year])

MaleYT1@pheno = as.matrix(calcGCA(MaleHybridYT1)$GCAm[,2])
FemaleYT1@pheno = as.matrix(calcGCA(FemaleHybridYT1)$GCAf[,2])


####>>>>>-------------- Year 1
MaleF1 = randCross2(MaleYT3, c(MaleYT4,MaleYT5), nCrosses)
FemaleF1 = randCross2(FemaleYT3, c(FemaleYT4,FemaleYT5), nCrosses)

MaleDH = makeDH(MaleF1,nDH)
FemaleDH = makeDH(FemaleF1,nDH)

