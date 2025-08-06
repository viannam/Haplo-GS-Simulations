###>>>---------------------------
###> 1. Creating the initial genome
###>>>---------------------------

# Generate initial haplotypes/genome
founderpop = runMacs(nInd=nParents*2,
                     nChr=10,
                     segSites=nQtl+nSnp,
                     inbred=TRUE,
                     split=nGenSplit,
                     species="MAIZE")

###>>>---------------------------
###> 2. Setting traits characteristics
###>>>---------------------------

SP = SimParam$new(founderpop)
SP$restrSegSites(nQtl,nSnp)

if(nSnp>0){
  SP$addSnpChip(nSnp)

}

SP$addTraitADG(nQtl,
               mean=initMeanG,var=initVarG,
               varGxE=initVarGE,
               meanDD=ddMean,varDD=ddVar)

SP$setVarE(varE=varE)

SP$setTrackPed(TRUE)

SP$setTrackRec(TRUE)

###>>>---------------------------
###> 3. Creating the base population
###>>>---------------------------

# Split heterotic pools to form initial parents
FemaleParents = newPop(founderpop[1:nParents])
MaleParents = newPop(founderpop[(nParents+1):(nParents*2)])

#Set hybrid parents for later yield trials
MaleElite = selectInd(MaleParents,nElite,use="gv")
FemaleElite = selectInd(FemaleParents,nElite,use="gv")

#Reverse order to keep best parent in longer
MaleElite = MaleElite[nElite:1]
FemaleElite = FemaleElite[nElite:1]

#Set initial testers for YT1 and YT2
#Requires nTesters to be smaller than nElite
MaleTester1 = MaleElite[1:nTester1]
FemaleTester1 = FemaleElite[1:nTester1]
MaleTester2 = MaleElite[1:nTester2]
FemaleTester2 = FemaleElite[1:nTester2]
MaleTester3 = MaleElite[1:nTester3]
FemaleTester3 = FemaleElite[1:nTester3]

