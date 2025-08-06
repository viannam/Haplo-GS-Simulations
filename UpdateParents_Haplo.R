# Male parents

MaleDH@ebv = as.matrix(runGBLUP(trainPop = trainPop, targetPop = MaleDH, 
                                  Model = "Haplo", nTraits = 1))                                  #***
MaleParents = selectInd(MaleDH, nParents, use = 'ebv')

# Female parents
FemaleDH@ebv = as.matrix(runGBLUP(trainPop = trainPopF, targetPop = FemaleDH, 
                                            Model = "Haplo", nTraits = 1))                       #***
FemaleParents = selectInd(FemaleDH, nParents, use = 'ebv')


