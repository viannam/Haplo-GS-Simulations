#Advance breeding program by 1 year
#Works backwards through pipeline to avoid copying data

########-----------------  Year 4  -------------########
# fall
PublicHybridYT3 = hybridCross(Tester3_private, PublicInbredYT3) # Fall

# Winter - target environment
PublicHybridYT3 = setPheno(PublicHybridYT3,varE=VarE, reps=repTC3, p=pgy) # Winter
PublicInbredYT3@pheno = as.matrix(calcGCA(PublicHybridYT3)$GCAm[,2])
PublicHybridYT4 = selectInd(PublicHybridYT3, 2, use="pheno") # Winter

HybridMean = PublicHybridYT4
DHTrain3 = PublicInbredYT3 # For the GS model


########-----------------  Year 3  -------------########
# Fall
PublicHybridYT2 = hybridCross(Tester2_private, PublicInbredYT2) #Fall

# Winter - target environment
PublicHybridYT2 = setPheno(PublicHybridYT2,varE=VarE, reps=repTC2,p=pgy) # Winter
PublicInbredYT2@pheno = as.matrix(calcGCA(PublicHybridYT2)$GCAm[,2])
PublicInbredYT3 = selectInd(PublicInbredYT2, 20, use="pheno") # Winter

PublicInbredYT2Phen = PublicInbredYT2 # Parental Selection
DHTrain2 = PublicInbredYT2 # For the GS model

########-----------------  Year 2  -------------########
# Fall
PublicInbredYT1 = setPheno(PublicInbredYT1,reps=1,p=pgye, varE = VarE) # Fall
PublicInbredYT1Sel = selectInd(PublicInbredYT1, nInd=800) # Fall
PublicHybridYT1 = hybridCross(Tester1_private, PublicInbredYT1Sel) # Fall

# Winter - target environment
PublicHybridYT1 = setPheno(PublicHybridYT1,varE=VarE, reps=repTC1, p=pgy) # Winter
PublicInbredYT1Sel@pheno = as.matrix(calcGCA(PublicHybridYT1)$GCAm[,2]) # Winter
PublicInbredYT2 = selectInd(PublicInbredYT1Sel,  150, use = "pheno")
DHTrain1 = PublicInbredYT1Sel # For the GS model

########-----------------  Year 1  -------------########
# Spring
F1_public = randCross(Parents_public_update, nCrosses_public, ignoreSexes = TRUE)

# Fall and Winter = Make DH
PublicDH = makeDH(F1_public, nDH = pub_DH, keepParents = FALSE) 
PublicInbredYT1 = PublicDH
