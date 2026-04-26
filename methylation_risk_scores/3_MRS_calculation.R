#################################################################################################
# Purpose of script: Calculate Methylation Scores
# This script is adapted from a script available from https://pubmed.ncbi.nlm.nih.gov/36908043/
# https://github.com/jche453/Pruning-Thresholding-MRS.git
################################################################################################

#####################################
##############STEP ONE###############
#####################################

#####################################################################################################
# Format the Co-methylation regions dataframe for later analysis - TOP123
#####################################################################################################
# Load residualised data TOP123
load("data/Residualised_Mvalues_TOP123_TenYearFollowupParticipants.Robj")
res = t(res)
dim(res)
load("analyses/mrs/CMR_TOP123.RData")
## Format the Co-methylation regions dataframe for later analysis, 
## you can either input cmrs calculated from your own dataset,
#CoMeRegion = GenCoMeRegion(cmrs = cmrs, beta = res, Overlap = F)
## or a cmrs calculated from a reference dataset,
CoMeRegion = GenCoMeRegion(beta = res, reference = cmrs, Overlap = F)
## or both and then overlap the cmrs
#CoMeRegion = GenCoMeRegion(cmrs = cmrs, beta = res, reference = cmrs_ref, Overlap = T)

#CoMeRegion is a matrix that assigned a unique number to each co-methylation region
save(CoMeRegion, file = "analyses/mrs/CoMeRegion_TenYearFolloup_refTOP123.rda")

#####################################################################################################
# Format the Co-methylation regions dataframe for later analysis - TOP4
#####################################################################################################
# Load residualised data TOP4. These were extracted from the residualisation of the full TOP4 dataset
# https://github.com/codedbyanne/BIP-meta/blob/main/meta_analysis_pipeline/9B_Residualise_TOP4.R
load("data/Residualised_Mvalues_TOP4_TenYearFollowupParticipants.Robj")
res = t(res)
dim(res)
# Load Reference CMRs - previously calculated for the BIP project: https://github.com/codedbyanne/BIP-meta/blob/main/meta_analysis_pipeline/9C_MRS_pipeline.R
load("Bipolar_meta_EWAS/analyses/mrs/CMR_EPICv1v2.TOP4.RData")
## Format the Co-methylation regions dataframe for later analysis, 
## you can either input cmrs calculated from your own dataset,
#CoMeRegion = GenCoMeRegion(cmrs = cmrs, beta = res, Overlap = F)
## or a cmrs calculated from a reference dataset,
CoMeRegion = GenCoMeRegion(beta = res, reference = cmrs, Overlap = F)
## or both and then overlap the cmrs
#CoMeRegion = GenCoMeRegion(cmrs = cmrs, beta = res, reference = cmrs_ref, Overlap = T)

#CoMeRegion is a matrix that assigned a unique number to each co-methylation region
save(CoMeRegion, file = "analyses/mrs/CoMeRegion_TenYearFolloup_refTOP4.rda")



#####################################
############## STEP TWO #############
#####################################

### Calculate MRS
library(dplyr) 
source("scripts/MRS_func.R")

#######################################################################################
# Calculation of MS for FEP-10y-followup samples extracted from TOP123 
#######################################################################################
# Load real phenotype and methylation data
DNAm = get(load("data/Residualised_Mvalues_TOP123_TenYearFollowupParticipants.Robj"))
rm(res);gc()
DNAm = t(DNAm)
DNAm[1:6,1:6] #check DNAm file

###Load SCZ summary statistics from discovery dataset
ss <- get(load("analyses/meta_analysis/Meta_SCZ_excl_10yfollowup.RData"))
ssn = as.data.frame(apply(ss,2,as.numeric))
rownames(ssn) = rownames(ss)
ssn = ssn[which(!is.na(ssn$All_Effect_Fixed)),]
ssn$Marker = rownames(ssn)
rownames(ssn) = 1:length(ssn$Marker)
SS = ssn[,c(12,5,6,7)]
colnames(SS) = c("Marker", "BETA", "SE", "Pvalue")
head(SS)
rm(ss,ssn);gc()
#Get the smallest p-value
minpvalue = min(SS$Pvalue[SS$Pvalue != 0])
minpvalue = sapply(strsplit(as.character(minpvalue), "-"), "[", 2)
###Load Co-methylation regions for newborns -> CoMeRegion
load("analyses/mrs/CoMeRegion_TenYearFolloup_refTOP123.rda")
#Specify how p-value threshold, for example, if you want 5 * 10 ^ (-2), specify pthread to be 2
Pthred = 2:minpvalue
MRS = GenMRS(DNAm, SS, Pthred, CoMeRegion, CoMeBack = T, weightSE = F)
#if weightSE = T, weights = BETA/SE, where BETA is the effect size
#Basic information of MRS
write.csv(MRS$pvalueinfo, "analyses/mrs/MRS_TenYearFollowupSamples_TOP123_pvalueinfo.csv", row.names = F)
write.csv(MRS$MRS, "analyses/mrs/MRS_TenYearFollowupSamples_TOP123.csv", row.names = F)

#######################################################################################
# Calculation of MS for FEP-10y-followup samples extracted from TOP4
#######################################################################################
# Load real phenotype and methylation data
DNAm = get(load("data/Residualised_Mvalues_TOP4_TenYearFollowupParticipants.Robj"))
rm(res);gc()
DNAm = t(DNAm)
DNAm[1:6,1:6] #check DNAm file

###Load SCZ summary statistics from discovery dataset
ss <- get(load("analyses/meta_analysis/Meta_SCZ_excl_10yfollowup.RData"))
ssn = as.data.frame(apply(ss,2,as.numeric))
rownames(ssn) = rownames(ss)
ssn = ssn[which(!is.na(ssn$All_Effect_Fixed)),]
ssn$Marker = rownames(ssn)
rownames(ssn) = 1:length(ssn$Marker)
SS = ssn[,c(12,5,6,7)]
colnames(SS) = c("Marker", "BETA", "SE", "Pvalue")
head(SS)
rm(ss,ssn);gc()

#Get the smallest p-value
minpvalue = min(SS$Pvalue[SS$Pvalue != 0])
minpvalue = sapply(strsplit(as.character(minpvalue), "-"), "[", 2)
###Load Co-methylation regions for newborns -> CoMeRegion
load("analyses/mrs/CoMeRegion_TenYearFolloup_refTOP4.rda")
#Specify how p-value threshold, for example, if you want 5 * 10 ^ (-2), specify pthread to be 2
Pthred = 2:minpvalue
MRS = GenMRS(DNAm, SS, Pthred, CoMeRegion, CoMeBack = T, weightSE = F)
#if weightSE = T, weights = BETA/SE, where BETA is the effect size
#Basic information of MRS
write.csv(MRS$pvalueinfo, "analyses/mrs/MRS_TenYearFollowupSamples_TOP4_pvalueinfo.csv", row.names = F)
write.csv(MRS$MRS, "analyses/mrs/MRS_TenYearFollowupSamples_TOP4.csv", row.names = F)
