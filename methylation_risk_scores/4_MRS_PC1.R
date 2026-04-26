###############################################################################################################
# Purpose of script: Calculate the PC1 of MRS calculated across different p-value thresholds
# Method: PRS-PCA used here is from the supplementary material of : https://onlinelibrary.wiley.com/doi/10.1002/gepi.22339
# The code was copied to the Util_functions.R file
# @author: Anne-Kristin Stavrum
###############################################################################################################

#####################################
# PC calculation
#####################################
source("scripts/Util_functions.R")

#########################################################################################################
# PCA on TOP123 and TOP4 combined
#########################################################################################################
MS123 = read.csv("analyses/mrs/MRS_TenYearFollowupSamples_TOP123.csv")
MS123$ID = sapply(MS123$ID, function(x) gsub("\\.","-",x))
rownames(MS123) = MS123$ID

MS4 = read.csv("analyses/mrs/MRS_TenYearFollowupSamples_TOP4.csv")
MS4$ID = gsub("X","",MS4$ID)
rownames(MS4) = MS4$ID

MS123 = MS123[,c(-2)]
MS = rbind(MS123,MS4)

PCA = prs.pc(as.data.frame(MS),"SCZ")
PC1 = PCA$data[,which(colnames(PCA$data) %in% c("ID","SCZ.prs.pc"))]
save(PC1, file="analyses/mrs/PC1_MS-SCZ_TOP1234.Robj")
