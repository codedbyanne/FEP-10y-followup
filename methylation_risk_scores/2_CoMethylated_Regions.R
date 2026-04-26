#################################################################################################
# Purpose of script: Find co-methylation regions for TOP123
# This script is adapted from a script available from https://pubmed.ncbi.nlm.nih.gov/36908043/
# https://github.com/jche453/Pruning-Thresholding-MRS.git
################################################################################################

##########################################################
###Run cmr function to generate Co-methylation regions####
##########################################################
### Generate co-methylated regions ###
# install.packages("comeback_0.1.0.tar.gz", repos = NULL, type="source")
#comeback_0.1.0.tar.gz can be found at https://bitbucket.org/flopflip/comeback.
#library(comeback)
source("scripts/comback_mod.R")

#####################################################################################################
# Finding comethylated regions from reference dataset TOP 123
#####################################################################################################

load("/tsd/p33/data/durable/groups/biostat/genetics/methylation/postQC/NORMENT/european/data/Residualised_Mvalues_AllSamples_residSexAgeSmokingscoreEcells_2geneticPCs.RData")

# Find co-methylated regions
res = t(top123_resid) #column are CpG sites, rows are for samples
rm(top123_resid);gc()
cmrs <- cmr(Mdata = res, corlo = 0.3, Iarray = "EPIC")
save(cmrs, file="analyses/mrs/CMR_TOP123.RData")
