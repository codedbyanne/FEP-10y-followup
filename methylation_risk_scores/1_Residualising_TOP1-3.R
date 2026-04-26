#!/usr/bin/Rscript

#################################################################################################
# Purpose of script: Residualisation of TOP1-3.
# This script was run through slurm om the cluster
# @author: Anne-Kristin Stavrum
################################################################################################

#####################################################################
# Loading Libraries
#####################################################################
library(doParallel)
library(foreach)
library(meta)

#####################################################################
# Setting working directory and paths
#####################################################################

args = commandArgs(trailingOnly=TRUE)
workdir = args[1]
input1 = args[2] #file with betas
covariate_file = args[3] # data_frame_with_covariates
no.clusters = as.integer(args[4]) # number of threads for parallisation
#####################################################################
# Loading Data
#####################################################################
print("Loading data")
betas = get(load(input1))
rm(list=ls()[which(!ls() %in% c("betas","covariate_file","no.clusters","workdir"))])
variables_df = get(load(covariate_file))
print(ls())
print(paste("No of clusters: ",no.clusters))

#####################################################################
# Preparing 
#####################################################################
print("Setting up parallell")
cl = makeCluster(no.clusters)
registerDoParallel(cl)

#####################################################################
# Residualisation
#####################################################################
residualise = function(CpG, samplesheet){
  
  y = lm(CpG ~ Sex + AgeAtSampling + smokingscore +
           CD8T + CD4T + NK + Bcell + Mono + Neu +
           gPC1 + gPC2, data = samplesheet)
  
  return(resid(y))
  
}

print("Running the loop")
top123_resid = foreach(i=1:length(betas[,1]), .combine=rbind) %dopar% {
  residualise(CpG=betas[i,],samplesheet=variables_df)
}

## create matrix to store results
print("setting rownames")
rownames(top123_resid)<-rownames(betas)
print("setting colnames")
print(length(colnames(betas)))
print(length(colnames(top123_resid)))
head(colnames(betas))
head(colnames(top123_resid))
colnames(top123_resid)<-colnames(betas)


# Saving result
print("Saving result")
save(top123_resid,file=file.path(workdir,"TOP123_residualised_Mvals_SexAgeSmokingEcells_2geneticPCs.RData"))

print("Done")


