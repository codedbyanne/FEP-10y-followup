#!/usr/bin/Rscript

#################################################################################################
# Purpose of script: Meta-analysis. This script was run through slurm on the cluster
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
input1 = args[2] #file with summary stats for samples to be metaanalysed
outfile = args[3] # name of output file for meta-analysis
nProbes = as.integer(args[4]) # number cohorts that needs to have result for a particular probe
no.clusters = as.integer(args[5]) # number of threads for parallelisation
#####################################################################
# Loading Data
#####################################################################
print("Loading data")
load(input1)
print(ls())
print(paste("No of clusters: ",no.clusters))
print(paste("nProbes: ",nProbes))
gc()
#####################################################################
# Preparing 
#####################################################################
print("Setting up parallell")
cl = makeCluster(no.clusters)
registerDoParallel(cl)

#####################################################################
# Meta analysis
#####################################################################
print("Setting up analysis")
n_cohorts = length(res)
print(n_cohorts)
print("Getting total no of probes")
probes = unlist(sapply(res,rownames))
print(length(probes))
print("Getting no of each probe")
probes<-table(probes)
print("Getting the no of probes over nProbes")
probes<-probes[which(probes >= nProbes)]
length_probes = length(probes)
print(length_probes)
print(paste("No of probes to be analysed: ",length_probes))
## run meta-analysis

metaCpG = function(row,probeID){
  #errind = c()
  ## first run meta analysis with dublin included
  meanDiff <- sapply(res, function(x) x[[1]][match(probeID,rownames(x))])
  seDiff <- sapply(res, function(x) x[[2]][match(probeID,rownames(x))])
  
  
  tryCatch({
    out<- meta::metagen(meanDiff, seDiff)
    return(c(sum(!is.na(meanDiff)), out$TE.fixed,out$seTE.fixed,out$pval.fixed, out$TE.random, out$seTE.random,out$pval.random, out$tau, out$I2, out$Q,1-pchisq(out$Q, out$df.Q)))
  },
  error = function(e) print(paste("error",i,sep=": ")))
  
}

print("Running the loop")
res.meta = foreach(i=1:length(probes), .combine=rbind) %dopar% {
  probeID<-names(probes)[i]
  metaCpG(row=i,probeID=probeID)
}

## create matrix to store results
rownames(res.meta)<-names(probes)
colnames(res.meta)<-c("N_cohorts", "All_Effect_Fixed", "All_Effect_SE_Fixed", "All_P_Fixed", "All_Effect_Random", "All_Effect_SE_Random","All_P_Random", "All_tau", "All_I2", "All_Q", "All_Het P")


# Saving result
print("Saving result")
save(res.meta,file=file.path(workdir,outfile))

print("Done")
