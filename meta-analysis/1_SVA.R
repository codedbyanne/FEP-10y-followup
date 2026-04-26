#################################################################################################
# Purpose of script: Surrogate Variable Analysis
# @author: Anne-Kristin Stavrum
################################################################################################

library(sva)

########################################################################
# SVA on the Female sample
########################################################################

load("data/Betas_SCZ_CC_excl_10y_followup_Females_incl_samplesheet.RData")

mod = model.matrix(~ Case_Control + smokingscore + AgeAtSampling + Bcell+Neu+CD4T+CD8T+NK+Mono+
                     cpPC1+cpPC2+cpPC3+cpPC4+cpPC5+gPC1+gPC2+gPC3+gPC4+gPC5, data=samplesheet_F)
mod0 = model.matrix(~ smokingscore + AgeAtSampling + Bcell+Neu+CD4T+CD8T+NK+Mono+
                      cpPC1+cpPC2+cpPC3+cpPC4+cpPC5+gPC1+gPC2+gPC3+gPC4+gPC5, data=samplesheet_F)

n.sv = num.sv(betaF,mod,method="be")

svobj = sva(betaF,mod,mod0,n.sv=n.sv,method="irw")
rownames(svobj$sv) = colnames(betaF)

save(svobj,n.sv, file="analyses/ewas_TOP123_scz_excluding_10y_followup/SVsOnBetas_methodIRW_SCZfemales_CCsmokingscoreAgeatsamplingEcells_cpPC1-5_gPC1-5.Robj")

rm(list=ls());gc()

########################################################################
# SVA on the male sample
########################################################################

load("data/Betas_SCZ_CC_excl_10y_followup_Males_incl_samplesheet.RData")

mod = model.matrix(~ Case_Control + smokingscore + AgeAtSampling + Bcell+Neu+CD4T+CD8T+NK+Mono+
                     cpPC1+cpPC2+cpPC3+cpPC4+cpPC5+gPC1+gPC2+gPC3+gPC4+gPC5, data=samplesheet_M)
mod0 = model.matrix(~ smokingscore + AgeAtSampling +Bcell+Neu+CD4T+CD8T+NK+Mono+cpPC1+cpPC2+cpPC3+cpPC4+cpPC5+
                      gPC1+gPC2+gPC3+gPC4+gPC5, data=samplesheet_M)

n.sv = num.sv(betaM,mod,method="be") # without mPCs: 113 | with 3 mPCs: 195

svobj = sva(betaM,mod,mod0,n.sv=n.sv,method="irw")
rownames(svobj$sv) = colnames(betaM)

save(svobj,n.sv, file="analyses/ewas_TOP123_scz_excluding_10y_followup/SVsOnBetas_methodIRW_SCZmales_CCsmokingscoreAgeatsamplingEcells_cpPC1-5_gPC1-5.Robj")

rm(list=ls());gc()