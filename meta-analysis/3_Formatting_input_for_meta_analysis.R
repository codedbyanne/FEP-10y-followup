################################################################################################################################
# Purpose of script: Formatting summary statisstics for meta-analysis
# @author: Anne-Kristin Stavrum
################################################################################################################################

load("analyses/ewas_TOP123_scz_excluding_10y_followup/sumstat/Ewas_CC_TOP1-3_Females_Excl_FEP10y_smokingscoreAgeatsamplingEcells_cpPC1-5_gPC1-5_irwSV1-10.Robj")
TOP123_F = data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=(DMPs$CI.R-DMPs$CI.L)/3.92 )
load("analyses/ewas_TOP123_scz_excluding_10y_followup/sumstat/Ewas_CC_TOP1-3_Males_Excl_FEP10y_smokingscoreAgeatsamplingEcells_cpPC1-5_gPC1-5_irwSV1-5.Robj")
TOP123_M = data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=(DMPs$CI.R-DMPs$CI.L)/3.92 )

DMPs = read.csv("analyses/meta_analysis/sumstat/SCZ_ABR_FEM_MPC3_SVA10_Imput.all.cpg.bacon.csv",row.names = 1)
ABR_F =  data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=(DMPs$CI.R-DMPs$CI.L)/3.92)
DMPs = read.csv("data/sumstat/SCZ_ABR_MALE_MPC3_SVA10_Imput.all.cpg.bacon.cs",row.names = 1)
ABR_M =  data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=(DMPs$CI.R-DMPs$CI.L)/3.92)

DMPs = read.csv("analyses/meta_analysis/sumstat/SCZ_IOP_FEM_SVA10_MPC3_Imput.all.cpg.bacon.csv",row.names = 1)
IOP_F =  data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=(DMPs$CI.R-DMPs$CI.L)/3.92)
DMPs = read.csv("analyses/meta_analysis/sumstat/SCZ_IOP_MAL_SVA10_MPC3_Imput.all.cpg.bacon.csv",row.names = 1)
IOP_M =  data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=(DMPs$CI.R-DMPs$CI.L)/3.92)

DMPs = read.csv("analyses/meta_analysis/sumstat/SCZ_UCL_FEM_MPC3_SVA10_Imput.all.cpg.bacon.csv",row.names = 1)
UCL_F =  data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=(DMPs$CI.R-DMPs$CI.L)/3.92)
DMPs = read.csv("analyses/meta_analysis/sumstat/SCZ_UCL_MALE_MPC3_SVA10_Imput.all.cpg.bacon.csv",row.names = 1)
UCL_M =  data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=(DMPs$CI.R-DMPs$CI.L)/3.92)

res = list(TOP123_F,TOP123_M,ABR_F,ABR_M,IOP_F,IOP_M,UCL_F,UCL_M)
names(res) = c("TOP123_F","TOP123_M","ABR_F", "ABR_M","IOP_F","IOP_M","UCL_F","UCL_M")


save(res, file="analyses/meta_analysis/sumstat/ReadyForMeta_SCZ_excluding_10yFollowup.RData")

