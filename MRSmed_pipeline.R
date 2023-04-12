setwd("/projects/huels_lab/MRS/02-Mediation/02-Simulations/01-Simulation1/02-E_MRS_O/simu_pheno_1")
### Load summary statistics, one for exposure EWAS and the other for outcome EWAS
load("SS.rda")
#SS_E denotes Exposure EWAS summary statistics
#SS_O denotes Outcomet EWAS summary statistics
#Summary statistics files need to include at least four coloumns named "Marker", "BETA", "SE", "Pvalue"
head(SS_E)
#  Marker       BETA          SE       Pvalue
#cg14327296 0.08903298 0.006227317 3.221431e-44
#cg23800778 0.10504546 0.008370557 7.971134e-35
#cg15137445 0.08265012 0.006900022 5.608122e-32
#cg13720022 0.12456846 0.010650265 1.296429e-30
#cg09975093 0.08201602 0.007022830 1.577188e-30
#cg04234631 0.08762171 0.007921750 1.212095e-27

head(SS_O)
#Marker       BETA          SE       Pvalue
#cg05752685 1.21150417 0.077415494 3.711070e-52
#cg13720022 0.12004906 0.010437600 1.091975e-29
#cg01387407 0.66939844 0.063716411 3.632927e-25
#cg14444099 0.63053177 0.063251524 7.094095e-23
#cg15137445 0.06723396 0.006829825 2.324909e-22
#cg06092310 0.09639765 0.010000103 1.588056e-21

#Mediation MRS
SS_merge = merge(SS_E, SS_O, by = "Marker") 
SS_merge$BETA = SS_merge$BETA.x*SS_merge$BETA.y
#SS_merge$BETA = SS_merge$BETA.x/SS_merge$SE.x*SS_merge$BETA.y/SS_merge$SE.y

#Get the smallest p-value
if (min(SS_merge$Pvalue.x) < min(SS_merge$Pvalue.y)){
  minpvalue = min(SS_merge$Pvalue.y)
}else{
  minpvalue = min(SS_merge$Pvalue.x)
}
minpvalue = sapply(strsplit(as.character(minpvalue), "-"), "[", 2)

Pthred = 2:minpvalue

#Load beta value
load("/projects/huels_lab/MRS/02-Mediation/01-Data/02-SimulationData/GSE55763_normalized_M.rda") #418418
load("Betas_causal.rda")
cpg = colnames(Betas_log2)
cpg_causal = colnames(Betas_causal)
beta = cbind(Betas_log2[,cpg[!(cpg%in%cpg_causal)]], Betas_causal)

MRS_list = GenMRS_mediation(beta, SS_merge, Pthred, method = "ranking")
pvalueinfo = MRS_list$pvalueinfo
MRS = MRS_list$MRS
MRS$ID = sapply(strsplit(MRS$ID, "X"), "[", 2)
write.csv(pvalueinfo, "pvalueinfo.csv", row.names = F)
write.csv(MRS, "MRS.mediation.csv", row.names = F)
