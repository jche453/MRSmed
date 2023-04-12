GenMRS_mediation = function(beta, SS_merge, Pthred, method = "Threshold"){
  if (method == "Threshold"){
    pvalueinfo = matrix(NA, length(Pthred),3)
    colnames(pvalueinfo) = c("Pvalue", "Number of CpG sites", 
                             "Numeber of CpG sites after matching with DNA methylation data")
    for (i in Pthred){
      pvalue = 5 * 10 ^ (-i)
      pvalueinfo[i-1,1] =  pvalue
      SS_sub = SS_merge[SS_merge$Pvalue.x < pvalue & SS_merge$Pvalue.y < pvalue, ]
      pvalueinfo[i-1,2] =  nrow(SS_sub)
      
      SS_final = data.frame(SS_sub)
      
      betas_final = data.frame(t(beta[,SS_final$Marker, drop = F]))
      betas_final$Marker = rownames(betas_final)
      mat_final = merge(betas_final, SS_final, by = "Marker")
      pvalueinfo[i-1,3] =  nrow(mat_final)
      if (nrow(mat_final) > 0){
        MRS = data.frame(t(as.matrix(mat_final[,2:(ncol(betas_final))])) %*% (mat_final$BETA))
        colnames(MRS) = c(paste0("P",pvalue))
        MRS$ID = rownames(MRS)
        if (i == Pthred[1]){
          MRS_final = MRS
        }else{
          MRS_final = merge(MRS, MRS_final, by = "ID", all = T)
        }
      }
    }
    
  }else if (method == "ranking"){
    SS_merge = SS_merge[order(SS_merge$Pvalue.x),]
    SS_merge$rank.x = c(1:nrow(SS_merge))
    
    SS_merge = SS_merge[order(SS_merge$Pvalue.y),]
    SS_merge$rank.y = c(1:nrow(SS_merge))
    
    SS_merge$rank = SS_merge$rank.x + SS_merge$rank.y
    SS_merge = SS_merge[order(SS_merge$rank),]
    
    pvalueinfo = matrix(NA, 100,3)
    colnames(pvalueinfo) = c("top  CpG sites", "Number of CpG sites", 
                             "Numeber of CpG sites after matching with DNA methylation data")
    n = 0
    for (i in seq(0.000005,0.01,0.0001)){
      n = n + 1
      top = round(i * nrow(SS_merge))
      pvalueinfo[n,1] =  top
      SS_sub = SS_merge[1:top, ]
      pvalueinfo[n,2] =  nrow(SS_sub)
      
      SS_final = data.frame(SS_sub)
      
      betas_final = data.frame(t(beta[,SS_final$Marker, drop = F]))
      betas_final$Marker = rownames(betas_final)
      mat_final = merge(betas_final, SS_final, by = "Marker")
      pvalueinfo[n,3] =  nrow(mat_final)
      if (nrow(mat_final) > 0){
        MRS = data.frame(t(as.matrix(mat_final[,2:(ncol(betas_final))])) %*% (mat_final$BETA))
        colnames(MRS) = c(paste0("top_",i))
        MRS$ID = rownames(MRS)
        if (i == 0.000005){
          MRS_final = MRS
        }else{
          MRS_final = merge(MRS, MRS_final, by = "ID", all = T)
        }
      }
    }
    
  }
  results = list(pvalueinfo = pvalueinfo, MRS = MRS_final)
  return(results)
}
