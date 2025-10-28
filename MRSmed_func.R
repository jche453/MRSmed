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
    
    pvalueinfo = matrix(NA, 200,3)
    colnames(pvalueinfo) = c("top  CpG sites", "Number of CpG sites", 
                             "Numeber of CpG sites after matching with DNA methylation data")
    n = 0
    for (i in seq(0.000005,0.001,0.000005)){
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
    
  }else if (method == "meta"){
    SS_merge$fisher.p = NA
    for (i in 1:nrow(SS_merge)){
      SS_merge$fisher.p[i] = fisher(c(SS_merge$Pvalue.x[i], SS_merge$Pvalue.y[i]), adjust = "none")$p
    }
    
    pvalueinfo = matrix(NA, length(Pthred),3)
    colnames(pvalueinfo) = c("Pvalue", "Number of CpG sites", 
                             "Numeber of CpG sites after matching with DNA methylation data")
    for (i in Pthred){
      pvalue = 5 * 10 ^ (-i)
      pvalueinfo[i-1,1] =  pvalue
      SS_sub = SS_merge[SS_merge$fisher.p < pvalue, ]
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
  }
  results = list(pvalueinfo = pvalueinfo, MRS = MRS_final)
  return(results)
}

GenMRS = function(beta, SS, Pthred, CoMeRegion = NULL, CoMeBack = T, weightSE = F){
  pvalueinfo = matrix(NA, length(Pthred),5)
  colnames(pvalueinfo) = c("Pvalue", "Number of CpG sites", "Numeber of CpG sites after matching with CoMeRegion", " Number of CpG sites after pruning", "Numeber of CpG sites after matching with DNA methylation data")
  if (weightSE == T){
    SS$BETA = SS$BETA/SS$SE
  }
  for (i in Pthred){
    pvalue = 5 * 10 ^ (-i)
    pvalueinfo[i-1,1] =  pvalue
    SS_sub = SS[SS$Pvalue < pvalue, ]
    pvalueinfo[i-1,2] =  nrow(SS_sub)
    
    if (CoMeBack){
      SS_final = merge(SS_sub, CoMeRegion, by.x = "Marker", by.y = "clustercpg")
      SS_final = data.frame(SS_final)
      pvalueinfo[i-1,3] =  nrow(SS_final)
      SS_final = SS_final %>% 
        group_by(CoMethylCluster) %>% 
        slice(which.min(Pvalue))
      pvalueinfo[i-1,4] =  nrow(SS_final)
    }else{
      #SS_final = merge(SS_sub, CoMeRegion, by.x = "Marker", by.y = "clustercpg")
      SS_final = data.frame(SS_sub)
    }
    betas_final = data.frame(t(beta[,SS_final$Marker, drop = F]))
    betas_final$Marker = rownames(betas_final)
    mat_final = merge(betas_final, SS_final, by = "Marker")
    pvalueinfo[i-1,5] =  nrow(mat_final)
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
  results = list(pvalueinfo = pvalueinfo, MRS = MRS_final)
  return(results)
}



