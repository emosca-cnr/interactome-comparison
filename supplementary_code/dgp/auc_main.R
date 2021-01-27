auc_main <- function(list_int, k_fold = 5, FP = 0.2){
  
  #AUC function
  AUC_calc <-function(list,j, perc ) {
    x <- list[[2]][j,]
    y <- list[[1]][j,]
    n <- length(x[x<=perc])
    x <- x[1:n]
    y <- y[1:n]
    auc_tot <- sum(diff(x)*zoo::rollmean(y,2))
    return(auc_tot)
  }
  
  #n_disease
  n_disease <- length(list_int[[1]])
  
  #create auc results list
  auc_list_res <- vector("list", length = length(list_int))
  
  #interactomes
  for(i in 1:length(list_int)){
    auc_tot <- matrix(0, nrow = n_disease, ncol = k_fold)
    
    #diseases
    for(j in 1: length(list_int[[i]])){
      list_disease <- list_int[[i]][[j]]
      #k-fold
      for(k in 1:k_fold){
        auc_tot[j,k] <- AUC_calc(list_disease, k, FP)
      }
    }
    
    #disease x k-fold matrix
    colnames(auc_tot) <- paste0("test",seq(1, k_fold, 1))
    rownames(auc_tot) <- names(list_int[[i]])
    
    #calculate mean
    auc_tot <- cbind(auc_tot, apply(auc_tot, 1, mean))
    auc_tot <- round(auc_tot, 4)
    colnames(auc_tot)[dim(auc_tot)[2]] <- "mean"
    auc_list_res[[i]] <- auc_tot
  }
  
  names(auc_list_res) <- names(list_int)
  return(auc_list_res)
}
