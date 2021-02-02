load("input_per_AUC_calc.RData") #obtained by means of 1_performance_assessment.R

#create list results - 1 list 1 interactome
list_int <- vector("list", length = length(int_results))
for(i in 1: length(int_results)){
  list_disease <- vector("list", length = length(int_results[[i]]))

    for(j in 1:length(int_results[[i]])){
    roc_TP <- vector("list", length = length(int_results[[i]][[j]]))
    roc_TN <- vector("list", length = length(int_results[[i]][[j]]))
    max_dim <- max(unlist(lapply(int_results[[i]][[j]], function(x) dim(x)[1])))
    #calculate roc curve for several points
    seq_n <- floor(seq(8, max_dim, length.out = floor(max_dim/10)))

      for(k in 1:length(int_results[[i]][[j]])){
      temp<- int_results[[i]][[j]][[k]]
      temp_roc <- roc_funct(names(int_results[[i]][[j]][k]), seq_n, temp, table(temp[,"db"])[1])
      roc_TP[[k]] <- temp_roc[[1]]
      roc_TN[[k]] <- temp_roc[[2]]
      }
    
    roc_TP <- do.call("rbind", roc_TP)
    roc_TN <- do.call("rbind", roc_TN)
    list_disease[[j]] <- list(roc_TP, roc_TN)
    }
  names(list_disease) <- names(int_results[[i]])
  list_int[[i]] <- list_disease
}
names(list_int) <- names(int_results)

save(list_int, file = "roc_curve_calc.RData")
