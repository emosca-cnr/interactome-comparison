roc_funct<-function(col_an, n_value, temp, TP){
	
  bpTP <- matrix(0,nrow = length(col_an), ncol = length(n_value)+1)
  bpTN <- matrix(0,nrow = length(col_an), ncol = length(n_value)+1)
  rownames(bpTP) <- col_an
  colnames(bpTP) <- c(0,n_value)
  rownames(bpTN) <- col_an
  colnames(bpTN) <- c(0,n_value)
  for(k in 1:length(col_an)){
    temp <- temp[order(temp[,col_an[k]], decreasing = T),]
    for(m in 1:length(n_value)){
      n_s <- n_value[m]
      tempTab <- table(temp[1:min(n_s, dim(temp)[1]),"db"])
      ind_T <- which(names(tempTab)=="0")
      ind_N <- which(names(tempTab)=="1")        
      if(length(ind_T)>0){
        bpTP[k,m+1] <- tempTab[ind_T]/TP
      } else {
        bpTP[k,m+1] <- 0
      }
      if(length(ind_N)>0){
        bpTN[k,m+1] <- tempTab[ind_N]/(dim(temp)[1]-TP)
      } else {
        bpTN[k,m+1] <- 0
      }
    }
  }
  
  bp <- list(bpTP, bpTN)
  return(bp)
}