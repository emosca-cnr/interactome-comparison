#pathway CROSS-talk
#require the network diffusion function ND
#gsl names list of gene sets
#W normalized adjacency matrix

#set up input matrix genes-by-pathways
X0 <- matrix(0, nrow=nrow(W), ncol = length(gsl), dimnames = list(rownames(W), names(gsl)))
cat("\tdefining X0\n")
for(j in 1:ncol(X0)){
	#rows of X0 that correspond to genes of the j-pathway
	X0[rownames(X0) %in% unlist(gsl[ names(gsl)==colnames(X0)[j] ]), j] <- 1
}

#network diffusion
cat("\tND")
if(!identical(rownames(X0), rownames(W))){
	stop("rownames of X0 and W must be identical")
}
Xs <- ND(X0, W)$Xs

#normalization of each pathway profile (by column)
#this will guarantee that PCT + t(PCT) are on the same scale
Xs_norm <- apply(Xs, 2, function(x) x / ifelse(any(x>0), sum(x), 1))

#calculation of pathway cross-talk
PCT <- matrix(0, ncol = ncol(X0), nrow = ncol(X0), dimnames = list(colnames(X0), colnames(X0)))
for(i in 1:ncol(X0)){
	cat(i)
	sum_X0_i <- sum(X0[, i]) #number of seeds
	if(sum_X0_i > 0){
		idx_X0_i <- which(X0[, i]==1) #index of seeds
		for(j in 1:ncol(X0)){
			#total fluid from pathway j in the genes of pathway i / #genes of i
			PCT[i, j] <- sum(Xs_norm[idx_X0_i, j]) / sum_X0_i
		}
	}
}
cat("\n")

#mean similarity
PCT <- (PCT + t(PCT)) / 2

write.table(PCT, file = "PCT.txt", sep="\t", row.names = T, col.names = NA)
