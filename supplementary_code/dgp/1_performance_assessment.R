########################################
######## PERFORMANCE ASSESSMENT ######## 
########################################

##### INPUT
#disease_gene_list <- 

library(igraph)

interattomi <- names(G)


## CREATE INPUT FOR DISEASE GENE PRIORITIZATION

disease_k_fold_int <- vector("list", length = length(interattomi))
df_list_int <- vector("list", length(interattomi))
X0_list_int <- vector("list", length(interattomi))

for(i in 1:length(G)){
	
	#gene id
	gene_id <- V(G[[i]])$name
	disease_k_fold <- vector("list",5)
	df_list_disease <- vector("list", 5)
	X0_disease <- vector("list", 5)
	
	#disease 
	for(j in 1:length(disease_gene_list)){
		gene_id_disease <- disease_gene_list[[j]]
		gene_com <- intersect(gene_id, gene_id_disease)
		
		#randomize input
		sampleIDX <- sample(1:length(gene_com), size = length(gene_com))
		gene_comSampl <- gene_com[sampleIDX]
		
		#divide size fold
		n_div <- floor(length(gene_comSampl)/5)
		largest_fold <- length(gene_comSampl)-n_div*4
		k_fold <- vector("list", 5)
		
		#k1
		k_fold[[1]] <- gene_comSampl[1:n_div]
		#k2
		k_fold[[2]] <- gene_comSampl[(n_div+1):(n_div*2)]
		#k3
		k_fold[[3]] <- gene_comSampl[(n_div*2+1):(n_div*3)]
		#k4
		k_fold[[4]] <- gene_comSampl[(n_div*3+1):(n_div*4)]
		#5
		k_fold[[5]] <- gene_comSampl[(n_div*4+1):length(gene_comSampl)]
		names(k_fold) <- c("k1","k2","k3","k4","k5")
		
		#check total number of genes
		dim_fold <- sum(unlist(lapply(k_fold, length))) - length(gene_comSampl)
		if(dim_fold>0){
			stop("wrong folds")
		}
		
		#check number fold
		df_list_disease[[j]] <- data.frame(int = interattomi[i],
																			 disease = names(disease_gene_list[j]), 
																			 k1 = length(k_fold[[1]]), 
																			 k2 = length(k_fold[[2]]),
																			 k3 = length(k_fold[[3]]), 
																			 k4 = length(k_fold[[4]]), 
																			 k5 = length(k_fold[[5]]), 
																			 tot = length(gene_comSampl), 
																			 stringsAsFactors = F)
		
		#create training and test sets
		seq1 <- c(2,3,4,5) #test 1
		seq2 <- c(1,3,4,5) #test 2
		seq3 <- c(1,2,4,5) #test 3
		seq4 <- c(1,2,3,5) #test 4
		seq5 <- c(1,2,3,4) #test 5
		
		training_list <- vector("list",5)
		training_list[[1]] <- unlist(k_fold[seq1])
		training_list[[2]] <- unlist(k_fold[seq2])
		training_list[[3]] <- unlist(k_fold[seq3])
		training_list[[4]] <- unlist(k_fold[seq4])
		training_list[[5]] <- unlist(k_fold[seq5])
		
		test_list <- vector("list",5)
		test_list[[1]] <- k_fold[[1]]
		test_list[[2]] <- k_fold[[2]]
		test_list[[3]] <- k_fold[[3]]
		test_list[[4]] <- k_fold[[4]]
		test_list[[5]] <- k_fold[[5]]
		
		names(training_list) <- seq(1:5)
		names(test_list) <- seq(1:5)
		
		#prepare X0 for j-th disease 
		mm <- matrix(0, nrow = length(gene_id), ncol = 5)
		rownames(mm) <- gene_id
		colnames(mm) <- c("test1","test2","test3","test4","test5")
		for(k in 1:5){
			mm[which(rownames(mm) %in% training_list[[k]]),k] <- 1
		}
		X0_disease[[j]] <- mm
		names(X0_disease)[j] <- names(disease_gene_list)[[j]]
		
		disease_k_fold[[j]]  <- list(training_list=training_list, test_list = test_list)
		names(disease_k_fold)[[j]] <- names(disease_gene_list)[[j]]
	}
	disease_k_fold_int[[i]] <- disease_k_fold
	names(disease_k_fold_int)[[i]] <- interattomi[i]
	
	X0_list_int[[i]] <- X0_disease
	names(X0_list_int)[[i]] <- interattomi[i]
	
	df_list_int[[i]] <- do.call("rbind", df_list_disease)
}

#save files
save(disease_k_fold_int, file="disease_k_fold.RData")
save(X0_list_int, file="X0_int.RData")



### DISEASE GENE PRIORITIZATION

load("X0_int.RData")

cores <- 4
Xs_int <- vector("list",19)

for(i in 1:length(G)){
	adj <- as.matrix(as_adjacency_matrix(G[[i]]))
	adj <- normalize_adj_mat(adj)
	X0_list <- X0_list_int[[i]]
	Xs_all <-  parallel::mclapply(X0_list, function(x) ND(x, adj)$Xs, mc.cores = cores)
	names(Xs_all) <- names(X0_list)
	Xs_int[[i]] <- Xs_all
}

names(Xs_int) <- names(X0_list_int)

save(Xs_int, file="ND_19int.RData") 



# MERGE DATA GENERATED ABOVE
load("disease_k_fold.RData")
load("ND_19int.RData")

'%!in%' <- function(x,y)!('%in%'(x,y))

int_results <- vector("list", length(disease_k_fold_int))
for(i in 1:length(disease_k_fold_int)){
	
	#setting (training/test) of i-th interactome
	setting <- disease_k_fold_int[[i]]
	#ND with the i-th interactome
	xs <- Xs_int[[i]]
	#list with pAUC for all disease, i-th interactome
	list_disease <- vector("list", length(setting))
	for(j in 1: length(xs)){
		#ND j-th disease
		xs_disease <- xs[[j]]
		
		#genes of training set
		setting_training <- setting[[j]]$training_list
		
		#genes to be predicted
		test_training <- setting[[j]]$test_list
		
		#matrix pAUC
		list_k <- vector("list", length = length(xs))
		for(k in 1:dim(xs_disease)[2]){
			
			#ND k-fold
			temp <- xs_disease[,k, drop=F]
			temp <- as.data.frame(temp[rownames(temp) %!in% setting_training[[k]],,drop=F])
			temp$gene <- rownames(temp)
			df <- data.frame(gene=test_training[[k]], db=0, stringsAsFactors = F)
			temp <- merge(temp, df, by="gene", all.x=T)
			temp[is.na(temp)] <- 1
			temp <- temp[order(temp[,2], decreasing = T),,drop=F]
			list_k[[k]] <- temp
		}
		names(list_k) <-  colnames(xs_disease)
		list_disease[[j]] <- list_k
	}
	names(list_disease) <- names(setting)
	int_results[[i]] <- list_disease
}
names(int_results) <- names(disease_k_fold_int)
save(int_results, file="input_per_AUC_calc.RData")


