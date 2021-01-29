# SIMILARITY NETWORK #

# INPUT
# C_mat_list: list of correlation matrices

# load required libraries
library(plotrix)
library(ggplot2)
library(igraph)


# define useful functions for matrix statistics computation
get_average <- function(mat){
  res <- rep(0, nrow(mat))
  names(res) <- colnames(mat)
  for (i in 1:nrow(mat)){
    idx_i <- (1:nrow(mat))!=i
    res[i] <- mean(mat[i,idx_i])
  }
  return(round(res, 3))
}

get_max <- function(mat){
  res <- rep(0, nrow(mat))
  names(res) <- colnames(mat)
  for (i in 1:nrow(mat)){
    idx_i <- (1:nrow(mat))!=i
    res[i] <- max(mat[i,idx_i])
  }
  return(round(res, 3))
}

get_top_int <- function(mat, n){
  res <- list()
  for (i in 1:nrow(mat)){
    idx_i <- (1:nrow(mat))!=i
    res[[i]] <- mat[i,idx_i]
    res[[i]] <- res[[i]][order(-res[[i]])]
    res[[i]] <- res[[i]][1:n]
    res[[i]] <- round(res[[i]],3)
  }
  return(res)
}

get_low_int <- function(mat, n){
  res <- list()
  for (i in 1:nrow(mat)){
    idx_i <- (1:nrow(mat))!=i
    res[[i]] <- mat[i,idx_i]
    res[[i]] <- res[[i]][order(res[[i]])]
    res[[i]] <- res[[i]][1:n]
    res[[i]] <- round(res[[i]],3)
  }
  return(res)
}

get_top_int_names <- function(mat, n){
  res <- list()
  for (i in 1:nrow(mat)){
    idx_i <- (1:nrow(mat))!=i
    res[[i]] <- mat[i,idx_i]
    res[[i]] <- res[[i]][order(-res[[i]])]
    res[[i]] <- res[[i]][1:n]
    res[[i]] <- round(res[[i]],3)
    res[[i]] <- names(res[[i]])
    res[[i]] <- paste0(res[[i]], collapse = ", ")
  }
  return(unlist(res))
}

get_low_int_names <- function(mat, n){
  res <- list()
  for (i in 1:nrow(mat)){
    idx_i <- (1:nrow(mat))!=i
    res[[i]] <- mat[i,idx_i]
    res[[i]] <- res[[i]][order(res[[i]])]
    res[[i]] <- res[[i]][1:n]
    res[[i]] <- round(res[[i]],3)
    res[[i]] <- names(res[[i]])
    res[[i]] <- paste0(res[[i]], collapse = ", ")
  }
  return(unlist(res))
}

# construction of network of interactomes

N_int <- ncol(C_mat_list[[1]])
	
# aggregate adjacency matrix
a_rri <- C_mat_list[[1]]
for(i in 2:length(C_mat_list)){
	a_rri <- a_rri + C_mat_list[[i]]
}
correct_names <- rownames(a_rri) # correct ordering of interactome names

# consider 4 closer neighbors for each interactome
nNeig <- 4 # 

# get the list of nNeigh closer neighbors for each interactome
g_rri_list <- list()
for (i in 1:N_int){
  g_rri_list[[i]] <- matrix(0, nNeig, 3)
  g_rri_list[[i]][, 1] <- rep(correct_names[i], nNeig) 
  g_rri_list[[i]][, 2] <- names(get_top_int(mat = a_rri, n = nNeig)[[i]])
  g_rri_list[[i]][, 3] <- get_top_int(mat = a_rri, n = nNeig)[[i]]
}

# initialize symmetric adyacency matrix
a_rri_adj <- matrix(0, N_int, N_int)
rownames(a_rri_adj) <- rownames(a_rri)
colnames(a_rri_adj) <- colnames(a_rri)

for (i in 1:N_int){
	a_rri_adj[i, match(g_rri_list[[i]][, 2], colnames(a_rri_adj))] <- as.numeric(g_rri_list[[i]][, 3])
}

a_rri_adj <- apply(a_rri_adj, MARGIN = 2, FUN = as.numeric)

a_rri <- a_rri_adj # a_rri: asymmetric adjacency matrix

a_rri_adj <- a_rri + t(a_rri) # a_rri_adj: symmetric adjacency matrix

# get directed and undirected graphs
g_rri <- graph.adjacency(a_rri, weighted = T, mode = "directed")
g_rri_adj <- graph.adjacency(a_rri_adj, weighted = T, mode = "plus")

# get communities based on undirected network modularity
rri_comm <- fastgreedy.community(g_rri_adj)

