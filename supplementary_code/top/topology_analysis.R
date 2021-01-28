###################################################################################
### Interactome Reproducibility: Topological Analysis #############################
###################################################################################

library(igraph)
library(rARPACK)
library(Hmisc)
library(gplots)
library(ggplot2) 
library(LaplacesDemon)
library(org.Hs.eg.db)


sym2eg <- as.data.frame(org.Hs.egSYMBOL2EG)

interactomes <- c('bioplex','cofrac15','huri','qubic','ghiassian','hint','hippie','irefindex','ncbi',
                'fp60','multinet','string400notm','string400','string700notm','string700',
                'biana','cpathdb','inbiowebmap','intact')

col_1 <- c('red','red','red','red','blue','blue','blue',
         'blue','blue','black','black','black','black','black','black','blue','black','black','blue')

labels <- c('BX','CF','HURI','QU','DMND','HN','HP','IR','NCBI',
            'FP60','MN','S04','S04T','S07','S07T','BN','CP','IBMP','INTC')     


G <- list()
Gnodes <- list()
N <- list()
E <- list()


# import interactomes
for (interactome in interactomes){
  filename <- paste("path_to_interactomes_edgelists",interactome,".txt",sep="")
  edgelist <- read.table(filename,sep = "\t")
  G[[interactome]] <- graph_from_data_frame(edgelist, directed = FALSE, vertices = NULL)
  Gnodes[[interactome]] <- vertex_attr(G[[interactome]])
  Gnodes[[interactome]] <- data.frame(Gnodes[[interactome]])
  N[[interactome]] <- length(V(G[[interactome]]))
  E[[interactome]] <- length(E(G[[interactome]]))
}

#####################################################################################
#global properties
for (interactome in interactomes){
  #max dist. between two nodes
  G[[interactome]]$diameter <- diameter(G[[interactome]],directed=F)
  #fraction of triangles (global clust.coeff.)
  G[[interactome]]$transitivity <- transitivity(G[[interactome]],type='global')
  #mean dist. between nodes
  G[[interactome]]$mean_distance <- mean_distance (G[[interactome]])
}


#network densities
for(interactome in interactomes){
  G[[interactome]]$prop_dens <- 2*length(E(G[[interactome]]))/length(V(G[[interactome]]))^2
  G[[interactome]]$density <- length(E(G[[interactome]]))/length(V(G[[interactome]]))
}


###################### SPECTRAL CENTRALITY ##############################################
# k1 -> k-centrality
# L_sparse: graph laplacian in sparse form
# A_sparse: graph laplacian in sparse form
spect_centrality<-function (L_sparse,A_sparse,k1){
  ev<-eigs(L_sparse, k1+1, which = "LM", sigma = 0)
  vec<-ev$vectors[,1]
  s<-rep(0, length(vec))
  for (i in 1:length(vec)) {
    vec1<-rep(1, length(vec))*vec[i]-vec 
    vec1<-vec1*vec1
    ss<-A_sparse[i,]
    ss<-ss*vec1
    s[i]<-sum(ss)
  }
  s<-data.frame(s)
  rownames(s)<-rownames(A_sparse)
  colnames(s)<-c('spectral_1')
  return(s)
}


############################ Overall CENTRALITIES  ####################################

i_cent <- list()
for (interactome in interactomes){
  #betweenness
  betw <- data.frame(betweenness(G[[interactome]], v = V(G[[interactome]]), directed = FALSE, weights = NULL,normalized = TRUE))
  colnames(betw)<-c('betweenness')
  #degree
  deg <- data.frame(degree(G[[interactome]], v = V(G[[interactome]]), normalized = TRUE))
  colnames(deg)<-c('degree')
  #page-rank
  p_rank <- page_rank(G[[interactome]], vids = V(G[[interactome]]),directed = FALSE, damping = 0.85)
  p_rank <- data.frame(p_rank$vector)
  colnames(p_rank)<-c('p_rank')
  #1-spectral
  L_sparse <- laplacian_matrix(G[[interactome]], normalized = FALSE, weights = NULL, sparse = igraph_opt("sparsematrices"))
  A_sparse <- get.adjacency(G[[interactome]],type='both',attr=NULL, names=TRUE,sparse=TRUE)
  spectral_1 <- spect_centrality(L_sparse,A_sparse,1)
  #closeness
  clns <-data.frame(closeness(G[[interactome]], v = V(G[[interactome]]), normalized = TRUE))
  colnames(clns)<-c('closeness')
  
  #built dataframe of interactomes centrality measures
  centrality<-data.frame(Gnodes[[interactome]]$name,deg$degree,p_rank$p_rank,betw$betweenness,spectral_1$spectral_1,clns$closeness)
  colnames(centrality) <- c('gene_entrez','degree','p_rank','betweenness','spectral_1','closeness')
  i_cent[[interactome]] <- centrality
}


########################## Centralities Spearman Correlation ########################

measures <- c('degree','p_rank','betweenness','spectral_1', 'closeness')
spearman_corr <- list()
spearman_pvalue <- list()

for (cent in measures){
  corr <- matrix(0, length(interactomes),length(interactomes))
  pvalue <- matrix(0, length(interactomes),length(interactomes))
  i=1
  for (interactome_1 in interactomes){
    j=1
    for (interactome_2 in interactomes){  
      nodes_1 <- data.frame(list(i_cent[[interactome_1]]['gene_entrez']))
      nodes_2 <- data.frame(list(i_cent[[interactome_2]]['gene_entrez']))
      a <- match(nodes_1$gene_entrez, nodes_2$gene_entrez)
      inter_centr_1<-i_cent[[interactome_1]][[cent]][which(!is.na(a))]
      inter_centr_2<-i_cent[[interactome_2]][[cent]][a[!is.na(a)]]
      inter_centr<-data.frame(inter_centr_1,inter_centr_2)
      spearman<-rcorr(as.matrix(inter_centr), type="spearman") 
      corr[i,j]<-spearman[[1]][1,2]
      pvalue[i,j]<-spearman[[3]][1,2]
      j=j+1
    }
    i=i+1
  }
  spearman_corr[[cent]] <- corr
  spearman_pvalue[[cent]]<-pvalue
}
central <- data.frame(spearman_corr$spectral_1)
colnames(central)[1:length(interactomes)] <- interactomes
rownames(central)[1:length(interactomes)] <- interactomes


############################## Scale-Freeness #######################################

library(igraph)
library(poweRlaw)

# FIT
deg <- list()
deg.dist <- list()
m_pl <- list()
m_exp <- list()
m_ln <- list()
est_pl <- list()
est_exp <- list()
est_ln <- list()
deg.max <- list()

plot.deg <- list()
fit.deg <- list()
plot.deg_exp <- list()
fit.deg_exp <- list()
plot.deg_ln <- list()
fit.deg_ln <- list()
d_est <- list()
p <- list()
p_exp <- list()
d_est_exp <- list()
p_ln<- list()
d_est_ln <- list()



for(interactome in interactomes){

  g <-  G[[interactome]]
  deg[[interactome]] <- degree(g)
  deg[[interactome]] <- deg[[interactome]][deg[[interactome]]>0]
  deg.dist[[interactome]] <- data.frame(k=0:max(deg[[interactome]]),p_k=degree_distribution(g))
  deg.dist[[interactome]] <- deg.dist[[interactome]][deg.dist[[interactome]]$p_k>0,]
  
  # powerlaw
  m_pl[[interactome]] <- displ$new(deg[[interactome]])
  est_pl[[interactome]] <- estimate_xmin(m_pl[[interactome]])
  m_pl[[interactome]]$setXmin(est_pl[[interactome]])
  #plot.deg[[interactome]] <- plot(m_pl[[interactome]], draw = F)
  #fit.deg[[interactome]]  <- lines(m_pl[[interactome]], draw = F)
  
  # exponential
  m_exp[[interactome]] <- disexp$new(deg[[interactome]])
  est_exp[[interactome]] <- estimate_xmin(m_exp[[interactome]])
  m_exp[[interactome]]$setXmin(est_exp[[interactome]])

  # lognormal 
  m_ln[[interactome]] <- dislnorm$new(deg[[interactome]])
  est_ln[[interactome]] <- estimate_xmin(m_ln[[interactome]])
  m_ln[[interactome]]$setXmin(est_ln[[interactome]])
}
  
  
# Semiparametric Bootstrap  -> GOF
gof_pl <- list()
gof_ln <- list()
gof_exp <- list()
  
nsims <- 500
i=1
for(interactome in interactomes){
  gof_pl[[interactome]] <- bootstrap_p(m_pl[[interactome]],no_of_sims=nsims, threads=12, seed = 123)
  gof_exp[[interactome]] <- bootstrap_p(m_exp[[interactome]],no_of_sims=nsims, threads=12, seed = 123)
  gof_ln[[interactome]] <- bootstrap_p(m_ln[[interactome]],no_of_sims=nsims, threads=12, seed = 123)
  i=i+1
}
  
#bootstrap -> overall fit uncertainty
err_pl <- list()
err_exp <- list()
  
for(interactome in interactomes){
  err_pl[[interactome]] <- bootstrap(m_pl[[interactome]], no_of_sims = 100, threads = 12)
  err_exp[[interactome]] = bootstrap(m_exp[[interactome]], no_of_sims = 100, threads = 12)    err_ln[[interactome]] = bootstrap(m_ln[[interactome]], no_of_sims = 100, threads = 12)
}
  
  
# ##### fit fixed PL range
# m_exp_pl <-list()
# est_exp_pl <-list()
# m_ln_pl <-list()
# est_ln_pl <-list()
#   
# for(interactome in interactomes){
#   # exponential 
#   m_exp_pl[[interactome]] <- disexp$new(deg[[interactome]])
#   m_exp_pl[[interactome]]$setXmin(est_pl[[interactome]])
#   est_exp_pl[[interactome]] <- estimate_pars(m_exp_pl[[interactome]])
#   #
#   # lognormal 
#   m_ln_pl[[interactome]] <- dislnorm$new(deg[[interactome]])
#   m_ln_pl[[interactome]]$setXmin(est_pl[[interactome]])
#   est_ln_pl[[interactome]] <- estimate_pars(m_ln_pl[[interactome]])
# }
#   
# ######   compare distributions (same range)
# comp_exp <- list()
# comp_ln <- list()
# for(interactome in interactomes){
#   #to exp
#   m1 <- disexp$new(deg[[interactome]])
#   m1$setXmin(m_pl[[interactome]]$getXmin())
#   m1$setPars(estimate_pars(m1))
#   comp_exp[[interactome]] <- compare_distributions(m_pl[[interactome]], m1)
#   #to lognorm
#   m2 <- dislnorm$new(deg[[interactome]])
#   m2$setXmin(m_pl[[interactome]]$getXmin())
#   m2$setPars(estimate_pars(m2))
#   comp_ln[[interactome]] <- compare_distributions(m_pl[[interactome]], m2)
# }

  
  
  
  
  
  





