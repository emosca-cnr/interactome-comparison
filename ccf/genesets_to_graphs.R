#CALCULATION OF THE CCF of a list of gene sets in an interactome
#INPUT
#G: igraph object of the interactome
#gsl: named list of gene sets; gene identifiers must match those in V(G)$name

library(igraph)

#gene set lenght
gsl_size <- unlist(lapply(gsl, length))

#initialize the output 
G_gs <- vector("list", length(gsl))
names(G_gs) <- names(gsl)
G_gs_cl <- G_gs
G_gs_stats <- G_gs

cat("induced subgraphs")
G_gs <- lapply(gsl, function(x) igraph::induced_subgraph(G, V(G)$name[ V(G)$name %in% x] ))
cat("\n")

cat("clusters")
G_gs_cl <- lapply(G_gs, function(x) igraph::clusters(x))
cat("\n")

cat("output table")
G_gs_stats <- data.frame(id=names(gsl_size), V0=gsl_size, V=unlist(lapply(G_gs, function(x) length(V(x)))), stringsAsFactors = F)
G_gs_stats$Vcc <- unlist(lapply(G_gs_cl, function(x) sum(x$csize[x$csize>1])))
G_gs_stats$d <- lapply(G_gs, graph.density)
G_gs_stats$CCF <-  ifelse(G_gs_stats$V0 > 0, G_gs_stats$Vcc / G_gs_stats$V0, 0)

write.table(G_gs_stats, file="ccf.txt", sep="\t", row.names = F)


