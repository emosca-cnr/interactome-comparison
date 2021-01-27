#pathway CROSS-talk
rotate90dx <- function(x) t(apply(x, 2, rev))

in_files <- list.files(pattern = "kegg_ND", full.names = T)
int_names <- sapply(in_files, function(x) gsub(".*kegg_ND_(.+)\\.RData$", "\\1", x))
pathND <- vector("list", length(int_names))
names(pathND) <- int_names
X0_at_least_10gg <- pathND
for(i in 1:length(in_files)){
	load(in_files[i])
	pathND[[i]] <- list(X0=X0, Xs=Xs)
	X0_at_least_10gg[[i]] <- colSums(X0)>=10 #pathways with at least 10 seeds
}
rm(X0, Xs)


S <- pathND
Sm <- S
for(k in 1:length(pathND)){
	cat(names(pathND)[k], "\n")
	temp_X0 <- pathND[[k]]$X0
	
	#normalization of each pathway profile (by column)
	#this will guarantee that S[[k]] + t(S[[k]] are on the same scale
	temp_Xs <- apply(pathND[[k]]$Xs, 2, function(x) x / ifelse(any(x>0), sum(x), 1))
	
	S[[k]] <- matrix(0, ncol = ncol(temp_X0), nrow = ncol(temp_X0), dimnames = list(colnames(temp_X0), colnames(temp_X0)))
	Sm[[k]] <- S[[k]]
	for(i in 1:ncol(temp_X0)){
		cat(i)
		sum_X0_i <- sum(temp_X0[, i]) #number of seeds
		if(sum_X0_i > 0){
			idx_X0_i <- which(temp_X0[, i]==1) #index of seeds
			for(j in 1:ncol(temp_X0)){
				#total fluid from pathway j in the genes of pathway i / #genes of i
				S[[k]][i, j] <- sum(temp_Xs[idx_X0_i, j]) / sum_X0_i
			}
		}
	}
	cat("\n")
	
	#mean similarity
	Sm[[k]] <- (S[[k]] + t(S[[k]])) / 2
}

#normalize Sm by the total of the column: fraction of PC over the total PC
Smn <- lapply(Sm, function(X) apply(X, 2, function(x) x / ifelse(any(x>0), sum(x), 1)))

#save similarity matrix to file
save(S, Sm, Smn, X0_at_least_10gg, file="PC.RData", compress = T) #modified 30 may; added S



all_val <- as.numeric(unlist(Smn))
min_val <- min(all_val[all_val>0]) / 2
max_val <- max(all_val)
rm(all_val)

jpeg("PC_path_dist.jpg", width = 200, height = 250, res=300, units="mm")
par(mfrow=c(5, 4))
par(mar=c(1,1,1,1))

for(i in 1:length(Smn)){
	cat(i)
	X <- Smn[[i]]
	X[X==0] <- min_val
	X <- log10(X)
	image(X, main=names(Smn)[i], col=heat.colors(32), zlim=c(log10(min_val), log10(max_val)))
}
plot.new()
legend("center", legend = c("low sim", "", "", "", "high sim"), col = heat.colors(5), pch=rep(15, 5), cex=0.8)
dev.off()

jpeg("PC_path_dist_density.jpg", width = 200, height = 250, res=300, units="mm")
par(mfrow=c(5, 4))
par(mar=c(1,1,1,1))
for(i in 1:length(Smn)){
	cat(i)
	X <- Smn[[i]]
	X[X==0] <- min_val
	X <- log10(X)
	plot(density(X), ylim=c(0, 1.5), xlab="log10(Smn)", col=1, main=names(Smn)[i], xlim=c(log10(min_val), log10(max_val)))
}
dev.off()

load("~/db/ncbi/2019_02_26/biosystems.RData")
kegg_cat <- read.csv("~/db/kegg/kegg_categories3.txt", stringsAsFactors = F, sep="\t")
kegg_cat$category <- factor(kegg_cat$category)
kegg_cat$category_number <- as.numeric(kegg_cat$category)
kegg_cat$bsid <- bsid$bs_id[match(kegg_cat$source_acc, bsid$source_acc)]

col <- c("green3", "dodgerblue", "hotpink", "darkgray", "brown", "black")

library(RColorBrewer)
jpeg("PC_path_dist_symmetric.jpg", width = 200, height = 250, res=300, units="mm")

lay_mat_temp <- matrix(1:4, nrow = 2, byrow = T)
lay_mat_temp <- cbind(lay_mat_temp, 5:6)
lay_mat <- lay_mat_temp
for(i in 1:4){
	lay_mat <- rbind(lay_mat, lay_mat_temp+6*i)
}
lay_mat_temp <- lay_mat
for(i in 1:3){
	lay_mat <- cbind(lay_mat, lay_mat_temp+30*i)
}
print(lay_mat)

lo <- layout(lay_mat, heights = c(rep(c(0.95, 0.05), time=5)/5), widths = rep(c(0.05, 0.90, 0.05), times=4)/4)
layout.show(lo)

overall_min <- unlist(Sm)
overall_min <- min(overall_min[overall_min>0]) / 2

X_to_plot <- Sm
for(i in 1:length(Sm)){
	
	X <- Sm[[i]]
	#cat(i, temp_min, temp_max, "\n")
	cat(i)
	
	temp_idx <- which(!X0_at_least_10gg[[i]])
	X[temp_idx, ] <- 0 #set to zero elements of pathway lt 10
	X[, temp_idx] <- 0 #set to zero elements of pathway lt 10
	
	temp_idx <- order(kegg_cat$category[match(rownames(X), kegg_cat$bs_id)])
	X <- X[temp_idx, temp_idx]
	X[X==0] <- overall_min
	X <- log10(X)
	
	X_to_plot[[i]] <- X
	
}


overall_min <- unlist(X_to_plot)
overall_max <- max(overall_min)
overall_min <- min(overall_min)

overall_min_rs <- unlist(lapply(X_to_plot, function(X) log10(rowSums(10^X))))
overall_max_rs <- max(overall_min_rs)
overall_min_rs <- min(overall_min_rs)

for(i in 1:length(Sm)){
	
	X <- X_to_plot[[i]]
	#1
	par(mar=c(0, 0.1, 1, 0))
	image(rbind(1:nrow(X)), axes=F, col=col[kegg_cat$category_number[match(rev(rownames(X)), kegg_cat$bsid)]])
	
	#2
	par(mar=c(0, 0, 1, 0))
	image(rotate90dx(X), main=names(Sm)[i], cex.main=0.7, col=(heat.colors(7)), xaxt="none", yaxt="none", zlim=c(overall_min, overall_max))
	
	#3
	par(mar=c(0, 0, 0, 0))
	plot.new()
	
	#4
	par(mar=c(0.1, 0, 0, 0))
	image(cbind(1:nrow(X)), axes=F, col=col[kegg_cat$category_number[match(colnames(X), kegg_cat$bsid)]])
	
	
	#5
	rS <- log10(rowSums(10^X))
	par(mar=c(0, 0, 1, 0.1))
	image(rotate90dx(matrix(rS, ncol=1)), axes=F, col=brewer.pal(5, "Greens"))# , zlim=c(overall_min_rs, overall_max_rs))
	
	#6
	par(mar=c(0, 0, 0, 0))
	plot.new()
	
	rm(rS)
}

#1
par(mar=c(0, 0, 0, 0))
plot.new()
#2
plot.new()
legend("topleft", legend = rev(round(seq(overall_min, overall_max, length.out = 7), 1)), col = rev(heat.colors(7)), pch=rep(15, 5), cex=0.8, title = "Similarity")
legend("topright", legend = rev(c("very low", "low", "mean", "high", "very high")), col = rev(brewer.pal(5, "Greens")), pch=rep(15, 5), cex=0.8, title = "Avg. Similarity")

legend("bottom", legend = levels(kegg_cat$category), col = col, pch=rep(15, max(kegg_cat$category_number)), cex=0.8, title = "KEGG Pathway")


#3-6
plot.new()
plot.new()
plot.new()
plot.new()

dev.off()
