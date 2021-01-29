#' aggregate_correlation
#' @param C_t_list list of correlation matrices

aggregate_correlation <- function(C_t_list){

	R <- numeric(ncol(C_t_list[[1]]))
	
	for(i in 1:length(C_t_list)){
		diag(C_t_list[[i]]) <- 0
	}
	
	Rt <- do.call(cbind, lapply(C_t_list, rowSums))
	
	R <- rowSums(Rt)
	
	return(list(R=R, Rt=Rt))
	
}