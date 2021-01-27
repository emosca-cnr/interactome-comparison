#' Network Diffusion using
#' @param X0 vector or matrix composed of column vectors with initial distribution of information
#' @param W symmetrically normalized adjacency matrix W = D^-1 A D^-1, see normalize_adj_mat function
#' @param alpha numeric, the smothing factor
#' @param nMax numeric, maximum number of iterations
#' @param eps numeric, the iteration will stop when the maximum difference between matrix Xs between two consecutive iteraction is smaller than \code{eps}
#' @param final.smooth TRUE/FALSE, whether to do the final step of smoothing
#' @param all.steps, TRUE/FALSE, whether to store all steps
#' @param verbose, TRUE/FALSE
#' @return a list with:
#' \itemize{
#' \item{\code{Xt}}{ the smoothed matrix;}
#' \item{\code{eps}}{ see above;}
#' \item{\code{max.abs.diff}}{ max(abs(F_t) - abs(F_{t-1}));}
#' \item{\code{Xs.all}}{ transient Xs matrices.}
#' }
#' @export
#'

ND <- function(X0, W, alpha=0.7, nMax=1e4, eps=1e-6, final.smooth=FALSE, all.steps=FALSE, verbose=FALSE){

  Xs <- X0
  Fprev <- X0

  if(all.steps){
    Xs.all <- list()
    Xs.all[[1]] <- X0
  }

  #X0 and A multiplied by their weight
  X0a <- (1 - alpha) * X0
  Wa <- alpha * W

  for(i in 2:nMax){

    if(i %% 5 == 0 & verbose)
      cat(i, " ")

    #current iteration
    Xs <- Wa %*%  Fprev + X0a

    if(all.steps){
      Xs.all[[i]] <- Xs
    }

    max.abs.diff <- max(abs(Xs-Fprev))

    #update Fprev for next iteration
    Fprev <- Xs

    if(max.abs.diff < eps){
      if(final.smooth){
        Xs <- Wa %*%  Fprev
      }
      break
    }

  }

  if(verbose)
    cat('\n')

  if(all.steps){
    return(list(Xs=Xs, eps=eps, max.abs.diff=max.abs.diff, Xs.all=Xs.all))
  }else{
    return(list(Xs=Xs, eps=eps, max.abs.diff=max.abs.diff))
  }
}



