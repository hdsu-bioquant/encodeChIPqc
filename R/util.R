#' Internal version of snow::clusterApplyLB
#'
#' \code{.clusterApplyLB} returns the result of snow::clusterApplyLB().
#'
#' The function works in very much the same way as \code{snow::clusterApplyLB()}.
#' The only difference is that it only uses a \code{snow} cluster, if
#' \code{cluster} is not \code{NULL}. Otherwise, it uses the function
#' \code{lapply()}.
#'
#' @param cluster \code{snow} cluster generated with \code{cluster <- snow::makeCluster(code)} or \code{NULL}
#' @param x List of objects which \code{clusterApplyLB()} shall be applied to
#' @param ... Additional parameters to \code{snow::clusterApplyLB()}
#' @return The result of \code{snow::clusterApplyLB()} or \code{lapply()}
#' @examples
#' # single-core invocation
#' .clusterApplyLB(NULL, c(1,2,3), identical, 1)
#'
#' # parallel invocation
#' cluster <- makeCluster(3)
#' .clusterApplyLB(cluster, c(1,2,3), identical, 1)
.clusterApplyLB <- function(cluster, x, ...) {
	# make sure x is subsettable (i.e. a list), so that we can use clusterApplyLB on it
	x <- tryCatch({ x[[1]]; x }, error=function(err) list(x))

	if (is.null(cluster)) {
		return(lapply(x, ...))
	} else {
		require(snow)
	        return(clusterApplyLB(cluster, x, ...))
	}
}

