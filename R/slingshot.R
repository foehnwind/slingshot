#' @rdname slingshot
#'   
#' @description Given a reduced-dimensional data matrix \code{n} by \code{p} and
#'   a vector of cluster labels (potentially including \code{-1}'s for 
#'   "unclustered"), this function performs lineage inference using a 
#'   cluster-based minimum spanning tree and constructing simulatenous principal
#'   curves for branching paths through the tree.
#'   
#' @description This wrapper function performs lineage inference in two steps:
#'   (1) identify lineage structure with a cluster-based minimum spanning tree
#'   with the \code{\link{getLineages}} function and (2) construct smooth
#'   representations of each lineage using simultaneous principal curves from
#'   the function \code{\link{getCurves}}.
#'   
#' @param reducedDim numeric matrix or \code{SlingshotDataSet} object containing
#'   low- dimensional representation of single cells.
#' @param clusterLabels character, a vector of length \code{n} denoting cluster
#'   labels, optionally including \code{-1}'s for "unclustered." If
#'   \code{reducedDim} is a \code{SlingshotDataSet}, cluster labels will be
#'   taken from it.
#' @param start.clus (optional) character, indicates the cluster(s) of origin.
#'   Lineages will be represented by paths coming out of this cluster.
#' @param end.clus (optional) character, indicates the cluster(s) which will be 
#'   forced leaf nodes. This introduces a constraint on the MST algorithm.
#' @param dist.fun (optional) function, method for calculating distances between
#'   clusters. Must take two matrices as input, corresponding to subsets of 
#'   \code{reducedDim}. If the minimum cluster size is larger than the number
#'   dimensions, the default is to use the joint covariance matrix to find 
#'   squared distance between cluster centers. If not, the default is to use the
#'   diagonal of the joint covariance matrix.
#' @param omega (optional) numeric, this granularity parameter determines the
#'   distance between every real cluster and the artificial cluster, OMEGA. It
#'   is parameterized such that this distance is \code{omega / 2}, making
#'   \code{omega} the maximum distance between two connected clusters. By
#'   default, \code{omega = Inf}.
#' @param lineages list, denotes which lineages each cluster is a part of and 
#'   contains the \code{K x K} connectivity matrix constructed on the clusters
#'   by \code{\link{getLineages}}.
#' @param thresh numeric, determines the convergence criterion. Percent change
#'   in the total distance from cells to their projections along curves must be
#'   less than \code{thresh}. Default is \code{0.001}, similar to 
#'   \code{\link{principal.curve}}.
#' @param maxit numeric, maximum number of iterations, see 
#'   \code{\link{principal.curve}}.
#' @param stretch numeric factor by which curves can be extrapolated beyond 
#'   endpoints. Default is \code{2}, see \code{\link{principal.curve}}.
#' @param smoother, choice of scatter plot smoother. Same as 
#'   \code{\link{principal.curve}}, but \code{"lowess"} option is replaced with 
#'   \code{"loess"} for additional flexibility.
#' @param shrink logical or numeric between 0 and 1, determines whether and how 
#'   much to shrink branching lineages toward their average prior to the split.
#' @param extend character, how to handle root and leaf clusters of lineages
#'   when constructing the initial, piece-wise linear curve. Accepted values are
#'   \code{'y'} (default), \code{'n'}, and \code{'pc1'}. See 'Details' for more.
#' @param reweight logical, whether to allow cells shared between lineages to be
#'   reweighted during curve-fitting. If \code{TRUE}, cells shared between 
#'   lineages will be weighted by: distance to nearest curve / distance to
#'   curve.
#' @param drop.multi logical, whether to drop shared cells from lineages which
#'   do not fit them well. If \code{TRUE}, shared cells with a distance to one 
#'   lineage above the 90th percentile and another below the 50th will be
#'   dropped from the further lineage.
#' @param shrink.method character denoting how to determine the appropriate
#'   amount of shrinkage for a branching lineage. Accepted values are the same
#'   as for \code{kernel} in \code{\link{density}} (default is \code{"cosine"}),
#'   as well as \code{"tricube"} and \code{"density"}. See 'Details' for more.
#' @param ... Additional parameters to pass to scatter plot smoothing function, 
#'   \code{smoother}.
#'   
#' @details The \code{connectivity} matrix is learned by fitting a (possibly
#'   constrained) minimum-spanning tree on the clusters and the artificial 
#'   cluster, OMEGA, which is a fixed distance away from every real cluster.
#'   This effectively limits the maximum branch length in the MST to twice the
#'   chosen distance, meaning that the output may contain multiple trees.
#'   
#' @details Once the \code{connectivity} is known, lineages are identified in
#'   any tree with at least two clusters. For a given tree, if there is an
#'   annotated starting cluster, every possible path out of a starting cluster
#'   and ending in a leaf that isn't another starting cluster will be returned.
#'   If no starting cluster is annotated, every leaf will be considered as a
#'   potential starting cluster and whichever configuration produces the longest
#'   average lineage length (in terms of number of clusters included) will be
#'   returned.
#'   
#' @details When there is only a single lineage, the curve-fitting algorithm is
#'  nearly identical to that of \code{\link{principal.curve}}. When there are 
#'  multiple lineages and \code{shrink == TRUE}, an additional step is added to
#'  the iterative procedure, forcing curves to be similar in the neighborhood
#' of shared points (ie., before they branch).
#' 
#' @details The \code{extend} argument determines how to construct the
#'   piece-wise linear curve used to initiate the recursive algorithm. The
#'   initial curve is always based on the lines between cluster centers and if
#'   \code{extend = 'n'}, this curve will terminate at the center of the
#'   endpoint clusters. Setting \code{extend = 'y'} will allow the first and
#'   last segments to extend beyond the cluster center to the orthogonal
#'   projection of the furthest point. Setting \code{extend = 'pc1'} is similar
#'   to \code{'y'}, but uses the first principal component of the cluster to
#'   determine the direction of the curve beyond the cluster center. These
#'   options typically have little to no impact on the final curve, but can
#'   occasionally help with stability issues.
#'   
#' @details When \code{shink == TRUE}, we compute a shrinkage curve,
#'   \eqn{w_l(t)}, for each lineage, a non-increasing function of pseudotime
#'   that determines how much that lineage should be shrunk toward a shared
#'   average curve. We set \eqn{w_l(0) = 1}, so that the curves will perfectly
#'   overlap the average curve at pseudotime \code{0}. The weighting curve
#'   decreases from \code{1} to \code{0} over the non-outlying pseudotime values
#'   of shared cells (where outliers are defined by the \code{1.5*IQR} rule).
#'   The exact shape of the curve in this region is controlled by
#'   \code{shrink.method}, and can follow the shape of any standard kernel
#'   function's cumulative density curve (or more precisely, survival curve,
#'   since we require a decreasing function). Different choices of
#'   \code{shrink.method} seem to have little impact on the final curves, in 
#'   most cases.
#'   
#' @references Hastie, T., and Stuetzle, W. (1989). "Principal Curves."
#'   \emph{Journal of the American Statistical Association}, 84:502–516.
#'   
#' @return An object of class \code{\link{SlingshotDataSet}} containing the 
#'   arguments provided to \code{slingshot} as well as the following output: 
#'   \itemize{ 
#'   \item{\code{lineages}}{ a list of \code{L} items, where \code{L} 
#'   is the number of lineages identified. Each lineage is represented by a 
#'   character vector with the names of the clusters included in that lineage, 
#'   in order.} 
#'   \item{\code{connectivity}}{ the inferred cluster connectivity 
#'   matrix.} 
#'   \item{\code{lineageControl}}{Additional parameters used for
#'   lineage inference or describing the process. This may include the elements
#'   \code{start.given} and \code{end.given}, logical values indicating whether
#'   the starting and ending clusters were specified a priori. Additionally,
#'   this will always include \code{dist}, the pairwise cluster distance
#'   matrix.} 
#'   \item{curves}{A list of \code{\link{principal.curve}} objects.}
#'   \item{curveControls}{Additional parameters used for fitting simultaneous
#'   principal curves.}}
#'   
#' @examples
#' data("slingshotExample")
#' sds <- slingshot(rd, cl, start.clus = '1')
#' 
#' plot(rd, col = cl, asp = 1)
#' lines(sds, lwd = 3)
#' 
#' @export
#' 
setMethod(f = "slingshot",
          signature = signature(reducedDim = "matrix", 
                                clusterLabels = "character"),
          definition = function(reducedDim, clusterLabels,
                                start.clus = NULL, end.clus = NULL,
                                dist.fun = NULL, omega = NULL,
                                lineages = list(),
                                shrink = TRUE,
                                extend = 'y',
                                reweight = TRUE,
                                drop.multi = TRUE,
                                thresh = 0.001, maxit = 15, stretch = 2,
                                smoother = 'smooth.spline',
                                shrink.method = 'cosine',
                                allow.breaks = TRUE, ...){
            sds <- getLineages(reducedDim, clusterLabels,
                               start.clus = start.clus, end.clus = end.clus,
                               dist.fun = dist.fun, omega = omega)
            sds <- getCurves(sds,
                             shrink = shrink, extend = extend,
                             reweight = reweight, drop.multi = drop.multi,
                             thresh = thresh, maxit = maxit,
                             stretch = stretch, smoother = smoother,
                             shrink.method = shrink.method,
                             allow.breaks = allow.breaks, ...)
            return(sds)
          }
)

#' @rdname slingshot
#' @export
setMethod(f = "slingshot",
          signature = signature(reducedDim = "SlingshotDataSet", 
                                clusterLabels = "ANY"),
          definition = function(reducedDim,
                                clusterLabels = reducedDim@clusterLabels,
                                start.clus = NULL, end.clus = NULL,
                                dist.fun = NULL, omega = NULL,
                                lineages = list(),
                                shrink = TRUE,
                                extend = 'y',
                                reweight = TRUE,
                                drop.multi = TRUE,
                                thresh = 0.001, maxit = 15, stretch = 2,
                                smoother = 'smooth.spline',
                                shrink.method = 'cosine', 
                                allow.breaks = TRUE, ...){
            return(slingshot(reducedDim = reducedDim@reducedDim,
                             clusterLabels = reducedDim@clusterLabels,
                             start.clus = start.clus, end.clus = end.clus,
                             dist.fun = dist.fun, omega = omega,
                             shrink = shrink, extend = extend,
                             reweight = reweight, drop.multi = drop.multi,
                             thresh = thresh, maxit = maxit,
                             stretch = stretch, smoother = smoother,
                             shrink.method = shrink.method,
                             allow.breaks = allow.breaks, ...))
          })

#' @rdname slingshot
#' @export
setMethod(f = "slingshot",
          signature = signature(reducedDim = "data.frame", 
                                clusterLabels = "ANY"),
          definition = function(reducedDim, clusterLabels,
                                start.clus = NULL, end.clus = NULL,
                                dist.fun = NULL, omega = NULL,
                                lineages = list(),
                                shrink = TRUE,
                                extend = 'y',
                                reweight = TRUE,
                                drop.multi = TRUE,
                                thresh = 0.001, maxit = 15, stretch = 2,
                                smoother = 'smooth.spline',
                                shrink.method = 'cosine',
                                allow.breaks = TRUE, ...){
            RD <- as.matrix(reducedDim)
            rownames(RD) <- rownames(reducedDim)
            return(slingshot(reducedDim = RD,
                             clusterLabels = clusterLabels,
                             start.clus = start.clus, end.clus = end.clus,
                             dist.fun = dist.fun, omega = omega,
                             shrink = shrink, extend = extend,
                             reweight = reweight, drop.multi = drop.multi,
                             thresh = thresh, maxit = maxit,
                             stretch = stretch, smoother = smoother,
                             shrink.method = shrink.method,
                             allow.breaks = allow.breaks, ...))
          })

#' @rdname slingshot
#' @export
setMethod(f = "slingshot",
          signature = signature(reducedDim = "matrix", 
                                clusterLabels = "numeric"),
          definition = function(reducedDim, clusterLabels,
                                start.clus = NULL, end.clus = NULL,
                                dist.fun = NULL, omega = NULL,
                                lineages = list(),
                                shrink = TRUE,
                                extend = 'y',
                                reweight = TRUE,
                                drop.multi = TRUE,
                                thresh = 0.001, maxit = 15, stretch = 2,
                                smoother = 'smooth.spline',
                                shrink.method = 'cosine', 
                                allow.breaks = TRUE, ...){
            return(slingshot(reducedDim = reducedDim,
                             clusterLabels = as.character(clusterLabels),
                             start.clus = start.clus, end.clus = end.clus,
                             dist.fun = dist.fun, omega = omega,
                             shrink = shrink, extend = extend,
                             reweight = reweight, drop.multi = drop.multi,
                             thresh = thresh, maxit = maxit,
                             stretch = stretch, smoother = smoother,
                             shrink.method = shrink.method, 
                             allow.breaks = allow.breaks, ...))
          })

#' @rdname slingshot
#' @export
setMethod(f = "slingshot",
          signature = signature(reducedDim = "matrix", 
                                clusterLabels = "factor"),
          definition = function(reducedDim, clusterLabels,
                                start.clus = NULL, end.clus = NULL,
                                dist.fun = NULL, omega = NULL,
                                lineages = list(),
                                shrink = TRUE,
                                extend = 'y',
                                reweight = TRUE,
                                drop.multi = TRUE,
                                thresh = 0.001, maxit = 15, stretch = 2,
                                smoother = 'smooth.spline',
                                shrink.method = 'cosine', 
                                allow.breaks = TRUE, ...){
            return(slingshot(reducedDim = reducedDim,
                             clusterLabels = as.character(clusterLabels),
                             start.clus = start.clus, end.clus = end.clus,
                             dist.fun = dist.fun, omega = omega,
                             shrink = shrink, extend = extend,
                             reweight = reweight, drop.multi = drop.multi,
                             thresh = thresh, maxit = maxit,
                             stretch = stretch, smoother = smoother,
                             shrink.method = shrink.method, 
                             allow.breaks = allow.breaks, ...))
          })


#' @rdname slingshot
#' @export
setMethod(f = "slingshot",
          signature = signature(reducedDim = "matrix", 
                                clusterLabels = "ANY"),
          definition = function(reducedDim, clusterLabels,
                                start.clus = NULL, end.clus = NULL,
                                dist.fun = NULL, omega = NULL,
                                lineages = list(),
                                shrink = TRUE,
                                extend = 'y',
                                reweight = TRUE,
                                drop.multi = TRUE,
                                thresh = 0.001, maxit = 15, stretch = 2,
                                smoother = 'smooth.spline',
                                shrink.method = 'cosine', 
                                allow.breaks = TRUE, ...){
            if(missing(clusterLabels)){
              message('Unclustered data detected.')
              clusterLabels <- rep('1', nrow(reducedDim))
            }
            return(slingshot(reducedDim = reducedDim,
                             clusterLabels = clusterLabels,
                             start.clus = start.clus, end.clus = end.clus,
                             dist.fun = dist.fun, omega = omega,
                             shrink = shrink, extend = extend,
                             reweight = reweight, drop.multi = drop.multi,
                             thresh = thresh, maxit = maxit,
                             stretch = stretch, smoother = smoother,
                             shrink.method = shrink.method, 
                             allow.breaks = allow.breaks, ...))
          })
