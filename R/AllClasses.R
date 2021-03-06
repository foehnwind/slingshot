#' @title Class \code{SlingshotDataSet}
#' @aliases SlingshotDataSet SlingshotDataSet-class
#'   
#' @description The \code{SlingshotDataSet} class holds data relevant for 
#'   performing lineage inference with the \code{slingshot} package, primarily a
#'   reduced dimensional representation of the data and a set of cluster labels.
#'   All \code{slingshot} methods can take an object of the class 
#'   \code{SlingshotDataSet} as input and will output the same.
#'   
#' @slot reducedDim matrix. An \code{n} by \code{p} numeric matrix or data frame
#'   giving the coordinates of the cells in a reduced dimensionality space.
#' @slot clusterLabels character. A character vector of length \code{n} denoting
#'   each cell's cluster label.
#' @slot lineages list. A list with each element a character vector of cluster 
#'   names representing a lineage as an ordered set of clusters.
#' @slot connectivity matrix. A binary matrix describing the connectivity 
#'   between clusters induced by the minimum spanning tree.
#' @slot lineageControl list. Additional parameters specifying how the minimum 
#'   spanning tree on clusters was constructed. \itemize{ 
#'   \item{\code{start.clus}}{character. The label of the root cluster.} 
#'   \item{\code{end.clus}}{character. Vector of cluster labels indicating the 
#'   terminal clusters.} \item{\code{start.given}}{logical. A logical value 
#'   indicating whether the initial state was pre-specified.} 
#'   \item{\code{end.given}}{logical. A vector of logical values indicating 
#'   whether each terminal state was pre-specified} \item{\code{dist}}{matrix. A
#'   numeric matrix of pairwise cluster distances.} }
#' @slot curves list. A list of \code{principal.curve} objects produced by 
#'   \code{\link{getCurves}}.
#' @slot curveControl list. Additional parameters specifying how the 
#'   simultaneous principal curves were constructed. \itemize{ 
#'   \item{\code{shrink}}{logical or numeric between 0 and 1. Determines whether
#'   and how much to shrink branching lineages toward their shared average 
#'   curve.} \item{\code{extend}}{character. Specifies the method for handling 
#'   root and leaf clusters of lineages when constructing the initial, 
#'   piece-wise linear curve. Accepted values are 'y' (default), 'n', and 'pc1'.
#'   See \code{\link{getCurves}} for details.} \item{\code{reweight}}{logical. 
#'   Indicates whether to reweight cells shared by multiple lineages during 
#'   curve-fitting. If \code{TRUE}, cells shared between lineages will have 
#'   lineage-specific weights determined by the ratio: (distance to nearest 
#'   curve) / (distance to specific curve).} \item{\code{drop.multi}}{logical. 
#'   Indicates whether to drop shared cells from lineages which do not fit them 
#'   well. If \code{TRUE}, shared cells with a distance to one lineage above the
#'   90th percentile and another lineage below the 50th percentile will be 
#'   dropped from the farther lineage.} \item{\code{shrink.method}}{character. 
#'   Denotes how to determine the amount of shrinkage for a branching lineage. 
#'   Accepted values are the same as for \code{kernel} in  the \code{density} 
#'   function (default is \code{"cosine"}), as well as \code{"tricube"} and 
#'   \code{"density"}. See \code{\link{getCurves}} for details.} \item{Other 
#'   parameters specified by \code{\link{principal.curve}}}. }
#'   
#' @return The accessor functions \code{reducedDim}, \code{clusterLabels}, 
#'   \code{lineages}, \code{connectivity}, \code{lineageControl}, \code{curves},
#'   and \code{curveControl} return the corresponding elements of a 
#'   \code{SlingshotDataSet}. The functions \code{pseudotime} and 
#'   \code{curveWeights} extract useful output elements of a 
#'   \code{SlingshotDataSet}, provided that curves have already been fit with 
#'   either \code{slingshot} or \code{getCurves}.
#'   
#' @import princurve
#' @import methods
#' @export
#' 
setClass(
    Class = "SlingshotDataSet",
    slots = list(
        reducedDim = "matrix",
        clusterLabels = "character",
        lineages = "list",
        connectivity = "matrix",
        lineageControl = "list",
        curves = "list",
        pseudotime = "matrix",
        curveWeights = "matrix",
        curveControl = "list"
    )
)

setValidity("SlingshotDataSet", function(object) {
    X <- reducedDim(object)
    n <- nrow(X)
    p <- ncol(X)
    if(!is.numeric(X)) {
        return("Reduced dimensional coordinates must be numeric.")
    }
    if(nrow(X)==0){
        return('reducedDim has zero rows.')
    }
    if(ncol(X)==0){
        return('reducedDim has zero columns.')
    }
    if(length(clusterLabels(object)) != n){
        return('nrow(reducedDim) must equal length(clusterLabels).')
    }
    # something requires row and column names. Princurve?
    if(is.null(rownames(reducedDim(object)))){
        rownames(reducedDim(object)) <- paste('Cell',
                                              seq_len(nrow(reducedDim(object))),
                                              sep='-')
    }
    if(is.null(colnames(reducedDim(object)))){
        colnames(reducedDim(object)) <- paste('Dim',
                                              seq_len(ncol(reducedDim(object))),
                                              sep='-')
    }
    
    # if lineages present
    if(length(lineages(object)) > 0){
        L <- length(lineages(object))
        clus.names <- unique(clusterLabels(object))
        K <- length(clus.names)
        if(any(sapply(lineages(object),class) != 'character')){
            return("lineages must be a list of character vectors.")
        }
        if(!all(sapply(lineages(object), 
                       function(lin){all(lin %in% clus.names)}))){
            return("lineages must be a list of character vectors composed of 
                   cluster names.")
        }
        if(!is.numeric(connectivity(object))) {
            return("Connectivity matrix must be numeric or logical.")
        }
        if(any(dim(connectivity(object)) != K)){
            return("Connectivity matrix must be square with number of dimensions
                    equal to number of clusters")
        }
        if(! is.null(lineageControl(object)$start.clus)){
            if(!all(lineageControl(object)$start.clus %in% clus.names)){
                return("Specified starting cluster not found in cluster labels")
            }
        }
        if(! is.null(lineageControl(object)$end.clus)){
            if(!all(lineageControl(object)$end.clus %in% clus.names)){
                return("Specified terminal cluster(s) not found in cluster 
                       labels")
            }
        }
        if(! is.null(lineageControl(object)$dist.fun)){
            if(!is.function(lineageControl(object)$dist.fun)){
                return("Pairwise cluster distance function must be a function.")
            }
        }
        if(! is.null(lineageControl(object)$omega)){
            if(lineageControl(object)$omega < 0 | 
               (lineageControl(object)$omega > 1 & 
                lineageControl(object)$omega != Inf)){
                return("Omega must be numeric element of [0,1] or Inf.")
            }
        }
    }
    
    # if curves present
    if(length(curves(object)) > 0){
        if(length(lineages(object)) > 0){
            L <- length(lineages(object))
            if(length(curves(object)) != L){
                return("Number of curves does not match number of lineages")
            }
        }
        L <- length(curves(object))
        if(any(sapply(curves(object),class) != 'principal.curve')){
            return("curves must be a list of principal.curve objects.")
        }
        if(dim(pseudotime(object))[1] > 0){
            if(any(dim(pseudotime(object)) != c(n,L))){
                return("Dimensions for pseudotime matrix are incorrect. Should 
                       be n (number of cells) by L (number of lineages).")
            }
        }
        if(dim(curveWeights(object))[1] > 0){
            if(any(dim(curveWeights(object)) != c(n,L))){
                return("Dimensions for curveWeights matrix are incorrect. 
                       Should be n (number of cells) by L 
                       (number of lineages).")
            }
        }
        if(!is.null(curveControl(object)$shrink)){
            if(curveControl(object)$shrink < 0 | 
               curveControl(object)$shrink > 1){
                stop("shrink argument must be logical or numeric between 
                     0 and 1.")
            }
        }
        if(!is.null(curveControl(object)$extend)){
            if(! curveControl(object)$extend %in% c('y','n','pc1')){
                stop("extend argument must be one of 'y', 'n', or 'pc1'.")
            }
        }
        if(!is.null(curveControl(object)$reweight)){
            if(!is.logical(curveControl(object)$reweight)){
                stop("reweight argument must be logical.")
            }
        }
        if(!is.null(curveControl(object)$drop.multi)){
            if(!is.logical(curveControl(object)$drop.multi)){
                stop("drop.multi argument must be logical.")
            }
        }
    }
    return(TRUE)
})

