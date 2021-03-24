#' @title Class \code{SlingshotDataSet}
#' @aliases SlingshotDataSet-class
#'   
#' @description This was the original class for storing \code{slingshot}
#'   results, but we now generally reommend using the
#'   \code{\link{PseudotimeOrdering}} class, instead. Most \code{slingshot}
#'   functions will still work with \code{SlingshotDataSet} objects, but will
#'   return \code{PseudotimeOrdering} objects, by default. To update old
#'   \code{SlingshotDataSet} objects, we have provided the
#'   \code{\link{as.PseudotimeOrdering}} conversion function. The only functions
#'   that require \code{SlingshotDataSet} objects are the plotting functions.
#'   
#' @description The \code{SlingshotDataSet} class holds data relevant for 
#'   performing lineage inference with the \code{slingshot} package, primarily a
#'   reduced dimensional representation of the data and a set of cluster labels.
#'   
#' @slot reducedDim matrix. An \code{n} by \code{p} numeric matrix or data frame
#'   giving the coordinates of the cells in a reduced dimensionality space.
#' @slot clusterLabels matrix or character. An \code{n} by \code{K} matrix of 
#'   weights indicating each cell's cluster assignment or a character vector of
#'   cluster assignments, which will be converted into a binary matrix.
#' @slot lineages list. A list with each element a character vector of cluster 
#'   names representing a lineage as an ordered set of clusters.
#' @slot adjacency matrix. A binary matrix describing the adjacency 
#'   between clusters induced by the minimum spanning tree.
#' @slot curves list. A list of \code{\link[princurve]{principal_curve}} objects
#'   produced by \code{\link{getCurves}}.
#' @slot slingParams list. Additional parameters used by Slingshot. These may 
#'   specify how the minimum spanning tree on clusters was constructed: 
#'   \itemize{ 
#'   \item{\code{start.clus}}{ character. The label of the root cluster, or a
#'   vector of cluster labels giving the root clusters of each disjoint
#'   component of the graph.}
#'   \item{\code{end.clus}}{ character. Vector of cluster labels indicating 
#'   terminal clusters.}
#'   \item{\code{start.given}}{ logical. A logical value 
#'   indicating whether the initial state was pre-specified.} 
#'   \item{\code{end.given}}{ logical. A vector of logical values indicating 
#'   whether each terminal state was pre-specified} 
#'   \item{\code{omega}}{ numeric or logical. Granularity parameter determining
#'   the maximum edge length for building the MST. See
#'   \code{\link{getLineages}}.}
#'   \item{\code{omega_scale}}{ numeric. Scaling factor used for setting maximum
#'   edge length when \code{omega = TRUE}. See \code{\link{getLineages}}.} }
#'   They may also specify how simultaneous principal curves were constructed
#'   (for a complete listing, see \code{\link{getCurves}}:
#'   \itemize{ 
#'   \item{\code{shrink}}{ logical or numeric between 0 and 1. Determines
#'   whether and how much to shrink branching lineages toward their shared
#'   average curve.}
#'   \item{\code{extend}}{ character. Specifies the method for handling 
#'   root and leaf clusters of lineages when constructing the initial, 
#'   piece-wise linear curve. Accepted values are 'y' (default), 'n', and 'pc1'.
#'   See \code{\link{getCurves}} for details.} 
#'   \item{\code{reweight}}{ logical. 
#'   Indicates whether to allow cells shared
#'   between lineages to be reweighted during curve-fitting. If \code{TRUE},
#'   cells shared between lineages will be iteratively reweighted based on the
#'   quantiles of their projection distances to each curve.} 
#'   \item{\code{reassign}}{ logical. 
#'   Indicates whether to reassign cells to lineages at each
#'   iteration. If \code{TRUE}, cells will be added to a lineage when their
#'   projection distance to the curve is less than the median distance for all
#'   cells currently assigned to the lineage. Additionally, shared cells will be
#'   removed from a lineage if their projection distance to the curve is above
#'   the 90th percentile and their weight along the curve is less than
#'   \code{0.1}.} 
#'   \item{\code{shrink.method}}{ character. 
#'   Denotes how to determine the amount of shrinkage for a branching lineage. 
#'   Accepted values are the same as for \code{kernel} in  the \code{density} 
#'   function (default is \code{"cosine"}), as well as \code{"tricube"} and 
#'   \code{"density"}. See \code{\link{getCurves}} for details.}
#'   \item{approx_points}{ numeric. Number of points to use in estimating
#'   curves. See \code{\link{getCurves}} for details.} \item{allow.breaks}{
#'   logical. Whether to allow curves that diverge very early on in a trajectory
#'   to have different starting points.}
#'   \item{Other parameters specified by 
#'   \code{\link[princurve]{principal_curve}}}. }
#'   
#' @return The accessor functions \code{reducedDim}, \code{clusterLabels}, 
#'   \code{lineages}, \code{adjacency}, \code{curves},
#'   and \code{slingParams} return the corresponding elements of a 
#'   \code{SlingshotDataSet}. The functions \code{slingPseudotime} and 
#'   \code{slingCurveWeights} extract useful output elements of a 
#'   \code{SlingshotDataSet}, provided that curves have already been fit with 
#'   either \code{slingshot} or \code{getCurves}.
#'   
#' @seealso \code{\link{PseudotimeOrdering}}
#'   
#' @import princurve
#' @import methods
#' @export
setClass(
    Class = "SlingshotDataSet",
    slots = list(
        reducedDim = "matrix",
        clusterLabels = "matrix",
        lineages = "list",
        adjacency = "matrix",
        curves = "list",
        slingParams = "list"
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
    if(nrow(slingClusterLabels(object)) != n){
        return(paste('Reduced dimensional coordinates and cluster labels', 
                     'contain different numbers of cells.'))
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
    if(length(slingLineages(object)) > 0){
        L <- length(slingLineages(object))
        clus.names <- colnames(slingClusterLabels(object))
        K <- length(clus.names)
        if(any(vapply(slingLineages(object),class,'') != 'character')){
            return("lineages must be a list of character vectors.")
        }
        if(!all(vapply(slingLineages(object), 
                       function(lin){all(lin %in% clus.names)}, TRUE))){
            return(paste0("lineages must be a list of character vectors ",
                          "composed of cluster names."))
        }
        if(!is.numeric(slingMST(object))) {
            return("adjacency matrix must be numeric or logical.")
        }
        if(any(dim(slingMST(object)) != K)){
            return(paste("adjacency matrix must be square with number of",
                         "dimensions equal to number of clusters"))
        }
        if(! is.null(slingParams(object)$start.clus)){
            if(!all(slingParams(object)$start.clus %in% clus.names)){
                return("Specified starting cluster not found in cluster labels")
            }
        }
        if(! is.null(slingParams(object)$end.clus)){
            if(!all(slingParams(object)$end.clus %in% clus.names)){
                return(paste0("Specified terminal cluster(s) not found in ",
                              "cluster labels"))
            }
        }
        if(! is.null(slingParams(object)$dist.fun)){
            if(!is.function(slingParams(object)$dist.fun)){
                return("Pairwise cluster distance function must be a function.")
            }
        }
        if(! is.null(slingParams(object)$omega)){
            if(length(slingParams(object)$omega) > 1){
                return('omega must be NULL or length 1')
            }
            if(is.na(as.numeric(slingParams(object)$omega))){
                stop('omega must be numeric, logical, or NULL')
            }
            if(slingParams(object)$omega < 0){
                stop('omega must be non-negative')
            }
        }
    }
    
    # if curves present
    if(length(slingCurves(object)) > 0){
        if(length(slingLineages(object)) > 0){
            L <- length(slingLineages(object))
            if(length(slingCurves(object)) != L){
                return("Number of curves does not match number of lineages")
            }
        }
        L <- length(slingCurves(object))
        if(any(vapply(slingCurves(object),class,'') != 'principal_curve')){
            return("curves must be a list of principal_curve objects.")
        }
        if(!is.null(slingParams(object)$shrink)){
            if(slingParams(object)$shrink < 0 | 
               slingParams(object)$shrink > 1){
                stop("shrink argument must be logical or numeric between ", 
                     "0 and 1.")
            }
        }
        if(!is.null(slingParams(object)$extend)){
            if(! slingParams(object)$extend %in% c('y','n','pc1')){
                stop("extend argument must be one of 'y', 'n', or 'pc1'.")
            }
        }
        if(!is.null(slingParams(object)$reweight)){
            if(!is.logical(slingParams(object)$reweight)){
                stop("reweight argument must be logical.")
            }
        }
        if(!is.null(slingParams(object)$reassign)){
            if(!is.logical(slingParams(object)$reassign)){
                stop("reassign argument must be logical.")
            }
        }
    }
    return(TRUE)
})
