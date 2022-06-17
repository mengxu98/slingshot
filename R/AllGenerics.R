#' @title Initialize an object of class \code{SlingshotDataSet}
#' @name newSlingshotDataSet
#' @docType methods
#'
#' @description Constructs a \code{SlingshotDataSet} object. Additional helper
#'   methods for manipulating \code{SlingshotDataSet} objects are  also
#'   described below. We now recommend using
#'   \code{\link[TrajectoryUtils]{PseudotimeOrdering}} objects, instead.
#'
#' @param reducedDim matrix. An \code{n} by \code{p} numeric matrix or data
#'   frame giving the coordinates of the cells in a reduced dimensionality
#'   space.
#' @param clusterLabels character. A character vector of length \code{n}
#'   denoting each cell's cluster label.
#' @param ... additional components of a \code{SlingshotDataSet} to specify.
#'   This may include any of the following:
#' @param lineages list. A list with each element a character vector of cluster
#'   names representing a lineage as an ordered set of clusters.
#' @param adjacency matrix. A binary matrix describing the connectivity
#'   between clusters induced by the minimum spanning tree.
#' @param slingParams list. Additional parameters used by Slingshot. These may
#'   specify how the minimum spanning tree on clusters was constructed:
#'   \itemize{
#'   \item{\code{start.clus}}{ character. The label of the root cluster.}
#'   \item{\code{end.clus}}{ character. Vector of cluster labels indicating the
#'   terminal clusters.}
#'   \item{\code{start.given}}{ logical. A logical value
#'   indicating whether the initial state was pre-specified.}
#'   \item{\code{end.given}}{ logical. A vector of logical values indicating
#'   whether each terminal state was pre-specified}
#'   \item{\code{dist}}{ matrix. A
#'   numeric matrix of pairwise cluster distances.} }
#'   They may also specify how simultaneous principal curves were constructed:
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
#'   Indicates whether to reassign cells to
#'   lineages at each iteration. If \code{TRUE}, cells will be added to a
#'   lineage when their projection distance to the curve is less than the median
#'   distance for all cells currently assigned to the lineage. Additionally,
#'   shared cells will be removed from a lineage if their projection distance to
#'   the curve is above the 90th percentile and their weight along the curve is
#'   less than \code{0.1}.}
#'   \item{\code{shrink.method}}{ character.
#'   Denotes how to determine the amount of shrinkage for a branching lineage.
#'   Accepted values are the same as for \code{kernel} in  the \code{density}
#'   function (default is \code{"cosine"}), as well as \code{"tricube"} and
#'   \code{"density"}. See \code{\link{getCurves}} for details.}
#'   \item{Other parameters specified by
#'   \code{\link[princurve]{principal_curve}}}. }
#' @param curves list. A list of \code{\link[princurve]{principal_curve}}
#'   objects produced by \code{\link{getCurves}}.
#'
#' @return A \code{SlingshotDataSet} object with all specified values.
#'
#' @seealso \code{\link[TrajectoryUtils]{PseudotimeOrdering}}
#'
#' @examples
#' rd <- matrix(data=rnorm(100), ncol=2)
#' cl <- sample(letters[seq_len(5)], 50, replace = TRUE)
#' sds <- newSlingshotDataSet(rd, cl)
#'
#' @import princurve
#' @import methods
#' @export
setGeneric(
    name = "newSlingshotDataSet",
    signature = c('reducedDim','clusterLabels'),
    def = function(reducedDim,  clusterLabels, ...) {
        standardGeneric("newSlingshotDataSet")
    }
)

#' @title Extract Slingshot output
#' @name SlingshotDataSet
#' @description This is a convenience function to extract a
#'   \code{SlingshotDataSet} from an object containing \code{\link{slingshot}}
#'   output. However, we now recommend using a
#'   \code{\link[TrajectoryUtils]{PseudotimeOrdering}} object, in most cases.
#'   The \code{SlingshotDataSet} is, however, still used for plotting purposes.
#' @param data an object containing \code{slingshot} output.
#' @param ... additional arguments to pass to object-specific methods.
#' @return A \code{SlingshotDataSet} object containing the output of
#' \code{slingshot}.
#' 
#' @seealso \code{\link[TrajectoryUtils]{PseudotimeOrdering}},
#'   \code{\link{as.SlingshotDataSet}}
#'
#' @examples 
#' data("slingshotExample")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' library(SingleCellExperiment)
#' u <- matrix(rpois(140*50, 5), nrow = 50)
#' sce <- SingleCellExperiment(assays = list(counts = u), 
#'                             reducedDims = SimpleList(PCA = rd),
#'                             colData = data.frame(clus = cl))
#' sce <- slingshot(sce, clusterLabels = 'clus', reducedDim = 'PCA')
#' SlingshotDataSet(sce)
#' 
#' @export
setGeneric(
    name = "SlingshotDataSet",
    signature = c('data'),
    def = function(data, ...) {
        standardGeneric("SlingshotDataSet")
    }
)

#' @title Infer Lineage Structure from Clustered Samples
#' @name getLineages
#' @param ... Additional arguments to specify how lineages are constructed from
#'   clusters.
#' @export
setGeneric(
    name = "getLineages",
    signature = c('data','clusterLabels'),
    def = function(data,
                   clusterLabels, ...) {
        standardGeneric("getLineages")
    }
)

#' @title Construct Simultaneous Principal Curves
#' @name getCurves
#' @export
setGeneric(
    name = "getCurves",
    signature = 'data',
    def = function(data, ...) {
        standardGeneric("getCurves")
    }
)

#' @title Perform trajectory inference with Slingshot
#' @description Perform trajectory inference with Slingshot
#' @name slingshot
#' @export
setGeneric(
    name = "slingshot",
    signature = c('data', 'clusterLabels'),
    def = function(data,
                   clusterLabels, ...) {
        standardGeneric("slingshot")
    }
)

#' @title Extract the Slingshot lineages
#' @name slingLineages
#'
#' @description Extract lineages (represented by ordered sets of clusters)
#'   identified by \code{\link{slingshot}}.
#'
#' @param x an object containing \code{\link{slingshot}} output.
#' @return A list of lineages, represented by ordered sets of clusters.
#' @examples
#' data("slingshotExample")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' pto <- slingshot(rd, cl, start.clus = '1')
#' slingLineages(pto)
#' @export
setGeneric(name = "slingLineages",
           signature = "x",
           def = function(x) standardGeneric("slingLineages"))

#' @title Extract dimensionality reduction used by Slingshot
#' @name slingReducedDim
#'
#' @description Extract the dimensionality reduction used by 
#' \code{\link{slingshot}}.
#'
#' @param x an object containing \code{\link{slingshot}} output.
#' @return A matrix of coordinates.
#' @examples
#' data("slingshotExample")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' pto <- slingshot(rd, cl, start.clus = '1')
#' slingReducedDim(pto)
#' @export
setGeneric(name = "slingReducedDim",
           signature = "x",
           def = function(x) standardGeneric("slingReducedDim"))

#' @title Extract cluster labels used by Slingshot
#' @name slingClusterLabels
#'
#' @description Extract the cluster labels used by \code{\link{slingshot}}.
#'
#' @param x an object containing \code{\link{slingshot}} output.
#' @return Typically returns a matrix of cluster assignment weights
#'   (\code{#cells} by \code{#clusters}). Rarely, a vector of cluster labels.
#' @examples
#' data("slingshotExample")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' pto <- slingshot(rd, cl, start.clus = '1')
#' slingClusterLabels(pto)
#' @export
setGeneric(name = "slingClusterLabels",
           signature = "x",
           def = function(x) standardGeneric("slingClusterLabels"))

#' @title Extract Slingshot minimum spanning tree
#' @name slingMST
#' @description Extract the minimum spanning tree from an object containing
#'   \code{\link{slingshot}} output.
#'
#' @param x an object containing \code{\link{slingshot}} output.
#' @param ... additional parameters to be passed to object-specific methods.
#' @return In most cases, output is an \code{\link[igraph]{igraph}} object
#'   representing the MST. If \code{x} is a \code{SlingshotDataSet}, then output
#'   is an adjacency matrix representing the MST.
#' @examples
#' data("slingshotExample")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' pto <- slingshot(rd, cl, start.clus = '1')
#' slingMST(pto)
#' @export
setGeneric(name = "slingMST",
           signature = "x",
           def = function(x, ...) standardGeneric("slingMST"))

#' @title Methods for parameters used by Slingshot
#' @name slingParams
#' @description Extracts additional control parameters used by Slingshot in
#' lineage inference and fitting simultaneous principal curves.
#'
#' @param x an object containing \code{\link{slingshot}} output.
#' @return The list of additional parameters used by Slingshot. These include
#'   parameters related to the cluster-based minimum spanning tree:
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
#' @examples
#' data("slingshotExample")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' pto <- slingshot(rd, cl, start.clus = '1')
#' slingParams(pto)
#' @export
setGeneric(name = "slingParams",
           signature = "x",
           def = function(x) standardGeneric("slingParams"))

#' @title Extract simultaneous principal curves
#' @name slingCurves
#' @description Extract the simultaneous principal curves from an object
#'   containing \code{\link{slingshot}} output.
#'
#' @param x an object containing \code{\link{slingshot}} output.
#' @param ... additional parameters to be passed to object-specific methods.
#' @return A list of smooth lineage curves, each of which is a
#'   \code{\link[princurve]{principal_curve}} object.
#' @examples
#' data("slingshotExample")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' pto <- slingshot(rd, cl, start.clus = '1')
#' slingCurves(pto)
#' @export
setGeneric(name = "slingCurves",
           signature = "x",
           def = function(x, ...) standardGeneric("slingCurves"))


#' @title Get Slingshot pseudotime values
#' @name slingPseudotime
#'
#' @description Extract the matrix of pseudotime values or cells' weights along
#'   each lineage.
#'
#' @param x an object containing \code{\link{slingshot}} output.
#' @param ... additional parameters to be passed to object-specific methods.
#' @return \code{slingPseudotime}: an \code{n} by \code{L} matrix representing
#'   each cell's pseudotime along each lineage.
#' @examples
#' data("slingshotExample")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' pto <- slingshot(rd, cl, start.clus = '1')
#' slingPseudotime(pto)
#' @export
setGeneric(name = "slingPseudotime",
           signature = "x",
           def = function(x, ...) standardGeneric("slingPseudotime"))

#' @rdname slingPseudotime
#' @return \code{slingCurveWeights}: an \code{n} by \code{L} matrix of cell
#'   weights along each lineage.
#' @examples
#' slingCurveWeights(pto)
#' @export
setGeneric(name = "slingCurveWeights",
           signature = "x",
           def = function(x, ...) standardGeneric("slingCurveWeights"))

#' @rdname slingPseudotime
#' @return \code{slingAvgPseudotime}: a length \code{n} vector of average cell
#'   pseudotimes, where the average is a weighted average across lineages,
#'   weighted by the assignment weights.
#' @examples
#' slingAvgPseudotime(pto)
#' @export
setGeneric(name = "slingAvgPseudotime",
           signature = "x",
           def = function(x, ...) standardGeneric("slingAvgPseudotime"))

#' @title Embed trajectory in new space
#' @rdname embedCurves
#' @export
setGeneric(name = "embedCurves",
           signature = c("x", "newDimRed"),
           def = function(x, newDimRed, ...) standardGeneric("embedCurves"))

#' @title Get slingshot branch labels
#' @rdname slingBranchID
#' @param ... additional arguments passed to object-specific methods.
#' @export
setGeneric(name = "slingBranchID",
           signature = c("x"),
           def = function(x, ...) standardGeneric("slingBranchID"))

#' @title Construct graph of slingshot branch labels
#' @rdname slingBranchGraph
#' @param ... additional arguments passed to object-specific methods.
#' @export
setGeneric(name = "slingBranchGraph",
           signature = c("x"),
           def = function(x, ...) standardGeneric("slingBranchGraph"))

#' @title Conversion to SlingshotDataSet
#' @rdname as.SlingshotDataSet
#' @param ... additional arguments passed to object-specific methods.
#' @export
setGeneric(name = "as.SlingshotDataSet",
           signature = c("x"),
           def = function(x, ...) standardGeneric("as.SlingshotDataSet"))

#' @title Conversion to PseudotimeOrdering
#' @rdname as.PseudotimeOrdering
#' @param ... additional arguments passed to object-specific methods.
#' @export
setGeneric(name = "as.PseudotimeOrdering",
           signature = c("x"),
           def = function(x, ...) standardGeneric("as.PseudotimeOrdering"))




