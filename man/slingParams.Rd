% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/AllHelperFunctions.R
\name{slingParams}
\alias{slingParams}
\alias{slingParams,PseudotimeOrdering-method}
\alias{slingParams,SingleCellExperiment-method}
\alias{slingParams,SlingshotDataSet-method}
\title{Methods for parameters used by Slingshot}
\usage{
slingParams(x)

\S4method{slingParams}{PseudotimeOrdering}(x)

\S4method{slingParams}{SingleCellExperiment}(x)

\S4method{slingParams}{SlingshotDataSet}(x)
}
\arguments{
\item{x}{an object containing \code{\link{slingshot}} output.}
}
\value{
The list of additional parameters used by Slingshot. These include
  parameters related to the cluster-based minimum spanning tree:
  \itemize{
  \item{\code{start.clus}}{ character. The label of the root cluster, or a
  vector of cluster labels giving the root clusters of each disjoint
  component of the graph.}
  \item{\code{end.clus}}{ character. Vector of cluster labels indicating
  terminal clusters.}
  \item{\code{start.given}}{ logical. A logical value
  indicating whether the initial state was pre-specified.}
  \item{\code{end.given}}{ logical. A vector of logical values indicating
  whether each terminal state was pre-specified}
  \item{\code{omega}}{ numeric or logical. Granularity parameter determining
  the maximum edge length for building the MST. See
  \code{\link{getLineages}}.}
  \item{\code{omega_scale}}{ numeric. Scaling factor used for setting maximum
  edge length when \code{omega = TRUE}. See \code{\link{getLineages}}.} }
  They may also specify how simultaneous principal curves were constructed
  (for a complete listing, see \code{\link{getCurves}}:
  \itemize{
  \item{\code{shrink}}{ logical or numeric between 0 and 1. Determines
  whether and how much to shrink branching lineages toward their shared
  average curve.}
  \item{\code{extend}}{ character. Specifies the method for handling
  root and leaf clusters of lineages when constructing the initial,
  piece-wise linear curve. Accepted values are 'y' (default), 'n', and 'pc1'.
  See \code{\link{getCurves}} for details.}
  \item{\code{reweight}}{ logical.
  Indicates whether to allow cells shared
  between lineages to be reweighted during curve-fitting. If \code{TRUE},
  cells shared between lineages will be iteratively reweighted based on the
  quantiles of their projection distances to each curve.}
  \item{\code{reassign}}{ logical.
  Indicates whether to reassign cells to lineages at each
  iteration. If \code{TRUE}, cells will be added to a lineage when their
  projection distance to the curve is less than the median distance for all
  cells currently assigned to the lineage. Additionally, shared cells will be
  removed from a lineage if their projection distance to the curve is above
  the 90th percentile and their weight along the curve is less than
  \code{0.1}.}
  \item{\code{shrink.method}}{ character.
  Denotes how to determine the amount of shrinkage for a branching lineage.
  Accepted values are the same as for \code{kernel} in  the \code{density}
  function (default is \code{"cosine"}), as well as \code{"tricube"} and
  \code{"density"}. See \code{\link{getCurves}} for details.}
  \item{approx_points}{ numeric. Number of points to use in estimating
  curves. See \code{\link{getCurves}} for details.} \item{allow.breaks}{
  logical. Whether to allow curves that diverge very early on in a trajectory
  to have different starting points.}
  \item{Other parameters specified by
  \code{\link[princurve]{principal_curve}}}. }
}
\description{
Extracts additional control parameters used by Slingshot in
lineage inference and fitting simultaneous principal curves.
}
\examples{
data("slingshotExample")
rd <- slingshotExample$rd
cl <- slingshotExample$cl
pto <- slingshot(rd, cl, start.clus = "1")
slingParams(pto)
}
