#' @rdname as.SlingshotDataSet
#' @description This function converts objects that contain \code{slingshot}
#'   results into a \code{SlingshotDataSet}.
#' @param x an object containing \code{slingshot} output.
#' @param ... additional arguments to pass to object-specific methods.
#' @return A \code{SlingshotDataSet} object containing the \code{slingshot}
#'   results from the original object, \code{x}.
#'
#' @seealso \code{\link[TrajectoryUtils]{PseudotimeOrdering}}
#'
#' @examples 
#' data("slingshotExample")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' pto <- slingshot(rd, cl, start.clus = '1')
#' as.SlingshotDataSet(pto)
#' 
#' @export
setMethod(
    f = "as.SlingshotDataSet",
    signature = "PseudotimeOrdering",
    definition = function(x){
        pto <- x
        if(length(slingCurves(pto)) > 0){
            crvs <- slingCurves(pto)
        }else{
            crvs <- list()
        }
        sds <- newSlingshotDataSet(reducedDim = slingReducedDim(pto),
                                   clusterLabels = slingClusterLabels(pto),
                                   lineages = slingLineages(pto),
                                   adjacency = as.matrix(
                                       igraph::as_adjacency_matrix(
                                           slingMST(pto))),
                                   curves = crvs,
                                   slingParams = slingParams(pto))
        return(sds)
    })

#' @rdname as.SlingshotDataSet
#' @export
setMethod(
    f = "as.SlingshotDataSet",
    signature = "SingleCellExperiment",
    definition = function(x){
        if("slingshot" %in% names(colData(x))){
            return(as.SlingshotDataSet(colData(x)$slingshot))
        }
        if("slingshot" %in% names(x@int_metadata)){
            return(x@int_metadata$slingshot)            
        }
        stop('No slingshot results found.')
    }
)

#' @rdname as.SlingshotDataSet
#' @export
setMethod(
    f = "as.SlingshotDataSet",
    signature = "SlingshotDataSet",
    definition = function(x){
        return(x)
    }
)


#' @rdname as.PseudotimeOrdering
#' @description This function converts objects that contain \code{slingshot}
#'   results into a \code{\link[TrajectoryUtils]{PseudotimeOrdering}}.
#' @param x an object containing \code{slingshot} output.
#' @param ... additional arguments to pass to object-specific methods.
#' @return A \code{PseudotimeOrdering} object containing the \code{slingshot}
#'   results from the original object, \code{x}.
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
#' as.PseudotimeOrdering(sce)
#' 
#' @import TrajectoryUtils
#' @export
setMethod(
    f = "as.PseudotimeOrdering",
    signature = "SlingshotDataSet",
    definition = function(x){
        sds <- x
        if(length(slingCurves(sds)) > 0){
            ps <- list(pseudotime = slingPseudotime(sds),
                       weights = slingCurveWeights(sds))
            meta <- list(lineages = slingLineages(sds),
                         mst = igraph::graph_from_adjacency_matrix(
                             slingMST(sds), mode = "undirected"),
                         curves = slingCurves(sds),
                         slingParams = slingParams(sds))
        }else if(length(slingLineages(sds)) > 0){
            cl <- slingClusterLabels(sds)
            wts <- vapply(slingLineages(sds), function(lin){
                rowSums(cl[,colnames(cl) %in% lin])
            }, rep(0,nrow(reducedDim(sds))))
            pst <- wts
            pst[,] <- NA
            ps <- list(pseudotime = pst,
                       weights = wts)
            meta <- list(lineages = slingLineages(sds),
                         mst = igraph::graph_from_adjacency_matrix(
                                        slingMST(sds), 
                                        mode = "undirected"),
                         slingParams = slingParams(sds))
        }else{
            stop("'as.PseudotimeOrdering' coversion failed: ",
                 "number of lineages could not be determined.")
        }
        
        # include cluster center coordinates
        centers <- rowmean(reducedDim(sds), slingClusterLabels(sds))
        coord.list <- vector("list", nrow(centers))
        names(coord.list) <- rownames(centers)
        for (r in rownames(centers)) {
            coord.list[[r]] <- centers[r,]
        }
        igraph::V(meta$mst)$coordinates <- 
            coord.list[names(igraph::V(meta$mst))]
        
        pto <- PseudotimeOrdering(pathStats = ps, metadata = meta)
        cellData(pto)$reducedDim <- reducedDim(sds)
        cellData(pto)$clusterLabels <- slingClusterLabels(sds)
        return(pto)
    })

#' @rdname as.PseudotimeOrdering
#' @export
setMethod(
    f = "as.PseudotimeOrdering",
    signature = "SingleCellExperiment",
    definition = function(x){
        if("slingshot" %in% names(colData(x))){
            return(colData(x)$slingshot)
        }
        if("slingshot" %in% names(x@int_metadata)){
            return(as.PseudotimeOrdering(x@int_metadata$slingshot))            
        }
        stop('No slingshot results found.')
    })

#' @rdname as.PseudotimeOrdering
#' @export
setMethod(
    f = "as.PseudotimeOrdering",
    signature = "PseudotimeOrdering",
    definition = function(x){ x })