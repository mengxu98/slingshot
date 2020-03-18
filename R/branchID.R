
#' @rdname slingBranchID
#' @title Get slingshot branch labels
#'
#' @description Extracts a categorical variable describing which lineage (or
#'   combination of lineages) each cell is assigned to.
#' @param x an object containing slingshot output.
#' @param thresh the minimum weight of assignment required to assign a cell to a
#'   lineage (default = 1/L)
#' @export
slingBranchID <- function(x, thresh = NULL){
    L <- length(slingLineages(x))
    if(is.null(thresh)){
        thresh <- 1/L
    }else{
        if(thresh < 0 | thresh > 1){
            stop("'thresh' value must be between 0 and 1.")
        }
    }
    return(factor(apply(slingCurveWeights(x) >= thresh, 1, function(bin){
        paste(which(bin), collapse = '+')
    })))
}

.under <- function(n, nodes){
    which.lin <- strsplit(nodes, split='[+]')
    nlins <- sapply(which.lin, length)
    out <- nodes[sapply(which.lin, function(wl){
        all(wl %in% unlist(strsplit(n, split='[+]')))
    })]
    return(out[out != n])
}

#' @rdname slingBranchGraph
#' @title Construct graph of slingshot branch labels
#'
#' @description Builds a graph describing the relationships between the
#'   different branch assignments
#' @param x an object containing slingshot output.
#' @param thresh the minimum weight of assignment required to assign a cell to a
#'   lineage (default = 1/L)
#' @importFrom igraph graph_from_literal graph_from_edgelist vertex_attr vertex_attr<-
#' @export
slingBranchGraph <- function(x, thresh = NULL){
    brID <- slingBranchID(x, thresh = thresh)
    nodes <- as.character(levels(brID))
    which.lin <- strsplit(nodes, split='[+]')
    nlins <- vapply(which.lin, length, 0)
    maxL <- max(nlins)
    if(maxL == 1){ # only one lineage
        g <- graph_from_literal(1)
        igraph::vertex_attr(g, 'size') <- length(brID)
        return(g)
    }
    
    el <- NULL
    # for each node n at level l
    for(l in seq(2,maxL)){
        for(n in nodes[nlins==l]){
            # find all descendants of n
            desc <- .under(n, nodes)
            for(d in desc){
                if(l - nlins[which(nodes==d)] >= 2){
                    # check for intermediates
                    granddesc <- unique(unlist(lapply(desc, .under, nodes)))
                    if(! d %in% granddesc){
                        # add edge
                        el <- rbind(el, c(n, d))
                    }
                }else{
                    # add edge
                    el <- rbind(el, c(n, d))
                }
            }
        }
    }
    g <- igraph::graph_from_edgelist(el)
    igraph::vertex_attr(g, 'size') <- table(brID)[igraph::vertex_attr(g)$name]
    return(g)
}






