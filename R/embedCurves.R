#' @rdname embedCurves
#'
#' @description This function takes the output of \code{\link{slingshot}} (or
#'   \code{\link{getCurves}}) and attempts to embed the curves in a different
#'   coordinate space than the one in which they were constructed. This should
#'   be considered a visualization tool, only.
#'
#' @param x an object containing \code{\link{slingshot}} output.
#' @param newDimRed a matrix representing the new coordinate space in which to
#'   embed the curves.
#' @param shrink logical or numeric between 0 and 1, determines whether and how
#'   much to shrink branching lineages toward their average prior to the split.
#' @param stretch numeric factor by which curves can be extrapolated beyond
#'   endpoints. Default is \code{2}, see
#'   \code{\link[princurve]{principal_curve}}.
#' @param approx_points numeric, whether curves should be approximated by a
#'   fixed number of points. If \code{FALSE} (or 0), no approximation will be
#'   performed and curves will contain as many points as the input data. If
#'   numeric, curves will be approximated by this number of points; preferably
#'   about 100 (see \code{\link[princurve]{principal_curve}}).
#' @param smoother, choice of scatter plot smoother. Same as
#'   \code{\link[princurve]{principal_curve}}, but \code{"lowess"} option is
#'   replaced with \code{"loess"} for additional flexibility.
#' @param shrink.method character denoting how to determine the appropriate
#'   amount of shrinkage for a branching lineage. Accepted values are the same
#'   as for \code{kernel} in \code{\link{density}} (default is \code{"cosine"}),
#'   as well as \code{"tricube"} and \code{"density"}. See 'Details' for more.
#' @param ... Additional parameters to pass to scatter plot smoothing function,
#'   \code{smoother}.
#'
#' @details Many of the same parameters are used here as in \code{getCurves}.
#'   This function attempts to translate curves from one reduced dimensional
#'   space to another by predicting each dimension as a function of pseudotime
#'   (ie. the new curve is determined by a series of scatterplot smoothers
#'   predicting the coordinates in the new space as a function of pseudotime).
#'   Because the pseudotime values are not changed, this amounts to a single
#'   iteration of the iterative curve-fitting process used by \code{getCurves}.
#'
#' @details Note that non-linear dimensionality reduction techniques (such as
#'   tSNE and UMAP) may produce discontinuities not observed in other spaces.
#'   Use caution when embedding curves in these spaces.
#'
#' @return a \code{\link{PseudotimeOrdering}} object containing curves in the
#'   new space.
#'
#' @examples
#' data("slingshotExample")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' pto <- slingshot(rd, cl, start.clus = '1')
#' rd2 <- cbind(rd[,2] + rnorm(nrow(rd)), -rd[,1] + rnorm(nrow(rd)))
#' pto.new <- embedCurves(pto, rd2)
#' pto.new
#'
#' plot(rd2, col = cl, asp = 1)
#' lines(SlingshotDataSet(pto.new), lwd = 3)
#'
#' @importFrom princurve project_to_curve
#' @export
setMethod(f = "embedCurves",
          signature = signature(x = "PseudotimeOrdering",
                                newDimRed = "matrix"),
          definition = function(x, newDimRed,
                                shrink = NULL, stretch = NULL,
                                approx_points = NULL, smoother = NULL,
                                shrink.method = NULL, ...){
              # SETUP for checks
              pto <- x
              X <- slingReducedDim(pto)
              newX <- newDimRed
              
              # if new arguments are not provided, use existing arguments
              if(is.null(shrink)){
                  shrink <- slingParams(pto)$shrink
              }
              # some were not previously included in slingParams output, so we
              # assume the default values were used
              if(is.null(stretch)){
                  stretch <- ifelse(is.null(slingParams(pto)$stretch), 2,
                                    slingParams(pto)$stretch)
              }
              if(is.null(approx_points)){
                  approx_points <- ifelse(
                      is.null(slingParams(pto)$approx_points),
                      FALSE, slingParams(pto)$approx_points)
              }
              if(is.null(smoother)){
                  smoother <- ifelse(
                      is.null(slingParams(pto)$smoother),
                      'smooth.spline', slingParams(pto)$smoother)
              }
              if(is.null(shrink.method)){
                  shrink.method <- slingParams(pto)$shrink.method
              }
              
              # CHECKS
              if(length(slingCurves(pto)) == 0){
                  stop("No slingshot curves found in original space.")
              }
              if(shrink < 0 | shrink > 1){
                  stop("'shrink' parameter must be logical or numeric between",
                       " 0 and 1")
              }
              if(nrow(X)!=nrow(newX)){
                  stop("'newX' must have same number of rows as original",
                       "'reducedDim'.")
              }
              if(any(is.na(newX))){
                  stop("'newX' cannot contain missing values.")
              }
              if(!all(apply(newX,2,is.numeric))){
                  stop("'newX' must only contain numeric values.")
              }
              if(is.null(rownames(newX))){
                  rownames(newX) <- rownames(X)
              }
              if(is.null(colnames(newX))){
                  colnames(newX) <- paste('Dim',seq_len(ncol(newX)),sep='-')
              }
              if(any(rownames(newX)=='')){
                  miss.ind <- which(rownames(newX) == '')
                  rownames(newX)[miss.ind] <- rownames(X)[miss.ind]
              }
              if(any(colnames(newX)=='')){
                  miss.ind <- which(colnames(newX) == '')
                  colnames(newX)[miss.ind] <- paste('Dim',miss.ind,sep='-')
              }
              
              # SETUP
              p.new <- ncol(newX)
              lineages <- slingLineages(pto)
              L <- length(grep("Lineage", names(lineages))) #number of lineages
              clusters <- colnames(slingClusterLabels(pto))
              d <- dim(X); n <- d[1]
              
              C <- as.matrix(vapply(lineages[seq_len(L)], function(lin) {
                  vapply(clusters, function(clID) {
                      as.numeric(clID %in% lin)
                  }, 0)
              }, rep(0,length(clusters))))
              rownames(C) <- clusters
              segmnts <- unique(C[rowSums(C)>1,,drop = FALSE])
              segmnts <- segmnts[order(rowSums(segmnts),decreasing = FALSE), ,
                                 drop = FALSE]
              avg.order <- list()
              for(i in seq_len(nrow(segmnts))){
                  idx <- segmnts[i,] == 1
                  avg.order[[i]] <- colnames(segmnts)[idx]
                  new.col <- rowMeans(segmnts[,idx, drop = FALSE])
                  segmnts <- cbind(segmnts[, !idx, drop = FALSE],new.col)
                  colnames(segmnts)[ncol(segmnts)] <- paste('average',i,sep='')
              }
              
              # DEFINE SMOOTHER FUNCTION
              smootherFcn <- switch(smoother, loess = function(lambda, xj,
                                                               w = NULL, ...){
                  loess(xj ~ lambda, weights = w, ...)$fitted
              }, smooth.spline = function(lambda, xj, w = NULL, ..., df = 5,
                                          tol = 1e-4){
                  # fit <- smooth.spline(lambda, xj, w = w, ..., df = df,
                  #                      tol = tol, keep.data = FALSE)
                  fit <- tryCatch({
                      smooth.spline(lambda, xj, w = w, ..., df = df,
                                    tol = tol, keep.data = FALSE)
                  }, error = function(e){
                      smooth.spline(lambda, xj, w = w, ..., df = df,
                                    tol = tol, keep.data = FALSE, spar = 1)
                  })
                  predict(fit, x = lambda)$y
              })
              
              pcurves <- slingCurves(pto)
              
              # for each curve,
              #   construct a new curve by predicting each (new) dimension as a
              #   function of pseudotime.
              for(l in seq_len(L)){
                  pcurve <- pcurves[[l]]
                  ordL <- order(pcurve$lambda)
                  s <- matrix(NA, nrow = n, ncol = p.new)
                  
                  if(approx_points > 0){
                      xout_lambda <- seq(min(pcurve$lambda), max(pcurve$lambda),
                                         length.out = approx_points)
                      s <- matrix(NA, nrow = approx_points, ncol = p.new)
                  }
                  for(jj in seq_len(p.new)){
                      yjj <- smootherFcn(pcurve$lambda, newX[,jj], w = pcurve$w,
                                         ...)[ordL]
                      if(approx_points > 0){
                          yjj <- approx(x = pcurve$lambda[ordL], y = yjj,
                                        xout = xout_lambda, ties = 'ordered')$y
                      }
                      s[, jj] <- yjj
                  }
                  new.pcurve <- project_to_curve(newX, s = s, stretch = stretch)
                  if(approx_points > 0){
                      xout_lambda <- seq(min(new.pcurve$lambda),
                                         max(new.pcurve$lambda),
                                         length.out = approx_points)
                      new.pcurve$s <- apply(new.pcurve$s, 2, function(sjj){
                          return(approx(x = new.pcurve$lambda[new.pcurve$ord],
                                        y = sjj[new.pcurve$ord],
                                        xout = xout_lambda, ties = 'ordered')$y)
                      })
                      new.pcurve$ord <- seq_len(approx_points)
                  }
                  new.pcurve$dist_ind <- abs(new.pcurve$dist_ind)
                  new.pcurve$lambda <- new.pcurve$lambda -
                      min(new.pcurve$lambda, na.rm = TRUE)
                  new.pcurve$w <- pcurve$w
                  pcurves[[l]] <- new.pcurve
              }
              
              # shrink together lineages near shared cells
              if(shrink > 0){
                  if(max(rowSums(C)) > 1){
                      
                      segmnts <- unique(C[rowSums(C)>1,,drop=FALSE])
                      segmnts <- segmnts[order(rowSums(segmnts),
                                               decreasing = FALSE),
                                         , drop = FALSE]
                      seg.mix <- segmnts
                      avg.lines <- list()
                      pct.shrink <- list()
                      
                      # determine average curves and amount of shrinkage
                      for(i in seq_along(avg.order)){
                          ns <- avg.order[[i]]
                          to.avg <- lapply(ns,function(n){
                              if(grepl('Lineage',n)){
                                  l.ind <- as.numeric(gsub('Lineage','',n))
                                  return(pcurves[[l.ind]])
                              }
                              if(grepl('average',n)){
                                  a.ind <- as.numeric(gsub('average','',n))
                                  return(avg.lines[[a.ind]])
                              }
                          })
                          avg <- .avg_curves(to.avg, newX, stretch = stretch,
                                             approx_points = approx_points)
                          avg.lines[[i]] <- avg
                          common.ind <- rowMeans(
                              vapply(to.avg, function(crv){ crv$w > 0 },
                                     rep(TRUE, n))) == 1
                          pct.shrink[[i]] <- lapply(to.avg,function(crv){
                              .percent_shrinkage(crv, common.ind,
                                                 approx_points = approx_points,
                                                 method = shrink.method)
                          })
                          # check for degenerate case (if one curve won't be
                          # shrunk, then the other curve shouldn't be,
                          # either)
                          all.zero <- vapply(pct.shrink[[i]], function(pij){
                              return(all(pij == 0))
                          }, TRUE)
                          if(any(all.zero)){
                              pct.shrink[[i]] <- lapply(pct.shrink[[i]],
                                                        function(pij){
                                                            pij[] <- 0
                                                            return(pij)
                                                        })
                          }
                      }
                      # do the shrinking in reverse order
                      for(j in rev(seq_along(avg.lines))){
                          ns <- avg.order[[j]]
                          avg <- avg.lines[[j]]
                          to.shrink <- lapply(ns,function(n){
                              if(grepl('Lineage',n)){
                                  l.ind <- as.numeric(gsub('Lineage','',n))
                                  return(pcurves[[l.ind]])
                              }
                              if(grepl('average',n)){
                                  a.ind <- as.numeric(gsub('average','',n))
                                  return(avg.lines[[a.ind]])
                              }
                          })
                          shrunk <- lapply(seq_along(ns),function(jj){
                              crv <- to.shrink[[jj]]
                              return(.shrink_to_avg(crv, avg,
                                            pct.shrink[[j]][[jj]] * shrink,
                                            newX, approx_points = approx_points,
                                            stretch = stretch))
                          })
                          for(jj in seq_along(ns)){
                              n <- ns[jj]
                              if(grepl('Lineage',n)){
                                  l.ind <- as.numeric(gsub('Lineage','',n))
                                  pcurves[[l.ind]] <- shrunk[[jj]]
                              }
                              if(grepl('average',n)){
                                  a.ind <- as.numeric(gsub('average','',n))
                                  avg.lines[[a.ind]] <- shrunk[[jj]]
                              }
                          }
                      }
                  }
              }
              
              # use the new curves, but keep existing pseudotime, weights, etc.
              newCurves <- lapply(seq_len(L), function(l){
                  crv <- list(s = pcurves[[l]]$s,
                              ord = pcurves[[l]]$ord,
                              lambda = slingCurves(pto)[[l]]$lambda,
                              dist_ind = slingCurves(pto)[[l]]$dist_ind,
                              dist = slingCurves(pto)[[l]]$dist,
                              w = slingCurves(pto)[[l]]$w)
                  class(crv) <- "principal_curve"
                  return(crv)
              })
              
              params <- slingParams(pto)
              params$shrink <- shrink
              params$stretch <- stretch
              params$approx_points <- approx_points
              params$smoother <- smoother
              params$shrink.method <- shrink.method
              params$embedding <- TRUE
              
              out <- pto
              metadata(out)$curves <- newCurves
              metadata(out)$slingParams <- params
              cellData(out)$reducedDim <- newX
              
              validObject(out)
              return(out)
          })


#' @rdname embedCurves
#' @export
setMethod(f = "embedCurves",
          signature = signature(x = "SingleCellExperiment",
                                newDimRed = "matrix"),
          definition = function(x, newDimRed,
                                shrink = NULL, stretch = NULL,
                                approx_points = NULL, smoother = NULL,
                                shrink.method = NULL, ...){
              embedCurves(x = as.PseudotimeOrdering(x), 
                          newDimRed = newDimRed,
                          shrink = shrink,
                          stretch = stretch,
                          approx_points = approx_points,
                          smoother = smoother,
                          shrink.method = shrink.method, ...)          
          })

#' @rdname embedCurves
#' @importFrom SingleCellExperiment reducedDim
#' @export
setMethod(f = "embedCurves",
          signature = signature(x = "SingleCellExperiment",
                                newDimRed = "character"),
          definition = function(x, newDimRed,
                                shrink = NULL, stretch = NULL,
                                approx_points = NULL, smoother = NULL,
                                shrink.method = NULL, ...){
              embedCurves(x = as.PseudotimeOrdering(x), 
                          newDimRed = reducedDim(x, newDimRed),
                          shrink = shrink,
                          stretch = stretch,
                          approx_points = approx_points,
                          smoother = smoother,
                          shrink.method = shrink.method, ...)
          })

