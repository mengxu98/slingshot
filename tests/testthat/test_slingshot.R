context("Test slingshot methods and SlingshotDataSet class.")

data("slingshotExample")
rd <- slingshotExample$rd
cl <- slingshotExample$cl
set.seed(1234)

# check for reordering

test_that("getLineages works for different input types", {
    reducedDim <- matrix(rnorm(100), ncol = 2)
    clusterLabels <- rep(seq_len(5), each = 10)

    # matrix / integer
    mi <- getLineages(reducedDim, clusterLabels)
    expect_is(mi, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(mi)))), 5)
    # 1-column matrix / integer
    m1i <- getLineages(reducedDim[,1,drop = FALSE], clusterLabels)
    expect_is(m1i, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(m1i)))), 5)
    # matrix / character
    mc <- getLineages(reducedDim, as.character(clusterLabels))
    expect_is(mc, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(mc)))), 5)
    # matrix / factor
    mf <- getLineages(reducedDim, as.factor(clusterLabels))
    expect_is(mf, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(mf)))), 5)
    # matrix / matrix
    cl.imb <- cbind(clusterLabels, sample(5,50, replace = TRUE))
    mm <- getLineages(reducedDim, cl.imb)
    expect_is(mm, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(mm)))), 2)



    df <- data.frame(reducedDim)
    # data frame / integer
    dfi <- getLineages(df, clusterLabels)
    expect_is(dfi, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(dfi)))), 5)
    # data frame / character
    dfc <- getLineages(df, as.character(clusterLabels))
    expect_is(dfc, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(dfc)))), 5)
    # data frame / factor
    dff <- getLineages(df, as.factor(clusterLabels))
    expect_is(dff, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(dff)))), 5)

    # SlingshotDataSet
    sds <- newSlingshotDataSet(reducedDim, clusterLabels)
    s <- getLineages(sds)
    expect_is(s, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(dff)))), 5)
    expect_equal(dim(slingMST(as.SlingshotDataSet(s))), c(5,5))

    # one cluster
    clus1 <- rep(1,50)
    c1 <- getLineages(reducedDim, clus1)
    expect_is(c1, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(c1)))), 1)

    # no clusters (default = make one cluster)
    c0 <- getLineages(reducedDim)
    expect_is(c0, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(c0)))), 1)
    
    # unclustered
    clmiss <- clusterLabels
    clmiss[1] <- -1
    unclus <- getLineages(reducedDim, clmiss)
    expect_equal(length(unique(unlist(slingLineages(unclus)))), 5)
    expect_true(all(slingClusterLabels(unclus)[1,] == rep(0,5)))
    clmiss[1] <- NA
    unclus <- getLineages(reducedDim, clmiss)
    expect_equal(length(unique(unlist(slingLineages(unclus)))), 5)
    expect_true(all(slingClusterLabels(unclus)[1,] == rep(0,5)))
    
    # with SingleCellExperiment objects
    require(SingleCellExperiment)
    u <- matrix(rpois(140*50, 5), nrow=50)
    sce <- SingleCellExperiment(assays=list(counts=u))
    expect_error(getLineages(sce), 'No dimensionality reduction found')

    reducedDims(sce) <- SimpleList(PCA = rd,
                                   tSNE = matrix(rnorm(140*2),ncol=2))
    # implicit reducedDim
    c0 <- getLineages(sce)
    expect_equal(length(unique(unlist(slingLineages(c0)))), 1)
    # reducedDim provided by name
    c0 <- getLineages(sce, reducedDim='tSNE')
    expect_equal(length(unique(unlist(slingLineages(c0)))), 1)
    # reducedDim provided as matrix
    c0 <- getLineages(sce, reducedDim = matrix(rnorm(140*2),ncol=2))
    expect_equal(length(unique(unlist(slingLineages(c0)))), 1)
    # cluster labels provided separately
    c0 <- getLineages(sce, clusterLabels = cl)
    expect_equal(length(unique(unlist(slingLineages(c0)))), 5)
    expect_true('slingshot' %in% names(colData(c0)))
    # accessor functions
    SlingshotDataSet(c0)
    expect_equal(length(slingLineages(c0)),2)
    expect_equal(length(slingCurves(c0)),0)
    expect_true(all(c('start.clus','end.clus','start.given','end.given',
                      'omega','omega_scale') %in% names(slingParams(c0)) ))

    # invalid inputs
    expect_error(getLineages(reducedDim[,-(seq_len(ncol(reducedDim)))],
                             clusterLabels), 'has zero columns')
    expect_error(getLineages(reducedDim[-(seq_len(nrow(reducedDim))),],
                             character(0)), 'has zero rows')
    expect_error(getLineages(reducedDim, clusterLabels[seq_len(10)]),
                 'must equal')
    rdna <- reducedDim; rdna[1,1] <- NA
    expect_error(getLineages(rdna, clusterLabels),
                 'cannot contain missing values')
    rdc <- reducedDim; rdc[1,1] <- 'a'
    expect_error(getLineages(rdc, clusterLabels),
                 'must only contain numeric values')
})

test_that("getLineages works as expected", {
    sds0 <- getLineages(rd, cl)
    expect_true(all(slingLineages(sds0)$Lineage1 == as.character(c(1,2,3,4))) ||
                    all(slingLineages(sds0)$Lineage1 == as.character(c(1,2,3,5))))
    expect_true(all(slingLineages(sds0)$Lineage2 == as.character(c(1,2,3,4))) ||
                    all(slingLineages(sds0)$Lineage2 == as.character(c(1,2,3,5))))
    expect_false(all(slingLineages(sds0)$Lineage1 ==
                         slingLineages(sds0)$Lineage2))
    # set start cluster
    sds1 <- getLineages(rd, cl, start.clus = 2)
    expect_true(all(vapply(slingLineages(sds1),function(l){ l[1] == '2' },
                           TRUE)))
    # set end cluster
    sds2 <- getLineages(rd, cl, start.clus = 1, end.clus = 3)
    expect_true(any(vapply(slingLineages(sds2),function(l){ (l[1] == '1') &&
            (l[length(l)] == '3') }, TRUE)))
    
    # omega
    # no effect
    sdsO <- getLineages(rd, cl, omega = 10)
    expect_identical(slingLineages(sdsO), slingLineages(sds0))
    # separate the clusters after the branching point
    sdsO <- getLineages(rd, cl, omega = 2.7)
    expect_identical(slingLineages(sdsO)$Lineage1, as.character(1:3))
    expect_true(all(as.character(4:5) %in% c(slingLineages(sdsO)$Lineage2,
                                             slingLineages(sdsO)$Lineage3)))
    # every cluster is its own lineage
    sdsO <- getLineages(rd, cl, omega = 2)
    expect_equal(length(slingLineages(sdsO)), 5)
    expect_true(all(as.character(1:5) %in% unlist(slingLineages(sdsO))))
    # with omega = TRUE
    # no effect
    sdsO <- getLineages(rd, cl, omega = TRUE)
    expect_identical(slingLineages(sdsO), slingLineages(sds0))
    # same as omega = 7, above
    sdsO <- getLineages(rd, cl, omega = TRUE, omega_scale = 1)
    expect_identical(slingLineages(sdsO)$Lineage1, as.character(1:3))
    expect_true(all(as.character(4:5) %in% c(slingLineages(sdsO)$Lineage2,
                                             slingLineages(sdsO)$Lineage3)))
    
    # two separate trajectories
    rd2 <- rbind(rd, cbind(rd[,2]-12, rd[,1]-6))
    cl2 <- c(cl, cl + 10)
    sds2 <- getLineages(rd2, cl2, omega = TRUE, start.clus = 11)
    expect_identical(slingParams(sds2)$start.clus, c('1','11'))
    expect_identical(slingParams(sds2)$start.given, c(FALSE,TRUE))
})

test_that("getCurves works as expected", {
    # 2 dim, 5 clus
    mi <- getLineages(rd, cl)
    mi <- getCurves(mi)
    expect_equal(length(slingCurves(mi)),2)

    # 3 lineages
    mi3 <- getLineages(rd, cl, end.clus = '3')
    mi3 <- getCurves(mi3)
    expect_equal(length(slingCurves(mi3)),3)

    # with SingleCellExperiment objects
    require(SingleCellExperiment)
    u <- matrix(rpois(140*50, 5), nrow=50)
    sce <- SingleCellExperiment(assays=list(counts=u))
    reducedDims(sce) <- SimpleList(PCA = rd,
                                   tSNE = matrix(rnorm(140*2),ncol=2))
    expect_error(getCurves(sce), 'No lineage information found')
    sce <- getLineages(sce, cl, 'PCA')
    sce <- getCurves(sce)
    expect_equal(length(slingCurves(sce)),2)

    # using approx_points produces similar curves
    mi_ap <- getCurves(mi, approx_points = 100)
    expect_true(cor(slingPseudotime(mi)[,2], slingPseudotime(mi_ap)[,2],
                    use='complete.obs') > .99)
    mi_ap <- getCurves(mi, approx_points = 300)
    expect_true(cor(slingPseudotime(mi)[,2], slingPseudotime(mi_ap)[,2],
                    use='complete.obs') > .99)

    # one dimension
    m1i <- getLineages(rd[,1,drop = FALSE], cl)
    m1i <- getCurves(m1i)
    expect_true(abs(abs(cor(slingReducedDim(m1i)[,1], slingPseudotime(m1i)[,1],
                            use='complete.obs'))-1) < .001)
    m1i <- getCurves(m1i, extend = 'n')
    expect_true(abs(abs(cor(slingReducedDim(m1i)[,1], slingPseudotime(m1i)[,1],
                            use='complete.obs'))-1) < .001)
    m1i <- getCurves(m1i, extend = 'pc1')
    expect_true(abs(abs(cor(slingReducedDim(m1i)[,1], slingPseudotime(m1i)[,1],
                            use='complete.obs'))-1) < .001)

    # one cluster
    clus1 <- cl; clus1[] <- 1
    c1 <- getLineages(rd, clus1)
    c1 <- getCurves(c1)
    expect_equal(length(slingCurves(c1)), 1)
    c1 <- getCurves(c1, extend = 'n')
    expect_equal(length(slingCurves(c1)), 1)
    c1 <- getCurves(c1, extend = 'pc1')
    expect_equal(length(slingCurves(c1)), 1)
})

test_that("slingshot works for different input types", {
    reducedDim <- matrix(rnorm(100), ncol = 2)
    clusterLabels <- rep(seq_len(5), length.out = 50)

    # matrix / integer
    mi <- slingshot(reducedDim, clusterLabels)
    expect_is(mi, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(mi)))), 5)
    # 1-column matrix / integer
    m1i <- slingshot(reducedDim[,1,drop = FALSE], clusterLabels)
    expect_is(m1i, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(m1i)))), 5)
    # matrix / character
    mc <- slingshot(reducedDim, as.character(clusterLabels))
    expect_is(mc, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(mc)))), 5)
    # matrix / factor
    mf <- slingshot(reducedDim, as.factor(clusterLabels))
    expect_is(mf, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(mf)))), 5)
    # matrix / matrix
    cl.imb <- cbind(clusterLabels, sample(3,50, replace = TRUE))
    mm <- slingshot(reducedDim, cl.imb)
    expect_is(mm, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(mm)))), 2)


    df <- data.frame(reducedDim)
    # data frame / integer
    dfi <- slingshot(df, clusterLabels)
    expect_is(dfi, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(dfi)))), 5)
    # data frame / character
    dfc <- slingshot(df, as.character(clusterLabels))
    expect_is(dfc, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(dfc)))), 5)
    # data frame / factor
    dff <- slingshot(df, as.factor(clusterLabels))
    expect_is(dff, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(dff)))), 5)

    sds <- newSlingshotDataSet(reducedDim, clusterLabels)
    # SlingshotDataSet
    s <- slingshot(sds)
    expect_is(s, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(s)))), 5)

    # diagonal distance measure
    slingshot(rd, cl, dist.method = "scaled.diag")
    # different shrinkage methods
    slingshot(rd, cl, shrink.method = 'tricube')
    slingshot(rd, cl, shrink.method = 'density')

    # one cluster
    clus1 <- rep(1,50)
    c1 <- slingshot(reducedDim, clus1)
    expect_is(c1, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(c1)))), 1)

    # no clusters (default = make one cluster)
    c0 <- slingshot(reducedDim)
    expect_is(c1, "PseudotimeOrdering")
    expect_equal(length(unique(unlist(slingLineages(c1)))), 1)

    # using approx_points produces similar curves
    v1 <- slingshot(rd, cl)
    v2 <- slingshot(rd, cl, approx_points = 100) # a_p < n
    expect_true(cor(slingPseudotime(v1)[,2], slingPseudotime(v2)[,2],
                    use='complete.obs') > .99)
    v2 <- slingshot(rd, cl, approx_points = 300) # a_p > n
    expect_true(cor(slingPseudotime(v1)[,2], slingPseudotime(v2)[,2],
                    use='complete.obs') > .99)

    # invalid inputs
    expect_error(slingshot(reducedDim[,-(seq_len(ncol(reducedDim)))],
                           clusterLabels), 'has zero columns')
    expect_error(slingshot(reducedDim[-(seq_len(nrow(reducedDim))),],
                           character(0)), 'has zero rows')
    expect_error(slingshot(reducedDim, clusterLabels[seq_len(10)]), 'must equal')
    rdna <- reducedDim; rdna[1,1] <- NA
    expect_error(slingshot(rdna, clusterLabels),
                 'cannot contain missing values')
    rdc <- reducedDim; rdc[1,1] <- 'a'
    expect_error(slingshot(rdc, clusterLabels),
                 'must only contain numeric values')

    # with SingleCellExperiment objects
    require(SingleCellExperiment)
    u <- matrix(rpois(140*50, 5), nrow=50)
    sce <- SingleCellExperiment(assays=list(counts=u))
    expect_error(slingshot(sce), 'No dimensionality reduction found')

    reducedDims(sce) <- SimpleList(PCA = rd,
                                   tSNE = matrix(rnorm(140*2),ncol=2))
    # implicit reducedDim
    c0 <- slingshot(sce)
    expect_equal(length(unique(unlist(slingLineages(c0)))), 1)
    # reducedDim provided by name
    c0 <- slingshot(sce, reducedDim='tSNE')
    expect_equal(length(unique(unlist(slingLineages(c0)))), 1)
    # reducedDim provided as matrix
    c0 <- slingshot(sce, reducedDim = matrix(rnorm(140*2),ncol=2))
    expect_equal(length(unique(unlist(slingLineages(c0)))), 1)
    # cluster labels provided separately
    c0 <- slingshot(sce, clusterLabels = cl)
    expect_equal(length(unique(unlist(slingLineages(c0)))), 5)
    expect_true('slingshot' %in% names(colData(c0)))
    # accessor functions
    SlingshotDataSet(c0)
    expect_equal(length(slingLineages(c0)),2)
    expect_equal(length(slingCurves(c0)),2)
    expect_true(all(c('start.clus','end.clus','start.given','end.given',
                      'shrink','extend','reweight','reassign',
                      'shrink.method') %in% names(slingParams(c0)) ))
    expect_equal(dim(slingPseudotime(c0)), c(140,2))
    expect_equal(dim(slingCurveWeights(c0)), c(140,2))
})

test_that("slingshot works with ClusterExperiment objects", {
    if(! requireNamespace('clusterExperiment', quietly = TRUE)){
        skip('clusterExperiment package not available.')
    }
    require(SingleCellExperiment)
    require(clusterExperiment)
    u <- matrix(rpois(140*50, 5), nrow=50)
    sce <- SingleCellExperiment(assays=list(counts=u))
    reducedDims(sce) <- SimpleList(PCA = rd,
                                   tSNE = matrix(rnorm(140*2),ncol=2))
    ce <- clusterExperiment::ClusterExperiment(sce, factor(cl),
                                               transformation = function(x){x})
    ce.sling <- slingshot(ce)
    expect_is(ce.sling, "ClusterExperiment")
    ce.sling <- slingshot(ce, reducedDim = 'tSNE')
    expect_is(ce.sling, "ClusterExperiment")
    ce.sling <- slingshot(ce, reducedDim = matrix(rnorm(140*2),ncol=2))
    expect_is(ce.sling, "ClusterExperiment")

    colData(ce) <- cbind(colData(ce), cl2 = sample(2,140, replace=TRUE))
    ce.sling <- slingshot(ce, 'cl2')
    expect_is(ce.sling, "ClusterExperiment")
    ce.sling <- slingshot(ce, sample(2,140, replace=TRUE))
    expect_is(ce.sling, "ClusterExperiment")
    expect_true('slingshot' %in% names(colData(ce.sling)))
})

test_that("2D Plotting functions don't give errors", {
    pto <- slingshot(rd,cl, start.clus = '1', end.clus = c('4','5'))
    sds <- as.SlingshotDataSet(pto)
    
    u <- matrix(rpois(140*50, 5), nrow=50)
    sce <- SingleCellExperiment(assays=list(counts=u))
    reducedDims(sce) <- SimpleList(PCA = rd)
    sce <- slingshot(sce, clusterLabels = cl)

    plot(sds)
    expect_error(plot(sds, linInd = 3:5),
                 'None of the provided lineage indices')
    plot(sds, type = "lineages", show.constraints = TRUE, linInd = 2:3)
    plot(sds, type = "lineages", show.constraints = TRUE)
    lines(sds, linInd = 2)
    lines(sds, type = "lineages", show.constraints = TRUE)
    lines(sds, type = "lineages", show.constraints = TRUE, linInd = c(1,3))
    pairs(sds, lower.panel = TRUE)
    pairs(sds, lower.panel = TRUE, type = "lineages", show.constraints = TRUE)
})

test_that("3D Plotting functions don't give errors", {
    if(! requireNamespace('rgl', quietly = TRUE)){
        skip('rgl package not available.')
    }
    rd3 <- cbind(rd, rnorm(140))
    sds3 <- as.SlingshotDataSet(slingshot(rd3, cl))
    plot3d.SlingshotDataSet(sds3)
    plot3d.SlingshotDataSet(sds3, type = 'lineages')
    plot3d.SlingshotDataSet(sds3, linInd = 1)
    plot3d.SlingshotDataSet(sds3, type = 'lineages', linInd = 2:3)
    expect_error(plot3d.SlingshotDataSet(sds3, linInd = 3:5),
                 'None of the provided lineage indices')
    rgl::plot3d(slingReducedDim(sds3))
    plot3d.SlingshotDataSet(sds3, add = TRUE)
})

test_that("predict works as expected", {
    sds <- slingshot(rd, cl)
    pred <- predict(sds)
    expect_identical(sds, pred)

    x.mat <- cbind(runif(100, min = -5, max = 10),
                   runif(100, min = -4, max = 4))
    pred <- predict(sds, x.mat)
    expect_true(all(slingClusterLabels(pred)==0))
    expect_equal(length(slingLineages(pred)), 2)
    expect_equal(length(slingCurves(pred)), 2)

    x.df <- as.data.frame(x.mat)
    pred <- predict(sds, x.df)
    expect_equal(length(slingCurves(pred)), 2)

    # invalid inputs
    x.text <- x.mat
    x.text[1,1] <- 'text'
    expect_error(predict(sds, x.text), 'must only contain numeric values')

    x.na <- x.mat
    x.na[1,1] <- NA
    expect_error(predict(sds, x.na), 'cannot contain missing values')

    x.big <- cbind(x.mat, rnorm(100))
    expect_error(predict(sds, x.big),
                 'does not match original number of dimensions')
})

test_that("Helper functions work as expected", {
    data("slingshotExample")
    rd <- slingshotExample$rd
    cl <- slingshotExample$cl
    pto <- slingshot(rd,cl, start.clus = '1', end.clus = c('4','5'))
    sds <- SlingshotDataSet(pto)
    show(sds)

    expect_equal(dim(reducedDim(sds)), c(140,2))
    expect_equal(dim(slingReducedDim(sds)), c(140,2))
    expect_equal(dim(slingClusterLabels(sds)), c(140,5))
    expect_equal(length(slingLineages(sds)),2)
    
    expect_equal(dim(slingMST(sds)), c(5,5))
    expect_s3_class(slingMST(pto), 'igraph')
    expect_s3_class(slingMST(sds, as.df = TRUE), 'data.frame')
    expect_s3_class(slingMST(pto, as.df = TRUE), 'data.frame')    
    expect_equal(ncol(slingMST(pto, as.df = TRUE)), 
                 ncol(slingReducedDim(pto)) + 3)
    expect_equal(sort(unique(slingMST(pto, as.df = TRUE)$Cluster)),
                 as.character(1:5))
    expect_equal(length(slingCurves(sds)),2)
    expect_s3_class(slingCurves(sds, as.df = TRUE), 'data.frame')
    expect_s3_class(slingCurves(pto, as.df = TRUE), 'data.frame')    
    expect_equal(ncol(slingCurves(pto, as.df = TRUE)), 
                 ncol(slingReducedDim(pto)) + 2)
    expect_equal(unique(slingCurves(pto, as.df = TRUE)$Lineage), 1:2)
    expect_true(all(c('start.clus','end.clus','start.given','end.given',
                      'shrink','extend','reweight','reassign',
                      'shrink.method') %in% names(slingParams(sds)) ))
    expect_equal(dim(slingPseudotime(sds)), c(140,2))
    expect_equal(sum(is.na(slingPseudotime(sds, na = FALSE))), 0)
    expect_equal(dim(slingCurveWeights(sds)), c(140,2))
    expect_true(all(
        abs(rowSums(slingCurveWeights(sds, as.probs = TRUE))-1) < .001))

    # newSlingshotDataSet
    # matrix / factor
    mf <- newSlingshotDataSet(rd, factor(cl))
    expect_is(mf, "SlingshotDataSet")
    # matrix / missing
    expect_message({m0 <- newSlingshotDataSet(rd)},
                   "Unclustered data detected.")
    expect_is(m0, "SlingshotDataSet")

    # data frame / character
    dfc <- newSlingshotDataSet(data.frame(rd))
    expect_is(dfc, "SlingshotDataSet")
    # data frame / missing
    expect_message({df0 <- newSlingshotDataSet(data.frame(rd))},
                   "Unclustered data detected.")
    expect_is(df0, "SlingshotDataSet")

    # matrix / matrix
    cl.mat <- outer(cl, unique(cl), '==') + 0.0
    rownames(cl.mat) <- NULL
    colnames(rd) <- NULL
    expect_error(newSlingshotDataSet(rd, cl.mat[-1,]), 'must equal')
})

test_that("embedCurves works as expected", {
    data("slingshotExample")
    rd <- slingshotExample$rd
    cl <- slingshotExample$cl
    sds <- getLineages(rd, cl)
    tsne <- rd + rnorm(nrow(rd)*2)

    # before running getCurves
    expect_error(embedCurves(sds, tsne), 'No slingshot curves found')
    sds <- getCurves(sds)
    # shrink argument out of bounds
    expect_error(embedCurves(sds, tsne, shrink = 2), 'numeric between 0 and 1')
    # wrong number of cells in newX
    expect_error(embedCurves(sds, tsne[-1, ]), 'must have same number of rows')
    # NAs in newX
    tsne2 <- tsne; tsne2[2,2] <- NA
    expect_error(embedCurves(sds, tsne2), 'cannot contain missing values')
    # non-numeric values in newX
    tsne2[2,2] <- 'a'
    expect_error(embedCurves(sds, tsne2), 'must only contain numeric values')
    # missing row/col names
    rownames(tsne)[1] <- ''
    colnames(tsne)[1] <- ''
    emb <- embedCurves(sds, tsne)
    rownames(tsne) <- NULL
    colnames(tsne) <- NULL
    emb <- embedCurves(sds, tsne)
    # approx_points
    emb <- embedCurves(sds, tsne, approx_points = 50)
    # loess
    emb <- embedCurves(sds, tsne, smoother = 'loess')

    # with SingleCellExperiment objects
    require(SingleCellExperiment)
    u <- matrix(rpois(140*50, 5), nrow=50)
    sce <- SingleCellExperiment(assays=list(counts=u))
    reducedDims(sce) <- SimpleList(PCA = rd,
                                   tSNE = tsne)
    # before running slingshot
    expect_error(embedCurves(sce, tsne), 'No slingshot results found')
    expect_error(embedCurves(sce, 'tSNE'), 'No slingshot results found')
    # after
    sce <- slingshot(sce, cl, 'PCA')
    emb <- embedCurves(sce, tsne)
    emb <- embedCurves(sce, 'tSNE')
    
    expect_is(emb, 'PseudotimeOrdering')
    expect_equal(length(slingCurves(emb)), 2)
})

test_that("branchID functions work as expected", {
    sds <- slingshot(rd, cl)
    
    # bad thresh
    expect_error(slingBranchID(sds, thresh = 5), 'between 0 and 1')
    expect_error(slingBranchID(sds, thresh = -1), 'between 0 and 1')
    # odd thresh
    expect_identical(levels(slingBranchID(sds, thresh = 0)), "1,2")
    g <- slingBranchGraph(sds, thresh = 0)
    expect_identical(names(g[[1]]), "1,2")
    
    id <- slingBranchID(sds)
    expect_equal(levels(id), c('1','1,2','2'))
    
    g <- slingBranchGraph(sds)
    expect_true(all(c('name','cells','size') %in% 
                        names(igraph::vertex_attr(g))))
    
    require(SingleCellExperiment)
    u <- matrix(rpois(140*50, 5), nrow=50)
    sce <- SingleCellExperiment(assays=list(counts=u))
    reducedDims(sce) <- SimpleList(PCA = rd)
    sce <- slingshot(sce, cl, 'PCA')
    
    id <- slingBranchID(sce)
    expect_equal(levels(id), c('1','1,2','2'))
    
    g <- slingBranchGraph(sce)
    expect_true(all(c('name','cells','size') %in% 
                        names(igraph::vertex_attr(g))))
    
    # one cluster
    sds1 <- slingshot(rd)
    id1 <- slingBranchID(sds1)
    expect_true(all(id1 == 1))
    g <- slingBranchGraph(sds1)
    expect_identical(igraph::vertex_attr(g)$name, '1')
    
    # case with missing intermediates (ie. '1,2,3' goes directly to '1')
    rd2 <- rbind(rd, c(-8.1,.1), c(-8,-.1))
    cl2 <- c(cl, 6,6)
    g <- slingBranchGraph(slingshot(rd2, cl2, start.clus = '1'))
    expect_true(all(c("1,2","1","2","1,2,3","3") %in% names(g[[1:5]])))
    expect_true(all(names(g[[1:5]]) %in% c("1,2","1","2","1,2,3","3")))
})

test_that("conversion functions work as expected", {
    require(SingleCellExperiment)
    u <- matrix(rpois(140*50, 5), nrow=50)
    sce <- SingleCellExperiment(assays=list(counts=u))
    reducedDims(sce) <- SimpleList(PCA = rd,
                                   tSNE = matrix(rnorm(140*2),ncol=2))
    sce <- slingshot(sce, reducedDim = 'PCA', clusterLabels = cl)
    pto <- as.PseudotimeOrdering(sce)
    sds <- as.SlingshotDataSet(pto)
    
    expect_identical(slingPseudotime(sce), slingPseudotime(pto))
    expect_identical(slingPseudotime(sce), slingPseudotime(sds))
    expect_identical(slingCurveWeights(sce), slingCurveWeights(pto))
    expect_identical(slingCurveWeights(sce), slingCurveWeights(sds))
    
    expect_identical(sds, as.SlingshotDataSet(sce))
    
    pto2 <- as.PseudotimeOrdering(sds)
    expect_identical(slingPseudotime(pto), slingPseudotime(pto2))
    expect_identical(slingCurveWeights(pto), slingCurveWeights(pto2))    
    
    expect_identical(sds, as.SlingshotDataSet(sds))
    expect_identical(pto, as.PseudotimeOrdering(pto))
    
    expect_false(all(assay(pto, 'weights') %in% c(0,1)))
    pto2 <- slingshot(rd, cl)
    expect_false(all(assay(pto2, 'weights') %in% c(0,1)))
    
    # old SCE integration
    colData(sce)$slingshot <- NULL
    int_metadata(sce)$slingshot <- sds
    expect_identical(sds, as.SlingshotDataSet(sce))
    pto2 <- as.PseudotimeOrdering(sce)
    expect_identical(slingPseudotime(pto), slingPseudotime(pto2))
    expect_identical(slingCurveWeights(pto), slingCurveWeights(pto2))    
    
    # no Slingshot results
    int_metadata(sce)$slingshot <- NULL
    expect_error(as.SlingshotDataSet(sce), 'No slingshot results found')
    
    # partial results
    sds <- getLineages(rd, cl)
    pto <- as.PseudotimeOrdering(sds)
    expect_identical(slingLineages(sds), slingLineages(pto))
    expect_true(all(is.na(assay(pto, 'pseudotime'))))
    expect_true(all(assay(pto, 'weights') %in% c(0,1)))
    
    # no results
    sds <- newSlingshotDataSet(rd, cl)
    expect_error(as.PseudotimeOrdering(sds), 
                 'number of lineages could not be determined')
    
    })

