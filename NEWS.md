# changes in version 2.0.0

* Added a `NEWS.md` file to track changes to the package.

* Changed default output of most functions from `SlingshotDataSet` to `PseudotimeOrdering` and added conversion functions between them and `SingleCellExperiment`.

* `getLineages` now relies on `createClusterMST`

* Removed `plotGenePseudotime`

* added `as.df` option to `slingCurves` and `slingMST`, which provide relevant
information for plotting as `data.frame` objects. This should help those
plotting Slingshot results with `ggplot`, especially for our `traffic` package.

* updated all documentation.

