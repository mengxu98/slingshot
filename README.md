# R package: slingshot
[![R build status](https://github.com/kstreet13/slingshot/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/kstreet13/slingshot/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/kstreet13/slingshot/master.svg)](https://codecov.io/github/kstreet13/slingshot?branch=master)

<img src=inst/slingshot_sticker.png height="200">

Provides functions for inferring continuous, branching lineage structures in low-dimensional data. Slingshot was designed to model developmental trajectories in single-cell RNA sequencing data and serve as a component in an analysis pipeline after dimensionality reduction and clustering. It is flexible enough to handle arbitrarily many branching events and allows for the incorporation of prior knowledge through supervised graph construction.

## Installation

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("kstreet13/slingshot")
```

## Python Version
An implementation of the Slingshot algorithm in `python` can be found here: https://github.com/mossjacob/pyslingshot

## Issues and bug reports

Please use https://github.com/kstreet13/slingshot/issues to submit issues, bug reports, and comments.
