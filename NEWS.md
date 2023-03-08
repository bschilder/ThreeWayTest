# ThreeWayTest 1.0.2

## New features

* Split functions into separate files (make easier to find them).
* Increase coverage to 100%
* Update README with icicle plot.
* Add `get_full_genotype` func to get full 1KG genotype data.
    - Added support function `get_data` 

## Bug fixes

* Add missing *tests/testthat.R* script.

# ThreeWayTest 1.0.1

## New features

* Add hex sticker!
* Add Docker vignette.
* Add link to Docker vignette in README.
* New functions:
    - `clustermap`
    - `postprocess_data`
    - `networkmap`

## Bug fixes

* Shorten function widths
* Shorten Roxygen note widths
* Remove `templateR` instructions from README
* Recompress built-in data with `tools::resaveRdaFiles`
* Fix `Dependence on R version ‘4.0.2’ not with patchlevel 0`
* Rename vignette *getting_started.Rmd* to *ThreeWayTest.Rmd* convention.
* Organize *getting_started.Rmd* chunks.


# ThreeWayTest 0.99.0

## New features
 
* Added a `NEWS.md` file to track changes to the package.

## Bug fixes

* Removed Dockerfile and switched to using `rworkflows`.