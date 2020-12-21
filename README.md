CyTOF Pipeline
================

**This README was made to check the data and make an example workflow.
To run the shiny app use** `shiny::runApp()`.

We are using renv with this project. You can load the library version
from the lockfile using the `restore` function:

``` r
# install.packages('renv')
renv::restore()
```

### Mouse Data

Download the mouse data set and save it at in extdata. *Sorry no idea
how to download from google drive using R.*

Check out the data.

``` r
suppressPackageStartupMessages(library(CATALYST))
suppressPackageStartupMessages(library(data.table))
pop_assign <- data.table::fread("extdata/Panorama BM 1-10/population_assignments.txt", header = FALSE)
sce <- CATALYST::prepData("extdata/Panorama BM 1-10/")
```
