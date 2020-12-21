CyTOF Pipeline
================

We are using renv with this project. You can load the library version
from the lockfile using the `restore` function:

``` r
# install.packages('renv')
renv::restore()
```

### Mouse Data

Download the mouse data set and save it at under extdata. *Sorry no idea
how to download from google drive using R.*

Check out the data.

``` r
suppressPackageStartupMessages(library(CATALYST))
suppressPackageStartupMessages(library(data.table))
pop_assign <- data.table::fread("~/Downloads/SystemsBioMedicine/Panorama BM 1-10/population_assignments.txt", header = FALSE)
sce <- CATALYST::prepData("~/Downloads/SystemsBioMedicine/Panorama BM 1-10/")
```

### Patient Data

``` r
pathToDir <- "~/Downloads/SystemsBioMedicine/systems_biomedicine_spill_applied_fcs_files/"
exp2 <- list.files(pathToDir, pattern = "*.fcs", full.names = T)
names <- list.files(pathToDir, pattern = "*.fcs")
exp2 <- exp2[c(1:6)]
names <- names[c(1:6)]
info <- sapply(names, function(x){strsplit(x, split = ".*\\d\\d\\d")[[1]][2]})
a_b <- sapply(info, function(x){factor(strsplit(x, split = "_")[[1]][1])})
dual_triple <- sapply(info, function(x){
  factor(
    strsplit(strsplit(x, split = "_")[[1]][2], split = "\\.")[[1]][1]
    )
  })
library(data.table)
metadata <- data.table("file_name" = names, 
                       "sample_id" = names, "a_b" = a_b, "dual_triple" = dual_triple)

sce2 <- prepData(exp2, md = metadata, md_cols = list(file = "file_name", 
                                                     id = "sample_id", 
                                                     factors = c("a_b", "dual_triple")))
```

``` r
set.seed(1601)
sce2 <- runDR(sce2, dr = "UMAP", features = "state", cells = 500)
sce2 <- runDR(sce2, dr = "TSNE", features = "state", cells = 500)
chs <- grep("^CD", rownames(sce2), value = TRUE)[c(1:5)]

plotDR(sce2, dr = "UMAP", color_by = "a_b", facet_by = "dual_triple")
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
plotDR(sce2, dr = "UMAP", color_by = chs, facet_by = "dual_triple")
```

![](README_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
plotDR(sce2, dr = "TSNE", color_by = "a_b", facet_by = "dual_triple")
```

![](README_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

``` r
plotDR(sce2, dr = "TSNE", color_by = chs, facet_by = "dual_triple")
```

![](README_files/figure-gfm/unnamed-chunk-3-4.png)<!-- -->
