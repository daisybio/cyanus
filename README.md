# A Systematic Comparison of Differential Analysis Methods for CyTOF Data

Systematic Evaluation of Existing and Novel Approaches for Differential Expression Anaylsis on Real and Simulated CyTOF Data Including the Shiny App CYANUS



## Table of Contents
* [Dataset Simulation](#dataset-simulation)
* [Methods Evaluation](#methods-evaluation)
* [CYANUS Shiny App](#cyanus-shiny-app)
* [Cite](#cite)
* [Contact](#contact)
* [Licenses](#licenses)


## Dataset Simulation
The PBMC dataset, the simulated CytoGLMM dataset (and its downsampled versions) and the semi-simulated COVID-19 dataset (and its downsampled versions) can be reproduced using the scripts in DataGeneration. The dual antiplatelet therapy dataset is available upon request from Melissa.Klug@tum.de.


## Methods Evaluation
The evaluation of the methods on the different datasets can be found in the DEComparison directory. Each dataset has a subdirectory containing the benchmarking scripts for running the methods, the resulting .rds objects and an evaluation.R script. Using the evaluation scripts, the plots from the paper can be reproduced.


## CYANUS Shiny App

User-friendly R Shiny App Analyzing Cytometry Data

The public web interace is available at https://exbio.wzw.tum.de/cyanus/.

The source code can be found in ui/, server/, functions/, server.R and ui.R. 

Additionally, a docker image for our shiny app is publicly available on [Docker
Hub](https://hub.docker.com/repository/docker/quirinmanz/cyanus).
First, make sure docker is installed on your machine. Then run

``` bash
docker pull quirinmanz/cyanus:latest
```

You can run the app in a single container using

``` bash
docker run --rm -p 3838:3838 quirinmanz/cyanus:latest
```

If you go to localhost:3838, you can see our Shiny app.

You can also run the Shiny App in your R session with `shiny::runApp()`.
In that case, you have to restore the project library using
[renv](https://rstudio.github.io/renv/articles/renv.html):

``` r
# install.packages('renv')
renv::restore()
```

To use multiple instances of the app, we configured nginx and [shinyproxy](https://www.shinyproxy.io/) with docker.
The command to use this production-ready version of CYANUS is
``` bash
docker-compose up -d
```

## Cite
Lis Arend, Judith Bernett, Quirin Manz, Melissa Klug, Olga Lazareva, Jan Baumbach, Dario Bongiovanni, Markus List, A systematic comparison of novel and existing differential analysis methods for CyTOF data, Briefings in Bioinformatics, 2021;, bbab471, [https://doi.org/10.1093/bib/bbab471](https://doi.org/10.1093/bib/bbab471)

## Contact
If you have difficulties using CYANUS, please open an issue. If issues on the simulation of data or evaluation of methods arise, you can either open an issue or contact the authors of the paper.

## Licenses
|Package              |License                        |
|:--------------------|:------------------------------|
|CATALYST             |GPL (>=2)                      |
|HDCytoData           |MIT + file LICENSE             |
|readxl               |GPL-3                          |
|data.table           |MPL-2.0 &#124; file LICENSE    |
|ggplot2              |MIT + file LICENSE             |
|BiocParallel         |GPL-2 &#124; GPL-3             |
|fs                   |GPL-3                          |
|SummarizedExperiment |Artistic-2.0                   |
|dplyr                |MIT + file LICENSE             |
|reshape2             |MIT + file LICENSE             |
|ggside               |MIT + file LICENSE             |
|cowplot              |GPL-2                          |
|ggpubr               |GPL-2                          |
|ComplexHeatmap       |MIT + file LICENSE             |
|ConsensusClusterPlus |GPL version 2                  |
|flowCore             |Artistic-2.0                   |
|FlowSOM              |GPL (>= 2)                     |
|ggridges             |GPL-2 &#124; file LICENSE      |
|graphics             |Part of R 4.0.5                |
|grDevices            |Part of R 4.0.5                |
|grid                 |Part of R 4.0.5                |
|igraph               |GPL (>= 2)                     |
|purrr                |GPL-3 &#124; file LICENSE      |
|RColorBrewer         |Apache License 2.0             |
|S4Vectors            |Artistic-2.0                   |
|Rcpp                 |GPL (>= 2)                     |
|RcppAlgos            |GPL (>= 2)                     |
|CytoGLMM             |LGPL-3                         |
|MASS                 |GPL-2 &#124; GPL-3             |
|diffcyt              |MIT + file LICENSE             |
|rstatix              |GPL-2                          |
|multcomp             |GPL-2                          |
|gamlss               |GPL-2 &#124; GPL-3             |
|gamlss.dist          |GPL-2 &#124; GPL-3             |
|lme4                 |GPL (>= 2)                     |
|circlize             |MIT + file LICENSE             |
|magrittr             |MIT + file LICENSE             |
|scales               |MIT + file LICENSE             |
|dimRed               |GPL-3 &#124; file LICENSE      |
|ggvenn               |MIT + file LICENSE             |
|RANN                 |GPL (>= 3)                     |
|uwot                 |GPL (>= 3)                     |
|renv                 |MIT + file LICENSE             |
|shinydashboard       |GPL (>= 2) &#124; file LICENSE |
|shinyjs              |MIT + file LICENSE             |
|shinyWidgets         |GPL-3                          |
|stats                |Part of R 4.0.5                |
|DT                   |GPL-3 &#124; file LICENSE      |
|shiny                |GPL-3 &#124; file LICENSE      |
|shinycssloaders      |MIT + file LICENSE             |
|shinyBS              |GPL-3                          |
|sortable             |MIT + file LICENSE             |
|xlsx                 |GPL-3                          |
|plotly               |MIT + file LICENSE             |
|waiter               |MIT + file LICENSE             |



