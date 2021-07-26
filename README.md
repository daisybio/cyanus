# A Systematic Comparison of Differential Analysis Methods for CyTOF Data

Systematic Evaluation of Existing and Novel Approaches for Differential Expression Anaylsis on Real and Simulated CyTOF Data Including the Shiny App CYANUS



## Table of Contents
* [Dataset Simulation](#dataset-simulation)
* [Methods Evaluation](#methods-evaluation)
* [CYANUS Shiny App](#cyanus-shiny-app)
* [Cite](#cite)
* [Contact](#contact)


## Dataset Simulation
The PBMC dataset, the simulated CytoGLMM dataset (and its downsampled versions) and the semi-simulated COVID-19 dataset (and its downsampled versions) can be reproduced using the scripts in DataGeneration. The dual antiplatelet therapy dataset is available upon request from Melissa.Klug@tum.de.


## Methods Evaluation
The evaluation of the methods on the different datasets can be found in the DEComparison directory. Each dataset has a subdirectory containing the benchmarking scripts for running the methods, the resulting .rds objects and an evaluation.R script. Using the evaluation scripts, the plots from the paper can be reproduced.


## CYANUS Shiny App

User-friendly R Shiny App Analyzing Cytometry Data

The public web interace is available at https://exbio.wzw.tum.de/cyanus/.

The source code can be found in ui/, server/, functions/, server.R and ui.R. 

Additionally, a docker image for our shiny app is publicly available at [Docker
Hub](https://hub.docker.com/repository/docker/quirinmanz/cytof_pipeline).
First, make sure docker is installed on your machine. Then run

``` bash
docker pull quirinmanz/cytof_pipeline
```

If you want to use the [PBMC example
data](https://www.nature.com/articles/nbt.2317), you can download it
from our [Google
Drive](https://drive.google.com/drive/folders/19hM51eoLLEJDQ_Oz4xqMu2t9bAY9Qcyf?usp=sharing).
Then run

``` bash
mkdir data
mv -r path/to/pbmcData data
docker run --rm -p 3838:3838 -v /absolute/path/to/data:/srv/cytof_pipeline/data quirinmanz/cytof_pipeline
```

in order to connect the data folder with the downloaded data to the
repository and run the app. Make sure that all parent directories of
data have at least read access for group (chmod 711).

If you go to localhost:3838, you can see our Shiny app.

You can also run the Shiny App in your R session with `shiny::runApp()`.
In that case, you have to restore the project library using
[renv](https://rstudio.github.io/renv/articles/renv.html):

``` r
# install.packages('renv')
renv::restore()
```

## Cite
TODO

## Contact
If you have difficulties using CYANUS, please open an issue. If issues on the simulation of data or evaluation of methods arise, you can either open an issue or contact the authors of the paper.


