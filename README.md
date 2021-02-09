CyTOF Pipeline
================

**This README has two parts: The Shiny App and an example analysis**

# Part 1: Shiny App

The docker image for our shiny app is publicly available at \[…\].
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
mkdir myData
mv -r path/to/pbmcData myData
docker run --rm -p 3838:3838 -v /absolute/path/to/myData:/srv/cytof_pipeline/data cytof_pipeline
```

in order to connect the data folder with the downloaded data to the
repository.

If you go to localhost:3838, you can see our Shiny app.

You can also run the Shiny App in your R session with `shiny::runApp()`.
In that case, you have to restore the project library using
[renv](https://rstudio.github.io/renv/articles/renv.html):

``` r
# install.packages('renv')
renv::restore()
```

# Part 2: Example analysis

Data from human platelets of patients with chronic coronary syndrome
undergoing different therapy: dual antiplatelet therapy versus triple
antiplatelet therapy, before and after platelet activation with 10µm
TRAP. There are files of 7 patients with triple therapy and 12 patients
with dual therapy (each in two conditions). For information about how to
create a SCE object, please refer to the [CATALYST
Vignette](https://www.bioconductor.org/packages/release/bioc/html/CATALYST.html)

# Read in SCE object

This SCE object contains 10 000 cells per sample and therefore 20 000
cells per patient.

``` r
sce <- readRDS("data/platelets_small/sce.rds")
metadata(sce)$experiment_info
```

    ##          sample_id activated_baseline dual_triple patient_id n_cells
    ## 1    RPS 096A_dual                  A        dual    RPS 096   10000
    ## 2    RPS 096B_dual                  B        dual    RPS 096   10000
    ## 3    RPS 097A_dual                  A        dual    RPS 097   10000
    ## 4    RPS 097B_dual                  B        dual    RPS 097   10000
    ## 5  RPS 098A_triple                  A      triple    RPS 098   10000
    ## 6  RPS 098B_triple                  B      triple    RPS 098   10000
    ## 7  RPS 099A_triple                  A      triple    RPS 099   10000
    ## 8  RPS 099B_triple                  B      triple    RPS 099   10000
    ## 9    RPS 101A_dual                  A        dual    RPS 101   10000
    ## 10   RPS 101B_dual                  B        dual    RPS 101   10000
    ## 11   RPS 103A_dual                  A        dual    RPS 103   10000
    ## 12   RPS 103B_dual                  B        dual    RPS 103   10000
    ## 13 RPS 111A_triple                  A      triple    RPS 111   10000
    ## 14 RPS 111B_triple                  B      triple    RPS 111   10000
    ## 15 RPS 127A_triple                  A      triple    RPS 127   10000
    ## 16 RPS 127B_triple                  B      triple    RPS 127   10000
    ## 17 RPS 136A_triple                  A      triple    RPS 136   10000
    ## 18 RPS 136B_triple                  B      triple    RPS 136   10000
    ## 19   RPS 137A_dual                  A        dual    RPS 137   10000
    ## 20   RPS 137B_dual                  B        dual    RPS 137   10000
    ## 21 RPS 138A_triple                  A      triple    RPS 138   10000
    ## 22 RPS 138B_triple                  B      triple    RPS 138   10000
    ## 23   RPS 139A_dual                  A        dual    RPS 139   10000
    ## 24   RPS 139B_dual                  B        dual    RPS 139   10000
    ## 25   RPS 144A_dual                  A        dual    RPS 144   10000
    ## 26   RPS 144B_dual                  B        dual    RPS 144   10000
    ## 27   RPS 146A_dual                  A        dual    RPS 146   10000
    ## 28   RPS 146B_dual                  B        dual    RPS 146   10000
    ## 29   RPS 150A_dual                  A        dual    RPS 150   10000
    ## 30   RPS 150B_dual                  B        dual    RPS 150   10000
    ## 31   RPS 154A_dual                  A        dual    RPS 154   10000
    ## 32   RPS 154B_dual                  B        dual    RPS 154   10000

# Quality Control

``` r
CATALYST::pbMDS(
      sce,
      label_by = "patient_id", 
      color_by = "activated_baseline",
      features = NULL, # possible inputs: unique(rowData(sce)$marker_class) or NULL (to use all features)
      assay = "exprs", # possible inputs: assayNames(sce)
    ) + theme(text = element_text(size=18))
```

![](README_files/figure-gfm/diagnostic%20plots-1.png)<!-- -->

``` r
CATALYST::plotNRS(
      sce,
      color_by = "activated_baseline",
      features = NULL,
      assay = "exprs"
    )
```

![](README_files/figure-gfm/diagnostic%20plots-2.png)<!-- -->

``` r
CATALYST::plotExprs(
      sce,
      color_by = "activated_baseline",
      features = NULL,
      assay = "exprs"
    )  + theme(text = element_text(size=18))
```

![](README_files/figure-gfm/diagnostic%20plots-3.png)<!-- -->

``` r
plotExprHeatmap(
      sce,
      features = NULL,
      assay = "exprs"
    )
```

![](README_files/figure-gfm/diagnostic%20plots-4.png)<!-- -->

``` r
CATALYST::plotCounts(
    sce,
    group_by = "patient_id", # possible inputs: names(colData(sce))
    color_by = "dual_triple", # possible inputs: names(colData(sce))
    prop = FALSE # TRUE for stacked relative abundances, FALSE for total cell counts
  ) 
```

![](README_files/figure-gfm/diagnostic%20plots-5.png)<!-- --> We can see
that sample 99 does not cluster with the other activated samples in the
MDS plot. We can also see in the expression plot, that the expressions
for 99 baseline and 99 activated are nearly the same. We therefore
exclude the patient from our analysis.

``` r
sce <- makePatientSelection(sce = sce, deselected_patients = c("RPS 099"))
```

The next step of quality control is clustering:

``` r
sce <- clusterSCE(sce)
```

    ## o running FlowSOM clustering...

    ## Building MST

    ## o running ConsensusClusterPlus metaclustering...

``` r
plotly::ggplotly(delta_area(sce))
```

    ## PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.

<div id="htmlwidget-90a947d9940ac62f7bf7" style="width:1152px;height:768px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-90a947d9940ac62f7bf7">{"x":{"data":[{"x":[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],"y":[0.0419737373737374,3.51860230062088,0.882406825516867,0.259758043525005,0.204193467246397,0.151220476552825,0.0995259742153302,0.0558096360689553,0.0424768831864781,0.0335441275688873,0.0311343213215314,0.0301565608859054,0.0276764375944985,0.0185520480103652,0.0149198442782025,0.0111998457410056,0.012183605366507,0.0113282362624669,0.0102943263669588],"text":["k:  2<br />y: 0.04197374","k:  3<br />y: 3.51860230","k:  4<br />y: 0.88240683","k:  5<br />y: 0.25975804","k:  6<br />y: 0.20419347","k:  7<br />y: 0.15122048","k:  8<br />y: 0.09952597","k:  9<br />y: 0.05580964","k: 10<br />y: 0.04247688","k: 11<br />y: 0.03354413","k: 12<br />y: 0.03113432","k: 13<br />y: 0.03015656","k: 14<br />y: 0.02767644","k: 15<br />y: 0.01855205","k: 16<br />y: 0.01491984","k: 17<br />y: 0.01119985","k: 18<br />y: 0.01218361","k: 19<br />y: 0.01132824","k: 20<br />y: 0.01029433"],"type":"scatter","mode":"lines+markers","line":{"width":1.88976377952756,"color":"rgba(70,130,180,1)","dash":"dash"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","marker":{"autocolorscale":false,"color":"rgba(0,0,128,1)","opacity":1,"size":9.4488188976378,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,0,128,1)"}},"frame":null}],"layout":{"margin":{"t":26.958904109589,"r":7.30593607305936,"b":40.9132420091324,"l":43.1050228310502},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[1.5,20.5],"tickmode":"array","ticktext":["2","4","6","8","10","12","14","16","18","20"],"tickvals":[2,4,6,8,10,12,14,16,18,20],"categoryorder":"array","categoryarray":["2","4","6","8","10","12","14","16","18","20"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":true,"linecolor":"rgba(0,0,0,1)","linewidth":0.66417600664176,"showgrid":true,"gridcolor":"rgba(190,190,190,1)","gridwidth":0.265670402656704,"zeroline":false,"anchor":"y","title":{"text":"k","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"scaleanchor":"y","scaleratio":0.25,"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.125,4.125],"tickmode":"array","ticktext":["0.0","0.5","1.0","1.5","2.0","2.5","3.0","3.5","4.0"],"tickvals":[0,0.5,1,1.5,2,2.5,3,3.5,4],"categoryorder":"array","categoryarray":["0.0","0.5","1.0","1.5","2.0","2.5","3.0","3.5","4.0"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":true,"linecolor":"rgba(0,0,0,1)","linewidth":0.66417600664176,"showgrid":true,"gridcolor":"rgba(190,190,190,1)","gridwidth":0.265670402656704,"zeroline":false,"anchor":"x","title":{"text":"Relative change<br />in area under CDF curve","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"scaleanchor":"x","scaleratio":4,"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":false,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.88976377952756,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","showSendToCloud":false},"source":"A","attrs":{"1f8fa966873":{"x":{},"y":{},"type":"scatter"},"1f8f18a35b73":{"x":{},"y":{}}},"cur_data":"1f8fa966873","visdat":{"1f8fa966873":["function (y) ","x"],"1f8f18a35b73":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>

Let’s have a look at meta 7:

``` r
CATALYST::plotAbundances(sce, "meta7", group_by = "activated_baseline") +  theme(text = element_text(size=18), axis.text.x = element_text(size=12)) 
```

![](README_files/figure-gfm/plot%20abundances%20QC-1.png)<!-- --> As we
can see, clusters 5,6, and 7 were mainly made because of sample 96.
Therefore, we also exclude this sample:

``` r
sce <- makePatientSelection(sce = sce, deselected_patients = c("RPS 096"))
```

We cluster again in order to see if the delta area changed and if the
cluster abundances now look good:

``` r
sce <- clusterSCE(sce)
```

    ## o running FlowSOM clustering...

    ## Building MST

    ## o running ConsensusClusterPlus metaclustering...

``` r
plotly::ggplotly(delta_area(sce))
```

<div id="htmlwidget-58621f2315e47bb9171c" style="width:1152px;height:768px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-58621f2315e47bb9171c">{"x":{"data":[{"x":[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],"y":[0.0783191919191919,3.38062319438712,0.262182914478178,0.29577707903748,0.112091591510504,0.0891906153990593,0.0809478092611141,0.0621932340316714,0.0409225169626007,0.0333151786002951,0.021398572786676,0.0159642401021712,0.0149642476940168,0.0129315171167219,0.00950432240249782,0.00759802004536457,0.00707941049269778,0.00718160135835995,0.00640882673208872],"text":["k:  2<br />y: 0.078319192","k:  3<br />y: 3.380623194","k:  4<br />y: 0.262182914","k:  5<br />y: 0.295777079","k:  6<br />y: 0.112091592","k:  7<br />y: 0.089190615","k:  8<br />y: 0.080947809","k:  9<br />y: 0.062193234","k: 10<br />y: 0.040922517","k: 11<br />y: 0.033315179","k: 12<br />y: 0.021398573","k: 13<br />y: 0.015964240","k: 14<br />y: 0.014964248","k: 15<br />y: 0.012931517","k: 16<br />y: 0.009504322","k: 17<br />y: 0.007598020","k: 18<br />y: 0.007079410","k: 19<br />y: 0.007181601","k: 20<br />y: 0.006408827"],"type":"scatter","mode":"lines+markers","line":{"width":1.88976377952756,"color":"rgba(70,130,180,1)","dash":"dash"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","marker":{"autocolorscale":false,"color":"rgba(0,0,128,1)","opacity":1,"size":9.4488188976378,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,0,128,1)"}},"frame":null}],"layout":{"margin":{"t":26.958904109589,"r":7.30593607305936,"b":40.9132420091324,"l":43.1050228310502},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[1.5,20.5],"tickmode":"array","ticktext":["2","4","6","8","10","12","14","16","18","20"],"tickvals":[2,4,6,8,10,12,14,16,18,20],"categoryorder":"array","categoryarray":["2","4","6","8","10","12","14","16","18","20"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":true,"linecolor":"rgba(0,0,0,1)","linewidth":0.66417600664176,"showgrid":true,"gridcolor":"rgba(190,190,190,1)","gridwidth":0.265670402656704,"zeroline":false,"anchor":"y","title":{"text":"k","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"scaleanchor":"y","scaleratio":0.25,"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.125,3.625],"tickmode":"array","ticktext":["0.0","0.5","1.0","1.5","2.0","2.5","3.0","3.5"],"tickvals":[0,0.5,1,1.5,2,2.5,3,3.5],"categoryorder":"array","categoryarray":["0.0","0.5","1.0","1.5","2.0","2.5","3.0","3.5"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":true,"linecolor":"rgba(0,0,0,1)","linewidth":0.66417600664176,"showgrid":true,"gridcolor":"rgba(190,190,190,1)","gridwidth":0.265670402656704,"zeroline":false,"anchor":"x","title":{"text":"Relative change<br />in area under CDF curve","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"scaleanchor":"x","scaleratio":4,"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":false,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.88976377952756,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","showSendToCloud":false},"source":"A","attrs":{"1f8f30e49eb9":{"x":{},"y":{},"type":"scatter"},"1f8f71709cc2":{"x":{},"y":{}}},"cur_data":"1f8f30e49eb9","visdat":{"1f8f30e49eb9":["function (y) ","x"],"1f8f71709cc2":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>

We do not seem to win much by looking at more than 7 clusters. Let’s
take a look at meta7 again. We can also make a star plot to compare the
marker abundances:

``` r
CATALYST::plotAbundances(sce, "meta7", group_by = "activated_baseline") +  theme(text = element_text(size=18), axis.text.x = element_text(size=12)) 
```

![](README_files/figure-gfm/clustering%20plots-1.png)<!-- -->

``` r
invisible(plotStarsCustom(sce, backgroundValues = cluster_codes(sce)[["meta7"]]))
```

![](README_files/figure-gfm/clustering%20plots-2.png)<!-- --> The
clustering now looks fine. We can therefore continue to do our analysis:
\# Data Visualization You can now run various dimensionality reduction
methods on your data. The options are (parameter dr\_chosen): “UMAP”,
TSNE“,”PCA“,”MDS“,”DiffusionMap" and “Isomap”.

Other run options are: cells\_chosen: number of cells to sample from

feature\_chosen: markers or marker class. We recommend using “type”

assay\_chosen: “counts” or “exprs”. We recommend using “exprs”

scale: TRUE or FALSE. We recommend using TRUE

k: Isomap-specific parameter. Specifies the number of neighbours for the
graph constructed by isomap

We try UMAP and PCA on our data and color by the expression of our state
markers:

``` r
DR_methods <- c("UMAP","PCA")

for( method in DR_methods ){
  sce <- runDimRed(sce, method, cells_chosen = 1000, feature_chosen = "type", assay_chosen = "exprs", scale = T)
}
```

    ## Warning in check_numbers(k = k, nu = nu, nv = nv, limit = min(dim(x)) - : more
    ## singular values/vectors requested than available

    ## Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth =
    ## TRUE, : You're computing too large a percentage of total singular values, use a
    ## standard svd instead.

    ## Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth
    ## = TRUE, : did not converge--results might be invalid!; try increasing work or
    ## maxit

``` r
state_markers <- rowData(sce)[rowData(sce)$marker_class=="state",]$marker_name
for( method in DR_methods ){
  g <- CATALYST::plotDR(sce, dr = method, color_by = state_markers, facet_by="activated_baseline") +  theme(text = element_text(size=18))
  print(g)
}
```

![](README_files/figure-gfm/dimensionality%20reduction-1.png)<!-- -->![](README_files/figure-gfm/dimensionality%20reduction-2.png)<!-- -->

After clustering and visualization, we can move to the differential
expression analysis.

# Differential Expression Analysis

We are especially interested in the markers that are differentially
expressed between baseline and activation. Because we have two treatment
groups, we make the differential expression analysis in these subgroups
of patients:

## Dual Activated vs. Dual Baseline

First, we have to filter the data to only dual patients.

``` r
sce_d <- filterSCE(sce, dual_triple == "dual")
```

Before performing the differential marker expression analysis, we plot
the median marker expressions.

``` r
CATALYST::plotPbExprs(sce_d, k = "all", features = NULL , color_by = "activated_baseline",  ncol=8, facet_by = "antigen") +  theme(text = element_text(size=16))
```

![](README_files/figure-gfm/dual%20expressions-1.png)<!-- -->

From this plot, we can see that our state markers CD63 and CD62P are
probably differentially expressed but the state marker CD107a has a
median marker expression of zero and the state marker CD154 just has one
outlier sample that does not have a median marker expression of zero.

We now perform out differential expression analysis with the LMM. We
include patient\_id as random effect such that the LMM can compare the
marker expressions sample-wise

``` r
res_d <- runDS(sce = sce_d,
      condition = "activated_baseline",
      de_methods = c("LMM"),
      k = "all",
      features = "all",
      markers_to_test = "all",
      random_effect = "patient_id"
      )
```

    ## Using LMM

    ## Warning in any(lapply(contrastVars, function(y) {: wandle Argument des Typs
    ## 'list' nach boolesch
    
    ## Warning in any(lapply(contrastVars, function(y) {: wandle Argument des Typs
    ## 'list' nach boolesch

    ##      [,1]
    ## [1,]    0
    ## [2,]    1

    ## using SingleCellExperiment object from CATALYST as input

    ## using cluster IDs from clustering stored in column 'all' of 'cluster_codes' data frame in 'metadata' of SingleCellExperiment object from CATALYST

    ## calculating features...

    ## calculating DS tests using method 'diffcyt-DS-LMM'...

    ## boundary (singular) fit: see ?isSingular

    ## Warning in vcov.merMod(model): Computed variance-covariance matrix problem: not a positive definite matrix;
    ## returning NA matrix

    ## Error in asMethod(object) : not a positive definite matrix

    ## boundary (singular) fit: see ?isSingular
    ## boundary (singular) fit: see ?isSingular

    ## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, :
    ## Model failed to converge with max|grad| = 0.00284889 (tol = 0.002, component 1)

Now, we can visualize that with a heatmap:

``` r
CATALYST::plotDiffHeatmap(sce_d, rowData(res_d$LMM$res), all=T, col_anno = c("activated_baseline", "patient_id"), normalize = TRUE)
```

![](README_files/figure-gfm/heatmap%20dual-1.png)<!-- --> As we
expected, CD62P and CD63 are differentially expressed. Additionally, the
LMM found the type markers PEAR, CD69, PAR1 and CD42a to be
differentially expressed as well.

## Triple Activated vs. Triple Baseline

We now subset our SCE object such that only the triple anticoagulation
therapy patients are included:

``` r
sce_t <- filterSCE(sce, dual_triple == "triple")
```

``` r
CATALYST::plotPbExprs(sce_t, k = "all", features = NULL , color_by = "activated_baseline",  ncol=8, facet_by = "antigen") +  theme(text = element_text(size=16))
```

![](README_files/figure-gfm/triple%20expressions-1.png)<!-- --> Again,
we expect CD62P and CD63 to be differentially expressed. CD107a and
CD154 again have a median marker expression in every sample except for
one.

``` r
res_t <- runDS(sce = sce_t,
      condition = "activated_baseline",
      de_methods = c("LMM"),
      k = "all",
      features = "all",
      random_effect = "patient_id"
      )
```

    ## Using LMM

    ## Warning in any(lapply(contrastVars, function(y) {: wandle Argument des Typs
    ## 'list' nach boolesch
    
    ## Warning in any(lapply(contrastVars, function(y) {: wandle Argument des Typs
    ## 'list' nach boolesch

    ##      [,1]
    ## [1,]    0
    ## [2,]    1

    ## using SingleCellExperiment object from CATALYST as input

    ## using cluster IDs from clustering stored in column 'all' of 'cluster_codes' data frame in 'metadata' of SingleCellExperiment object from CATALYST

    ## calculating features...

    ## calculating DS tests using method 'diffcyt-DS-LMM'...

    ## boundary (singular) fit: see ?isSingular
    ## boundary (singular) fit: see ?isSingular

``` r
CATALYST::plotDiffHeatmap(sce_t, rowData(res_t$LMM$res), all=T, col_anno = c("activated_baseline", "patient_id"), normalize = TRUE)
```

![](README_files/figure-gfm/heatmap%20triple-1.png)<!-- --> Here, we
also find two type markers to be differentially expressed: PAR1 and
CD42a.

\#Comparison of limma and LMM for Triple A vs. B

Lastly, let’s compare which markers are found by limma and LMM:

``` r
res_triple <- runDS(sce = sce_t,
      condition = "activated_baseline",
      random_effect = "patient_id",
      de_methods = c("limma","LMM"),
      k = "all",
      features = "all")
```

    ## Using limma

    ## Warning in any(lapply(contrastVars, function(y) {: wandle Argument des Typs
    ## 'list' nach boolesch
    
    ## Warning in any(lapply(contrastVars, function(y) {: wandle Argument des Typs
    ## 'list' nach boolesch

    ##      [,1]
    ## [1,]    0
    ## [2,]    1

    ## using SingleCellExperiment object from CATALYST as input

    ## using cluster IDs from clustering stored in column 'all' of 'cluster_codes' data frame in 'metadata' of SingleCellExperiment object from CATALYST

    ## calculating features...

    ## calculating DS tests using method 'diffcyt-DS-limma'...

    ## Using LMM

    ## Warning in any(lapply(contrastVars, function(y) {: wandle Argument des Typs
    ## 'list' nach boolesch
    
    ## Warning in any(lapply(contrastVars, function(y) {: wandle Argument des Typs
    ## 'list' nach boolesch

    ##      [,1]
    ## [1,]    0
    ## [2,]    1

    ## using SingleCellExperiment object from CATALYST as input

    ## using cluster IDs from clustering stored in column 'all' of 'cluster_codes' data frame in 'metadata' of SingleCellExperiment object from CATALYST

    ## calculating features...

    ## calculating DS tests using method 'diffcyt-DS-LMM'...

    ## boundary (singular) fit: see ?isSingular
    ## boundary (singular) fit: see ?isSingular

``` r
# significant markers found by all methods:
createVennDiagram(res_triple)
```

![](README_files/figure-gfm/venn%20triple-1.png)<!-- --> limma only
finds the two state markers. When we look into the heatmap from the LMM
and the boxplots of the median marker expression, we can see why: The
difference between the median marker expressions does not look so
different overall in the boxplots. This is the only thing that limma can
evaluate. Because we tell the LMM to take the patient\_id into account
as an random effect, it compares the median marker expressions
patient-wise. In the heatmap, the patient samples are shown next to each
other. We can see that PAR1 and CD42a are downregulated for the
activated platelets in each patient but the expressions per patient
vary. This is why limma cannot find these two markers
