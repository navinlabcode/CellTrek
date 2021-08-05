# A Quick Tour of GSimp

# 1. Introduction

SChart is a computational framework that can directly map single cells back to their spatial coordinates in tissue sections based on scRNA-seq and ST data. This method provides a new paradigm that is distinct from ST deconvolution, enabling a more flexible and direct investigation of single cell data with spatial topography. The SChart toolkit also provides two downstream analysis modules, including SColoc for spatial colocalization analysis and SCoexp for spatial co-expression analysis.

In this tutorial, we will demonstrate the cell charting workflow based on the mouse brain data as part of our paper Figure 2

# 2. Loading the packages and datasets (scRNA-seq and ST data)
We start by loading the packages needed for the analyses.
``` r
options(stringsAsFactors = F)
library("dplyr")
library("magrittr")
library("dbscan")
library("pheatmap")
library("spatstat")
library("Seurat")
library("SeuratData")
library("reshape2")
library("visNetwork")
library("shiny")
library("plotly")
library("viridis")
```

We then load mouse brain scRNA-seq and ST data, respectively. For ST data, we only used the frontal cortex region for this study. 
``` r
brain_st_cortex <- readRDS("brain_st_cortex.rds")
brain_sc <- readRDS("brain_sc.rds")
## Visualize the ST data
SpatialDimPlot(brain_st_cortex)
```
![](vignette_files/F1_ST_spatialdimplot.png)

``` r
## Visualize the scRNA-seq data
DimPlot(brain_sc, label = T, label.size = 4.5)
```
![](vignette_files/F2_SC_dimplot.png)

# 3. Cell charting using SChart
We first co-embed ST and scRNA-seq datasets using *traint*
``` r
brain_traint <- SChartPack::traint(st_data=brain_st_cortex, sc_data=brain_sc, sc_assay='RNA', cell_names='cell_type')
```
``` r
## Finding transfer anchors... 
## Using 2000 features for integration... 
## Data transfering... 
## Creating new Seurat object... 
## Scaling -> PCA -> UMAP...
```
``` r
## We can check the co-embedding result to see if there is overlap between these two data modalities
DimPlot(brain_traint, group.by = "type") 
```
![](vignette_files/F3_scst_coembed.png)
After coembedding, we can chart single cells to their spatial locations. Here, we use the non-linear interpolation (intp = T, intp_lin=F) approach to augment the ST spots.
``` r
brain_schart <- SChartPack::schart(st_sc_int=brain_traint, int_assay='traint', sc_data=brain_sc, sc_assay = 'RNA', 
                                   reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000, 
                                   dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)$schart
```
``` r
## Distance between spots is: 138 
## Interpolating...
## Random Forest training... 
## Random Forest prediction...  
## Making distance matrix... 
## Making graph... 
## Pruning graph...
## Spatial Charting SC data...
## Repelling points...
## Creating Seurat Object... 
## sc data...
```

Mouse brain data tutorial: [Here](vignette_files/MouseBrainSChartExample.pdf)
