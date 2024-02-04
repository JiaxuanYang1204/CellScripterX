# CellScripter -- the automated analysis workflow for single-cell transcriptome data

This automated analysis package is a summary of my learning experience for single-cell transcriptome data. It contains a large number of built-in functions for analysis, allowing for streamlined completion of the entire single-cell transcriptome analysis in just 1-2 lines of commands. Additionally, it includes visually appealing visualization plots and can generate reports summarizing key findings. This current version, version 1, is implemented entirely in the R programming language. I will continue to update and add more analysis applications in future versions. Your feedback and suggestions are welcome and appreciated.

## Dependencies

autoscRNA requires the following R packages installed: 
* dplyr (1.1.3) 
* patchwork (1.1.3) 
* reshape2 (1.4.4)
* ggplot2 (3.4.4)
* ggcorrplot (0.1.4.1)
* corrplot (0.92)
* pheatmap (1.0.12)
* cols4all (0.6)
* Seurat (<= 4.3.0)(Ver5.0.0 is not workable for this pipeline)
* decontX (1.0.0)
* DoubletFinder (2.0.3)
* AnnotationHub (3.10.0)
* org.Hs.eg.db (3.18.0)
* clusterProfiler (4.10.0)
* Rgraphviz (2.46.0)
* SingleR (2.4.0)
* NOTE: These package versions were used in the my workflow, but other versions may also work as well.

If you can load all the following packages successfully, you are all good to drive autoscRNA!
``` r
## basic packages:
library(dplyr)
library(patchwork)
library(reshape2)
## visualizing packages:
library(ggcorrplot)
library(ggplot2)
library(pheatmap)
library(cols4all)
## Seurat based
library(Seurat)
## contamination predicting
library(decontX)
## Doublet predicting
library(DoubletFinder)
## pathway analysis
options(connectionObserver = NULL)
library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)
library(Rgraphviz)
## celltype predicting
library(SingleR)
```


## Installation (in R/RStudio)

autoscRNA is available on GitHub, please run the following command:
``` r
remotes::install_github('JiaxuanYang1204/CellScripterX/CellScripterX')
```
