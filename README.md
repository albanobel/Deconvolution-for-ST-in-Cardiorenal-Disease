<img src="Figures/CVD2_macrophages.png" width="800" align="center" /> <br>


# Single Cell Deconvolution of Cardiorenal Spatial Transcriptomics Data

High-throughput RNA-sequencing technologies that provide spatial resolution of transcripts, popularly known as spatial transcriptomics, are on the rise. This technology seems highly promising and has an untapped potential for expression-driven discovery in development and disease. Nonetheless, it also faces the central challenge of mixed cell type signals due to limitations in resolution. This is apparent in sequencing-based 10X Visium where slides have larger spots of 55 um. This mixed transcriptional signal can pose inferential problems; however, it can theoretically be deconvoluted into underlying cell types. To this end, we developed a systematic deconvolution framework and performed benchmarking in previously unvalidated healthy and disease samples from human coronary arterial and kidney disease. We used: Cell2location, RCTD and spatialDWLS that have previously been shown to perform well in mouse brain and simulated data (1). We show that all three methods are capable of deconvoluting verifiable cell types when benchmarked against expert provided ground truth based on accuracy scores (0.7-0.73). Kidney podocyte cells and major populations of macrophages, smooth muscle cells and fibroblasts in arteries are all deconvoluted with a high level of agreement. Bayesian Cell2location is more computationally demanding, however it provides quality solutions, when less reference data is available. 

Authors: Alban Obel Slabowska, Charles Pyke, Henning Hvid, Leon Eyrich Jessen, Simon Baumgart, Vivek Das

## Data

Single cell RNA reference data were obtained from public sources.


Atherosclerosis scRNA from [Wirka et al. 2019](https://doi.org/10.1038/s41591-019-0512-5), [Pan et al. 2020](https://pubmed.ncbi.nlm.nih.gov/32962412/), and [Alsaigh et al. 2020](https://doi.org/10.1038/s42003-022-04056-7). All of which can be obtained from the [PlaqView Portal](https://www.plaqview.com/).

Kidney scRNA data was obtained from The Kidney Precision Medicine Project [KPMP](https://www.kpmp.org/available-data).

Spatial transcriptomics data was generated internally at Novo Nordisk using 10X Visium technology for FFPE samples.

## Software and Packages
RCTD was run using R 4.1.2 and the spacexr package (2.0.1). SpatialDWLS was run
using the Giotto package (1.1.2) and R 4.1.2 with python 3.10 dependencies.
Cell2location was run using the cell2location package (0.1) and python 3.10 and
libraries pyro-ppl 1.8.0; scvi-tools 1.0.0; torch 1.9.0; numpy 1.23.4; scanpy 1.9.1;
anndata 0.8.0.
Packages Seurat 4.3.0, SeuratDisk 0.0.0.9020, SeuratData 0.2.2, tidyverse 1.3.2,
ggplot 3.4.0 were used for data wrangling and visualizations.


## Poster presented at ISMB-ECCB 2023 in Lyon

<img src="Figures/Poster.png" width="1000" align="center" /> <br>
