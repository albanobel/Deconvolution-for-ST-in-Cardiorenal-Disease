<img src="Figures/CVD2_macrophages.png" width="800" align="center" /> <br>


# Single Cell Deconvolution of Cardiorenal Spatial Transcriptomics Data

High-throughput RNA-sequencing technologies that provide spatial resolution of transcripts, popularly known as spatial transcriptomics, are on the rise. This technology seems highly promising and has an untapped potential for expression-driven discovery in development and disease. Nonetheless, it also faces the central challenge of mixed cell type signals due to limitations in resolution. This is apparent in sequencing-based 10X Visium where slides have larger spots of 55 um. This mixed transcriptional signal can pose inferential problems; however, it can theoretically be deconvoluted into underlying cell types. To this end, we developed a systematic deconvolution framework and performed benchmarking in previously unvalidated healthy and disease samples from human coronary arterial and kidney disease. We used: Cell2location, RCTD and spatialDWLS that have previously been shown to perform well in mouse brain and simulated data (1). We show that all three methods are capable of deconvoluting verifiable cell types when benchmarked against expert provided ground truth based on accuracy scores (0.7-0.73). Kidney podocyte cells and major populations of macrophages, smooth muscle cells and fibroblasts in arteries are all deconvoluted with a high level of agreement. Bayesian Cell2location is more computationally demanding, however it provides quality solutions, when less reference data is available. 

Authors: Alban Obel Slabowska, Charles Pyke, Henning Hvid, Leon Eyrich Jessen, Simon Baumgart, Vivek Das

## Data

Single cell RNA reference data were obtained from public sources.


Atherosclerosis scRNA from [Wirka et al. 2019](https://doi.org/10.1038/s41591-019-0512-5), [Pan et al. 2020](https://pubmed.ncbi.nlm.nih.gov/32962412/), and [Alsaigh et al. 2020](https://doi.org/10.1038/s42003-022-04056-7). All of which can be obtained from the [PlaqView Portal](https://www.plaqview.com/).

Kidney scRNA data was obtained from The Kidney Precision Medicine Project [KPMP](https://www.kpmp.org/available-data).

Spatial transcriptomics data was generated internally at Novo Nordisk using 10X Visium technology for FFPE samples.


## Poster

<img src="Figures/Poster.png" width="1000" align="center" /> <br>
