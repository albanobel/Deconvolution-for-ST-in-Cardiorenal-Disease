library(Giotto)
library(SeuratDisk)

my_python_path = "C:/Users/alban/miniconda3/python.exe"

output_path = "./DWLS_output"

instrs = createGiottoInstructions(python_path = my_python_path)
sc <- LoadH5Seurat("athero_subset_4000.h5seurat")

datadir <- "C:/Users/alban/Data/CVD"

#Initialize giotto single cell object
sc_data <- createGiottoObject(
    raw_exprs = sc@assays$RNA@counts,
    instructions = instrs
)

sc_data <- normalizeGiotto(gobject = sc_data)
sc_data <- calculateHVG(gobject = sc_data)
gene_metadata = fDataDT(sc_data)
featgenes = gene_metadata[hvg == 'yes']$gene_ID
sc_data <- runPCA(gobject = sc_data, genes_to_use = featgenes, scale_unit = F)
signPCA(sc_data, genes_to_use = featgenes, scale_unit = F)
sc_data@cell_metadata$leiden_clus <- as.character(sc@meta.data[,'Final_annotation'])
scran_markers_subclusters = findMarkers_one_vs_all(gobject = sc_data,
                                                   method = 'scran',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_clus')

Sig_scran <- unique(scran_markers_subclusters$genes[which(scran_markers_subclusters$ranking <= 100)])
norm_exp<-2^(sc_data@norm_expr)-1
id<-sc_data@cell_metadata$leiden_clus
ExprSubset<-norm_exp[Sig_scran,]
Sig_exp<-NULL

for (i in unique(id)){
  Sig_exp<-cbind(Sig_exp,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
}

colnames(Sig_exp)<-unique(id)




samples = c('CVD1', 'CVD2', 'CVD3', 'CVD4', 'CVD5',
            'CVD6', 'CVD7', 'CVD8', 'CVD9', 'CVD10')

for (sample in samples) {
  
  path = file.path(datadir, sample)
  print(path)
  
  resultsdir <- paste0('DWLS_fullRun2/', sample, '_results')
  dir.create(resultsdir)
  
  instrs = createGiottoInstructions(save_plot = TRUE,
                                    show_plot = TRUE,
                                    return_plot = FALSE,
                                    save_dir = resultsdir,
                                    python_path = my_python_path)
  
  
  h5matrix <- get10Xmatrix_h5(paste0(path, "/", sample, " filtered_feature_bc_matrix.h5"))
  cells <- h5matrix$`Gene Expression`@Dimnames[[2]]
  locs <- read.csv(paste0(path, "/spatial/", "tissue_positions_list.csv"), header = FALSE)
  locs <- locs[locs$V1 %in% cells,]
  
  locs <- subset(locs, select = -c(V2, V5, V6))
  
  locs <- locs[, c("V3", "V4", "V1")]
  
  
  
  st_data <- createGiottoObject(
    raw_exprs = h5matrix$`Gene Expression`,
    spatial_locs = locs,
    instructions = instrs
  )
  
  st_data <- filterGiotto(gobject = st_data,
                          expression_threshold = 1,
                          gene_det_in_min_cells = 10,
                          min_det_genes_per_cell = 50,
                          expression_values = c('raw'),
                          verbose = T)
  
  st_data <- normalizeGiotto(gobject = st_data)
  st_data <- calculateHVG(gobject = st_data,
                          method = "cov_loess")
  gene_metadata = fDataDT(st_data)
  featgenes = gene_metadata[hvg == 'yes']$gene_ID
  st_data <- runPCA(gobject = st_data, genes_to_use = featgenes, scale_unit = F)
  signPCA(st_data, genes_to_use = featgenes, scale_unit = F)
  st_data <- runUMAP(st_data, dimensions_to_use = 1:10)
  st_data <- createNearestNetwork(gobject = st_data, dimensions_to_use = 1:10, k = 15)
  st_data <- doLeidenCluster(gobject = st_data, resolution = 0.4, n_iterations = 1000)
  
  
  
  st_data <- runDWLSDeconv(st_data,sign_matrix = Sig_exp, n_cell = 15)
  write.csv(st_data@spatial_enrichment$DWLS, paste0(resultsdir, '/SpatialDWLS_result.txt'))
}