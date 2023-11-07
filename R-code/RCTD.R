library(spacexr)
library(SeuratDisk)
library(Matrix)
library(Giotto)

samples = c('CVD1', 'CVD2', 'CVD3', 'CVD4', 'CVD5',
            'CVD6', 'CVD7', 'CVD8', 'CVD9', 'CVD10')

#Name of annotation in seurat object
Anno_col <- "Final_annotation"

sc_obj <- LoadH5Seurat("athero_subset3000.h5seurat")

datadir <- "C:/Users/alban/Data/CVD"


sc_obj@meta.data[,Anno_col] <- as.factor(as.character(sc_obj@meta.data[,Anno_col]))

# Create table with reference data counts
counts <- data.frame(sc_obj@assays$RNA@counts)
colnames(counts) <- colnames(sc_obj)

# Extract annotations to vector
meta_data <- data.frame(sc_obj@meta.data[,Anno_col])
cell_types <- meta_data[,1]
names(cell_types) <- rownames(sc_obj@meta.data)
cell_types <- as.factor(cell_types)

# Sum counts to nUMI for each spot into another vector
nUMI_df <- data.frame(colSums(sc_obj@assays$RNA@counts))
nUMI <- nUMI_df$colSums.sc_obj.assays.RNA
names(nUMI) <- rownames(nUMI_df)


### Create the sc Reference object
reference <- Reference(counts, cell_types, nUMI)


for (sample in samples) {
  path = file.path(datadir, sample)
  print(path)
  
  h5matrix <- get10Xmatrix_h5(paste0(path, "/", sample, " filtered_feature_bc_matrix.h5"))
  cells <- h5matrix$`Gene Expression`@Dimnames[[2]]
  locs <- read.csv(paste0(path, "/spatial/", "tissue_positions_list.csv"), header = FALSE)
  locs <- locs[locs$V1 %in% cells,]
  
  locs <- subset(locs, select = -c(V2, V5, V6))
  
  locs <- locs[, c("V3", "V4", "V1")]
  
  counts <- h5matrix$`Gene Expression`
  counts <- as.data.frame(as.matrix(counts))
  #counts <- t(counts)
  
  coords <- locs[, -3]
  rownames(coords) <- locs[, 3]
  
  #Create spacexr spatial object
  puck <- SpatialRNA(coords, counts)
  
  barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 
  plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                       title ='plot of nUMI')
  
  myRCTD <- create.RCTD(puck, reference, max_cores = 8)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  results <- myRCTD@results
  
  resultsdir <- paste0('RCTD_results/', sample, '_results')
  dir.create(resultsdir)
  
  # normalize the cell type proportions to sum to 1.
  norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
  cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
  spatialRNA <- myRCTD@spatialRNA
  write.csv(norm_weights, paste0(resultsdir, '/norm_weights.txt'))
  
  
  # Uncomment to save plots of cell type weights
  
  #plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
  
  # Plots all weights for each cell type as in full_mode. (saved as 
  # 'results/cell_type_weights.pdf')
  
  
}
