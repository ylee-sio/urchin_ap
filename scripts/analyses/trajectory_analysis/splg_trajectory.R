library(tidyverse)
library(tictoc)
library(reticulate)
library(Seurat)
library(monocle3)

# Step 1: Normalizing and pre-processing the data
cds = load_mm_data(
  mat_path = "data/unmodified/sp8_lg/matrix.mtx.gz",
  feature_anno_path = "data/unmodified/sp8_lg/features.tsv.gz",
  cell_anno_path = "data/unmodified/sp8_lg/barcodes.tsv.gz",
  
)

cds_matrix = Matrix::readMM("data/unmodified/sp8_lg/matrix.mtx.gz")
cds_cell_metadata = read_tsv("data/unmodified/sp8_lg/barcodes.tsv.gz", col_names = FALSE) %>% as.data.frame()
cds_feature_metadata = read_tsv("data/unmodified/sp8_lg/features.tsv.gz", col_names = FALSE) %>% as.data.frame()

# removing duplicated features and setting rownames for feature file
cds_duplicates_index = cds_feature_metadata$X2 %>% duplicated() %>% which()
cds_feature_metadata_duplicates_removed = cds_feature_metadata[-cds_duplicates_index,]
rownames(cds_feature_metadata_duplicates_removed) = cds_feature_metadata_duplicates_removed$X2
colnames(cds_feature_metadata_duplicates_removed) = c("gene_num", "gene_short_name", "assay")

# setting rownames to cell metadata
rownames(cds_cell_metadata) = cds_cell_metadata$X1

# removing duplicate genes from expression matrix
cds_matrix = cds_matrix[-cds_duplicates_index,]

# setting rownames and colnames for expression matrix
rownames(cds_matrix) = cds_feature_metadata_duplicates_removed$X2
colnames(cds_matrix) = cds_cell_metadata$X1

cds = new_cell_data_set(
  cds_matrix,
  cell_metadata = cds_cell_metadata,
  gene_metadata = cds_feature_metadata_duplicates_removed
  )

cds = preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds = reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")

# refer to this post https://github.com/cole-trapnell-lab/monocle-release/issues/262
# all_stages = readRDS("data/working_data/all_stages.clustered_seurat_obj_list.rds")

cds = estimate_size_factors(cds)
cds = cluster_cells(cds, resolution=1e-4)
cds = learn_graph(cds)
plot_cells(cds, genes = c("foxq2", "gcm"))
