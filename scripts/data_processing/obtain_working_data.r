tic()
# Run this script in terminal using 'Rscript obtain_working_data.R'
# This script only needs to be run one-time to obtain working data that has been processed with the parameters below. Adjust accordingly.

# setting workers for parallel operations
future::plan("multicore", workers = 8)

# make a list of paths for each directory (developmental stage specific) containing the matrix, feature, and barcode files 
indv_stage_dirs = list.files("data/unmodified", full.names = T)

# creates a list of matrices (developmental stage specific) of expression data to be used to create seurat_objects

all_stages.data = future_map(
  .x = indv_stage_dirs, 
  .f = function(x)(
    Read10X(
      data.dir = x,
    )
  )
)

dir.create("~/Projects/urchin_ap/data/working_data")

saveRDS(
  object = all_stages.data, 
  file = "~/Projects/urchin_ap/data/working_data/all_stages.data.rds"
)
# creates a list of seurat_objects for each stage using list of scrna matrices from all_stages.data
all_stages.seurat_obj_list = future_map(
  .x = all_stages.data,
  .f = function(x)(
    CreateSeuratObject(
      counts = x
    )
  )
)
saveRDS(
  object = all_stages.seurat_obj_list,
  file = "~/Projects/urchin_ap/data/working_data/all_stages.seurat_obj_list.rds"
  )

# creates a list of normalized, clustered seurat_objects using a list of unclustered seurat_objects from all_stages.seurat_obj_list
all_stages.clustered_seurat_obj_list = future_map(
  .x = all_stages.seurat_obj_list,
  .f = function(x)(
    subcluster_suite(
      seurat_object = x,
      res = 0.8,
      ndims = 1:15,
      strat = "tsne",
      jackstraw = F,
      coord_strat = "pca"
    )
  )
)
saveRDS(
  object = all_stages.clustered_seurat_obj_list,
  file = "~/Projects/urchin_ap/data/unmodified/all_stages.clustered_seurat_obj_list.rds"
  )

toc()
