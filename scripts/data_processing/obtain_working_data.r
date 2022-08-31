tic()
# Run this script in terminal using 'Rscript obtain_working_data.R'
# This script only needs to be run one-time to obtain working data that has been processed with the parameters below. Adjust accordingly.

# setting workers for parallel operations
future::plan("multicore", workers = availableCores()-1)

# make a list of paths for each directory (developmental stage specific) containing the matrix, feature, and barcode files 
indv_stage_dirs = list.files("data/unmodified", full.names = T)

# creates a list of matrices (developmental stage specific) of expression data to be used to create seurat_objects
print("Reading in data to be processed... ")
all_stages.data = future_map(
  .x = indv_stage_dirs, 
  .f = function(x)(
    Read10X(
      data.dir = x,
    )
  )
)
print("Done reading in data.")

# creates a list of seurat_objects for each stage using list of scrna matrices from all_stages.data
print("Creating seurat objects... ")
all_stages.seurat_obj_list = future_map(
  .x = all_stages.data,
  .f = function(x)(
    CreateSeuratObject(
      counts = x
    )
  )
)
print("Done creating seurat objects. ")

#switching to sequential for cluster. future_map does not seem to effectively parallelize clustering
future::plan("sequential")
# creates a list of normalized, clustered seurat_objects using a list of unclustered seurat_objects from all_stages.seurat_obj_list
print("Clustering seurat objects... ")
all_stages.clustered_seurat_obj_list = map(
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
print("Done clustering seurat objects. ")

# creates separate directory outside of raw data directory containing working data
dir.create("~/Projects/urchin_ap/data/working_data")

# creates list of intermediate and final working data objects to be saved
working_data_object_list = c(
  all_stages.data,
  all_stages.seurat_obj_list,
  all_stages.clustered_seurat_obj_list
)

# creates list of intermediate and final working data file names to be saved
working_data_filename_list = c(
  "~/Projects/urchin_ap/data/working_data/all_stages.data.rds",
  "~/Projects/urchin_ap/data/working_data/all_stages.seurat_obj_list.rds",
  "~/Projects/urchin_ap/data/unmodified/all_stages.clustered_seurat_obj_list.rds"
)

future::plan("multicore", workers = 3)
# saves working data in parallel
print("Saving intermediary and final working data files... ")
future_map2(
  .x = working_data_object_list,
  .y = working_data_filename_list,
  .f = function(x,y)(
    saveRDS(
      object = x,
      file = y
      )
  )
)
print("Done saving files. ")
toc()