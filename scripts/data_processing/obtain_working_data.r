library(tictoc)
# Run this script in terminal using 'Rscript obtain_working_data.R'
# This script only needs to be run one-time to obtain working data that has been processed with the parameters below. Adjust accordingly.

# make a list of paths for each directory (developmental stage specific) containing the matrix, feature, and barcode files 
indv_stage_dirs = list.files("data/unmodified", full.names = T)

# creates a list of matrices (developmental stage specific) of expression data to be used to create seurat_objects
print("Reading in data to be processed... ")
all_stages.data = map(
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
all_stages.seurat_obj_list = map(
  .x = all_stages.data,
  .f = function(x)(
    CreateSeuratObject(
      counts = x
    )
  )
)
print("Done creating seurat objects. ")

# remove features with 0 counts
assay_data_counts = map(
  .x = all_stages.seurat_obj_list,
  .f = function(x)(
    GetAssayData(
      object = x,
      slot = "counts"
    )
  )
)
positive_features_each_stage = map(
  .x = assay_data_counts,
  .f = function(x)(
    which(rowSums(x) > 0) %>% 
      names()
  )
)
stages_filtered_positive_features = map2(
  .x = all_stages,
  .y = positive_features_each_stage,
  .f = function(x,y)(
    subset(
      x = x,
      features = y
    )
  )
)


# creates a list of normalized, clustered seurat objects using a list of unclustered seurat_objects from all_stages.seurat_obj_list
print("Clustering seurat objects... ")
all_stages.clustered_seurat_obj_list = map(
  .x = stages_filtered_positive_features,
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

print("Saving intermediary and final working data files... ")

saveRDS(all_stages.data,  "~/Projects/urchin_ap/data/working_data/all_stages.data.rds")
saveRDS(all_stages.seurat_obj_list,  "~/Projects/urchin_ap/data/working_data/all_stages.seurat_obj_list.rds")

#this will not contain features with 0 counts
saveRDS(all_stages.clustered_seurat_obj_list, "~/Projects/urchin_ap/data/all_stages.clustered_seurat_obj_list.rds")

print("Done saving files. ")
toc()
