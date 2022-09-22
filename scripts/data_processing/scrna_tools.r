library(tidyverse)
library(Seurat)
library(tictoc)

# returns which GeneIDs ("LOC...") from a list exist in the Foster 2020 S. purpuratus dataset
in_urchin = function(seurat_obj, gene_id_list){
  res = which(gene_id_list %in% rownames(seurat_obj))
  return(res)
  }

# returns cells which have counts of at least x transcripts of a list of specific GeneIDs
get_cells_with_features = function(seurat_obj, features, thresh){

  list_of_lists_of_cells_with_feature = map(
	.x = features,
	.f = function(x)(names(which(GetAssayData(seurat_obj)[which(GetAssayData(seurat_obj) %>% rownames() %in% x),] > thresh))))
  return(list_of_lists_of_cells_with_feature)

}

# wrapper function for preparing for clustering cells with options for clustering strategies PCA, ICA, and UMAP
subcluster_suite = function(seurat_object, res, ndims, strat, jackstraw, coord_strat){

  tic()
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object, nfeatures=2000)
  seurat_object <- ScaleData(seurat_object)

  if (coord_strat == "pca"){

    seurat_object <- RunPCA(seurat_object, verbose = FALSE)

  }

  if (coord_strat == "ica"){

    seurat_object <- RunICA(seurat_object, verbose = FALSE, ica.function = "icaimax")

  }
  seurat_object <- FindNeighbors(seurat_object, dims = ndims)
  seurat_object <- FindClusters(seurat_object, resolution = res, verbose = FALSE, algorithm = 3)

  if (jackstraw == T){
  seurat_object = JackStraw(seurat_object, num.replicate = 200)
  seurat_object = ScoreJackStraw(seurat_object, dims = 1:20)
  }

  if (strat == "umap"){
    seurat_object <- RunUMAP(seurat_object, dims = ndims, n.epochs = 500, n.neighbors = 5L, min.dist = 0.05)
  } else {
    seurat_object <- RunTSNE(seurat_object, dims = ndims)
  }

  toc()
  print("****** Job finished. ******")

  return(seurat_object)

}

# gene set quick search
# from a tibble of a gene set where columns = c(GeneID,Name)
# quickly search name matches and obtain tibble of search by name
gsqs = function(gene_set, gene_name){
  res = gene_set[gene_set$Name %>% str_which(gene_name),]
  return(res)
}

# returns numbers of cells from each Ident in a seurat_obj
get_cell_stats = function(seurat_obj){

  cluster_ids = unique(Idents(seurat_obj))
  cell_nums = map(.x = unique(Idents(seurat_obj)), .f = function(x)(seurat_obj %>%
                                                                      subset(idents = x) %>%
                                                                      GetAssayData() %>%
                                                                      ncol())
  ) %>%
    unlist()
  cell_stats_tibble = tibble(cluster_num = cluster_ids, num_cells = cell_nums)
  return(cell_stats_tibble)
}

# returns a list of seurat_obj objects that have been clustered with a user provided 1) single ndim range and 2) list of resolutions
# depends on subcluster_suite
cluster_multiple_res = function(seurat_obj, ndims, res_list){
  tic()
  plan(strategy = multiprocess, workers = 8)
  cluster_list = future_map(.x = res_list, .f = function(x)(subcluster_suite(seurat_obj, x, ndims)))
  return(cluster_list)
  toc()
}

# returns a list of seurat_obj objects that ahve been clustered with a user provided 1) single resolution and 2) list of ndim ranges
# depends on subcluster_suite
cluster_multiple_ndims = function(seurat_obj, ndims_list, res){
  tic()
  plan(strategy = multiprocess, workers = 16)
  cluster_list = future_map(.x = ndims_list, .f = function(x)(subcluster_suite(seurat_obj, res, x)))
  return(cluster_list)
  toc()
}

# added get_pos_cells, thresholded extractions of cells using a single target
get_pos_cells = function(seurat_obj, gene_id, thresh, rev = NULL){

  if (rev == T){

    cells = which(GetAssayData(seurat_obj, slot = "counts")[gene_id,] < thresh) %>% names()

  } else {
    cells = which(GetAssayData(seurat_obj, slot = "counts")[gene_id,] > thresh) %>% names()

  }

  if (length(cells) > 0){

    pos_cells = subset(seurat_obj, cells = cells)
    pos_cells_raw_count_values = GetAssayData(pos_cells, slot = "counts")[gene_id,]
    pos_cells_normalized_cell_type_relative_avg_log2FC_values = GetAssayData(pos_cells)[gene_id,]
    pos_cells_orig_cluster_id = map(.x = cells, .f = function(x)(seurat_obj@meta.data[x,])) %>% bind_rows() %>% select(seurat_clusters)

    pos_cells_tibble = tibble(cells = cells,
                              raw_counts = pos_cells_raw_count_values,
                              normalized_cell_type_relative_avg_log2FC = pos_cells_normalized_cell_type_relative_avg_log2FC_values,
                              orig_cluster_id = pos_cells_orig_cluster_id$seurat_clusters,
                              gene = gene_id,
                              stage = seurat_obj@meta.data$orig.ident %>% levels(),
                              num_cells = length(cells)
    )
  } else {

    pos_cells_tibble = tibble(
      cells = c(NA) ,
      raw_counts = 0,
      normalized_cell_type_relative_avg_log2FC = 0,
      orig_cluster_id = c(NA) %>% as.factor(),
      gene = gene_id,
      stage = seurat_obj@meta.data$orig.ident %>% levels(),
      num_cells = 0 %>% as.integer()
    )
  }
  return(pos_cells_tibble)
}
