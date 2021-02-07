library(tidyverse)

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
