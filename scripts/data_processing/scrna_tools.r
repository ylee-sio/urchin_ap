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
subcluster_suite = function(seurat_object, res, ndims){

  tic()
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object, nfeatures=2000)
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- RunPCA(seurat_object, verbose = FALSE)
  seurat_object <- RunICA(seurat_object, verbose = FALSE, ica.function = "icaimax")
  seurat_object <- FindNeighbors(seurat_object, dims = ndims)
  seurat_object <- FindClusters(seurat_object, resolution = res, verbose = FALSE, algorithm = 3)
  seurat_object <- RunUMAP(seurat_object, dims = ndims, n.epochs = 500, n.neighbors = 5L, min.dist = 0.05)
  seurat_object <- RunTSNE(seurat_object, dims = ndims)
  
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

get_percent_cells_represented = function(new_seurat_obj_cell_stats, old_seurat_obj_cell_stats){
  
  new_seurat_obj_cell_stats = new_seurat_obj_cell_stats %>% mutate(num_cells_new = num_cells)
  old_seurat_obj_cell_stats = old_seurat_obj_cell_stats %>% mutate(num_cells_old = num_cells)
  new_old_joined = left_join(new_seurat_obj_cell_stats, old_seurat_obj_cell_stats, by = "cluster_num")
  percent_repped = new_old_joined$num_cells_new/new_old_joined$num_cells_old
  cell_stats = mutate(new_seurat_obj_cell_stats, percent_representation = percent_repped, total_num_cells_in_cell_type = old_seurat_obj_cell_stats$num_cells)
  return(cell_stats)
  
}

# requires get_cell_stats
find_new_cluster_cell_idents_in_old_cluster_cell_idents = function(new_seurat_obj, old_seurat_obj){
  
  
  # get cluster labels for each cluster in new_seurat_obj
  new_seurat_obj_idents = levels(Idents(new_seurat_obj))
  
  # get a list of subsets of new_seurat_obj, which are delineated by cluster label
  new_seurat_obj_idents_subsetted_list = map(.x = new_seurat_obj_idents,
                                             .f = function(x)(subset(new_seurat_obj, 
                                                                     idents = x)
                                             )
  )
  
  # for each of the seurat_obj in new_seurat_obj, obtain a subsets of the old_seurat_obj where the old_seurat_obj contains only cells from new_seurat_obj.
  # each seurat_obj in this list will then show cells from new_seurat_obj present in old_seurat_obj
  cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj = map(.x = new_seurat_obj_idents_subsetted_list, 
                                                                 .f = function(x)(subset(old_seurat_obj,
                                                                                         cells = Cells(x)
                                                                 )
                                                                 )
  )
  
  # get cell stats for each item in cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj. 
  # the output here can be understood in the following way: 
  # 1. the first item in this list contains only cells from the first cluster in new_seurat_obj.
  # 2. cluster_num denotes the cluster ids of clusters in an old_seurat_obj.
  # 3. num_cells denotes the number of cells from the first cluster from a new_seurat_obj that belong to each
  # cluster in the old_seurat_obj
  cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj = map(.x = cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj, 
                                                                               .f = function(x)(
                                                                                 get_cell_stats(x)
                                                                               )
  )
  
  cells_stats_of_old_seurat_obj = get_cell_stats(old_seurat_obj)
  cells_stats_of_old_seurat_obj$cluster_num = levels(cells_stats_of_old_seurat_obj$cluster_num)
  
  cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column = map2(.x = cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj,
                                                                                                                  .y = new_seurat_obj_idents,
                                                                                                                  .f = function(x,y)(
                                                                                                                    mutate(x,
                                                                                                                           cluster_origin = as.character(y)
                                                                                                                    )
                                                                                                                  )
  )
  
  for (j in 1:length(cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column)){
    missing_cell_identities = which(cells_stats_of_old_seurat_obj$cluster_num %in% 
                                      cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column[[j]]$cluster_num == F)
    
    if (length(missing_cell_identities) >= 1){
      for (i in 1:length(missing_cell_identities)){
        cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column[[j]] =
          add_row(cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column[[j]],
                  cluster_num = cells_stats_of_old_seurat_obj$cluster_num[missing_cell_identities[i]],
                  num_cells = 0,
                  cluster_origin = cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column[[j]]$cluster_origin %>% unique())
        
      }
      
    } else {print(paste0("All cross representative clusters present for cluster ", j))}
    
    
    
  }
  res = map(.x = cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column,
            .f = function(x)(arrange(x, desc(as.character(cluster_num)))))

  new_seurat_obj_cell_stats = map(
    res,
    .f = function(x)(
      arrange(
        x,
        desc(cluster_num)
      )
    )
  )
  
  old_seurat_obj_cell_stats = get_cell_stats(old_seurat_obj)
  old_seurat_obj_cell_stats$cluster_num = old_seurat_obj_cell_stats$cluster_num %>% as.character()
  old_seurat_obj_cell_stats = old_seurat_obj_cell_stats %>% arrange(desc(cluster_num))
  
   rep_data = map(
    .x = res,
    .f = function(x)(
      get_percent_cells_represented(
        new_seurat_obj_cell_stats = x,
        old_seurat_obj_cell_stats = old_seurat_obj_cell_stats
      )
    )
  ) %>%
    bind_rows()
  
  return(rep_data)
  
}
