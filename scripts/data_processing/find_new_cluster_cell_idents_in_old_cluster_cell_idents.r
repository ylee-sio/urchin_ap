find_new_cluster_cell_idents_in_old_cluster_cell_idents = function(new_seurat_obj, old_seurat_obj){
  
  new_seurat_obj_idents = levels(Idents(new_seurat_obj))
  new_seurat_obj_idents_subsetted_list = map(.x = new_seurat_obj_idents,
                                             .f = function(x)(subset(new_seurat_obj, 
                                                                     idents = x)
                                             )
  )
  cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj = map(.x = new_seurat_obj_idents_subsetted_list, 
                                                                 .f = function(x)(subset(old_seurat_obj,
                                                                                         cells = Cells(x)
                                                                 )
                                                                 )
  )
  
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
    missing_cell_identities = which(cells_stats_of_old_seurat_obj$cluster_num %in% cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column[[j]]$cluster_num == F)
    
    if (length(missing_cell_identities >= 1)){
      for (i in 1:length(missing_cell_identities)){
        cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column[[j]] =
          add_row(cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column[[j]],
                  cluster_num = cells_stats_of_old_seurat_obj$cluster_num[missing_cell_identities[i]],
                  num_cells = 0,
                  cluster_origin = cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column[[j]]$cluster_origin %>% unique())
        
      }
      
    } else {print("mop")}
    
    
    
  }
  res = map(.x = cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column,
            .f = function(x)(arrange(x, desc(cluster_num))))
  # res = bind_rows(cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column)
  # res$cluster_origin = as.factor(res$cluster_origin)
  return(res)
  
}
