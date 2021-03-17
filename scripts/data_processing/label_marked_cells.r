label_marked_cells = function(seurat_obj, marker_list, cell_type_labels, thresh){
  
  marked_cells_df = map(
    .x = marker_list,
    .f = function(x)(
      get_pos_cells(
        seurat_obj, 
        gene_id = x, 
        thresh = thresh, 
        rev = F)
    )
  ) %>% 
    bind_rows() 
  
  marked_cells = marked_cells_df$cells %>% unique()
  print(marked_cells_df)
  Idents(seurat_obj, cells = marked_cells) = cell_type_labels
  
  return(seurat_obj)
}
