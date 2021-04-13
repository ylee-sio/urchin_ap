get_num_pos_cells_multiple_genes = function(seurat_obj, gene_list, thresh, rev = NULL){
  
  res1 = map(
    .x = gene_list, 
    .f = function(x)(
      get_pos_cells(seurat_obj, 
                    x, 
                    thresh = thresh, 
                    rev = rev)
    )
  )  
  
  res2 = map(
    .x = res1,
    .f = function(x)(
      tibble(
        gene = x$gene %>% unique(),
        stage = x$stage %>% unique(),
        num_cells = x$num_cells %>% unique()
      )
    )
  )
  return(res2)
  
}
