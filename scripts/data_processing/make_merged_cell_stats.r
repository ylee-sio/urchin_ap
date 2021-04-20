make_merged_cell_stats = function(mappable_list_of_seurat_objects, gene_set, gene_set_name, thresh){
  res = map(
    .x = mappable_list_of_seurat_objects,
    .f = function(x)(get_basic_stats_2(x, gene_set, thresh = thresh))
  ) %>%
    bind_rows()
  
  util_rate = res$num_genes/(nrow(gene_set))
  
  res = mutate(
    res, 
    gene_set = gene_set_name, 
    cell_type_grouped_state = "unmerged", 
    gene_util_rate = util_rate)
  
  print(
    paste0(
      "Completed merging statistics for gene set: ", 
      gene_set_name)
  )
  return(res)
}
