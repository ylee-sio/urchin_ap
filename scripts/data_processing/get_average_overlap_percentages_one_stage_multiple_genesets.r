get_average_overlap_percentages_one_stage_multiple_genesets = function(seurat_obj, gene_set_list){
  
  res = map(.x = gene_set_list, 
            .f = function(x)(
              get_overlap_percentages_one_stage_one_geneset(seurat_obj, x)
            )
  ) %>% bind_rows()
  
  res = res %>% mutate(Stage = seurat_obj@meta.data$orig.ident %>% levels())
  print(paste0(seurat_obj@meta.data$orig.ident %>% levels(), " done."))
  return(res)
  
}
get_p
