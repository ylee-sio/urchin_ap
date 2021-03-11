get_expression_data_from_cells_with_gradients_of_marker_expression = function(seurat_obj, gene){
  
  gene_pos_cells = get_pos_cells(seurat_obj, gene_id = gene, thresh = 0, rev = F)
  gene_transcript_count_quantiles = quantile(gene_pos_cells$raw_counts, c(0.95, 0.90, 0.85, 0.80, 0.75))
  gene_transcript_count_quantiles_tibble = tibble(
    
    quantile = names(gene_transcript_count_quantiles),
    transcript_count = gene_transcript_count_quantiles
    
  )
  
  mapped_gene_transcript_count_quantiles_tibble = 
    map(
      .x = names(gene_transcript_count_quantiles), 
      .f = function(x)(
        subset(
          gene_transcript_count_quantiles_tibble,
          quantile == x)
      )
    )
  
  gene_pos_cells_95_100_quantile = gene_pos_cells %>% subset(raw_counts >= mapped_gene_transcript_count_quantiles_tibble[[1]]$transcript_count)
  gene_pos_cells_90_95_quantile = gene_pos_cells %>% subset(raw_counts >= mapped_gene_transcript_count_quantiles_tibble[[2]]$transcript_count & raw_counts < mapped_gene_transcript_count_quantiles_tibble[[1]]$transcript_count) 
  gene_pos_cells_85_90_quantile = gene_pos_cells %>% subset(raw_counts >= mapped_gene_transcript_count_quantiles_tibble[[3]]$transcript_count & raw_counts < mapped_gene_transcript_count_quantiles_tibble[[2]]$transcript_count)
  gene_pos_cells_80_85_quantile = gene_pos_cells %>% subset(raw_counts >= mapped_gene_transcript_count_quantiles_tibble[[4]]$transcript_count & raw_counts < mapped_gene_transcript_count_quantiles_tibble[[3]]$transcript_count)
  gene_pos_cells_75_80_quantile = gene_pos_cells %>% subset(raw_counts >= 1 & raw_counts < mapped_gene_transcript_count_quantiles_tibble[[4]]$transcript_count) 
  
  gene_pos_cells_quantile_list = list(
    
    gene_pos_cells_95_100_quantile,
    gene_pos_cells_90_95_quantile,
    gene_pos_cells_85_90_quantile,
    gene_pos_cells_80_85_quantile,
    gene_pos_cells_75_80_quantile
    
  )
  
  gene_pos_cells_quantile_names_list = list(
    "_95_100_quantile",
    "_90_95_quantile",
    "_85_90_quantile",
    "_80_85_quantile",
    "_75_80_quantile"
  )
  
  gene_pos_cells_nrows = map(gene_pos_cells_quantile_list, nrow)
  
  for (i in 1:length(gene_pos_cells_quantile_list)){
    
    if(gene_pos_cells_nrows[[i]] > 0){
      Idents(seurat_obj, cells = gene_pos_cells_quantile_list[[i]]$cells) = paste0(gene, gene_pos_cells_quantile_names_list[[i]])
    }
  }
  
  return(seurat_obj)
