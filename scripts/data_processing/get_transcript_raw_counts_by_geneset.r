# returns raw transcript counts per gene set 
get_transcript_raw_counts_by_geneset = function(genesets_key, seurat_obj){
  
  stage = 
    seurat_obj@meta.data$orig.ident %>%
    levels()
  
  geneset_names = 
    genesets_key$GeneSet %>% 
    unique() %>% 
    sort()
  
  total_transcript_count = 
    GetAssayData(seurat_obj, slot = "counts") %>% 
    sum()
  
  geneset_key_grouped_by_geneset = 
    genesets_key %>% 
    group_by(GeneSet) %>% 
    group_map(~ (.x),
              .keep = TRUE)
  
  seurat_obj_subsets_by_geneset = map(
    .x = geneset_key_grouped_by_geneset, 
    .f = function(x)(
      subset(
        seurat_obj,
        features = x$GeneID
      )
    )
  )
  
  transcript_counts_by_geneset = map(
    .x = seurat_obj_subsets_by_geneset,
    .f = function(x)(
      GetAssayData(
        x,
        slot = "counts") %>% 
        sum()
    )
  )
  
  res = tibble(
    stage = stage,
    GeneSet = geneset_names,
    geneset_raw_transcript_counts = transcript_counts_by_geneset %>% unlist(),
    total_raw_transcript_counts = total_transcript_count,
    geneset_expression_ratio = geneset_raw_transcript_counts/total_raw_transcript_counts
    
  )
  
  return(res)
}
