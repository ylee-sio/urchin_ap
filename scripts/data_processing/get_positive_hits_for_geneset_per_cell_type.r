get_positive_hits_for_geneset_per_cell_type = function(seurat_obj, gene_set, gene_set_name, thresh, relative, from_kegg = NULL, missing_germline_smts = NULL){
  
  # seurat_obj = urchin_5_spmb_cell_type_labelled
  # gene_set = urchin_5_smts_final_cleaned
  # gene_set_name = "smts"
  # thresh = 0
  # relative = F
  # from_kegg = F

  
  if (from_kegg == T){
    geneset_subset = subset(urchin_5_keys[[1]], GeneSet == gene_set_name)
  } else {
    geneset_subset = gene_set
  }
  
  seurat_obj = subset(
    seurat_obj, 
    features = geneset_subset$GeneID
  )
  
  cell_types = Idents(seurat_obj) %>%
    levels()
  
  stage = seurat_obj@meta.data$orig.ident %>%
    levels()
  
  num_all_genes_in_gene_set = seurat_obj %>% 
    rownames() %>% 
    length()
  
  # for one stage, create separate seurat_objs by cell type
  t1.b = map(
    .x = cell_types,
    .f = function(x)(
      subset(seurat_obj, idents = x)
    )
  )
  
  if(relative == T){
    assay_data_choice = "scale.data"
  } else{
    assay_data_choice = "counts"
  }
  # for one stage, for each cell type, get expression count matrix data
  t1.c = map(
    .x = t1.b,
    .f = function(x)(
      GetAssayData(
        x,
        slot = assay_data_choice
      )
    )
  )
  
  # for one stage, for each cell type, get cells for genes with any positive hits
  t2 = map(
    .x = t1.c,
    .f = function(x)(
      get_positive_hits(x, thresh = thresh)
    )
  )
  
  # t2_germline_test = map(t2, names)
  # t2_non_germline = t2_germline_test[-c(8)]
  # t2_germline = t2_germline_test[[8]]
  # 
  # stage_germline_absent_smts = map(.x = t2_non_germline, .f = function(x)(setdiff(x, t2_germline))) %>% 
  #   unlist() %>% 
  #   unique()
  # 
  # stage_germline_absent_smts_tibble = tibble(
  #   germline_missing_smts = stage_germline_absent_smts,
  #   stage = stage
  # )
  # 
  # get the number of smts present in each cell type with at least one transcript count
  t3 = map(
    .x = t2,
    .f = function(x)(
      length(x)
    )
  ) %>% 
    unlist()
  
  res = tibble(
    Stage = stage,
    "Cell type" = cell_types,
    "Gene set" = gene_set_name,
    "Number present" = t3,
    "Total existing" = num_all_genes_in_gene_set
  )
  
  
  
  if (missing_germline_smts == T){
    return(stage_germline_absent_smts_tibble)
  } else {
    return(res)
  }
  
}
