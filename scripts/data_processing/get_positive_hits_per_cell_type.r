# returns a dataframe with proportions of a given gene set that is expressed by each cell type 
get_positive_hits_per_cell_type = function(seurat_obj, gene_set, thresh, relative){
  
  cell_types = Idents(seurat_obj) %>% levels()
  stage = seurat_obj@meta.data$orig.ident %>% levels()
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
    "Gene set" = gene_set,
    "Number present" = t3,
    "Total existing" = num_all_genes_in_gene_set
  )
  
  return(res)
}
