# returns a data frame with percentages of a geneset shared by various cell types in a seurat_obj

get_overlap_percentages_one_stage_one_geneset = function(seurat_obj, gene_set){
  
  # seurat_obj = urchin_5_clustered_seurat_objs[[1]]
  # gene_set = urchin_5_smts_final_cleaned_with_mito_smts %>% mutate(GeneSet = "smts")
  
  b = Idents(seurat_obj) %>% levels()
  
  
  # subset a single seurat obj by ident or "cell type"
  c = map(
    .x = b,
    .f = function(x)(
      subset(
        seurat_obj,
        ident = x)
    )
  )
  
  # subset each seurat obj in list of seurat objs by GeneIDs of a single gene set
  d = map(
    .x = c,
    .f = function(x)(
      subset(
        x,
        features = gene_set$GeneID
      )
    )
  )
  
  # get transcript count data from each seurat obj in list of seurat objs that have been subsetted by GeneIDs in a single gene set
  e = map(
    .x = d,
    .f = function(x)(
      GetAssayData(x, slot = "counts")
    )
  )
  
  # two things achieved here: get cell IDs for cells expressing each gene in gene set. 
  # not only will this tell us #1 which and how many cells express a gene, but also #2 how many genes have any expression
  # at all. remember, this is still being done for each cell type.
  f = map(
    .x = e,
    .f = function(x)(
      get_positive_hits(
        x,
        0
      )
    )
  )
  
  # get the names of the genes that have any number of cells expressing each gene for each cell type.
  g = map(
    .x = f,
    .f = function(x)(
      names(x)
    )
  )
  
  # create Venn object- requires a list of vectors to be compared
  g = Venn(g)
  
  # find pairwise differences
  h = discern_pairs(g)
  
  i = map(
    .x = h,
    .f = function(x)(
      length(x)/nrow(gene_set)
    )
  )
  
  j = names(i)
  k = unlist(i, use.names = F)
    
  l = tibble(
    gene_set = gene_set$GeneSet %>% unique(),
    comparison_pair = j,
    difference_percentage = k,
    mean_pairwise_difference = mean(k),
    stage = seurat_obj@meta.data$orig.ident %>% levels()
  )
  return(l)
  
}
