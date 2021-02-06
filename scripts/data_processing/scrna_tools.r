library(tidyverse)

# returns which GeneIDs ("LOC...") from a list exist in the Foster 2020 S. purpuratus dataset 
in_urchin = function(seurat_obj, gene_id_list){
  res = which(gene_id_list %in% rownames(seurat_obj))
  return(res)
  }

# returns cells which have counts of at least x transcripts of a list of specific GeneIDs
get_cells_with_features = function(seurat_obj, features, thresh){

  list_of_lists_of_cells_with_feature = map(
	.x = features, 
	.f = function(x)(names(which(GetAssayData(seurat_obj)[which(GetAssayData(seurat_obj) %>% rownames() %in% x),] > thresh))))
  return(list_of_lists_of_cells_with_feature)
  
}
