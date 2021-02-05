library(tidyverse)

# returns which GeneIDs ("LOC...") exist in the Foster 2020 S. purpuratus dataset 

in_urchin = function(seurat_obj, gene_id_list){
  res = which(gene_id_list %in% rownames(seurat_obj))
  return(res)
  }
