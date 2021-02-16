merge_isoform_feature_names_gene_label = function(features_path, features_df){
  
  features_source = read_tsv(features_path, col_names = F)
  features_df = genes_to_be_merged
  
  list_of_paralog_dfs = features_df %>%
    mutate(source_name = Name) %>% 
    group_by(source_name) %>% 
    group_map(~ (.x))
  
  plan("multisession", workers = 8)
  list_of_paralog_indices = future_map(list_of_paralog_dfs, ~(features_source$X2 %in% .x$GeneID %>% which()))
  first_paralog_index_list = map(.x = list_of_paralog_indices, .f = function(x)(x[1]))
  merged_gene_key = bind_cols(features_source[unlist(first_paralog_index_list),]$X2, paste0("merged_gene_",1:length(features_source[unlist(first_paralog_index_list),]$X2)))
  names(merged_gene_key) = c("GeneID", "Name")
  features_source[unlist(first_paralog_index_list),]$X2 = paste0("merged_gene_",1:nrow(merged_gene_key))
  
  return(features_source)
  
  }
