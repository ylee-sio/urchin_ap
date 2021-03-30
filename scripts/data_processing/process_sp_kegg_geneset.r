process_sp_kegg_geneset = function(parsed_and_filtered_sp_kegg_geneset, unmerged_feature_df){
  
  res1 = parsed_and_filtered_sp_kegg_geneset[colSums(!is.na(parsed_and_filtered_sp_kegg_geneset)) > 0]
  res2 = res1[[length(res1)]] %>% str_split("; ", simplify = T) %>% as_tibble()
  res3 = res2 %>% mutate(GeneID = paste0("LOC", str_extract(res2[[1]], "[:digit:]*")))
  res4 = tibble(GeneID = res3$GeneID, Name = res3$V2, additional_info = res3$V1)
  res5 = which(res4$GeneID %in% rownames(unmerged_feature_df))
  res6 = res4[res5,]
  
  return(res6)
  
}
