search_purp = function(term){
  
  spuranno_res_1 = 
    spuranno[spuranno$product %>% str_which(term),] %>% 
    as_tibble() %>% 
    select(gene, product) %>% 
    transmute(GeneID = gene, Name = product) %>%
    unique()
  
  unique_gene_names = spuranno_res_1 %>% 
    group_by(GeneID) %>% 
    group_map(~.x[1,], .keep = T) %>% 
    bind_rows()
  
  spuranno_res_2 = tibble(
    GeneID = spuranno_res_1$GeneID %>% unique(),
    Name = spuranno_res_1$GeneID %>% unique()
  ) %>% 
    right_join(unique_gene_names, by = "GeneID") %>% 
    transmute(
      GeneID = GeneID,
      Name = Name.y
      )

    return(spuranno_res_2)
}
