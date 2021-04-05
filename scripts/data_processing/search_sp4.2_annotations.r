search_sp4.2_annotations = function(term){
  
  spuranno_v2 = spuranno %>% as_tibble(spuranno)
  spuranno_v3 = spuranno_v2$product
  spuranno_v4 = str_which(spuranno_v3, term)
  spuranno_v5 = spuranno_v2[spuranno_v4,] %>% select(gene, product) %>% unique()
  spuranno_v6 = spuranno_v2[spuranno_v4,] %>% select(gene) %>% unique()
  spuranno_v7 = tibble(GeneID = spuranno_v6$gene, Name = spuranno_v6$gene)
  spuranno_v8 = map(.x = map(.x = spuranno_v7$GeneID, .f = function(x)(subset(spuranno_v5, gene == x))), .f = function(x)(x[1,])) %>% 
    bind_rows() %>%
    mutate(GeneID = gene, Name = product)
  
  return(spuranno_v8)
}
