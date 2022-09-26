library(jsonlite)

# read in kegg brite pathway categorized genes from KEGG
sp_kegg_all_genes = read_json("data/source/spu00001_all.json") %>% 
  unlist(use.names = F)

# filter for items starting with digits by getting an index for sp_kegg_transporters. these digits are entrez GeneIDs
present_geneid = sp_kegg_all_genes %>% 
  str_starts("[:digit:]") %>% 
  which()

# filtered list of items starting with digits
sp_kegg_all_genes_cleaned = sp_kegg_all_genes[present_geneid]

sp_kegg_all_genes_geneid_num = sp_kegg_all_genes_cleaned %>% 
  str_extract("[:digit:]{6,}") 

sp_kegg_slc_geneid_isolated = paste0("LOC", sp_kegg_all_genes_geneid_num)

sp_kegg_all = tibble(
  GeneID = c(sp_kegg_slc_geneid_isolated),
  Name = c(sp_kegg_all_genes_cleaned)
)

sp_kegg_all_final = subset(sp_kegg_all, GeneID != "LOCNA")
write_csv(sp_kegg_all_final, "data/working_data/sp_kegg_all_final.csv")

# # obtain index for those annotated as "SLC". after manual checking, this does get all of them since they are annotated as SLC by standard
# sp_kegg_slc_index = sp_kegg_transporters_cleaned %>% 
#   str_detect("SLC") %>% 
#   which()
# 
# # obtain index for those annotated as "ABC". after manual checking, this does get all of them since they are annotated as ABC by standard
# sp_kegg_abc_index = sp_kegg_transporters_cleaned %>% 
#   str_detect("ABC") %>% 
#   which()

# use the indices obtained above to filter from a cleaned df 
sp_kegg_slc_digit_start_extracted = sp_kegg_transporters_cleaned[sp_kegg_slc_index]
sp_kegg_abc_digit_start_extracted = sp_kegg_transporters_cleaned[sp_kegg_abc_index]

#SLC
#create usable tibble with GeneIDs and annotations for each final filtered list of genes of interest
sp_kegg_slc_geneid_num = sp_kegg_slc_digit_start_extracted %>% 
  str_extract("[:digit:]{6,}") 

sp_kegg_slc_geneid = paste0("LOC", sp_kegg_slc_geneid_num)

sp_kegg_slc_gene_anno = sp_kegg_slc_digit_start_extracted %>% 
  str_extract("SLC[:digit:]*[:alpha:]*[:digit:]*") %>% 
  unlist()

sp_kegg_slc = tibble(
  GeneID = sp_kegg_slc_geneid,
  Name = sp_kegg_slc_gene_anno
)

#ABC
sp_kegg_abc_geneid_num = sp_kegg_abc_digit_start_extracted %>% 
  str_extract("[:digit:]{6,}") 

sp_kegg_abc_geneid = paste0("LOC", sp_kegg_abc_geneid_num)

sp_kegg_abc_gene_anno = sp_kegg_abc_digit_start_extracted %>% 
  str_extract("ABC[:digit:]*[:alpha:]*[:digit:]*") %>% 
  unlist()

sp_kegg_abc = tibble(
  GeneID = sp_kegg_abc_geneid,
  Name = sp_kegg_abc_gene_anno
)

sp_kegg_smts = bind_rows(sp_kegg_slc, sp_kegg_abc)

write_csv(sp_kegg_smts, "data/working_data/sp_kegg_smts.csv")

# parses raw kegg brite files into manipulatable dataframes
parse_kegg_brite = function(kegg_json_file_path){
  awful_kegg_json = tibble(L2_l1 = "L2_l1",
                           L2_l2 = "L2_l2",
                           L3_l1 = "L3_l1",
                           L3_l2 = "L3_l2",
                           L3_l3 = "L3_l3",
                           L4_l1 = "L4_l1",
                           L4_l2 = "L4_l2",
                           L4_l3 = "L4_l3",
                           L4_l4 = "L4_l4",
                           L5_l1 = "L5_l1",
                           L5_l2 = "L5_l2",
                           L5_l3 = "L5_l3",
                           L5_l4 = "L5_l4",
                           L5_l5 = "L5_l5")
  
  kegg_sp_gene_set = read_json(kegg_json_file_path)
  kegg_sp_gene_set_tibble = kegg_sp_gene_set %>% as_tibble()
  kegg_sp_gene_set_tibble_processed = kegg_sp_gene_set_tibble[2]
  
  l1_children_names = map(.x = kegg_sp_gene_set_tibble_processed$children, .f = function(x)(x$name))
  l1_children_lengths = map(.x = kegg_sp_gene_set_tibble_processed$children, .f = function(x)(x$children %>% length()))
  l1_children = map(.x = kegg_sp_gene_set_tibble_processed$children, .f = function(x)(x$children))
  
  for (g in 1:length(l1_children_names)){
    l2_children_names = map(.x = l1_children[[g]], .f = function(x)(x$name))
    l2_children_lengths = map(.x = l1_children[[g]], .f = function(x)(x$children %>% length()))
    
    awful_kegg_json = add_row(awful_kegg_json, 
                              L2_l2 = l2_children_names %>% unlist(),
                              L2_l1 = l1_children_names[[g]]
    )
    
    l2_children = map(.x = l1_children[[g]], .f = function(x)(x$children))
    l2_non_null = which(map(l2_children, is_null) == F)
    
    # print(paste0("g:", g))
    

    for (j in l2_non_null){
      l3_children_names = map(.x = l2_children[[j]], .f = function(x)(x$name))
      l3_children_lengths = map(.x = l2_children[[j]], .f = function(x)(x$children %>% length()))
      
      awful_kegg_json = add_row(awful_kegg_json, 
                                L3_l3 = l3_children_names %>% unlist(),
                                L3_l2 = l2_children_names[[j]],
                                L3_l1 = l1_children_names[[g]]
      )
      
      l3_children = map(.x = l2_children[[j]], .f = function(x)(x$children))
      l3_non_null = which(map(l3_children, is_null) == F)
      
      # print(paste0("j:", j))
      
      for (i in l3_non_null){
        l4_children_names = map(.x = l3_children[[i]], .f = function(x)(x$name))
        awful_kegg_json = add_row(awful_kegg_json, 
                                  L4_l4 = l4_children_names %>% unlist(), 
                                  L4_l3 = l3_children_names[[i]],
                                  L4_l2 = l2_children_names[[j]],
                                  L4_l1 = l1_children_names[[g]]
        )
        
        l4_children = map(.x = l3_children[[i]], .f = function(x)(x$children))
        l4_non_null = which(map(l4_children, is_null) == F)
        
        # print(paste0("i:", i))
        
        for (k in l4_non_null){
          l5_children_names = map(.x = l4_children[[k]], .f = function(x)(x$name))
          awful_kegg_json = add_row(awful_kegg_json,
                                    L5_l5 = l5_children_names %>% unlist(),
                                    L5_l4 = l4_children_names[[k]], 
                                    L5_l3 = l3_children_names[[i]],
                                    L5_l2 = l2_children_names[[j]],
                                    L5_l1 = l1_children_names[[g]]
          )
        
          # print(paste0("k:", k))
        
      }
    }
    }
  }
  
    awful_kegg_json = awful_kegg_json[-1,]
  return(awful_kegg_json)
  
}
