# using the output from merge_isoform_names.r (list of merged isoform feature names), this function merges the associated counts for the new isoforms 
merge_paralog_feature_counts = function(matrix_path, features_path, features_df, key_name, path){

  orig_matrix = readMM(matrix_path)
  features_source = read_tsv(features_path, col_names = F)
  list_of_paralog_dfs = features_df %>%
    mutate(source_name = Name) %>% 
    group_by(source_name) %>% 
    group_map(~ (.x))

  list_of_paralog_indices = map(.x = list_of_paralog_dfs, function(x)(features_source$X2 %in% x$GeneID %>% which()))
  # list_of_paralog_indices = list_of_paralog_indices[-c(which((map(list_of_paralog_indices, length) %>% unlist()) == 0))]
  
  new_matrix = Matrix(ncol = ncol(orig_matrix), sparse = T)
  print("Performing matrix addition operations...")

  for (i in 1:length(list_of_paralog_indices)){
    if (length(list_of_paralog_indices[[i]])>1){
      summed_appended_row = Matrix(orig_matrix[list_of_paralog_indices[[i]],] %>% colSums(), nrow = 1, sparse = F)
      new_matrix = rbind(new_matrix, summed_appended_row)
    } else {
      unmodified_appended_row = orig_matrix[list_of_paralog_indices[[i]],]
      new_matrix = rbind(new_matrix, unmodified_appended_row)
    }
    
    if (i%%25 == 0){
      percent_complete = 100*(i/length(list_of_paralog_indices)) %>% round(digits = 4)
      print(paste0(percent_complete, "% of matrix operations for ", key_name, " completed..."))
    }
  }
  print(paste0("100% of matrix operations for ", key_name, " completed."))
  new_matrix = as(new_matrix, "dgCMatrix")
  rownames(new_matrix) = NULL
  new_matrix = new_matrix[-1,]
  
  first_paralog_index_list = map(.x = list_of_paralog_indices, .f = function(x)(x[1])) %>% unlist()
  other_paralog_index_list = map(.x = list_of_paralog_indices, .f = function(x)(x[-c(1)])) %>% unlist()
  print("Copying summed matrix data...")

  orig_matrix[first_paralog_index_list[1:1000],] = new_matrix[1:1000,]
  orig_matrix[first_paralog_index_list[1001:2000],] = new_matrix[1001:2000,]
  orig_matrix[first_paralog_index_list[2001:length(first_paralog_index_list)],] = new_matrix[2001:length(first_paralog_index_list),]

  gene_names_for_key = map(.x = list_of_paralog_dfs, .f = function(x)(x[1,3])) %>% unlist(use.names = F)
  gene_set_names_for_key = map(.x = list_of_paralog_dfs, .f = function(x)(x[1,4])) %>% unlist(use.names = F)
  key = features_source[first_paralog_index_list,] %>% 
    mutate(
      GeneID = paste0("merged-gene-", 
                  as.character(1:length(first_paralog_index_list))),
      Name = gene_names_for_key,
      GeneSet = gene_set_names_for_key
      )
  
  features_source[first_paralog_index_list,]$X2 = key$GeneID
  features_source = features_source[-other_paralog_index_list,]
  write_tsv(features_source, paste0(path, "features.tsv.gz"), col_names = F)
  
  orig_matrix = orig_matrix[-other_paralog_index_list,]
  writeMM(orig_matrix, paste0(path,"matrix.mtx.gz"))
  write_csv(key, paste0(key_name,"_key.csv"))
  
}
}
