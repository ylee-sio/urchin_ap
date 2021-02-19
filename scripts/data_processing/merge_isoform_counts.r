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
  
}
