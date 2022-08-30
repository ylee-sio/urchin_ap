all_mean_pairwise_differences = 
  map(
  .x = urchin_5_clustered_labelled_seurat_objs,
  .f = function(x)(
    get_average_overlap_percentages_one_stage_multiple_genesets(
      seurat_obj = x,
      gene_set_list = grouped_genesets
    )
  )
)
all_mean_pairwise_differences = 
  all_pairwise_differences_bound %>% 
  group_by(stage) %>% 
  group_map(~.x)

all_pairwise_differences_bound = all_mean_pairwise_differences %>% bind_rows()  
# write_csv(all_pairwise_differences_bound, "~/sp_transportomics/out.dump/all_mean_pairwise_differences_bound_labelled_cell_clusters.csv")
all_pairwise_differences_bound = read_csv("~/sp_transportomics/out.dump/all_mean_pairwise_differences_bound_labelled_cell_clusters.csv")
all_pairwise_differences_bound$Stage = factor(all_pairwise_differences_bound$Stage,
                                                   levels = c(all_pairwise_differences_bound$Stage %>% unique()))



all_pairwise_differences_bound_final_1 = all_pairwise_differences_bound %>% 
  subset(stage == "Late Gastrula")

all_pairwise_differences_bound_final_2 = 
  all_pairwise_differences_bound_final_1 %>% 
  group_by(stage) %>% 
  group_map(~.x) 

all_pairwise_differences_bound_final_3 = 
  all_pairwise_differences_bound_final_2[[1]] %>% 
  as_tibble()
  

renamed_legend_1 = case_when(all_pairwise_differences_bound_final_3$gene_set == "enzymes" ~ "Enzymes",
          all_pairwise_differences_bound_final_3$gene_set == "membrane_receptors_gpcrs" ~ "GPCRs",
          all_pairwise_differences_bound_final_3$gene_set == "membrane_trafficking_exocytosis" ~ "Membrane exocytosis proteins",
          all_pairwise_differences_bound_final_3$gene_set == "pattern_recognition_receptors" ~ "Pattern recognition receptors",
          all_pairwise_differences_bound_final_3$gene_set == "smts" ~ "SMTs",
          all_pairwise_differences_bound_final_3$gene_set == "transcription_factors" ~ "Transcription factors",
          all_pairwise_differences_bound_final_3$gene_set == "translation_factors" ~ "Translation factors",
          all_pairwise_differences_bound_final_3$gene_set == "All genes" ~ "All genes"
          )
all_pairwise_differences_bound_final_3$gene_set = renamed_legend_1

fig = plot_ly(
  all_pairwise_differences_bound_final_3 %>% bind_rows() %>% subset(gene_set != "pattern_recognition_receptors"),
  x = ~Stage,
  y = ~difference_percentage*100,
  color = ~gene_set,
  type = "box"
) %>%
  layout(boxmode = "group",
         yaxis = list(
           title = "Distribution of pairwise difference percentages\n of expression proportions of gene sets among cell types",
           titlefont = list(size = 25), 
           tickfont = list(size = 25)
           ),
         xaxis = list(
           title = "Stage",
           titlefont = list(size = 25), 
           tickfont = list(size = 25)
         ),
         legend = list(font = list(size = 25))
         )

fig

pairwise_t_test_all_stages = map(
  .x = all_pairwise_differences_bound_final_2,
  .f = function(x)(
    pairwise.t.test(
      x$difference_percentage, 
      x$gene_set, 
      p.adjust.method = "bonf"
      )
  )
)


