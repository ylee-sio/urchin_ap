all_genes = tibble(
  X1 = named_locus_list,
  X2 = named_locus_list,
  X3 = named_locus_list,
  GeneID = named_locus_list,
  Name = named_locus_list,
  GeneSet = "All genes"
)

merged_geneset_names =
  urchin_5_keys[[1]]$GeneSet %>%
  unique() %>%
  sort()

merged_geneset_names[8] = "All genes"

grouped_genesets =
  urchin_5_keys[[1]] %>%
  group_by(GeneSet) %>%
  group_map(~ (.x),
            .keep = TRUE)

grouped_genesets[[8]] = all_genes


tic()
merged_stats_v1.1.1 = map2(
  .x = grouped_genesets, 
  .y = merged_geneset_names, 
  .f = function(x,y)(
    make_merged_cell_stats(
      mappable_list_of_seurat_objects = urchin_5_clustered_seurat_objs,
      gene_set = x, 
      gene_set_name = y, 
      thresh = 0)
  )
)
toc()

merged_stats_v1.1.1_bound = merged_stats_v1.1.1 %>% bind_rows
merged_stats_v1.1.1_bound = read_csv("~/sp_transportomics/out.dump/merged_stats_v1.1.1.csv")

merged_stats_v1.1.1_bound$stage = factor(merged_stats_v1.1.1_bound$stage, levels = c(merged_stats_v1.1.1_bound$stage %>% unique()))

merged_stats_v1.1.1_bound = subset(
  merged_stats_v1.1.1_bound,
  stage != "8 Cell" & stage != "64 Cell" & stage != "Morula")

renamed_legend_2 = case_when(merged_stats_v1.1.1_bound$gene_set == "enzymes" ~ "Enzymes",
                             merged_stats_v1.1.1_bound$gene_set == "membrane_receptors_gpcrs" ~ "GPCRs",
                             merged_stats_v1.1.1_bound$gene_set == "membrane_trafficking_exocytosis" ~ "Membrane exocytosis proteins",
                             merged_stats_v1.1.1_bound$gene_set == "pattern_recognition_receptors" ~ "Pattern recognition receptors",
                             merged_stats_v1.1.1_bound$gene_set == "smts" ~ "SMTs",
                             merged_stats_v1.1.1_bound$gene_set == "transcription_factors" ~ "Transcription factors",
                             merged_stats_v1.1.1_bound$gene_set == "translation_factors" ~ "Translation factors",
                             merged_stats_v1.1.1_bound$gene_set == "All genes" ~ "All genes"
)
merged_stats_v1.1.1_bound$gene_set = renamed_legend_2

fig = plot_ly(
  merged_stats_v1.1.1_bound %>% subset(gene_set != "Pattern recognition receptors"),
  x = ~stage,
  y = ~gene_util_rate*100,
  color = ~gene_set,
  type = "box"
) %>%
  layout(boxmode = "group",
         yaxis = list(
           title = "Distribution of pairwise difference percentages\n of expression proportions of gene sets among cell types",
           titlefont = list(size = 30), 
           tickfont = list(size = 30)
         ),
         xaxis = list(
           title = "Stage",
           titlefont = list(size = 30), 
           tickfont = list(size = 30)
         ),
         legend = list(font = list(size = 20))
  )

fig
