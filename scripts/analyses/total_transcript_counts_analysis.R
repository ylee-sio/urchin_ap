urchin_5_transcript_stats = map2(
  .x = urchin_5_keys, 
  .y = urchin_5_clustered_seurat_objs,
  .f = function(x,y)(
    get_transcript_raw_counts_by_geneset(
      genesets_key = x,
      seurat_obj = y
      )
    )
  )

urchin_5_transcript_stats_bound = bind_rows(urchin_5_transcript_stats)

urchin_5_transcript_stats_bound_names_cleaned = mutate(
  urchin_5_transcript_stats_bound,
  "Gene set" = GeneSet %>% str_replace_all("_", " "),
  geneset_expression_ratio = geneset_expression_ratio*100,
  Stage = stage
)

urchin_5_transcript_stats_bound_names_cleaned$Stage = factor(
  urchin_5_transcript_stats_bound_names_cleaned$Stage,
  levels = c(
    urchin_5_transcript_stats_bound_names_cleaned$Stage %>% 
      unique()
  )
)

urchin_5_transcript_stats_bound_names_cleaned$"Gene set" = urchin_5_transcript_stats_bound_names_cleaned$"Gene set" %>% str_replace_all("membrane receptors gpcrs", "Membrane receptors: GPCRs")
urchin_5_transcript_stats_bound_names_cleaned$"Gene set" = urchin_5_transcript_stats_bound_names_cleaned$"Gene set" %>% str_replace_all("membrane trafficking exocytosis", "Membrane trafficking/exocytosis proteins")
urchin_5_transcript_stats_bound_names_cleaned$"Gene set" = urchin_5_transcript_stats_bound_names_cleaned$"Gene set" %>% str_replace_all("pattern recognition receptors", "PRRs")
urchin_5_transcript_stats_bound_names_cleaned$"Gene set" = urchin_5_transcript_stats_bound_names_cleaned$"Gene set" %>% str_replace_all("smts", "SMTs")
urchin_5_transcript_stats_bound_names_cleaned$"Gene set" = urchin_5_transcript_stats_bound_names_cleaned$"Gene set" %>% str_replace_all("transcription factors", "Transcription factors")
urchin_5_transcript_stats_bound_names_cleaned$"Gene set" = urchin_5_transcript_stats_bound_names_cleaned$"Gene set" %>% str_replace_all("translation factors", "Translation factors")
urchin_5_transcript_stats_bound_names_cleaned$"Gene set" = urchin_5_transcript_stats_bound_names_cleaned$"Gene set" %>% str_replace_all("enzymes", "Enzymes")

urchin_5_transcript_stats_bound_names_cleaned_1 = subset(
  urchin_5_transcript_stats_bound_names_cleaned,
  stage != "8 Cell" & stage != "64 Cell" & stage != "Morula")

total_transcript_counts_analysis_plot = ggplot(
  data = urchin_5_transcript_stats_bound_names_cleaned_1 %>% 
    subset(
      `Gene set` != "Enzymes" &
      `Gene set` != "PRRs"
      )
  ) +
  geom_line(
    aes(
      group = `Gene set`,
      x = Stage,
      y = geneset_expression_ratio,
      col = `Gene set`
    )
  ) +
  geom_point(
    mapping = aes(
      x = stage,
      y = geneset_expression_ratio,
      col = `Gene set`
    )
  ) +
  RotatedAxis() +
  ylab("Percentage of raw transcript counts") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1)) 



ggplotly(total_transcript_counts_analysis_plot) %>% 
  layout( 
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
  legend = list(
    size = 20,
    font = list(
      size = 20
      )
    )
  )

#all genes analysis
# 
# urchin_5_all_genes = rownames(urchin_5_sp1)
# urchin_5_all_genes_key = tibble(
#   X1 = urchin_5_all_genes,
#   X2 = urchin_5_all_genes,
#   X3 = urchin_5_all_genes,
#   GeneID = urchin_5_all_genes,
#   Name = urchin_5_all_genes,
#   GeneSet = "All genes")
# 
# urchin_5_all_genes_key_list = list(
#   urchin_5_all_genes_key,
#   urchin_5_all_genes_key,
#   urchin_5_all_genes_key,
#   urchin_5_all_genes_key,
#   urchin_5_all_genes_key,
#   urchin_5_all_genes_key,
#   urchin_5_all_genes_key,
#   urchin_5_all_genes_key
# )
# 
# urchin_5_all_genes_transcript_stats = map2(
#   .x = urchin_5_all_genes_key_list, 
#   .y = urchin_5_clustered_seurat_objs,
#   .f = function(x,y)(
#     get_transcript_raw_counts_by_geneset(
#       genesets_key = x,
#       seurat_obj = y
#     )
#   )
# )
# 
# urchin_5_all_genes_transcript_stats_bound = bind_rows(urchin_5_all_genes_transcript_stats)
# 
# urchin_5_all_genes_transcript_stats_bound_names_cleaned = mutate(
#   urchin_5_all_genes_transcript_stats_bound,
#   "Gene set" = GeneSet %>% str_replace_all("_", " "),
#   geneset_expression_ratio = geneset_expression_ratio*100,
#   Stage = stage
# )
