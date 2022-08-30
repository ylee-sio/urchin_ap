urchin_5_splg_clustered_labelled_v2 = urchin_5_splg_clustered

stage_names = list("sp8", "sp64", "spmor", "speb", "sphb", "spmb", "speg", "splg")

splg_nanos2 = map(
  .x = urchin_5_clustered_seurat_objs[2:8],
  .f = function(x)(get_expression_data_from_cells_with_gradients_of_marker_expression(x, "Nanos2"))
)

splg_seawi = map(
  .x = urchin_5_clustered_seurat_objs,
  .f = function(x)(get_expression_data_from_cells_with_gradients_of_marker_expression(x, "LOC373434"))
)

splg_tdrd1 = map(
  .x = urchin_5_clustered_seurat_objs,
  .f = function(x)(get_expression_data_from_cells_with_gradients_of_marker_expression(x, "LOC100887863"))
)

splg_vasa = map(
  .x = urchin_5_clustered_seurat_objs,
  .f = function(x)(get_expression_data_from_cells_with_gradients_of_marker_expression(x, "vasa"))
)

map2(
  .x = splg_nanos2,
  .y = stage_names[2:8],
  .f = function(x,y)(
    
    lab_DP(x, 
           named_locus_df,
           plot_title = paste0(y, "_nanos2"),
           F, 
           6000,
           T,
           0,
           3, 
           0,
           cols = c("grey",
                    "blue")
    )
    
  )
)

map2(
  .x = splg_seawi,
  .y = stage_names,
  .f = function(x,y)(
    
    lab_DP(x, 
           named_locus_df,
           plot_title = paste0(y, "_seawi"),
           name = paste0(y, "_seawi"),
           F, 
           6000,
           T,
           0,
           3, 
           0,
           cols = c("grey",
                    "blue")
    )
    
  )
)


map2(
  .x = splg_tdrd1,
  .y = stage_names,
  .f = function(x,y)(
    
    lab_DP(x, 
           named_locus_df,
           plot_title = paste0(y, "_tdrd1"),
           name = paste0(y, "_tdrd1"),
           F, 
           6000,
           T,
           0,
           3, 
           0,
           cols = c("grey",
                    "blue")
    )
    
  )
)

map2(
  .x = splg_vasa,
  .y = stage_names,
  .f = function(x,y)(
    
    lab_DP(x, 
           named_locus_df,
           plot_title = paste0(y, "_vasa"),
           name = paste0(y, "_vasa"),
           F, 
           6000,
           T,
           0,
           3, 
           0,
           cols = c("grey",
                    "blue")
    )
    
  )
)

#*************************************
get_cell_stats(urchin_5_sp3_clustered_labelled_v2_germline) %>% View()

DimPlot(marked_cells_test)
lab_DP(marked_cells_test, 
       standard_markers_for_display,
       plot_title = "SMT expression profiles of cell types at the 8 cell stage",
       F, 
       1200,
       T,
       0,
       3, 
       0,
       cols = c("grey",
                "blue")
)




# 1. find all markers for each cluster in each stage 
urchin_5_sp2_clustered_markers = FindAllMarkers(urchin_5_sp2_clustered)
# 2. merge clusters which contain at least X number of markers in common
urchin_5_sp2_clustered_markers_grouped =
  urchin_5_sp2_clustered_markers %>% 
  group_by(cluster) %>% 
  slice_min(order_by = p_val, n = 50) %>% 
  select(gene, cluster) %>% 
  group_map(~.x %>% unlist(use.names = F), .keep = F)

urchin_5_sp2_clustered_venn = Venn(urchin_5_sp2_clustered_markers_grouped)
overlapping_markers = overlap_pairs(urchin_5_sp2_clustered_venn) 
comparison_set_names = names(overlapping_markers)
overlapping_marker_numbers = overlap_pairs(urchin_5_sp2_clustered_venn) %>% map(length) %>% unlist(use.names = F)
marker_numbers_comparison_stats = tibble(
  comparison = comparison_set_names,
  overlapping_numbers_of_genes = overlapping_marker_numbers
)
quantile_distribution_of_number_of_overlapping_markers = 
  quantile(marker_numbers_comparison_stats$overlapping_numbers_of_genes, prob = c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.96, 0.97, 0.98, 0.99, 1))
marker_numbers_comparison_stats %>% subset(overlapping_numbers_of_genes >= quantile_distribution_of_number_of_overlapping_markers["95%"])

urchin_5_sp2_dp = lab_DP(urchin_5_sp1_clustered_2, 
                         standard_markers_for_display,
                         plot_title = "SMT expression profiles of cell types at the 8 cell stage",
                         F, 
                         4000,
                         T,
                         0,
                         3, 
                         0,
                         cols = c("grey",
                                  "blue")
)





pks_gradient_expression_cells = map(
  .x = urchin_5_clustered_labelled_seurat_objs,
  .f = function(x)(
    get_expression_data_from_cells_with_gradients_of_marker_expression(
      seurat_obj = x,
      gene = "Endo16"
      )
  )
)
    
gcm_gradient_expression_cells = map(
  .x = urchin_5_clustered_labelled_seurat_objs,
  .f = function(x)(
    get_expression_data_from_cells_with_gradients_of_marker_expression(
      seurat_obj = x,
      gene = "gcm"
    )
  )
)

pks_cell_stats = map(pks_gradient_expression_cells, get_cell_stats)
gcm_cell_stats = map(gcm_gradient_expression_cells, get_cell_stats)

# erg_alx1_specific_smts = tibble(
#   
#   GeneID = c("ABCA3", "ABCB1", "ABCC4", "SLC4A3", "SLC4A10", "SLC4A11", "SLC13A2", "SLC22A4", "SLC29A1", "SLC39A13"),
#   Name = c("ABCA3", "ABCB1", "ABCC4", "SLC4A3", "SLC4A10", "SLC4A11", "SLC13A2", "SLC22A4", "SLC29A1", "SLC39A13")
# )
# 
# erg_alx1_specific_smts_cleaned = urchin_5_smts_final_cleaned %>% right_join(erg_alx1_specific_smts, by = "Name") %>% transmute(GeneID = GeneID.x, Name = Name)
# gene_specific_cell_stats = splg_alx1_gradient_expression_cells %>% get_cell_stats
# gene_specific_cell_stats %>% View()

pigment_specific_smts = tibble(
  GeneID = c("ABCB1", "ABCC5", "ABCG11", "SLC5A9", "SLC7A11", "SLC26A11", "SLC35A1", "SLC35A5", "SLC35B2", "SLC35F3", "SLC45A1", "SLC45A3", "SLC52A3"),
  Name = c("ABCB1", "ABCC5", "ABCG11", "SLC5A9", "SLC7A11", "SLC26A11", "SLC35A1", "SLC35A5", "SLC35B2", "SLC35F3", "SLC45A1", "SLC45A3", "SLC52A3")
)

pigment_specific_smts_cleaned = urchin_5_smts_final_cleaned %>% right_join(pigment_specific_smts, by = "Name") %>% transmute(GeneID = GeneID.x, Name = Name)
pigment_specific_smts_cleaned = add_row(
  pigment_specific_smts_cleaned,
  GeneID = c("gcm", "LOC588806"),
  Name = c("gcm", "LOC588806")
)

pks_plots = map(
  .x = pks_gradient_expression_cells,
  .f = function(x)(
    lab_DP(
      x, 
      pigment_specific_smts_cleaned,
      plot_title = "PKS expression co-gradients",
      F, 
      600,
      F,
      0,
      3, 
      0,
      cols = c(
        "grey",
        "blue"
        )
      )
    )
  )

gcm_plots = map(
  .x = gcm_gradient_expression_cells,
  .f = function(x)(
    lab_DP(
      x, 
      pigment_specific_smts_cleaned,
      plot_title = "GCM expression co-gradients",
      F, 
      600,
      F,
      0,
      3, 
      0,
      cols = c(
        "grey",
        "blue"
      )
    )
  )
)

dev.off()
pdf(onefile = T, file = "~/sp_transportomics/out.dump/pks_plots.pdf", height = 8, width = 12)
pks_plots[1:5]
dev.off()

dev.off()
pdf(onefile = T, file = "~/sp_transportomics/out.dump/gcm_plots.pdf", height = 8, width = 12)
gcm_plots[1:5]
dev.off()


germline_pos_cell_stats =
  map(
    .x = urchin_5_clustered_seurat_objs,
    .f = function(x)(
      get_num_pos_cells_multiple_genes(x,
                                       germline_markers,
                                       thresh = 0,
                                       rev = F
      )
    )
  ) %>% bind_rows()

get_positive_hits_for_geneset_specific_cells_for_cells_pos_for_specific_gene = function(seurat_obj, gene_set, specific_gene, specific_gene_name, thresh, relative){
  
  
  # seurat_obj = labelled_urchin_5_seurat_objs[[1]]
  # gene_set = available_genesets[[5]]
  # specific_gene = "Nanos2"
  # specific_gene_name = "Nanos2"
  # thresh = 0
  # relative = F
  
  specific_gene_filtered_cells = get_pos_cells(
    seurat_obj = seurat_obj,
    gene_id = specific_gene,
    thresh = thresh,
    rev = F
    )
  
  geneset_subset = subset(urchin_5_keys[[1]], GeneSet == gene_set)
  
  seurat_obj = subset(
    seurat_obj, 
    features = geneset_subset$GeneID,
    cells = specific_gene_filtered_cells$cells
  )
  
  cell_types = Idents(seurat_obj) %>%
    levels()
  
  stage = seurat_obj@meta.data$orig.ident %>%
    levels()
  
  num_all_genes_in_gene_set = seurat_obj %>% 
    rownames() %>% 
    length()
  
  # for one stage, create separate seurat_objs by cell type
  t1.b = map(
    .x = cell_types,
    .f = function(x)(
      subset(seurat_obj, idents = x)
    )
  )
  
  if(relative == T){
    assay_data_choice = "scale.data"
  } else{
    assay_data_choice = "counts"
  }
  # for one stage, for each cell type, get expression count matrix data
  t1.c = map(
    .x = t1.b,
    .f = function(x)(
      GetAssayData(
        x,
        slot = assay_data_choice
      )
    )
  )
  
  # for one stage, for each cell type, get cells for genes with any positive hits
  t2 = map(
    .x = t1.c,
    .f = function(x)(
      get_positive_hits(x, thresh = thresh)
    )
  )
  
  # get the number of smts present in each cell type with at least one transcript count
  t3 = map(
    .x = t2,
    .f = function(x)(
      length(x)
    )
  ) %>% 
    unlist()
  
  res = tibble(
    Stage = stage,
    "Cell type" = cell_types,
    "Gene set" = gene_set,
    "Number present" = t3,
    "Total existing" = num_all_genes_in_gene_set,
    "Gene filter condition" = specific_gene_name
  )
  
  return(res)
}


germline_number_smts_per_cell_type_each_stage = map(
  .x = urchin_5_clustered_labelled_seurat_objs[7:8],
  .f = function(x)(
    get_positive_hits_for_geneset_specific_cells_for_cells_pos_for_specific_gene(
      seurat_obj = x,
      gene_set = available_genesets[[5]],
      specific_gene = "Nanos2",
      specific_gene_name = "Nanos2",
      thresh = 1,
      relative = F)
  )
) %>% 
  bind_rows()

germline_number_smts_per_cell_type_each_stage %>% 
  group_by(Stage) %>% 
  summarise(sum(`Number present`))

ggplot(data = germline_pos_cell_stats) +
  geom_bar(aes(x = stage, y = num_cells, col = gene), stat = "identity", position = "dodge")

germline_markers_marked_cells = intersect(germline_pos_nanos2$cells, germline_pos_seawi$cells)
urchin_5_spmb_clustered_germline = subset(urchin_5_speb_clustered, cells = germline_markers_marked_cells)
