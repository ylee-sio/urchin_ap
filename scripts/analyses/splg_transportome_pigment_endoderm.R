urchin_5_splg_clustered_smts_genes_present =
  rownames(urchin_5_splg_clustered)[which((urchin_5_splg_clustered %>% rownames()) %in% urchin_5_smts_final_cleaned_with_mito_smts$GeneID)]

urchin_5_splg_clustered_smts_only =
  subset(
    urchin_5_splg_clustered,
    features = urchin_5_splg_clustered_smts_genes_present
  )

urchin_5_splg_clustered_smts_only_clustered =
  subcluster_suite(
    seurat_object = urchin_5_splg_clustered_smts_only,
    ndims = 1:20,
    strat = "tsne",
    res = 0.8,
    jackstraw = F,
    coord_strat = "ica"
  )

lab_DP(urchin_5_splg_clustered_smts_only_clustered,
       urchin_5_smts_final_cleaned_with_mito_smts,
       plot_title = "",
       F,
       3000,
       T,
       0,
       3,
       0,
       cols = c("grey",
                "blue")
)

get_cell_stats(urchin_5_splg_clustered_smts_only_clustered)

urchin_5_splg_smt_cluster_rep_data = find_new_cluster_cell_idents_in_old_cluster_cell_idents(urchin_5_splg_clustered_smts_only_clustered, urchin_5_splg_cell_type_labelled)
urchin_5_splg_smt_cluster_rep_data_sorted = map(
  urchin_5_splg_smt_cluster_rep_data,
  .f = function(x)(
    arrange(
      x,
      desc(cluster_num)
    )
  )
)
urchin_5_splg_cell_type_labelled_cell_stats = get_cell_stats(urchin_5_splg_cell_type_labelled)
urchin_5_splg_cell_type_labelled_cell_stats$cluster_num = urchin_5_splg_cell_type_labelled_cell_stats$cluster_num %>% as.character()
cells_stats_of_old_seurat_obj = urchin_5_splg_cell_type_labelled_cell_stats %>% arrange(desc(cluster_num))

get_percent_cells_represented = function(new_seurat_obj_cell_stats, old_seurat_obj_cell_stats){
  new_seurat_obj_cell_stats = new_seurat_obj_cell_stats %>% mutate(num_cells_new = num_cells)
  old_seurat_obj_cell_stats = old_seurat_obj_cell_stats %>% mutate(num_cells_old = num_cells)
  new_old_joined = left_join(new_seurat_obj_cell_stats, old_seurat_obj_cell_stats, by = "cluster_num")
  percent_repped = new_old_joined$num_cells_new/new_old_joined$num_cells_old
  cell_stats = mutate(new_seurat_obj_cell_stats, percent_representation = percent_repped, total_num_cells_in_cell_type = old_seurat_obj_cell_stats$num_cells)
  return(cell_stats)

}

rep_data = map(
  .x = urchin_5_splg_smt_cluster_rep_data_sorted,
  .f = function(x)(
    get_percent_cells_represented(
      new_seurat_obj_cell_stats = x,
      old_seurat_obj_cell_stats = cells_stats_of_old_seurat_obj
    )
  )
) %>%
  bind_rows()

fig = rep_data
fig = fig %>% plot_ly(x = ~cluster_num, y = ~percent_representation, color = ~cluster_origin)

fig %>%
  layout(
    yaxis = list(
      title = list(
        text = '% of Cells in Cell Identities from Cells in SMT activity clusters',
        font = list(size = 30)
      ),
      tickfont = list(size = 30)
    ),
    xaxis = list(
      title = list(
        text = 'Cell Identities',
        font = list(size = 30)
      ),
      tickfont = list(size = 30)
    ),
    legend = list(
      title = list(
        text = "<b> SMT Activity Clusters</b>",
        font = list(size = 25)
      ),
      font = list(
        size = 25
      )
    )
  )

urchin_5_splg_clustered_smts_only_clustered_pigment = subset(urchin_5_splg_clustered_smts_only_clustered, idents = c("10", "13"))
# urchin_5_splg_clustered_smts_only_clustered_pigment_common = subset(urchin_5_splg_clustered_smts_only_clustered, idents = c("7"))
urchin_5_splg_cell_type_labelled_pigment = subset(urchin_5_splg_cell_type_labelled, idents = "Pigment", cells = Cells(urchin_5_splg_clustered_smts_only_clustered_pigment))
# urchin_5_splg_cell_type_labelled_pigment_common = subset(urchin_5_splg_cell_type_labelled, idents = "Pigment", cells = Cells(urchin_5_splg_clustered_smts_only_clustered_pigment_common))

urchin_5_splg_clustered_smts_only_clustered_endoderm = subset(urchin_5_splg_clustered_smts_only_clustered, idents = c("3","10","11"))
# urchin_5_splg_clustered_smts_only_clustered_common_1 = subset(urchin_5_splg_clustered_smts_only_clustered, idents = c("3"))
# urchin_5_splg_clustered_smts_only_clustered_common_2 = subset(urchin_5_splg_clustered_smts_only_clustered, idents = c("7"))

urchin_5_splg_cell_type_labelled_endoderm = subset(urchin_5_splg_cell_type_labelled, idents = "Endoderm", cells = Cells(urchin_5_splg_clustered_smts_only_clustered_endoderm))
# urchin_5_splg_cell_type_labelled_endoderm_common_1 = subset(urchin_5_splg_cell_type_labelled, idents = "Endoderm", cells = Cells(urchin_5_splg_clustered_smts_only_clustered_common_1))
# urchin_5_splg_cell_type_labelled_endoderm_common_2 = subset(urchin_5_splg_cell_type_labelled, idents = "Endoderm", cells = Cells(urchin_5_splg_clustered_smts_only_clustered_common_2))

# b = get_pos_cells(urchin_5_splg_clustered, "merged-gene-175", thresh = 0, rev = F)
# urchin_5_splg_cell_type_labelled_endoderm = subset(urchin_5_splg_cell_type_labelled, cells = b$cells)
# urchin_5_splg_cell_type_labelled_endoderm_smts_only = 
#   subset(
#     urchin_5_splg_cell_type_labelled_endoderm,
#     features = urchin_5_splg_clustered_smts_genes_present
#   )

#find which cluster contains cells with SMT AP cluster 10 cell barcodes

#1 find SMT AP cluster 10 barcodes
smtap_cluster3_barcodes = 
  subset(urchin_5_splg_clustered_smts_only_clustered, idents = c("3")) %>% 
  Cells()

smtap_cluster10_barcodes = 
  subset(urchin_5_splg_clustered_smts_only_clustered, idents = c("10")) %>% 
  Cells()

smtap_cluster11_barcodes = 
  subset(urchin_5_splg_clustered_smts_only_clustered, idents = c("11")) %>% 
  Cells()

smtap_cluster13_barcodes = 
  subset(urchin_5_splg_clustered_smts_only_clustered, idents = c("13")) %>% 
  Cells()

#2 find all barcodes in both endoderm and pigment datasets

urchin_5_splg_clustered_endoderm_with_AP_cells = 
  subcluster_suite(
    seurat_object = urchin_5_splg_cell_type_labelled_endoderm, 
    ndims = 1:10, 
    strat = "tsne", 
    res = 0.4, 
    jackstraw = F,
    coord_strat = "ica"
  )
endoderm_barcodes = Cells(urchin_5_splg_clustered_endoderm_with_AP_cells) 

urchin_5_splg_clustered_pigment_with_AP_cells = 
  subcluster_suite(
    seurat_object = urchin_5_splg_cell_type_labelled_pigment, 
    ndims = 1:15, 
    strat = "tsne", 
    res = 0.4, 
    jackstraw = F,
    coord_strat = "ica"
  )

pigment_barcodes = Cells(urchin_5_splg_clustered_pigment_with_AP_cells)

#3 get cells from both datasets which match SMT AP cluster 10 barcodes
smtap_cluster3_barcodes_in_endoderm_barcodes = endoderm_barcodes[endoderm_barcodes %in% smtap_cluster3_barcodes %>% which()]
smtap_cluster3_barcodes_in_pigment_barcodes = pigment_barcodes[pigment_barcodes %in% smtap_cluster3_barcodes %>% which()]

smtap_cluster10_barcodes_in_endoderm_barcodes = endoderm_barcodes[endoderm_barcodes %in% smtap_cluster10_barcodes %>% which()]
smtap_cluster10_barcodes_in_pigment_barcodes = pigment_barcodes[pigment_barcodes %in% smtap_cluster10_barcodes %>% which()]

smtap_cluster11_barcodes_in_endoderm_barcodes = endoderm_barcodes[endoderm_barcodes %in% smtap_cluster11_barcodes %>% which()]
smtap_cluster11_barcodes_in_pigment_barcodes = pigment_barcodes[pigment_barcodes %in% smtap_cluster11_barcodes %>% which()]

smtap_cluster13_barcodes_in_endoderm_barcodes = endoderm_barcodes[endoderm_barcodes %in% smtap_cluster13_barcodes %>% which()]
smtap_cluster13_barcodes_in_pigment_barcodes = pigment_barcodes[pigment_barcodes %in% smtap_cluster13_barcodes %>% which()]
#4 get cluster number(s) of the cells containing those barcodes in each dataset
endoderm_cluster_nums_with_smtapcluster3_cells = Idents(urchin_5_splg_clustered_endoderm_with_AP_cells)[urchin_5_splg_clustered_endoderm_with_AP_cells %>% Idents() %>% names() %in% smtap_cluster3_barcodes_in_endoderm_barcodes %>% which()]
endoderm_cluster_nums_with_smtapcluster3_cells_cluster_percentages = map(.x = levels(endoderm_cluster_nums_with_smtapcluster3_cells), .f = function(x)(sum(endoderm_cluster_nums_with_smtapcluster3_cells == x)))
pigment_cluster_nums_with_smtapcluster3_cells = Idents(urchin_5_splg_clustered_pigment_with_AP_cells)[urchin_5_splg_clustered_pigment_with_AP_cells %>% Idents() %>% names() %in% smtap_cluster3_barcodes_in_pigment_barcodes %>% which()]
pigment_cluster_nums_with_smtapcluster3_cells_cluster_percentages = map(.x = levels(pigment_cluster_nums_with_smtapcluster3_cells), .f = function(x)(sum(pigment_cluster_nums_with_smtapcluster3_cells == x)))

endoderm_cluster_nums_with_smtapcluster10_cells = Idents(urchin_5_splg_clustered_endoderm_with_AP_cells)[urchin_5_splg_clustered_endoderm_with_AP_cells %>% Idents() %>% names() %in% smtap_cluster10_barcodes_in_endoderm_barcodes %>% which()]
endoderm_cluster_nums_with_smtapcluster10_cells_cluster_percentages = map(.x = levels(endoderm_cluster_nums_with_smtapcluster10_cells), .f = function(x)(sum(endoderm_cluster_nums_with_smtapcluster10_cells == x)))
pigment_cluster_nums_with_smtapcluster10_cells = Idents(urchin_5_splg_clustered_pigment_with_AP_cells)[urchin_5_splg_clustered_pigment_with_AP_cells %>% Idents() %>% names() %in% smtap_cluster10_barcodes_in_pigment_barcodes %>% which()]
pigment_cluster_nums_with_smtapcluster10_cells_cluster_percentages = map(.x = levels(pigment_cluster_nums_with_smtapcluster10_cells), .f = function(x)(sum(pigment_cluster_nums_with_smtapcluster10_cells == x)))

endoderm_cluster_nums_with_smtapcluster11_cells = Idents(urchin_5_splg_clustered_endoderm_with_AP_cells)[urchin_5_splg_clustered_endoderm_with_AP_cells %>% Idents() %>% names() %in% smtap_cluster11_barcodes_in_endoderm_barcodes %>% which()]
endoderm_cluster_nums_with_smtapcluster11_cells_cluster_percentages = map(.x = levels(endoderm_cluster_nums_with_smtapcluster11_cells), .f = function(x)(sum(endoderm_cluster_nums_with_smtapcluster11_cells == x)))
pigment_cluster_nums_with_smtapcluster11_cells = Idents(urchin_5_splg_clustered_pigment_with_AP_cells)[urchin_5_splg_clustered_pigment_with_AP_cells %>% Idents() %>% names() %in% smtap_cluster11_barcodes_in_pigment_barcodes %>% which()]
pigment_cluster_nums_with_smtapcluster11_cells_cluster_percentages = map(.x = levels(pigment_cluster_nums_with_smtapcluster11_cells), .f = function(x)(sum(pigment_cluster_nums_with_smtapcluster11_cells == x)))

endoderm_cluster_nums_with_smtapcluster13_cells = Idents(urchin_5_splg_clustered_endoderm_with_AP_cells)[urchin_5_splg_clustered_endoderm_with_AP_cells %>% Idents() %>% names() %in% smtap_cluster13_barcodes_in_endoderm_barcodes %>% which()]
endoderm_cluster_nums_with_smtapcluster13_cells_cluster_percentages = map(.x = levels(endoderm_cluster_nums_with_smtapcluster13_cells), .f = function(x)(sum(endoderm_cluster_nums_with_smtapcluster13_cells == x)))
pigment_cluster_nums_with_smtapcluster13_cells = Idents(urchin_5_splg_clustered_pigment_with_AP_cells)[urchin_5_splg_clustered_pigment_with_AP_cells %>% Idents() %>% names() %in% smtap_cluster13_barcodes_in_pigment_barcodes %>% which()]
pigment_cluster_nums_with_smtapcluster13_cells_cluster_percentages = map(.x = levels(pigment_cluster_nums_with_smtapcluster13_cells), .f = function(x)(sum(pigment_cluster_nums_with_smtapcluster13_cells == x)))

get_cell_stats(urchin_5_splg_clustered_endoderm_with_AP_cells)
get_cell_stats(urchin_5_splg_clustered_pigment_with_AP_cells)

urchin_5_splg_clustered_endoderm_with_AP_cells_renamed = RenameIdents(
  urchin_5_splg_clustered_endoderm_with_AP_cells,
  '1' = "Putative neural-endodermal cluster:\nConsists of 98.01% cells\n from SMT AP cluster 10",
  '0' = "1",
  '2' = "2",
  '3' = "3",
  '4' = "4"
)


urchin_5_splg_clustered_pigment_with_AP_cells_renamed = RenameIdents(
  urchin_5_splg_clustered_pigment_with_AP_cells,
  '3' = "Putative neural-pigment cluster:\nConsists of 98.32% cells\n from SMT AP cluster 10",
  '0' = "1",
  '1' = "2",
  '2' = "3"
)


t1 = lab_DP(urchin_5_splg_clustered_endoderm_with_AP_cells_renamed, 
            urchin_5_smts_final_cleaned %>% bind_rows(standard_markers_for_display),
            plot_title = "",
            F, 
            4000,
            T,
            0,
            3, 
            0,
            cols = c(
              "grey",
              "blue"
            )
)

t2 = lab_DP(urchin_5_splg_clustered_pigment_with_AP_cells_renamed, 
            urchin_5_smts_final_cleaned %>% bind_rows(standard_markers_for_display),
            plot_title = "",
            F, 
            4000,
            T,
            0,
            3, 
            0,
            cols = c(
              "grey",
              "blue"
            )
)

t3 = lab_DP(urchin_5_splg_clustered_endoderm_with_AP_cells_renamed, 
            common_genes_v2 %>% unique(),
            plot_title = "",
            F, 
            1100,
            F,
            0,
            3, 
            0,
            cols = c(
              "blue",
              "gold"
            )
) 

fig3a = t3 +
  ylab("Endoderm identity subclusters") +
  xlab("Genes") + 
  theme_bw(base_size = 22) + 
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))+
  theme(panel.background = element_rect(fill = "black"),
        panel.grid.major = element_line(color = 'grey', linetype = "dotted", size = 0.15)
  )

t3 %>% 
  layout(
    yaxis = list(
      title = list(
        text = 'Genes', 
        font = list(size = 25)
      ),
      tickfont = list(size = 25)
    ),
    xaxis = list(
      title = list(
        text = 'Endoderm subclusters', 
        font = list(size = 25)
      ),
      tickfont = list(size = 25)
    ),
    legend = list(
      title = list(
        text = "<b> SMT Activity Clusters</b>",
        font = list(size = 25)
      ),
      font = list(
        size = 25
      )
    )
  )

t4 = lab_DP(urchin_5_splg_clustered_pigment_with_AP_cells_renamed, 
            common_genes_v2 %>% unique(),
            plot_title = "",
            F, 
            1100,
            F,
            0,
            3, 
            0,
            cols = c(
              "blue",
              "gold"
            )
)

fig3b = t4 +
  ylab("Pigment identity subclusters") +
  xlab("Genes") + 
  theme_bw(base_size = 22) + 
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))+
  theme(panel.background = element_rect(fill = "black"),
        panel.grid.major = element_line(color = 'grey', linetype = "dotted", size = 0.15)
  )

t4 %>% 
  layout(
    yaxis = list(
      title = list(
        text = 'Genes', 
        font = list(size = 25)
      ),
      tickfont = list(size = 25)
    ),
    xaxis = list(
      title = list(
        text = 'Pigment subclusters', 
        font = list(size = 25)
      ),
      tickfont = list(size = 25)
    ),
    legend = list(
      title = list(
        text = "<b> SMT Activity Clusters</b>",
        font = list(size = 25)
      ),
      font = list(
        size = 25
      )
    )
  )

common_genes_v1 = tibble(
  
  GeneID = c("LOC588806","gcm", "Endo16", "GATAe", "blimp1/krox", "merged-gene-165", "Six1", "vasa",
             "merged-gene-2312", "merged-gene-2261", "merged-gene-2389","merged-gene-2383", "merged-gene-2397",
             "merged-gene-2376", "merged-gene-164", "merged-gene-150",
             "LOC580931","COLP3alpha", "merged-gene-1075", "merged-gene-1073","hh"),
  Name = c("PKS1","gcm", "Endo16", "GATAe", "blimp1/krox", "ABCC5", "Six1", "vasa",
           "SLC28A1", "SLC16A14", "SLC6A5","SLC5A7", "SLC7A9",
           "SLC4A3", "ABCC4", "ABCA1",
           "collagen alpha-1(I) chain","COLP3alpha", "frizzled-5", "frizzled-1","hh")
  
)

common_genes_v2 = tibble(
  
  GeneID = c("LOC588806","gcm","merged-gene-165", "Endo16", "vasa", "GATAe",
             "merged-gene-2312", "merged-gene-2261","merged-gene-2383", "merged-gene-2397",
             "merged-gene-2312", "merged-gene-164", "merged-gene-150",
             "LOC580931","COLP3alpha", "merged-gene-1075", "merged-gene-1073","hh",
             "merged-gene-1669", "LOC577764", "FoxA", "LOC752279", "LOC580689",
             "merged-gene-593", "merged-gene-591", "merged-gene-1245"),
  
  Name = c("PKS1","gcm","ABCC5", "Endo16",  "vasa","GATAe",
           "SLC28A1", "SLC16A14","SLC5A7", "SLC7A9",
           "SLC28A1", "ABCC4", "ABCA1",
           "collagen alpha-1(I) chain","COLP3alpha", "frizzled-5", "frizzled-1","hh",
           "neural cell adhesion molecule 1 (NCAM1)", "neuronal acetylcholine receptor subunit alpha-2", "FoxA", "patched domain-containing protein 2", "dispatched homolog 1 like protein",
           "CHST13 carbohydrate sulfotransferase 13","CHST11 carbohydrate sulfotransferase 11", "heparan-sulfate 6-O-sulfotransferase 3")
  
)



# Idents(urchin_5_splg_clustered_endoderm_with_AP_cells, cells = Cells(urchin_5_splg_clustered_smts_only_clustered_common_1)) = "SMT AP cluster 3"
# Idents(urchin_5_splg_clustered_endoderm_with_AP_cells, cells = Cells(urchin_5_splg_clustered_smts_only_clustered_common_2)) = "SMT AP cluster 7"

# urchin_5_splg_clustered_endoderm_with_AP_cells_renamed = 
#   RenameIdents(urchin_5_splg_clustered_endoderm_with_AP_cells,
#                '0' = "Endoderm 1",
#                '1' = "Endoderm AP cluster 2",
#                '2' = "Endoderm 2",
#                '3' = "Endoderm 3",
#                '4' = "Endoderm 4"
#                )



#3 and 4
t2 = lab_DP(urchin_5_splg_clustered_endoderm_with_AP_cells, 
            endoderm_pigment_selected_smts_v2 %>% bind_rows(standard_markers_for_display),
            plot_title = "",
            F, 
            1200,
            T,
            0,
            3, 
            0,
            cols = c(
              "grey",
              "blue"
            )
)


get_cell_stats(urchin_5_splg_clustered_endoderm_with_AP_cells)
patched = search_purp("patched")
sulfotransferase = search_purp("sulfotransferase")

pdx = search_purp("PD")
enter = search_purp("enter")
pks = search_purp("polyketide synthase 2")
intest = search_purp("intest")
smooth = search_purp("smooth")
cyclic = search_purp("cyclic")
lipid = search_purp("lipid")
ficolin = search_purp("ficolin")
bmp = search_purp("BMP")
cholest = search_purp("cholest")
wnt = search_purp("wnt")
neur = search_purp("neur")
enzymes = urchin_5_keys[[1]] %>% subset(GeneSet == "enzymes") %>% select(GeneID, Name)
pattern_recognition_receptors = urchin_5_keys[[1]] %>% subset(GeneSet == "pattern_recognition_receptors") %>% select(GeneID, Name)

membrane_trafficking_exocytosis = urchin_5_keys[[1]] %>% subset(GeneSet == "membrane_trafficking_exocytosis") %>% select(GeneID, Name)
membrane_receptors_gpcrs = urchin_5_keys[[1]] %>% subset(GeneSet == "membrane_receptors_gpcrs") %>% select(GeneID, Name)

wnt = search_purp("wnt") %>% mutate(X2 = GeneID)
a = urchin_5_keys[[1]] %>% right_join(wnt, by = "X2") 
a$GeneID.x[which(is.na(a$GeneID.x))] = a$X2[which(is.na(a$GeneID.x))]

a = a %>% 
  transmute(
    GeneID = GeneID.x,
    Name = Name.y
  )

wnt = wnt %>% 
  left_join(a, by = "Name") %>% 
  select(GeneID.y, Name) %>% 
  transmute(GeneID = GeneID.y, Name = Name)

urchin_5_splg_clustered_endoderm_with_AP_cells_markers = 
  FindAllMarkers(urchin_5_splg_clustered_endoderm_with_AP_cells)

lab_DP(urchin_5_splg_cell_type_labelled, 
       urchin_5_smts_final_cleaned,
       plot_title = "",
       F, 
       2500,
       T,
       0,
       1.5, 
       0,
       cols = c(
         "grey",
         "blue"
       )
)

temp_endoderm_genes = tibble(
  GeneID = c(
    "merged-gene-593",
    "merged-gene-493",
    "merged-gene-684",
    "merged-gene-608",
    "merged-gene-1667",
    "merged-gene-1364",
    "merged-gene-1363",
    "merged-gene-2019",
    "merged-gene-1966",
    "merged-gene-1967",
    "merged-gene-1968",
    "merged-gene-1948",
    "merged-gene-2459",
    "merged-gene-1115"
  ),
  Name = c(
    
    "merged-gene-593",
    "merged-gene-493",
    "merged-gene-684",
    "merged-gene-608",
    "merged-gene-1667",
    "merged-gene-1364",
    "merged-gene-1363",
    "merged-gene-2019",
    "merged-gene-1966",
    "merged-gene-1967",
    "merged-gene-1968",
    "merged-gene-1948",
    "merged-gene-2459",
    "merged-gene-1115"
  )
)


urchin_5_splg_clustered_endoderm_with_AP_cells_markers_top_markers = 
  urchin_5_splg_clustered_endoderm_with_AP_cells_markers %>% 
  group_by(cluster) %>% 
  arrange(desc(avg_log2FC)) %>% 
  top_n(40) %>%
  select(gene) %>% 
  unique()

b = tibble(
  GeneID = urchin_5_splg_clustered_endoderm_with_AP_cells_markers_top_markers$gene,
  Name = urchin_5_splg_clustered_endoderm_with_AP_cells_markers_top_markers$gene
) %>% unique()
