library(tidyverse)
library(Seurat)
library(plotly)
library(gridExtra)
library(cowplot)

sp_kegg_smts = read_csv("data/working_data/sp_kegg_smts.csv")


splg_expression_matrix = Read10X("data/unmodified/sp8_lg/")
speg_expression_matrix = Read10X("data/unmodified/sp7_eg/")

all_stages = readRDS("data/working_data/all_stages.clustered_seurat_obj_list.rds")
splg = all_stages[[8]]

isolate_geneset_in_clusters = function(expression_matrix, geneset) {
  
  filtrate_index = rownames(expression_matrix) %in% geneset$GeneID %>% which()
  new_expression_matrix = expression_matrix[filtrate_index,] > 0
  seurat_obj = CreateSeuratObject(counts = new_expression_matrix)
  return(seurat_obj)
  
}
get_min_top10_markers_each_cluster = function(seurat_obj, geneset){
  
  markers = FindAllMarkers(seurat_obj, only.pos = T)
  markers_top_10 = 
    markers %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = T) %>%
    top_n(n = 10) %>% 
    select(gene) %>% 
    mutate(GeneID = gene) %>% 
    left_join(geneset, by = "GeneID") %>% 
    mutate(Name = paste0(Name,"-", GeneID))
  
  markers_top_10_duplicates = markers_top_10$GeneID %>% duplicated() %>% which()
  markers_top_10_duplicates_removed = markers_top_10[-c(markers_top_10_duplicates),]
  return(markers_top_10_duplicates_removed)
  
}
get_cell_subsets_all_idents = function(seurat_obj){
  
  idents_list = Idents(seurat_obj) %>% levels()
  
  cell_subsets_by_idents = map(
    .x = idents_list,
    .f = function(x)(
      seurat_obj %>% 
        subset(idents = x) %>% 
        Cells()
    )
  )
  
  return(cell_subsets_by_idents)
  
}
cross_rep_dimplot_overlay = function(control_seurat_obj, geneset_specific_seurat_obj, cells_for_highlight_by_idents_list, base_title){
  
  idents_list = Idents(control_seurat_obj) %>% levels()
  
  cell_subsets_by_idents = map(
    .x = idents_list,
    .f = function(x)(
      control_seurat_obj %>% 
        subset(idents = x) %>% 
        Cells()
    )
  )
  
  geneset_dimplot_overlays = map2(
    .x = cell_subsets_by_idents,
    .y = idents_list,
    .f = function(x,y)(
      DimPlot(
        object = geneset_specific_seurat_obj,
        cells.highlight = x,
        sizes.highlight = 0.2,
        shuffle = T,
        label = T
        # label.box = T,
        # label.size = 3
      ) + 
        ggtitle(
          paste0(
            "Control cluster ",
            y,
            " cells/\n",
            base_title,
            " clusters overlay"
          )
        ) +
        NoLegend()
    )
  )
  
  # control_dimplot = DimPlot(control_seurat_obj, label = T) + ggtitle("Control")
  # treatment_dimplot = DimPlot(geneset_specific_seurat_obj, label = T) + ggtitle(paste0("Treatment: ", base_title))
  # 
  # base_plots = list(control_dimplot, treatment_dimplot)
  # final_plot_list = c(base_plots, geneset_dimplot_overlays)
  # return(final_plot_list)
  
  return(geneset_dimplot_overlays)
  
}

splg_smts_only = isolate_geneset_in_clusters(splg_expression_matrix, smts_in_scrnaseq_with_common_names_cleaned)

splg = subcluster_suite(
  seurat_object = splg,
  res = 1.0,
  ndims = 1:25
)

splg_smts_only = subcluster_suite(
  seurat_object = splg_smts_only,
  res = 1.0,
  ndims = 1:25
)

celltype_markers = tibble(
  GeneID = c(
    "Klf7", "NK2.2", "spec2c",
    "AnkAT-1", "foxq2", "NK2.1",
    "hnf6",
    "CHRD", "FoxG", "Lim1","LOC575591",
    "LOC577601","LOC576365", "LOC100889993",
    
    "blimp1/krox", "FoxA", "Endo16","LOC575684",  "LOC764728",
   
    
    "gcm", "Six1", "ABCC5a", "LOC584852", 

    "GATAc", "Alx1", "erg", "LOC752471","SPARC", "LOC593609", "LOC583799", "LOC576148",
    "Nanos2", "ago1",
    
     "vasa"
    ),
  Name = c(
    "Klf7", "NK2.2", "spec2c",
    "AnkAT-1", "foxq2", "NK2.1",
    "onecut",
    "CHRD", "FoxG", "Lim1","dmbx",
    "brn1/2/4","islet1", "dmrt3",
    
    "blimp1/krox", "FoxA", "Endo16","a1cf", "hnf1", 
    
    
    "gcm", "Six1", "ABCC5a", "strp4", 

    "GATAc", "Alx1", "erg", "calumenin", "SPARC/osteonectin", "tetraspanin-4", "MMP18","prox1",
    "Nanos2", "ago1",
    
    "vasa"
    
    )
)

labeled_dotplot(
  scrna_df = splg,
  feature_df = celltype_markers,
  col.min = 0,
  col.max = 3,
  interactive = T,
  plot_height = 1000,
  dot.min = 0,
  cols = c("blue", "gold"),
  plot_title = "",
  by_stage = F
)

# splg_labeled = RenameIdents(
#   splg,
#   '0' = "Aboral ectoderm/neural",
#   '1' = "Aboral ectoderm/neural",
#   '2' = "Hindgut",
#   '3' = "Cilliary band",
#   '4' = "Oral ectoderm",
#   '5' = "Aboral ectoderm/neural",
#   '6' = "Cilliary band",
#   '7' = "Oral ectoderm",
#   '8' = "Neural",
#   '9' = "Cilliary band",
#   '10' = "Mid-gut",
#   '11' = "Oral ectoderm",
#   '12' = "Aboral ectoderm/neural",
#   '13' = "Foregut",
#   '14' = "Foregut",
#   '15' = "Pigment",
#   '16' = "Aboral ectoderm/neural",
#   '17' = "Oral ectoderm",
#   '18' = "Foregut",
#   '19' = "Primary mesenchyme",
#   '20' = "Pigment",
#   '21' = "Secondary mesenchyme",
#   '22' = "Germline"
# )

splg_labeled = RenameIdents(
  splg,
  '0' = "Aboral ectoderm/neural",
  '1' = "Aboral ectoderm/neural",
  '8' = "Neural",
  '3' = "Cilliary band",
  '4' = "Oral ectoderm",
  '5' = "Aboral ectoderm/neural",
  '6' = "Cilliary band",
  '7' = "Oral ectoderm",
  '9' = "Cilliary band",
  '11' = "Oral ectoderm",
  '12' = "Aboral ectoderm/neural",
  '18' = "Foregut",
  '13' = "Foregut",
  '14' = "Foregut",
  '10' = "Mid-gut",
  '2' = "Hindgut",
  '15' = "Pigment",
  '16' = "Aboral ectoderm/neural",
  '17' = "Oral ectoderm",
  '19' = "Skeletal",
  '20' = "Pigment",
  '21' = "Skeletal",
  '22' = "Germline"
)

labeled_dotplot(
  scrna_df = splg_labeled,
  feature_df = celltype_markers,
  col.min = 0,
  col.max = 3,
  interactive = T,
  plot_height = 1000,
  dot.min = 0,
  cols = c("blue", "gold"),
  plot_title = "",
  by_stage = F
)

DimPlot(splg_labeled, pt.size = 0.25 , label = T, label.box = F, label.size = 3, order = T, repel = T)

splg_smts_overlay_dimplots_by_idents_list = cross_rep_dimplot_overlay(
  control_seurat_obj = splg_labeled,
  geneset_specific_seurat_obj = splg_smts_only,
  cells_for_highlight_by_idents_list = splg_cells_by_idents_list,
  base_title = "SMTs"
)

pdf(file = "output/geneset_overlays/smts/base_splg_smts_plots.pdf", width = 9, height = 4, onefile = T)
control_dimplot + smts_dimplot
dev.off()

pdf(file = "output/geneset_overlays/smts/splg_smts_overlay.pdf", width = 18, height = 18, onefile = T)
grid.arrange(grobs = splg_smts_overlay_dimplots_by_idents_list)
dev.off()


# labeling cell types

# celltype_markers = tibble(
#   GeneID = c(
#     "Klf7", "NK2.2", "spec2c",
#     "blimp1/krox", "FoxA", "Endo16",
#     "CHRD", "FoxG", "Lim1", "hnf6",
#     "AnkAT-1", "foxq2", "NK2.1",
#     "gcm", "Six1", "ABCC5a", 
#     "delta", "snail", "ago1",
#     "GATAc", "Alx1", "LOC592979", "LOC752471", "LOC584852", "LOC575684", "LOC580597", 
#     "LOC576365", "LOC100889993", "LOC590650", "LOC581729", "LOC100892692", "LOC575591", "LOC764728", "LOC577601"),
#   Name = c(
#     "Klf7", "NK2.2", "spec2c",
#     "blimp1/krox", "FoxA", "Endo16",
#     "CHRD", "FoxG", "Lim1","hnf6",
#     "AnkAT-1", "foxq2", "NK2.1",
#     "gcm", "Six1", "ABCC5a", 
#     "delta", "snail", "ago1",
#     "GATAc", "Alx1", "SLC5A11", "calumenin", "strp4", "a1cf", "SLC6A5",
#     "islet1", "dmrt3", "SLC16A14", "SoxE", "FoxF", "mbx", "hnf1", "brn1/2/4")
# )

labeled_dotplot(
  scrna_df = splg_labeled,
  feature_df = celltype_markers,
  col.min = 0,
  col.max = 3,
  interactive = F,
  plot_height = 800,
  dot.min = 0,
  cols = c("blue", "gold"),
  plot_title = "",
  by_stage = F
) +
  ylab("IP Clusters") +
  xlab("IP Genes") +
  ggtitle("Control: IP Clusters") +
  theme(axis.title = element_text(face="bold", size=20),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))+
  theme(plot.title = element_text(size = 22, face = "bold"))


#obtaining genes for top most expressed smts from each smt ap cluster
all_top_markers_smtap = FindAllMarkers(splg_smts_only, logfc.threshold = 1.0, only.pos = T, min.pct = 0.20)
all_top_markers_smtap_df_cleaned =
  all_top_markers_smtap %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC))

all_top_markers_smtap_df_cleaned_for_dp = 
  all_top_markers_smtap_df_cleaned %>% 
  transmute(GeneID = gene, Name = gene) %>% 
  select(GeneID, Name, cluster) %>% 
  as_tibble() %>% 
  unique()
all_top_markers_smtap_df_cleaned_for_dp$cluster = unfactor(all_top_markers_smtap_df_cleaned_for_dp$cluster)

filtered_top_markers_smtap = 
  merge(all_top_markers_smtap_df_cleaned_for_dp, sp_kegg_smts, by = "GeneID") %>% 
  select(GeneID, Name.y, cluster) %>% 
  mutate(Name = Name.y)

filtered_top_markers_smtap_ordered = 
  filtered_top_markers_smtap %>% 
  arrange(-desc(cluster))

filtered_top_markers_smtap_ordered_duplicated_index = filtered_top_markers_smtap_ordered$GeneID %>% duplicated() %>% which()

filtered_top_markers_smtap_ordered_final = 
  filtered_top_markers_smtap_ordered[-c(filtered_top_markers_smtap_ordered_duplicated_index),] %>% 
  mutate(cluster=as.numeric(cluster)) %>% 
  arrange(-desc(cluster)) %>% 
  left_join(smts_in_scrnaseq_with_common_names_cleaned, by="GeneID") %>% 
  unique() %>% 
  mutate(Name=Name.y.y)

write.csv(filtered_top_markers_smtap_ordered_final, "~/Projects/urchin_ap/data/working_data/filtered_top_markers_smtap_ordered_final.csv")
smtap_final_fig_input = read.csv("~/Projects/urchin_ap/data/working_data/filtered_top_markers_smtap_ordered_final.csv")
# smtap_final_fig_input = smtap_final_fig_input %>% mutate(Name = paste0(description, "-", "(",GeneID.x,") ", Name.y))
smt_ap_dp = 
  labeled_dotplot(
    scrna_df = splg_smts_only,
    feature_df = smtap_final_fig_input,
    col.min = 0,
    col.max = 3,
    interactive = F,
    plot_height = 4000,
    dot.min = 0,
    cols = c("blue", "gold"),
    plot_title = "",
    by_stage = F
  ) +
    ylab("SMT AP Clusters") +
    xlab("SMT AP Genes") +
    ggtitle("Treatment: SMTAP Clusters") +
    theme(
      axis.title = element_text(face="bold", size=20),
      axis.text.x = element_text(size=16),
      axis.text.y = element_text(size=16)) +
  theme(plot.title = element_text(size = 22, face = "bold"))

ggsave2(smt_ap_dp,width = 24, height = 20, filename = "~/Projects/urchin_ap/output/smt_ap_dp.pdf", )


#*** cross cluster type representation analysis
splg_smts_cross_representation = find_new_cluster_cell_idents_in_old_cluster_cell_idents(splg_smts_only, splg_labeled)

splg_smts_cross_representation$cluster_origin = factor(splg_smts_cross_representation$cluster_origin, 
                                                       levels = splg_smts_cross_representation$cluster_origin %>% unique())
splg_smts_cross_representation = splg_smts_cross_representation %>% mutate("Cell Identity"=cluster_num, percent_representation=percent_representation*100)
cell_ident_names=Idents(splg_labeled) %>% levels()

#FIG2 free scale 2800X2000
ggplot(data=splg_smts_cross_representation) +
  geom_bar(mapping = aes(x=cluster_origin, y=percent_representation), stat = "identity") +
  # geom_text(aes(label = "hello"), vjust = 0) +
  facet_wrap(~cluster_num, nrow = 4, ncol = 4, scales = "free") +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  
  ylab("Percentages of IP cluster cells \n mapped to SMT AP cluster cells") +
  xlab("SMT AP Cluster") +
  theme(strip.text = element_text(face="bold", size=22),
        strip.background = element_rect(colour="black",size=1),
        axis.title = element_text(face="bold", size=22),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18)) +
  geom_text(aes(x=cluster_origin, y=percent_representation, label=paste0(num_cells)), angle=90, size = 6, hjust=-0.1)

cross_representation_plot(splg_smts_cross_representation)


# cross_representation_plot(splg_cross_representation_with_labeled_clusters)

# skeletal: post cross representation analysis


skeletal_smt_ap_cells = subset(splg_smts_only, idents = c("8","12")) %>% Cells()
splg_mapped_skeletal_smt_ap_cells = subset(splg, cells = skeletal_smt_ap_cells)

splg_mapped_skeletal_smt_ap_cells = subcluster_suite(
  seurat_object = splg_mapped_skeletal_smt_ap_cells,
  res = 0.1,
  ndims = 1:4
)

skeletal_geneset = tibble(
  
  GeneID = c(
    # canonical PMC markers set 1
    "SM50", "SM30A", "msp130", "LOC579101", "LOC579173", "Alx1", 
    
    # highly enriched markers found 
    "LOC576560", "LOC576018", "LOC583799","LOC586737", "LOC586753", "LOC577128",
    # canonical PMC markers set 2
     "LOC575628", "LOC100892657", "LOC582336",
    # PMC and SMC markers
    "SFK7","LOC574611", "SPARC", "LOC752471",
    
    # differentially expressed collagens
    "LOC589794", "COLP1alpha", "COLP2alpha","COLP3alpha", "COLP4alpha", 
    # IP Markers
    "LOC100887885", "LOC100892373", "LOC575826", "LOC580069", "GATAc", "erg","LOC586389","LOC592057",
    # SMTs which have defined these two primary subpopulation of skeletal cells
    "LOC576148", "LOC591586", "LOC585510", "LOC592979", "LOC587674",
    #additional
    "calret", "LOC575365","LOC594798", "LOC580975", "adv", "LOC577597", "LOC752782", 
    "LOC593355"),
  
  Name = c(
    # canonical PMC markers set 1
    "SM50", "SM30A", "msp130", "carbonic anhydrase 2", "otop2l", "Alx1",
    
    # highly enriched markers found 
    "neurexin-4", "Protein identified and sequenced from purified spicules in Mann et al. 2010: matrilysin (MMP18/19L3)",  "Protein identified and sequenced from purified spicules in Mann et al. 2010: MMP18 (MMP18/19L5)", "Protein identified and sequenced from purified spicules in Mann et al. 2010: macrophage metalloelastase (Sp-Mt1-4/MmpL5)","Protein identified and sequenced from purified spicules in Mann et al. 2010: macrophage metalloelastase (Sp-Mt1-4/MmpL6)", "Protein identified and sequenced from purified spicules in Mann et al. 2010: matrix metalloproteinase-24 (Sp-Mt5/MmpL2)", 
    # canonical PMC markers set 2
     "cyclophilin 1", "cyclophilin Livingston et al. 2006", "cyclophilin Livingston et al. 2006",
    # PMC and SMC markers
    "SFK7", "cofilin", "SPARC/osteonectin",  "calumenin", 
    
    # differentially expressed collagens
    "collagen alpha-1(V) chain", "COLP1alpha", "COLP2alpha", "COLP3alpha", "COLP4alpha", 
    # IP Markers
    "Kirrel2L_4", "KirreL2L_6", "hex",  "Alx4","GATAc","erg","tbr","HesC",
    # SMTs which have defined these two primary subpopulation of skeletal cells
    "prox1", "SLC9A2", "SLC13A2", "SLC5A11", "SLC26A5",
    #addtional
    "calret", "calmodulin", "alpha-2-macroglobulin", "neprilysin-4", "adv",
    "zinc metalloproteinase nas-13", "actin", "Calumenin-A")
  
)

labeled_dotplot(
  scrna_df = splg_mapped_skeletal_smt_ap_cells,
  # scrna_df = splg_labeled,
  feature_df = skeletal_geneset,
  col.min = 0,
  col.max = 4,
  interactive = T,
  plot_height = 1000,
  dot.min = 0,
  cols = c("blue", "gold"),
  plot_title = "",
  by_stage = F
)

 
labeled_dotplot(
  scrna_df = splg_mapped_skeletal_smt_ap_cells,
  # scrna_df = splg,
  feature_df = skeletal_geneset,
  col.min = 0,
  col.max = 2,
  interactive = T,
  plot_height = 1000,
  dot.min = 0,
  cols = c("blue", "gold"),
  plot_title = "",
  by_stage = F
) +
  ylab("Skeletogenic-SMT AP clusters") +
  xlab("Skeletogenic IP and AP markers\nand effector genes") +
  ggtitle("Treatment: Subclusters of SMT AP cells mapped to the skeletal IP cluster") +
  theme(axis.title = element_text(face="bold", size=20),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))+
  theme(plot.title = element_text(size = 22, face = "bold"))

labeled_dotplot(
  scrna_df = splg,
  feature_df = skeletal_geneset,
  col.min = 0,
  col.max = 3,
  interactive = F,
  plot_height = 3000,
  dot.min = 0,
  cols = c("blue", "gold"),
  plot_title = "",
  by_stage = F
) +
  ylab("Skeletogenic-SMT AP clusters") +
  xlab("Skeletogenic IP and AP markers\nand effector genes") +
  ggtitle("Treatment: Subclusters of SMT AP cells mapped to the skeletal IP cluster") +
  theme(axis.title = element_text(face="bold", size=20),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))+
  theme(plot.title = element_text(size = 22, face = "bold"))





top_smtap_sk_markers = FindAllMarkers(splg_mapped_skeletal_smt_ap_cells, only.pos = T)
get_cell_stats(splg_mapped_skeletal_smt_ap_cells)

top_smtap_sk_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>% 
  slice(1:50) %>% 
  View()

#DEAF1, LOC100893739
FeaturePlot(splg, features="LOC763226")

