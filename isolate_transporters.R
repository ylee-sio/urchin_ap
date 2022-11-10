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
speg = all_stages[[7]]

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

# splg_smts_only = isolate_geneset_in_clusters(splg_expression_matrix, sp_kegg_smts)
splg_smts_only = isolate_geneset_in_clusters(splg_expression_matrix, smts_in_scrnaseq)
speg_smts_only = isolate_geneset_in_clusters(speg_expression_matrix, sp_kegg_smts)
splg_gpcrs_only = isolate_geneset_in_clusters(splg_expression_matrix, sp_kegg_gpcr)
speg_gpcrs_only = isolate_geneset_in_clusters(speg_expression_matrix, sp_kegg_gpcr)
splg_tfs_only = isolate_geneset_in_clusters(splg_expression_matrix, sp_kegg_tfs)
speg_tfs_only = isolate_geneset_in_clusters(speg_expression_matrix, sp_kegg_tfs)
splg_exo_only = isolate_geneset_in_clusters(splg_expression_matrix, sp_kegg_exo)
speg_exo_only = isolate_geneset_in_clusters(speg_expression_matrix, sp_kegg_exo)
splg_mt_only = isolate_geneset_in_clusters(splg_expression_matrix, sp_kegg_mt)
speg_mt_only = isolate_geneset_in_clusters(speg_expression_matrix, sp_kegg_mt)
splg_enzymes_only = isolate_geneset_in_clusters(splg_expression_matrix, sp_kegg_enzymes)
speg_enzymes_only = isolate_geneset_in_clusters(speg_expression_matrix, sp_kegg_enzymes)
splg_splice_only = isolate_geneset_in_clusters(splg_expression_matrix, sp_kegg_splice)
speg_splice_only = isolate_geneset_in_clusters(speg_expression_matrix, sp_kegg_splice)
splg_pk_only = isolate_geneset_in_clusters(splg_expression_matrix, sp_kegg_pk)
speg_pk_only = isolate_geneset_in_clusters(speg_expression_matrix, sp_kegg_pk)

splg = subcluster_suite(
  seurat_object = splg,
  res = 1.0,
  ndims = 1:25
)

speg = subcluster_suite(
  seurat_object = speg,
  res = 0.5,
  ndims = 1:15
)

splg_pk_only = subcluster_suite(
  seurat_object = splg_pk_only,
  res = 0.5,
  ndims = 1:15
)

speg_pk_only = subcluster_suite(
  seurat_object = speg_pk_only,
  res = 0.5,
  ndims = 1:15
)

splg_splice_only = subcluster_suite(
  seurat_object = splg_splice_only,
  res = 0.5,
  ndims = 1:15
)

speg_splice_only = subcluster_suite(
  seurat_object = speg_splice_only,
  res = 0.5,
  ndims = 1:15
)

splg_enzymes_only = subcluster_suite(
  seurat_object = splg_enzymes_only,
  res = 0.5,
  ndims = 1:15
)

speg_enzymes_only = subcluster_suite(
  seurat_object = speg_enzymes_only,
  res = 0.5,
  ndims = 1:15
)

splg_mt_only = subcluster_suite(
  seurat_object = splg_mt_only,
  res = 0.5,
  ndims = 1:15
)

speg_mt_only = subcluster_suite(
  seurat_object = speg_mt_only,
  res = 0.5,
  ndims = 1:15
)

splg_exo_only = subcluster_suite(
  seurat_object = splg_exo_only,
  res = 0.5,
  ndims = 1:15
)

speg_exo_only = subcluster_suite(
  seurat_object = speg_exo_only,
  res = 0.5,
  ndims = 1:15
)

# splg_ics_only = subcluster_suite(
#   seurat_object = splg_ics_only,
#   res = 0.5,
#   ndims = 1:10
# )
# 
# speg_ics_only = subcluster_suite(
#   seurat_object = speg_ics_only,
#   res = 0.5,
#   ndims = 1:10
# )

splg_tfs_only = subcluster_suite(
  seurat_object = splg_tfs_only,
  res = 0.5,
  ndims = 1:15
)

speg_tfs_only = subcluster_suite(
  seurat_object = speg_tfs_only,
  res = 0.5,
  ndims = 1:15
)

# splg_gpcrs_only = subcluster_suite(
#   seurat_object = splg_gpcrs_only,
#   res = 0.5,
#   ndims = 1:15
# )
# 
# speg_gpcrs_only = subcluster_suite(
#   seurat_object = speg_gpcrs_only,
#   res = 0.5,
#   ndims = 1:15
# )

splg_smts_only = subcluster_suite(
  seurat_object = splg_smts_only,
  res = 1.0,
  ndims = 1:25
)

speg_smts_only = subcluster_suite(
  seurat_object = speg_smts_only,
  res = 0.5,
  ndims = 1:15
)

control_dimplot = DimPlot(splg_labeled, label = T, pt.size = 0.1, reduction="umap") + ggtitle("Control")
exo_dimplot = DimPlot(splg_exo_only, label = T, pt.size = 0.1) + ggtitle("Exosome related proteins")
tfs_dimplot = DimPlot(splg_tfs_only, label = T, pt.size = 0.1) + ggtitle("Transcription factors")
smts_dimplot = DimPlot(splg_smts_only, label = T, pt.size = 0.1, reduction="tsne") + ggtitle("SMTs")
mt_dimplot = DimPlot(splg_mt_only, label = T, pt.size = 0.1) + ggtitle("Membrane trafficking related proteins")
enzymes_dimplot = DimPlot(splg_enzymes_only, label = T, pt.size = 0.1) + ggtitle("Enzymes")
splice_dimplot = DimPlot(splg_splice_only, label = T, pt.size = 0.1) + ggtitle("Spliceosome proteins")
pk_dimplot = DimPlot(splg_pk_only, label = T, pt.size = 0.1) + ggtitle("Protein kinases")


#***ENZYMES
splg_enzymes_overlay_dimplots_by_idents_list = cross_rep_dimplot_overlay(
  control_seurat_obj = splg,
  geneset_specific_seurat_obj = splg_enzymes_only,
  cells_for_highlight_by_idents_list = splg_cells_by_idents_list,
  base_title = "enzyme"
  )

pdf(file = "output/geneset_overlays/enzymes/base_splg_enzymes_plots.pdf", width = 9, height = 4, onefile = T)
control_dimplot + enzymes_dimplot
dev.off()

pdf(file = "output/geneset_overlays/enzymes/splg_enzyme_overlay.pdf", width = 18, height = 18, onefile = T)
grid.arrange(grobs = splg_enzymes_overlay_dimplots_by_idents_list)
dev.off()

#***EXOSOME
splg_exo_overlay_dimplots_by_idents_list = cross_rep_dimplot_overlay(
  control_seurat_obj = splg,
  geneset_specific_seurat_obj = splg_exo_only,
  cells_for_highlight_by_idents_list = splg_cells_by_idents_list,
  base_title = "Exosome"
)

pdf(file = "output/geneset_overlays/exosome/base_splg_exo_plots.pdf", width = 9, height = 4, onefile = T)
control_dimplot + exo_dimplot
dev.off()

pdf(file = "output/geneset_overlays/exosome/splg_enzyme_overlay.pdf", width = 18, height = 18, onefile = T)
grid.arrange(grobs = splg_exo_overlay_dimplots_by_idents_list)
dev.off()

#***TFs
splg_tfs_overlay_dimplots_by_idents_list = cross_rep_dimplot_overlay(
  control_seurat_obj = splg,
  geneset_specific_seurat_obj = splg_tfs_only,
  cells_for_highlight_by_idents_list = splg_cells_by_idents_list,
  base_title = "transcription factors"
)

pdf(file = "output/geneset_overlays/tfs/base_splg_tfs_plots.pdf", width = 9, height = 4, onefile = T)
control_dimplot + tfs_dimplot
dev.off()

pdf(file = "output/geneset_overlays/tfs/splg_tfs_overlay.pdf", width = 18, height = 18, onefile = T)
grid.arrange(grobs = splg_tfs_overlay_dimplots_by_idents_list)
dev.off()

#***SMTs
splg_smts_overlay_dimplots_by_idents_list = cross_rep_dimplot_overlay(
  control_seurat_obj = splg_labeled,
  geneset_specific_seurat_obj = splg_smts_only,
  cells_for_highlight_by_idents_list = splg_cells_by_idents_list,
  base_title = "SMTs"
)

pdf(file = "output/geneset_overlays/smts/base_splg_smts_plots.pdf", width = 9, height = 4, onefile = T)
control_dimplot + smts_dimplot
dev.off()

library(gridExtra)
pdf(file = "output/geneset_overlays/smts/splg_smts_overlay.pdf", width = 18, height = 18, onefile = T)
grid.arrange(grobs = splg_smts_overlay_dimplots_by_idents_list)
dev.off()

#***MTs
splg_mt_overlay_dimplots_by_idents_list = cross_rep_dimplot_overlay(
  control_seurat_obj = splg,
  geneset_specific_seurat_obj = splg_mt_only,
  cells_for_highlight_by_idents_list = splg_cells_by_idents_list,
  base_title = "membrane trafficking"
)

pdf(file = "output/geneset_overlays/membrane_trafficking/base_splg_mt_plots.pdf", width = 9, height = 4, onefile = T)
control_dimplot + mt_dimplot
dev.off()

pdf(file = "output/geneset_overlays/membrane_trafficking/splg_mt_overlay.pdf", width = 18, height = 18, onefile = T)
grid.arrange(grobs = splg_mt_overlay_dimplots_by_idents_list)
dev.off()

#***splice
splg_splice_overlay_dimplots_by_idents_list = cross_rep_dimplot_overlay(
  control_seurat_obj = splg,
  geneset_specific_seurat_obj = splg_splice_only,
  cells_for_highlight_by_idents_list = splg_cells_by_idents_list,
  base_title = "spliceosome"
)

pdf(file = "output/geneset_overlays/splice/base_splg_splice_plots.pdf", width = 9, height = 4, onefile = T)
control_dimplot + splice_dimplot
dev.off()

pdf(file = "output/geneset_overlays/splice/splg_splice_overlay.pdf", width = 18, height = 18, onefile = T)
grid.arrange(grobs = splg_splice_overlay_dimplots_by_idents_list)
dev.off()

#***pk
splg_pk_overlay_dimplots_by_idents_list = cross_rep_dimplot_overlay(
  control_seurat_obj = splg,
  geneset_specific_seurat_obj = splg_pk_only,
  cells_for_highlight_by_idents_list = splg_cells_by_idents_list,
  base_title = "pkosome"
)

pdf(file = "output/geneset_overlays/pk/base_splg_pk_plots.pdf", width = 9, height = 4, onefile = T)
control_dimplot + pk_dimplot
dev.off()

pdf(file = "output/geneset_overlays/pk/splg_pk_overlay.pdf", width = 18, height = 18, onefile = T)
grid.arrange(grobs = splg_pk_overlay_dimplots_by_idents_list)
dev.off()


# ****************************************** SPEG

control_dimplot = DimPlot(speg, label = T, pt.size = 0.1) + ggtitle("Control")
exo_dimplot = DimPlot(speg_exo_only, label = T, pt.size = 0.1) + ggtitle("Exosome related proteins")
tfs_dimplot = DimPlot(speg_tfs_only, label = T, pt.size = 0.1) + ggtitle("Transcription factors")
smts_dimplot = DimPlot(splg_smts_only, label = T, pt.size = 0.1) + ggtitle("SMTs")
mt_dimplot = DimPlot(speg_mt_only, label = T, pt.size = 0.1) + ggtitle("Membrane trafficking related proteins")
enzymes_dimplot = DimPlot(speg_enzymes_only, label = T, pt.size = 0.1) + ggtitle("Enzymes")
splice_dimplot = DimPlot(speg_splice_only, label = T, pt.size = 0.1) + ggtitle("Spliceosome proteins")
pk_dimplot = DimPlot(speg_pk_only, label = T, pt.size = 0.1) + ggtitle("Protein kinases")


#***ENZYMES
speg_enzymes_overlay_dimplots_by_idents_list = cross_rep_dimplot_overlay(
  control_seurat_obj = speg,
  geneset_specific_seurat_obj = speg_enzymes_only,
  cells_for_highlight_by_idents_list = speg_cells_by_idents_list,
  base_title = "enzyme"
)

pdf(file = "output/geneset_overlays/enzymes/base_speg_enzymes_plots.pdf", width = 9, height = 4, onefile = T)
control_dimplot + enzymes_dimplot
dev.off()

pdf(file = "output/geneset_overlays/enzymes/speg_enzyme_overlay.pdf", width = 18, height = 18, onefile = T)
grid.arrange(grobs = speg_enzymes_overlay_dimplots_by_idents_list)
dev.off()

#***EXOSOME
speg_exo_overlay_dimplots_by_idents_list = cross_rep_dimplot_overlay(
  control_seurat_obj = speg,
  geneset_specific_seurat_obj = speg_exo_only,
  cells_for_highlight_by_idents_list = speg_cells_by_idents_list,
  base_title = "Exosome"
)

pdf(file = "output/geneset_overlays/exosome/base_speg_exo_plots.pdf", width = 9, height = 4, onefile = T)
control_dimplot + exo_dimplot
dev.off()

pdf(file = "output/geneset_overlays/exosome/speg_enzyme_overlay.pdf", width = 18, height = 18, onefile = T)
grid.arrange(grobs = speg_exo_overlay_dimplots_by_idents_list)
dev.off()

#***TFs
speg_tfs_overlay_dimplots_by_idents_list = cross_rep_dimplot_overlay(
  control_seurat_obj = speg,
  geneset_specific_seurat_obj = speg_tfs_only,
  cells_for_highlight_by_idents_list = speg_cells_by_idents_list,
  base_title = "transcription factors"
)

pdf(file = "output/geneset_overlays/tfs/base_speg_tfs_plots.pdf", width = 9, height = 4, onefile = T)
control_dimplot + tfs_dimplot
dev.off()

pdf(file = "output/geneset_overlays/tfs/speg_tfs_overlay.pdf", width = 18, height = 18, onefile = T)
grid.arrange(grobs = speg_tfs_overlay_dimplots_by_idents_list)
dev.off()

#***SMTs
speg_smts_overlay_dimplots_by_idents_list = cross_rep_dimplot_overlay(
  control_seurat_obj = speg,
  geneset_specific_seurat_obj = speg_smts_only,
  cells_for_highlight_by_idents_list = speg_cells_by_idents_list,
  base_title = "SMTs"
)

pdf(file = "output/geneset_overlays/smts/base_speg_smts_plots.pdf", width = 9, height = 4, onefile = T)
control_dimplot + smts_dimplot
dev.off()

pdf(file = "output/geneset_overlays/smts/speg_smts_overlay.pdf", width = 18, height = 18, onefile = T)
grid.arrange(grobs = speg_smts_overlay_dimplots_by_idents_list)
dev.off()

#***MTs
speg_mt_overlay_dimplots_by_idents_list = cross_rep_dimplot_overlay(
  control_seurat_obj = speg,
  geneset_specific_seurat_obj = speg_mt_only,
  cells_for_highlight_by_idents_list = speg_cells_by_idents_list,
  base_title = "membrane trafficking"
)

pdf(file = "output/geneset_overlays/membrane_trafficking/base_speg_mt_plots.pdf", width = 9, height = 4, onefile = T)
control_dimplot + mt_dimplot
dev.off()

pdf(file = "output/geneset_overlays/membrane_trafficking/speg_mt_overlay.pdf", width = 18, height = 18, onefile = T)
grid.arrange(grobs = speg_mt_overlay_dimplots_by_idents_list)
dev.off()

#***splice
speg_splice_overlay_dimplots_by_idents_list = cross_rep_dimplot_overlay(
  control_seurat_obj = speg,
  geneset_specific_seurat_obj = speg_splice_only,
  cells_for_highlight_by_idents_list = speg_cells_by_idents_list,
  base_title = "spliceosome"
)

pdf(file = "output/geneset_overlays/splice/base_speg_splice_plots.pdf", width = 9, height = 4, onefile = T)
control_dimplot + splice_dimplot
dev.off()

pdf(file = "output/geneset_overlays/splice/speg_splice_overlay.pdf", width = 18, height = 18, onefile = T)
grid.arrange(grobs = speg_splice_overlay_dimplots_by_idents_list)
dev.off()

#***pk
speg_pk_overlay_dimplots_by_idents_list = cross_rep_dimplot_overlay(
  control_seurat_obj = speg,
  geneset_specific_seurat_obj = speg_pk_only,
  cells_for_highlight_by_idents_list = speg_cells_by_idents_list,
  base_title = "protein kinases"
)

pdf(file = "output/geneset_overlays/pk/base_speg_pk_plots.pdf", width = 9, height = 4, onefile = T)
control_dimplot + pk_dimplot
dev.off()

pdf(file = "output/geneset_overlays/pk/speg_pk_overlay.pdf", width = 18, height = 18, onefile = T)
grid.arrange(grobs = speg_pk_overlay_dimplots_by_idents_list)
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

celltype_markers = tibble(
  GeneID = c(
    "Klf7", "NK2.2", "spec2c",
    "blimp1/krox", "FoxA", "Endo16",
    "CHRD", "FoxG", "Lim1", "hnf6",
    "AnkAT-1", "foxq2", "NK2.1",
    "gcm", "Six1", "ABCC5a", 
    "delta", "snail", "ago1",
    "GATAc", "Alx1", "erg"),
  Name = c(
    "Klf7", "NK2.2", "spec2c",
    "blimp1/krox", "FoxA", "Endo16",
    "CHRD", "FoxG", "Lim1","hnf6",
    "AnkAT-1", "foxq2", "NK2.1",
    "gcm", "Six1", "ABCC5a", 
    "delta", "snail", "ago1",
    "GATAc", "Alx1", "erg")
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

splg_labeled = RenameIdents(
  splg,
  '0' = "Aboral ectoderm/neural",
  '1' = "Aboral ectoderm/neural",
  '2' = "Endoderm",
  '3' = "Cilliary band",
  '4' = "Oral ectoderm",
  '5' = "Aboral ectoderm/neural",
  '6' = "Cilliary band",
  '7' = "Oral ectoderm",
  '8' = "Neural",
  '9' = "Cilliary band",
  '10' = "Endoderm",
  '11' = "Oral ectoderm",
  '12' = "Aboral ectoderm/neural",
  '13' = "Endoderm",
  '14' = "Endoderm",
  '15' = "Pigment",
  '16' = "Aboral ectoderm/neural",
  '17' = "Endoderm",
  '18' = "Secondary mesenchyme",
  '19' = "Skeletal",
  '20' = "Secondary mesenchyme",
  '21' = "Skeletal",
  '22' = "Secondary mesenchyme"
)

DimPlot(splg_labeled, label=T, label.size = 2.5, label.box = T, repel = T, pt.size = 0.05) + 
  ggtitle("Control: IP Clusters")

DimPlot(splg_smts_only, label=T, label.size = 2.5, label.box = T, repel = T, pt.size = 0.05) +
  ggtitle("Treatment: SMT AP Clusters")

#FIG2
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
all_top_markers_smtap = FindAllMarkers(splg_smts_only, logfc.threshold = 1.0, only.pos = T, min.pct = 0.25)
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
  left_join(smts_in_scrnaseq, by="GeneID") %>% 
  unique() %>% 
  mutate(Name=Name.y.y)

write.csv(filtered_top_markers_smtap_ordered_final, "~/Projects/urchin_ap/data/working_data/filtered_top_markers_smtap_ordered_final.csv")
smtap_final_fig_input = read.csv("~/Projects/urchin_ap/data/working_data/filtered_top_markers_smtap_ordered_final.csv")
smtap_final_fig_input = smtap_final_fig_input %>% mutate(Name = paste0(description, "-", "(",GeneID.x,") ", Name.y))
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
  theme(axis.title = element_text(face="bold", size=20),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))+
  theme(plot.title = element_text(size = 22, face = "bold"))



#*** cross cluster type representation analysis
splg_enzyme_cross_representation = find_new_cluster_cell_idents_in_old_cluster_cell_idents(splg_enzymes_only, splg_labeled)
splg_smts_cross_representation = find_new_cluster_cell_idents_in_old_cluster_cell_idents(splg_smts_only, splg_labeled)
splg_exo_cross_representation = find_new_cluster_cell_idents_in_old_cluster_cell_idents(splg_exo_only, splg_labeled)
splg_mt_cross_representation = find_new_cluster_cell_idents_in_old_cluster_cell_idents(splg_mt_only, splg_labeled)
splg_tfs_cross_representation = find_new_cluster_cell_idents_in_old_cluster_cell_idents(splg_tfs_only, splg_labeled)

# speg_enzyme_cross_representation = find_new_cluster_cell_idents_in_old_cluster_cell_idents(speg_enzymes_only, speg_labeled)
# speg_smts_cross_representation = find_new_cluster_cell_idents_in_old_cluster_cell_idents(speg_smts_only, speg_labeled)
# speg_exo_cross_representation = find_new_cluster_cell_idents_in_old_cluster_cell_idents(speg_exo_only, speg_labeled)
# speg_mt_cross_representation = find_new_cluster_cell_idents_in_old_cluster_cell_idents(speg_mt_only, speg_labeled)
# speg_tfs_cross_representation = find_new_cluster_cell_idents_in_old_cluster_cell_idents(speg_tfs_only, speg_labeled)


# splg_cross_representation_with_labeled_clusters = find_new_cluster_cell_idents_in_old_cluster_cell_idents(splg_enzymes_only, splg_labeled)
splg_smts_cross_representation$cluster_origin = factor(splg_smts_cross_representation$cluster_origin, 
                                                       levels = splg_smts_cross_representation$cluster_origin %>% unique())
splg_smts_cross_representation = splg_smts_cross_representation %>% mutate("Cell Identity"=cluster_num, percent_representation=percent_representation*100)
cell_ident_names=Idents(splg_labeled) %>% levels()

#FIG2 free scale 2800X2000
ggplot(data=splg_smts_cross_representation) +
  geom_bar(mapping = aes(x=cluster_origin, y=percent_representation), stat = "identity") +
  # geom_text(aes(label = "hello"), vjust = 0) +
  facet_wrap(~cluster_num, nrow = 2, ncol = 4, scales = "free") +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  
  ylab("Percentages of IP cluster cells \n mapped to SMT AP cluster cells") +
  xlab("SMT AP Cluster") +
  theme(strip.text = element_text(face="bold", size=22),
        strip.background = element_rect(colour="black",size=1),
        axis.title = element_text(face="bold", size=22),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18)) +
  geom_text(aes(x=cluster_origin, y=percent_representation, label=paste0(num_cells)), angle=90, size = 6, hjust=-0.1)

# #FIG2 fixed scale
# ggplot(data=splg_smts_cross_representation) +
#   geom_bar(mapping = aes(x=cluster_origin, y=percent_representation), stat = "identity") +
#   # geom_text(aes(label = "hello"), vjust = 0) +
#   facet_wrap(~cluster_num, nrow = 2, ncol = 4) +
#   scale_y_continuous(expand = expansion(mult = c(0, .1))) +
#   
#   ylab("Percentages of IP cluster cells \n mapped to SMT AP cluster cells") +
#   xlab("SMT AP Cluster") +
#   theme(strip.text = element_text(face="bold", size=18),
#         strip.background = element_rect(colour="black",size=1),
#         axis.title = element_text(face="bold", size=18),
#         axis.text.x = element_text(size=10),
#         axis.text.y = element_text(size=10)) +
#   geom_text(aes(x=cluster_origin, y=percent_representation, label=paste0(num_cells, " cell(s)")), angle=90, size = 4, hjust=-0.05)
# 

cross_representation_plot(splg_enzyme_cross_representation)
cross_representation_plot(splg_smts_cross_representation)
cross_representation_plot(splg_exo_cross_representation)
cross_representation_plot(splg_mt_cross_representation)
cross_representation_plot(splg_tfs_cross_representation)

# cross_representation_plot(splg_cross_representation_with_labeled_clusters)

# skeletal: post cross representation analysis


skeletal_smt_ap_cells = subset(splg_smts_only, idents = c("9","10")) %>% Cells()
splg_mapped_skeletal_smt_ap_cells = subset(splg, cells = skeletal_smt_ap_cells)

splg_mapped_skeletal_smt_ap_cells = subcluster_suite(
  seurat_object = splg_mapped_skeletal_smt_ap_cells,
  res = 0.1,
  ndims = 1:4
)

skeletal_geneset = tibble(
  
  GeneID = c(
    # canonical PMC markers set 1
    "SM50", "SM30A", "msp130", "LOC579101", "LOC579173", 
    
    # highly enriched markers found 
    "LOC576560", "LOC576018",
    # canonical PMC markers set 2
            "LOC583799", 
    # PMC and SMC markers
            "SFK7","LOC574611", "SPARC", "LOC752471",

    # differentially expressed collagens
            "LOC589794", "COLP1alpha", "COLP2alpha","COLP3alpha", "COLP4alpha", 
    # IP Markers
    "LOC100887885", "LOC100892373", "LOC575826", "Alx1", "LOC580069", "GATAc", "erg",
    # SMTs which have defined these two primary subpopulation of skeletal cells
            "LOC591586", "LOC585510", "LOC592979", "LOC587674"),
  
  Name = c(
    # canonical PMC markers set 1
     "SM50", "SM30A", "msp130", "carbonic anhydrase 2", "otop2l", 
     
     # highly enriched markers found 
      "neurexin-4", "matrilysin (MMP8)",
    # canonical PMC markers set 2
            "MMP18",
    # PMC and SMC markers
           "SFK7", "cofilin", "SPARC",  "calumenin", 
    
    # differentially expressed collagens
           "collagen alpha-1(V) chain", "COLP1alpha", "COLP2alpha", "COLP3alpha", "COLP4alpha", 
    # IP Markers
    "Kirrel2L_4", "KirreL2L_6", "hex", "Alx1", "Alx4","GATAc","erg",
    # SMTs which have defined these two primary subpopulation of skeletal cells
           "SLC9A2", "SLC13A2", "SLC5A11", "SLC26A5")
  
)

labeled_dotplot(
  scrna_df = splg_mapped_skeletal_smt_ap_cells,
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


