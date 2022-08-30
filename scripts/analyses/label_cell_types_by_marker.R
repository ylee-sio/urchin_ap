urchin_5_splg_clustered_test = urchin_5_splg_clustered

#aboral-ectoderm
nk2.2_pos_cells = get_pos_cells(urchin_5_splg_clustered_test, gene_id = "NK2.2", thresh = 3, rev = F)
klf7_pos_cells = get_pos_cells(urchin_5_splg_clustered_test, gene_id = "Klf7", thresh = 2, rev = F)
ars_pos_cells = get_pos_cells(urchin_5_splg_clustered_test, gene_id = "ARS", thresh = 36, rev = F)
aboral_ectoderm_markers_pos_cells = bind_rows(nk2.2_pos_cells, klf7_pos_cells, ars_pos_cells)
Idents(urchin_5_splg_clustered_test, cells = aboral_ectoderm_markers_pos_cells$cells) = "Aboral ectoderm"

#endoderm
endo16_pos_cells = get_pos_cells(urchin_5_splg_clustered_test, gene_id = "Endo16", thresh = 15, rev = F)
blimp1_krox_pos_cells = get_pos_cells(urchin_5_splg_clustered_test, gene_id = "blimp1/krox", thresh = 10, rev = F)
endoderm_markers_pos_cells = bind_rows(endo16_pos_cells, blimp1_krox_pos_cells)
Idents(urchin_5_splg_clustered_test, cells = endoderm_markers_pos_cells$cells %>% unique()) = "Endoderm"

#neural
foxq2_pos_cells = get_pos_cells(urchin_5_splg_clustered_test, gene_id = "foxq2", thresh = 3, rev = F)
ankat1_pos_cells = get_pos_cells(urchin_5_splg_clustered_test, gene_id = "AnkAT-1", thresh = 1, rev = F)
neural_markers_pos_cells = bind_rows(foxq2_pos_cells, ankat1_pos_cells)
Idents(urchin_5_splg_clustered_test, cells = neural_markers_pos_cells$cells %>% unique()) = "Neural"

#pigment
gcm_pos_cells = get_pos_cells(urchin_5_splg_clustered_test, gene_id = "gcm", thresh = 1, rev = F)
six1_pos_cells = get_pos_cells(urchin_5_splg_clustered_test, gene_id = "Six1", thresh = 1, rev = F)
abcc5_pos_cells = get_pos_cells(urchin_5_splg_clustered_test, gene_id = "merged-gene-165", thresh = 1, rev = F)
pigment_markers_pos_cells = bind_rows(gcm_pos_cells, six1_pos_cells, abcc5_pos_cells)
Idents(urchin_5_splg_clustered_test, cells = pigment_markers_pos_cells$cells %>% unique) = "Pigment"

#skeletal
alx1_pos_cells = get_pos_cells(urchin_5_splg_clustered_test, gene_id = "Alx1", thresh = 1, rev = F)
erg_pos_cells = get_pos_cells(urchin_5_splg_clustered_test, gene_id = "erg", thresh = 1, rev = F)
delta_pos_cells = get_pos_cells(urchin_5_splg_clustered_test, gene_id = "delta", thresh = 1, rev = F)
sm50_pos_cells = get_pos_cells(urchin_5_splg_clustered_test, gene_id = "SM50", thresh = 1, rev = F)
skeletal_markers_pos_cells = bind_rows(alx1_pos_cells, erg_pos_cells, delta_pos_cells, sm50_pos_cells)
Idents(urchin_5_splg_clustered_test, cells = skeletal_markers_pos_cells$cells %>% unique) = "Skeletal"

#germline
nanos2_pos_cells = get_pos_cells(urchin_5_splg_clustered_test, gene_id = "Nanos2", thresh = 0, rev = F)
Idents(urchin_5_splg_clustered_test, cells = nanos2_pos_cells$cells) = "Germline"

# #mtb1 enriched
# mtb1_pos_cells = get_pos_cells(urchin_5_splg_clustered_test, gene_id = "MTB1", thresh = 1.5, rev = F)
# Idents(urchin_5_splg_clustered_test, cells = mtb1_pos_cells$cells) = "MTB1-enriched"
# 
# #mvp enriched
# mvp_pos_cells = get_pos_cells(urchin_5_splg_clustered_test, gene_id = "mvp", thresh = 1.5, rev = F)
# Idents(urchin_5_splg_clustered_test, cells = mvp_pos_cells$cells) = "mvp-enriched"
get_cell_stats(urchin_5_splg_clustered_test)
DimPlot(urchin_5_splg_clustered_test, label = T) %>% ggplotly()
