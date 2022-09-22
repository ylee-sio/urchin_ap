all_stages = readRDS("data/working_data/all_stages.clustered_seurat_obj_list_positive_features_filtered.rds")


all_stages[[1]]$orig.ident
sp_kegg_smts[sp_kegg_smts$Name %>% str_which("SLC26"),]
DimPlot(all_stages[[4]])
FeaturePlot(all_stages[[8]], features = c("delta", "COLP3alpha", "Alx1", "erg", "LOC576018", "COLP1alpha"), order = T)
FeaturePlot(all_stages[[4]], features = c("gcm", "FoxA", "Endo16", "blimp1/krox"))

FeaturePlot(all_stages[[5]], features = c("gcm", "blimp1/krox"), blend = T, order = T)
FeaturePlot(all_stages[[5]], features = c("gcm", "delta"), blend = T, order = T)
FeaturePlot(all_stages[[5]], features = c("gcm", "FoxA"), blend = T, order = T)

                                          
sp4_markers = FindAllMarkers(all_stages[[4]], only.pos = T)
sp4_markers %>% group_by(cluster) %>% top_n(10) %>% View()
DimPlot(all_stages[[4]])
