# Run this script in terminal using 'R CMD BATCH obtain_working_data.R'
# This script only needs to be run one-time to obtain working data that has been processed with the parameters below. Adjust accordingly.

# setting workers for parallel operations
plan(multicore, workers = 8)

# make a list of paths for each directory (developmental stage specific) containing the matrix, feature, and barcode files 
indv_stage_dirs = list.files("data/unmodified", full.names = T)

# creates a list of matrices (developmental stage specific) of expression data to be used to create seurat_objects
all_stages.data = future_map(
  .x = indv_stage_dirs, 
  .f = function(x)(
    Read10X(
      data.dir = x,
    )
  )
)
saveRDS(all_stages.data, "data/unmodified/all_stages.data.rds")  

# creates a list of seurat_objects for each stage using list of scrna matrices from all_stages.data
all_stages.seurat_obj_list = future_map(
  .x = all_stages.data,
  .f = function(x)(
    CreateSeuratObject(
      counts = x
      )
  )
)
saveRDS(all_stages.seurat_obj_list, "data/unmodified/all_stages.seurat_obj_list.rds")

# creates a list of normalized, clustered seurat_objects using a list of unclustered seurat_objects from all_stages.seurat_obj_list
all_stages.clustered_seurat_obj_list = future_map(
  .x = all_stages.seurat_obj_list,
  .f = function(x)(
    subcluster_suite(
      seurat_object = x,
      res = 0.8,
      ndims = 1:15,
      strat = "tsne",
      jackstraw = F,
      coord_strat = "pca"
    )
  )
)
saveRDS(all_stages.clustered_seurat_obj_list, "data/unmodified/all_stages.clustered_seurat_obj_list.rds")


unmod_integ = readRDS("sp_transportomics/data_sources/primary/scrna/GSE149221_integrated/GSE149221_SpInteg.rds")


standard_markers_for_display = 
  tibble(
    GeneID = c("ARS", "SPEC1", "spec2d", "NK2.2", "Klf7", 
               "egf1", "ebr1", "CHRD", "FoxG", "Lim1", "foxq2", 
               "NK2.1", "AnkAT-1", "Endo16", "FoxA", "GATAe", "blimp1/krox",
               "Six1", "gcm", "rab3", "LOC588806", "COLP3alpha","SpP19L", "SM50", 
               "SM37", "msp130",  "erg", "adv", "Alx1","LOC100887863", "LOC373434", "Nanos2", 
               "GCF1", "HDAC1", "NSF", "pabp", "vasa",
               "LOC592057", "LOC373476", "LOC577601", "delta", "VAMP"),
    
    Name = c("ARS", "SPEC1", "spec2d", "NK2.2", "Klf7", 
             "egf1", "ebr1", "CHRD", "FoxG", "Lim1", "foxq2", 
             "NK2.1", "AnkAT-1", "Endo16", "FoxA", "GATAe", "blimp1/krox",
             "Six1", "gcm", "rab3", "PKS1", "COLP3alpha","SpP19L" ,"SM50", 
             "SM37", "msp130", "erg", "adv", "Alx1", "tdrd1", "seawi", "Nanos2",
             "GCF1", "HDAC1", "NSF", "pabp", "vasa",
             "HesC", "syntaxin", "brn1/2/4", "delta", "VAMP")
  )

abc_slc_3_geneid_final_cleaned = read_csv("/Users/yl/sp_transportomics/data_sources/secondary/pm_abcs_slcs_merged.csv")

splg_unmodified.data = Read10X(data.dir = "~/sp_transportomics/data_sources/primary/scrna/GSE149221_RAW/SpLG")
splg_unmodified = CreateSeuratObject(counts = splg_unmodified.data, project = "Late\ngastrula")
urchin_5_splg_clustered_unmodified = subcluster_suite(seurat_object = splg_unmodified, res = 1.2, ndims = 1:20, strat = "tsne", jackstraw = F, coord_strat = "pca")

splg.data = Read10X(data.dir = "~/sp_transportomics/data_sources/primary/scrna/GSE149221_RAW_modified/SpLG")
splg = CreateSeuratObject(counts = splg.data, project = "Late\ngastrula")
urchin_5_splg_clustered = subcluster_suite(seurat_object = splg, res = 1.2, ndims = 1:20, strat = "tsne", jackstraw = F, coord_strat = "pca")

speg.data = Read10X(data.dir = "~/sp_transportomics/data_sources/primary/scrna/GSE149221_RAW_modified/SpEG")
speg = CreateSeuratObject(counts = speg.data, project = "Early\ngastrula")
urchin_5_speg_clustered = subcluster_suite(seurat_object = speg, res = 1.2, ndims = 1:20, strat = "tsne", jackstraw = F,coord_strat = "pca")

spmb.data = Read10X(data.dir = "~/sp_transportomics/data_sources/primary/scrna/GSE149221_RAW_modified/SpMB")
spmb = CreateSeuratObject(counts = spmb.data, project = "Mesenchymal\nblastula")
urchin_5_spmb_clustered = subcluster_suite(seurat_object = spmb, ndims = 1:20, strat = "tsne", res = 1.2, jackstraw = F,coord_strat = "pca")

sphb.data = Read10X(data.dir = "~/sp_transportomics/data_sources/primary/scrna/GSE149221_RAW_modified/SpHB")
sphb = CreateSeuratObject(counts = sphb.data, project = "Hatched\nblastula")
urchin_5_sphb_clustered = subcluster_suite(seurat_object = sphb, ndims = 1:20, strat = "tsne", res = 1.2, jackstraw = F, coord_strat = "pca")

speb.data = Read10X(data.dir = "~/sp_transportomics/data_sources/primary/scrna/GSE149221_RAW_modified/SpEB")
speb = CreateSeuratObject(counts = speb.data, project = "Early\nblastula")
urchin_5_speb_clustered = subcluster_suite(speb, res = 1.2, ndims = 1:20, strat = "tsne", jackstraw = F, coord_strat = "pca")

sp3.data = Read10X(data.dir = "~/sp_transportomics/data_sources/primary/scrna/GSE149221_RAW_modified/Sp3")
sp3 = CreateSeuratObject(counts = sp3.data, project = "Morula")
urchin_5_sp3_clustered = subcluster_suite(sp3, res = .4, ndims = 1:5, jackstraw = F, strat = "tsne", coord_strat = "pca")

sp2.data = Read10X(data.dir = "~/sp_transportomics/data_sources/primary/scrna/GSE149221_RAW_modified/Sp2")
sp2 = CreateSeuratObject(counts = sp2.data, project = "64 Cell")
urchin_5_sp2_clustered = subcluster_suite(sp2, res = .4, ndims = 1:5, jackstraw = F, strat = "tsne",coord_strat = "pca")

sp1.data = Read10X(data.dir = "~/sp_transportomics/data_sources/primary/scrna/GSE149221_RAW_modified/Sp1")
sp1 = CreateSeuratObject(counts = sp1.data, project = "8 Cell")
urchin_5_sp1_clustered = subcluster_suite(sp1, res = .4, ndims = 1:15, strat = "tsne", jackstraw = F,coord_strat = "pca")
# 
# sp48gcm_mo.data = Read10X(data.dir = "~/sp_transportomics/data_sources/primary/scrna/GSE155427_RAW_modified//sp48gcm_mo/")
# sp48gcm_mo = CreateSeuratObject(counts = sp48gcm_mo.data, project = "48HPF GCM MO ")
# urchin_5_sp48gcm_mo_clustered = subcluster_suite(sp48gcm_mo, res = 1.2, ndims = 1:20, strat = "tsne", jackstraw = F,coord_strat = "pca")
# DimPlot(urchin_5_sp48gcm_mo_clustered)
# 
# sp48gcm_con_mo.data = Read10X(data.dir = "~/sp_transportomics/data_sources/primary/scrna/GSE155427_RAW/sp48con_mo/")
# sp48gcm_con_mo = CreateSeuratObject(counts = sp48gcm_con_mo.data, project = "48HPF GCM MO ")
# urchin_5_sp48gcm_con_mo_clustered = subcluster_suite(sp48gcm_con_mo, res = 1.2, ndims = 1:20, strat = "tsne", jackstraw = F,coord_strat = "pca")
# DimPlot(urchin_5_sp48gcm_con_mo_clustered)

named_locus_list_a = rownames(urchin_5_splg_clustered_unmodified)[which(str_sub(rownames(urchin_5_splg_clustered_unmodified), 1,3) != "LOC")]
named_locus_list_b = named_locus_list_a[which(str_sub(named_locus_list_a, 1,3) != "trn")]
named_locus_list_c = named_locus_list_b[which(str_sub(named_locus_list_b, 1,3) != "Mir")]
named_locus_list_d = named_locus_list_c[which(str_sub(named_locus_list_c, 1,3) != "mer")] %>% sort()
named_locus_df = tibble(
  
  GeneID = named_locus_list_d,
  Name = named_locus_list_d
  
) 
spuranno = readGFF("~/sp_transportomics/data_sources/primary/spur4.2_annotations_echinobase/Spur4_2_models.gff3")
sp_kegg_smts = parse_kegg_brite("~/sp_transportomics/data_sources/primary/kegg_genesets/sp_smts.json")
sp_kegg_tfs = parse_kegg_brite("~/sp_transportomics/data_sources/primary/kegg_genesets/sp_tfs.json")
sp_kegg_mem_traf = parse_kegg_brite("~/sp_transportomics/data_sources/primary/kegg_genesets/sp_membrane_trafficking.json")
sp_kegg_gpcr = parse_kegg_brite("~/sp_transportomics/data_sources/primary/kegg_genesets/sp_gpcr.json")
sp_kegg_enzymes = parse_kegg_brite("~/sp_transportomics/data_sources/primary/kegg_genesets/sp_enzymes.json")
sp_kegg_prr = parse_kegg_brite("~/sp_transportomics/data_sources/primary/kegg_genesets/sp_pattern_recog_recptors.json")
sp_kegg_nuclear_receptors = parse_kegg_brite("~/sp_transportomics/data_sources/primary/kegg_genesets/sp_nuclear_receptors.json")
sp_kegg_translation_factors = parse_kegg_brite("~/sp_transportomics/data_sources/primary/kegg_genesets/sp_translations_factors.json")


urchin_5_keys_file_list = list.files("~/sp_transportomics/data_sources/secondary/merged_paralogs_key/",full.names = T)
urchin_5_keys = map(urchin_5_keys_file_list, read_csv)
# urchin_5_keys = urchin_5_keys[[1]]
urchin_5_clustered_seurat_objs = c(
  urchin_5_sp1_clustered,
  urchin_5_sp2_clustered,
  urchin_5_sp3_clustered,
  urchin_5_speb_clustered,
  urchin_5_sphb_clustered,
  urchin_5_spmb_clustered,
  urchin_5_speg_clustered,
  urchin_5_splg_clustered
)

# urchin_5_clustered_labelled_seurat_objs = c(
# 
#   urchin_5_sp1_clustered,
#   urchin_5_sp2_clustered,
#   urchin_5_sp3_clustered,
#   urchin_5_speb_cell_type_labelled,
#   urchin_5_sphb_cell_type_labelled,
#   urchin_5_spmb_cell_type_labelled,
#   urchin_5_speg_cell_type_labelled,
#   urchin_5_splg_cell_type_labelled
# 
# )

stage_names = c("8 cell stage",
                "64 cell stage",
                "Morula stage",
                "Early blastula",
                "Hatched blastula",
                "Mesenchymal blastula",
                "Early gastrula",
                "Late gastrula"
                )

available_genesets = urchin_5_keys[[1]]$GeneSet %>% unique()

urchin_5_smts = 
  urchin_5_keys[[1]] %>% 
  subset(
    GeneSet == "smts"
  ) %>% 
  select(
    GeneID,
    Name
  )

urchin_5_smts_final_cleaned = 
  left_join(
    abc_slc_3_geneid_final_cleaned,
    urchin_5_smts,
    by = "Name"
  ) %>% 
  transmute(
    GeneID = GeneID.y,
    Name = Name)

urchin_5_smts_final_cleaned_with_mito_smts = 
  bind_rows(
    urchin_5_smts[urchin_5_smts$Name %>% str_which("25"),],
    urchin_5_smts_final_cleaned
  )

germline_markers = tibble(
  GeneID = c("LOC373434", "LOC100887863", "Nanos2"),
  Name = c("seawi", "tdrd1", "Nanos2")
)
neural_markers = c("foxq2", "AnkAT-1", "NK2.1")
pigment_markers = c("gcm", "Six1")
skeletal_markers = c("Alx1", "erg", "SM50", "SM30A")
oral_ectoderm_markers = c("CHRD", "FoxG", "Lim1")
endoderm_markers = c("Endo16", "blimp1/krox")
aboral_ectoderm_markers = c("SPEC1", "spec2c", "NK2.2", "Tbx2/3", "Klf7")

markers = list(
  
  germline_markers,
  neural_markers,
  pigment_markers,
  skeletal_markers
  # oral_ectoderm_markers,
  # endoderm_markers,
  # aboral_ectoderm_markers
  
)

marker_names_list = list(
  
  "Germline", 
  "Neural",
  "Pigment",
  "Skeletal",
  "Oral ectoderm",
  "Endoderm",
  "Aboral ectoderm"
  
  
  )

