library(rtracklayer)
library(tidyverse)
library(Seurat)
# Run this script in terminal using 'R CMD BATCH obtain_working_data.R'
# This script only needs to be run one-time to obtain working data that has been processed with the parameters below. Adjust accordingly.

sp_kegg_gpcr = read_csv("data/working_data/sp_kegg_gpcr.csv")

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

named_locus_list_a = rownames(splg)[which(str_sub(rownames(splg), 1,3) != "LOC")]
named_locus_list_b = named_locus_list_a[which(str_sub(named_locus_list_a, 1,3) != "trn")]
named_locus_list_c = named_locus_list_b[which(str_sub(named_locus_list_b, 1,3) != "Mir")]
named_locus_list_d = named_locus_list_c[which(str_sub(named_locus_list_c, 1,3) != "mer")] %>% sort()
named_locus_df = tibble(
  
  GeneID = named_locus_list_d,
  Name = named_locus_list_d
  
) 
spuranno = readGFF("~/Projects/urchin_ap/data/source/sp5.0.gff3") %>% as_tibble()

# contains 32,948 records. however, does not include transcript isoforms which are nested under the same gene id
spuranno_type_gene = spuranno %>% subset(type=="gene")

# spuranno_type_gene will not contain give prodcut information
# subsetting "type" by mrna, will allow us to make a dataframe with product information
# which is the only feature that seems to contain transcript variant information
# however, subsetting "type" my mrna gives a longer df because mrnas will have transcript variant rows. 
# rows with transcript variants will have the same gene id, but will have different transcript ids

# contains 38,427 mrna records; does not mean 38,427 genes
spuranno_type_mrna = spuranno %>% subset(type=="mRNA")

# confirmation that any gene id that is in mrna df, is in gene df
spuranno_type_mrna$gene %in% spuranno_type_gene$gene %>% which() %>% length
spuranno_type_mrna$gene %>% length

# subset for rows containing "transcript variant" 
# in the "product" variable to get all genes which have any transcript and all their respective transcripts.
# to do this, first get an index of all the rows that contain this information.
# it appears this information is only available as text in the "product" variable
# no separate column/variable is used to contain transcript variant information
spuranno_type_mrna_tv_only_index = spuranno_type_mrna$product %>% str_which("transcript variant")
spuranno_type_mrna_nontv_only_index = spuranno_type_mrna$product %>% str_which("transcript variant", negate = T)

# quickly check lengths of resulting dfs. if using a search term to pick out transcript variants is the way to do this,
# which doesn't seem ideal, the sum of the lengths of both indices should add up to the total number of rows of 
# spuranno_type_mrna
spuranno_type_mrna_tv_only_index %>% length
spuranno_type_mrna_nontv_only_index %>% length
(spuranno_type_mrna_tv_only_index %>% length) + (spuranno_type_mrna_nontv_only_index %>% length)
spuranno_type_mrna %>% nrow

# the lengths seem to confirm that exclusive portions of the df are being captured
# 16,566/38,427 records are variants
# TVs have the same gene id. however, the language is a bit tricky... 
# might the product column for a gene NOT have "transcript variant" in the text if it is the "first" of the variants?
# to check for this, check if any of the gene ids in spuranno_type_mrna_nontv_only_index are in spuranno_type_mrna_tv_only_index
tv_genes = spuranno_type_mrna$gene[spuranno_type_mrna_nontv_only_index] 
nontv_genes = spuranno_type_mrna$gene[spuranno_type_mrna_tv_only_index]

# it seems like there are some matches. 
tv_nontv_overlap_genes = nontv_genes[non_tv_genes %in% tv_genes %>% which()] %>% unique()
tv_nontv_overlap_genes_full_info = filter(spuranno_type_mrna, gene %in% tv_nontv_overlap_genes) 
tv_nontv_overlap_genes_product = tv_nontv_overlap_genes_full_info %>% select(product)

# it seems that there are indeed some genes that have isoforms but there is a record of the gene whose 
# product is not annotated to be an isoform... despite having the same gene ID
# it seems this is the case with manually annotated/curated genes
# to solve this issue, we can simply add "transcript variant" into the product column for records that do not have this
# first, locate in spuranno_type_mrna the records for these genes, and obtain unique record identifier for these records.
# in this case, it would be "transcript_id"
tv_nontv_overlap_genes_records_transcript_id= subset(
  spuranno_type_mrna, 
  gene %in% 
    tv_nontv_overlap_genes & 
    str_detect(
      product, 
      "transcript variant", 
      negate = T
      )
  ) %>% 
  select("transcript_id")

# obtain indices of these transcript ids
tv_nontv_overlap_genes_records_transcript_id_index = spuranno_type_mrna$transcript_id %in% tv_nontv_overlap_genes_records_transcript_id$transcript_id %>% which()

# finally, add the product annotation to spuranno_type_mrna for the above found genes
# then, among these, locate genes whose "product" does not contain "transcript variant" and add "transcript variant X0"
spuranno_type_mrna[tv_nontv_overlap_genes_records_transcript_id_index,]$product = paste0(spuranno_type_mrna[tv_nontv_overlap_genes_records_transcript_id_index,]$product, " transcript variant 0")

# remember that we found 16,566/38,427 records to have variants. if everything has been processed as planned,
# we should see in spuranno_type_mrna that there are an additional 134, for a total of 167000 records with transcript variants.
spuranno_type_mrna$product %>% str_which("transcript variant") %>% length()

# all of the transcript variant containing genes have been accounted for.
# now we can count the things we need to count:

# number of unique gene IDs
spuranno_type_mrna$gene %>% unique() %>% length()

# number of unique genes with genes that have transcript isoforms counted as single genes
genes_with_tv_count = spuranno_type_mrna[spuranno_type_mrna$product %>% str_which("transcript variant"),] %>% select(gene) %>% unique() %>% nrow()
genes_without_tv_count = spuranno_type_mrna[spuranno_type_mrna$product %>% str_which("transcript variant", negate = T),] %>% select(gene) %>% unique() %>% nrow()
genes_with_tv_count + genes_without_tv_count

# number of unique transcripts
spuranno_type_mrna$ID %>% unique() %>% length()

spuranno_transporters = spuranno_type_mrna[spuranno_type_mrna$gene %in%  sp_kegg_smts$GeneID %>% which(),]

# number of unique transporter genes
spuranno_transporters$gene %>% unique() %>% length()

# number of total transporter transcripts and all transcript variants
spuranno_transporters$transcript_id %>% unique() %>% length()

# spuranno_transporters only gene ids
spuranno_transporters_geneid_only = spuranno_transporters$gene

# number of transporters each stage of scrnaseq dataset = 722
smts_in_scrnaseq = tibble(
  GeneID = rownames(all_stages[[1]])[rownames(all_stages[[1]]) %in% spuranno_transporters_geneid_only],
  Name = rownames(all_stages[[1]])[rownames(all_stages[[1]]) %in% spuranno_transporters_geneid_only]
)

smts_in_scrnaseq_with_common_names = left_join(smts_in_scrnaseq, sp_kegg_smts, by = "GeneID") %>% unique()

smts_in_scrnaseq_with_common_names_cleaned = 
  subset(
    smts_in_scrnaseq_with_common_names, 
    str_detect(Name.y, "transporter") |
      str_detect(Name.y, "SLC[:digit:]*[:alpha:]*[:digit:]*") |
      str_detect(Name.y,"ABC[:digit:]*[:alpha:]*[:digit:]*") |
      str_detect(Name.y,"solute") |
      str_detect(Name.y,"carrier") 
  ) %>% 
  transmute(GeneID = GeneID, Name = Name.y)
# cluster analysis on all stages separately with the same resolution and number of dimenions
all_stages_clustered = 
  map(
    .x = all_stages, 
    .f = function(x)(
      subcluster_suite(
        seurat_object = x,
        res = 1.0,
        ndims = 1:25
        )
      )
    )

top_markers_all_stages = 
  map(
    .x = all_stages_clustered,
    .f = function(x)(
      FindAllMarkers(
        object = x,
        logfc.threshold = 1,
        only.pos = T,
        min.pct = 0.25
        )
      )
    )

smt_only_all_stages_clustered = 
  map(
    .x = all_stages_clustered,
    .f = function(x)(
      isolate_geneset_in_clusters(
        GetAssayData(object = x, slot = "count"),
        smts_in_scrnaseq
        )
    )
  )

smt_only_all_stages_reclustered = 
  map(
    .x = smt_only_all_stages_clustered, 
    .f = function(x)(
      subcluster_suite(
        seurat_object = x,
        res = 1.0,
        ndims = 1:25
      )
    )
  )

smt_only_top_markers_all_stages = 
  map(
    .x = smt_only_all_stages_reclustered,
    .f = function(x)(
      FindAllMarkers(
        object = x,
        logfc.threshold = 1,
        only.pos = T,
        min.pct = 0.25
      )
    )
  )





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
sp_kegg_smts = sp_kegg_smts %>% unique()
slc4a = gsqs(sp_kegg_smts, "SLC4A") %>% unique() %>% add_row(GeneID = "COLP3alpha", Name = "COLP3alpha")
(all_stages[[7]] %>% GetAssayData(slot = "counts"))["erg",] %>% sort(decreasing = T) %>% head(20)
(all_stages[[7]] %>% GetAssayData(slot = "counts"))["erg",] %>% sum()
FeaturePlot(all_stages[[8]], features = c("erg"), order = T)

tt = tibble(GeneID = "LOC591481", Name = "SLC22A4/5")

c19_markers = FindMarkers(all_stages[[8]], ident.1 = "19")
c19_markers_2 = c19_markers %>% arrange(desc(avg_log2FC))
c19_markers_2 %>% head(15)

labeled_dotplot(all_stages[[8]], tt, plot_title = "", by_stage = F, plot_height = 500, interactive = T, col.min = -1, col.max = 5, dot.min = 0, cols = c("blue", "gold"))
